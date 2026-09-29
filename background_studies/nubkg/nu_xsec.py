"""
nu_xsec.py -- neutrino-nucleon cross sections and inelasticity, table-first.

Design
------
The cross section is a *table* of (E_nu, sigma) points plus a log-log
interpolator, exactly the workflow of digitising a published figure.  Nothing
in this module hard-codes a functional form as the primary source of truth.

Two built-in tables are provided (see `DEFAULT_CC_NU` / `DEFAULT_CC_NUBAR`) so
the code runs out of the box, but the intent is that you replace them with your
own digitisation:

    sigma = CrossSection.from_csv("data/xsec_cc_nu_mydigitisation.csv",
                                  label="Formaggio&Zeller Fig 9 + IceCube Fig 1")

Provenance of the built-in tables
---------------------------------
E_nu < 350 GeV
    Measured world average of the CC inclusive cross section per nucleon on an
    isoscalar target.  PDG "Neutrino Cross Section Measurements" review quotes
    sigma_CC/E = 0.677e-38 cm^2/GeV (nu) and 0.334e-38 (nubar) in the linear
    DIS-scaling regime; the mild turn-on below ~30 GeV follows the shape in
    Formaggio & Zeller, arXiv:1305.7513, Fig. 9.
        https://arxiv.org/abs/1305.7513
        https://pdg.lbl.gov/2025/reviews/rpp2025-rev-nu-cross-sections.pdf

E_nu = 1e4 GeV
    sigma_CC(nu) = 4.5-4.6e-35 cm^2 from the NLO QCD calculations
    (Gandhi-Quigg-Reno-Sarcevic 1998; Cooper-Sarkar-Mertsch-Sarkar 2011).
    These are the curves plotted in IceCube's Nature 551, 596 (2017)
    cross-section measurement, arXiv:1711.08119, Fig. 1.
        https://arxiv.org/abs/1106.3723   (CSMS)
        https://arxiv.org/abs/hep-ph/9807264  (GQRS)
        https://arxiv.org/abs/1711.08119  (IceCube measurement)

E_nu > 1e4 GeV
    sigma ~ E^0.363, the standard high-energy behaviour of those calculations,
    anchored on the 1e4 GeV value.

Entries flagged `interp` in the CSV comments are log-log bridges between the
measured linear regime and the 1e4 GeV anchor, i.e. they carry the largest
model dependence (~10-15%).  For EarthShine this region contributes little,
but see `sigma_uncertainty()`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np

Species = Literal["nu", "nubar"]

__all__ = ["CrossSection", "Inelasticity", "default_cc", "builtin_cc",
           "load_cc_pair", "set_default_xsec", "get_default_xsec",
           "sigma_uncertainty", "DEFAULT_XSEC_FILE"]


# ---------------------------------------------------------------------------
# Built-in tables:  (E_nu / GeV, sigma_CC per nucleon / cm^2)
# ---------------------------------------------------------------------------

def _linear_regime(e, slope):
    """slope * E, with the sub-30 GeV turn-on of Formaggio & Zeller Fig 9."""
    e = np.asarray(e, float)
    # empirical turn-on: sigma/E rises to the plateau by ~30 GeV
    return slope * e * (1.0 - 0.18 * np.exp(-e / 12.0))


def _build_default(slope, sigma_1e4, index_hi=0.363):
    lo = np.geomspace(1.0, 350.0, 25)
    hi = np.geomspace(1.0e4, 1.0e8, 25)
    mid = np.geomspace(350.0 * 1.2, 1.0e4 / 1.2, 12)

    s_lo = _linear_regime(lo, slope)
    s_hi = sigma_1e4 * (hi / 1.0e4) ** index_hi
    # log-log bridge
    x = np.log(mid)
    x0, x1 = np.log(350.0), np.log(1.0e4)
    y0, y1 = np.log(_linear_regime(350.0, slope)), np.log(sigma_1e4)
    s_mid = np.exp(y0 + (y1 - y0) * (x - x0) / (x1 - x0))

    e = np.concatenate([lo, mid, hi])
    s = np.concatenate([s_lo, s_mid, s_hi])
    return e, s


# nu: anchored at sigma_CC(1e4 GeV) = 4.60e-35 cm^2
DEFAULT_CC_NU = _build_default(0.677e-38, 4.60e-35)
# nubar: converges to the nu value at high energy (valence contribution dies)
DEFAULT_CC_NUBAR = _build_default(0.334e-38, 4.15e-35)


# ---------------------------------------------------------------------------
# Interpolator
# ---------------------------------------------------------------------------

@dataclass
class CrossSection:
    """Log-log interpolated cross section, cm^2 per nucleon.

    Parameters
    ----------
    energy : array of E_nu in GeV, strictly increasing
    sigma  : array of sigma in cm^2
    label  : provenance string, carried into output metadata
    extrapolate : if False, raise outside the table range.  Default True with a
        log-log linear extrapolation using the end-point slope.
    """

    energy: np.ndarray
    sigma: np.ndarray
    label: str = "unlabelled"
    extrapolate: bool = True

    def __post_init__(self):
        self.energy = np.asarray(self.energy, float)
        self.sigma = np.asarray(self.sigma, float)
        if np.any(np.diff(self.energy) <= 0):
            order = np.argsort(self.energy)
            self.energy, self.sigma = self.energy[order], self.sigma[order]
        if np.any(self.sigma <= 0):
            raise ValueError("cross section table contains non-positive values")
        self._lx = np.log(self.energy)
        self._ly = np.log(self.sigma)

    @classmethod
    def from_csv(cls, path, label=None, energy_col=0, sigma_col=1,
                 energy_unit: Literal["GeV", "log10GeV"] = "GeV",
                 sigma_unit: Literal["cm2", "cm2_per_GeV"] = "cm2"):
        """Load a two-column CSV (comments with '#').

        energy_unit='log10GeV' handles tables digitised straight off a plot
        with a log10(E/GeV) axis.  sigma_unit='cm2_per_GeV' handles the very
        common sigma/E presentation (e.g. Formaggio & Zeller Fig 9).
        """
        arr = np.genfromtxt(path, delimiter=",", comments="#")
        e = arr[:, energy_col]
        s = arr[:, sigma_col]
        if energy_unit == "log10GeV":
            e = 10.0 ** e
        if sigma_unit == "cm2_per_GeV":
            s = s * e
        return cls(e, s, label=label or str(Path(path).name))

    def to_csv(self, path, header=""):
        with open(path, "w") as fh:
            fh.write(f"# {self.label}\n# {header}\n# E_nu[GeV],sigma[cm^2]\n")
            for e, s in zip(self.energy, self.sigma):
                fh.write(f"{e:.6e},{s:.6e}\n")

    def __call__(self, e_nu):
        x = np.log(np.asarray(e_nu, float))
        y = np.interp(x, self._lx, self._ly)
        if self.extrapolate:
            # log-log slope extrapolation off both ends
            lo = x < self._lx[0]
            hi = x > self._lx[-1]
            if np.any(lo):
                m = (self._ly[1] - self._ly[0]) / (self._lx[1] - self._lx[0])
                y = np.where(lo, self._ly[0] + m * (x - self._lx[0]), y)
            if np.any(hi):
                m = (self._ly[-1] - self._ly[-2]) / (self._lx[-1] - self._lx[-2])
                y = np.where(hi, self._ly[-1] + m * (x - self._lx[-1]), y)
        elif np.any((x < self._lx[0]) | (x > self._lx[-1])):
            raise ValueError("energy outside cross-section table range")
        return np.exp(y)

    @property
    def range_GeV(self):
        return float(self.energy[0]), float(self.energy[-1])

    def as_metadata(self, prefix="xsec"):
        return {f"{prefix}_label": self.label,
                f"{prefix}_e_min_GeV": self.range_GeV[0],
                f"{prefix}_e_max_GeV": self.range_GeV[1],
                f"{prefix}_n_points": len(self.energy)}


def builtin_cc(species: Species = "nu") -> CrossSection:
    """The original hand-built table, kept only for comparison.

    SUPERSEDED.  It extrapolated upward from the 10 TeV anchor with E^0.363,
    which is the *asymptotic* index (valid above ~1e7 GeV).  Between 10 TeV and
    1 PeV the true curve is steeper, so this table is 2x too low at 100 TeV and
    3x too low at 1 PeV.  Its nubar is also ~25% high at 10 TeV.
    """
    e, s = DEFAULT_CC_NU if species == "nu" else DEFAULT_CC_NUBAR
    lab = f"built-in placeholder (SUPERSEDED) [{species}]"
    return CrossSection(e, s, label=lab)


# ---------------------------------------------------------------------------
# Table-driven default
# ---------------------------------------------------------------------------

#: file in nubkg/data/ used by default_cc().  Change it per-notebook with
#: set_default_xsec("other_file.csv").
DEFAULT_XSEC_FILE = "xsec_cc_digitised.csv"

_DATA_DIR = Path(__file__).resolve().parent / "data"
_current_file = DEFAULT_XSEC_FILE
_cache: dict = {}


def load_cc_pair(path, label=None) -> dict:
    """Load a combined table  E_GeV, sigma_nu_cm2, sigma_nubar_cm2  into
    {"nu": CrossSection, "nubar": CrossSection}.  Accepts a bare file name
    (looked up in nubkg/data/) or a path."""
    p = Path(path)
    if not p.is_absolute() and not p.exists():
        p = _DATA_DIR / p
    key = (str(p.resolve()), p.stat().st_mtime)
    if key not in _cache:
        arr = np.genfromtxt(p, delimiter=",", comments="#")
        lab = label or p.name
        _cache[key] = {"nu": CrossSection(arr[:, 0], arr[:, 1], label=f"{lab} [nu]"),
                       "nubar": CrossSection(arr[:, 0], arr[:, 2], label=f"{lab} [nubar]")}
    return _cache[key]


def set_default_xsec(name_or_path) -> None:
    """Choose the cross-section table every calculation uses by default.

    One line in a notebook's setup cell:
        nb.set_default_xsec("xsec_cc_digitised.csv")
    Pass "builtin" to fall back to the old hand-built table.
    """
    global _current_file
    if name_or_path != "builtin":
        load_cc_pair(name_or_path)        # fail now, not mid-calculation
    _current_file = name_or_path


def get_default_xsec() -> str:
    return _current_file


def default_cc(species: Species = "nu") -> CrossSection:
    """CC cross section used whenever none is passed explicitly."""
    if _current_file == "builtin":
        return builtin_cc(species)
    try:
        return load_cc_pair(_current_file)[species]
    except (OSError, FileNotFoundError):
        import warnings
        warnings.warn(f"cross-section table {_current_file!r} not found in "
                      f"{_DATA_DIR}; falling back to the SUPERSEDED built-in "
                      "table, which is 2-3x too low above 100 TeV.",
                      stacklevel=2)
        return builtin_cc(species)


def sigma_uncertainty(e_nu) -> np.ndarray:
    """Fractional 1-sigma uncertainty on sigma_CC.  These are choices, set to
    match the provenance of the default table (xsec_cc_digitised.csv):

      < 30 GeV       10%  QE/resonance/DIS transition; generator prediction
      30 - 340 GeV    3%  directly measured world average (PDG)
      340 - 600 GeV   8%  the PDG/CSMS junction; the two differ by -5.6% (nu)
                          and +7.9% (nubar) at 350 GeV
      600 GeV - 1 PeV 5%  CSMS SM calculation plus ~1% digitisation error
      > 1 PeV        10%  beyond the digitised range (extrapolated)

    Not the IceCube measurement's error: that measured 1.30 x SM with a +-~45%
    total uncertainty over 6.3-980 TeV, which is a test of the SM rather than a
    better determination of it.  Use xsec_icecube_measured.csv to apply the
    measured scale as an explicit variation if wanted.
    """
    e = np.asarray(e_nu, float)
    frac = np.full_like(e, 0.05)
    frac = np.where(e < 30, 0.10, frac)
    frac = np.where((e >= 30) & (e < 340), 0.03, frac)
    frac = np.where((e >= 340) & (e < 600), 0.08, frac)
    frac = np.where(e > 1e6, 0.10, frac)
    return frac


# ---------------------------------------------------------------------------
# Inelasticity  y = 1 - E_mu / E_nu
# ---------------------------------------------------------------------------

@dataclass
class Inelasticity:
    """Distribution of the muon energy fraction u = 1 - y = E_mu / E_nu.

    model
      "qpm"  : quark-parton model.  dsigma/dy flat for nu (<y>=0.5) and
               3(1-y)^2 for nubar (<y>=0.25).  Self-consistent with the
               sigma_nubar/sigma_nu -> 1/3 valence limit in the cross sections.
      "flat" : flat for both species.
      "none" : E_mu = E_nu exactly.  The naive ansatz, kept for comparison.

    Only two quantities are ever needed, and both are analytic -- which is why
    no Monte Carlo over y is required anywhere in this package:

      pdf_u(u)      : density of u
      survival(z)   : P(u > z), used for the arriving-muon spectrum
    """

    model: Literal["qpm", "flat", "none"] = "qpm"

    def pdf_u(self, u, species: Species = "nu") -> np.ndarray:
        u = np.asarray(u, float)
        inside = (u >= 0.0) & (u <= 1.0)
        if self.model == "none":
            raise ValueError("model='none' has a delta-function pdf; use "
                             "survival() or a different model")
        if self.model == "flat" or species == "nu":
            return np.where(inside, 1.0, 0.0)
        return np.where(inside, 3.0 * u**2, 0.0)   # from 3(1-y)^2

    def survival(self, z, species: Species = "nu") -> np.ndarray:
        """P(u > z).

        NB z must NOT be clipped into [0, 1] before the branch: for
        model='none' the whole distribution sits at u = 1, so S(z) is a step
        and clipping z_hi down to 1 would wrongly cancel the band
        S(z_lo) - S(z_hi) in `arriving_muon_spectrum`.
        """
        z = np.asarray(z, float)
        if self.model == "none":
            return np.where(z < 1.0, 1.0, 0.0)
        zc = np.clip(z, 0.0, 1.0)
        s = (1.0 - zc) if (self.model == "flat" or species == "nu") else (1.0 - zc**3)
        return np.where(z > 1.0, 0.0, np.where(z < 0.0, 1.0, s))

    def mean_u(self, species: Species = "nu") -> float:
        if self.model == "none":
            return 1.0
        return 0.5 if (self.model == "flat" or species == "nu") else 0.75

    def as_metadata(self):
        return {"inelasticity_model": self.model}
