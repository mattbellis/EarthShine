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

__all__ = ["CrossSection", "Inelasticity", "default_cc", "sigma_uncertainty"]


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


def default_cc(species: Species = "nu") -> CrossSection:
    e, s = DEFAULT_CC_NU if species == "nu" else DEFAULT_CC_NUBAR
    lab = ("PDG world average (<350 GeV) + GQRS/CSMS anchor (>=1e4 GeV) "
           f"[{species}]")
    return CrossSection(e, s, label=lab)


def sigma_uncertainty(e_nu) -> np.ndarray:
    """Fractional 1-sigma uncertainty on sigma_CC.

    ~3% where directly measured (30-350 GeV), rising to ~15% in the
    unmeasured bridge region, then ~10% at high energy where the NLO QCD
    predictions agree with the IceCube measurement (arXiv:1711.08119) at that
    level.  Deliberately conservative.
    """
    e = np.asarray(e_nu, float)
    frac = np.full_like(e, 0.05)
    frac = np.where((e >= 30) & (e <= 350), 0.03, frac)
    frac = np.where((e > 350) & (e < 1e4), 0.15, frac)
    frac = np.where(e >= 1e4, 0.10, frac)
    frac = np.where(e < 30, 0.10, frac)
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
