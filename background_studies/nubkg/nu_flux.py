"""
nu_flux.py -- atmospheric and astrophysical neutrino fluxes.

All fluxes return dPhi/dE in GeV^-1 cm^-2 s^-1 sr^-1 for a single species
("nu" or "nubar" muon neutrinos), given (E_nu / GeV, cos_zenith).

Two independent routes are provided, and the point of having both is that they
are genuinely independent:

TabulatedFlux  -- digitised from a *measurement*.  The reference case is
    IceCube's unfolded atmospheric nu_mu spectrum (IC-59),
        https://icecube.wisc.edu/news/research/2014/09/an-improved-measurement-of-atmospheric-neutrino-flux-in-icecube/
    published as Aartsen et al., Eur. Phys. J. C 75, 116 (2015),
        https://arxiv.org/abs/1409.4535
    The unfolded points are nu_mu + nubar_mu, averaged over the northern
    (up-going) sky, and span roughly 100 GeV - 1 PeV.  Below that, the same
    figure carries the Frejus measurement (~0.3-50 GeV),
        https://arxiv.org/abs/1409.4535  (Fig. 8 / Frejus points)

ChirkinAtmospheric -- an analytic parameterisation of the conventional pi/K
    flux, fitted by Chirkin to CORSIKA air-shower simulations,
        https://arxiv.org/abs/hep-ph/0407078

IS THE ANALYTIC MODEL CIRCULAR WITH THE ICECUBE DATA?  No.  CORSIKA is an
air-shower Monte Carlo; its inputs are the primary cosmic-ray spectrum and
composition (Hoerandel's poly-gonato model in Chirkin's fit) plus a hadronic
interaction model (QGSJET).  It takes no neutrino-telescope data as input.
The IceCube unfolding is a measurement of the same quantity.  Comparing the
two is therefore a real cross-check, not a tautology -- and where they
disagree, the data should win on normalisation.

The hybrid case is usually the best of both: take the *normalisation* from the
measured table and the *zenith shape* from the analytic model, via
    TabulatedFlux(..., zenith_shape_from=ChirkinAtmospheric())
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Sequence

import numpy as np

Species = Literal["nu", "nubar"]

__all__ = [
    "cos_theta_star", "ChirkinAtmospheric", "TabulatedFlux",
    "AstrophysicalPowerLaw", "SumFlux", "earth_column_depth",
    "earth_transmission",
]

N_AVOGADRO = 6.02214076e23


# ---------------------------------------------------------------------------
# Zenith geometry
# ---------------------------------------------------------------------------

# Chirkin hep-ph/0407078 Table 1: the cos(theta) -> cos(theta*) correction that
# extends the flat-atmosphere formula to the horizon.
_COSTHETA_STAR_P = (0.102573, -0.068287, 0.958633, 0.0407253, 0.817285)


def cos_theta_star(cos_zenith) -> np.ndarray:
    """Effective zenith cosine including atmospheric curvature.

    Uses |cos(zenith)|: an up-going neutrino was produced in the atmosphere on
    the far side of the Earth at the mirrored zenith angle.
    """
    x = np.abs(np.asarray(cos_zenith, float))
    p1, p2, p3, p4, p5 = _COSTHETA_STAR_P
    num = x**2 + p1**2 + p2 * x**p3 + p4 * x**p5
    den = 1.0 + p1**2 + p2 + p4
    return np.sqrt(np.maximum(num / den, 0.0))


# ---------------------------------------------------------------------------
# Analytic conventional atmospheric flux
# ---------------------------------------------------------------------------

@dataclass
class ChirkinAtmospheric:
    """Conventional (pi/K) atmospheric nu_mu + nubar_mu flux.

        dN/dE = 2.85e-2 A E^-gamma [ 1/(1 + 6 E c*/121)
                                     + 0.213/(1 + 1.44 E c*/897) ]

    A = 0.646, gamma = 2.684 (Chirkin, hep-ph/0407078, Table 2).

    Validity: fitted over 600 GeV - 60 TeV.  Below ~100 GeV it omits nu_mu from
    muon decay and the primary-spectrum fit is extrapolating; treat as +-30%
    there.  No oscillations (a <10% effect above 30 GeV over an Earth
    diameter).  No prompt/charm component (negligible below ~50 TeV).
    """

    norm: float = 2.85e-2
    A: float = 0.646
    gamma: float = 2.684
    nu_fraction: float = 0.635
    label: str = "Chirkin/Volkova conventional atmospheric (hep-ph/0407078)"

    def __call__(self, e_nu, cos_zenith, species: Species = "nu") -> np.ndarray:
        e = np.asarray(e_nu, float)
        cs = cos_theta_star(cos_zenith)
        pion = 1.0 / (1.0 + 6.0 * e * cs / 121.0)
        kaon = 0.213 / (1.0 + 1.44 * e * cs / 897.0)
        total = self.norm * self.A * e ** (-self.gamma) * (pion + kaon)
        frac = self.nu_fraction if species == "nu" else 1.0 - self.nu_fraction
        return frac * total

    def zenith_averaged(self, e_nu, species: Species = "nu",
                        hemisphere: Literal["up", "down", "all"] = "up",
                        n_cos: int = 200) -> np.ndarray:
        lo, hi = {"up": (-1.0, 0.0), "down": (0.0, 1.0), "all": (-1.0, 1.0)}[hemisphere]
        cz = np.linspace(lo, hi, n_cos)
        vals = self(np.atleast_1d(e_nu)[:, None], cz[None, :], species)
        return np.trapezoid(vals, cz, axis=1) / (hi - lo)

    def as_metadata(self):
        return {"flux_model": self.label, "flux_A": self.A,
                "flux_gamma": self.gamma, "flux_nu_fraction": self.nu_fraction}


# ---------------------------------------------------------------------------
# Tabulated / digitised flux
# ---------------------------------------------------------------------------

@dataclass
class TabulatedFlux:
    """Flux from a digitised curve, log-log interpolated.

    Parameters
    ----------
    energy : GeV
    flux   : dPhi/dE for nu_mu + nubar_mu combined, GeV^-1 cm^-2 s^-1 sr^-1
    flux_err : optional 1-sigma, same units.  Used by the uncertainty budget.
    nu_fraction : split of the tabulated total between nu and nubar.
    zenith_shape_from : optional model with a `zenith_averaged` method.  If
        given, the tabulated value is treated as the average over
        `shape_hemisphere` and the zenith dependence is taken from the model:

            Phi(E, cz) = Phi_table(E) * model(E, cz) / <model(E)>

        This preserves the measured normalisation while restoring the
        (substantial, factor ~3 horizon-to-vertical) zenith dependence.
    """

    energy: np.ndarray
    flux: np.ndarray
    flux_err: np.ndarray | None = None
    nu_fraction: float = 0.635
    label: str = "unlabelled tabulated flux"
    zenith_shape_from: object | None = None
    shape_hemisphere: Literal["up", "down", "all"] = "up"
    extrapolate: bool = True
    outside: Literal["model", "extrapolate", "zero", "error"] = "model"

    def __post_init__(self):
        self.energy = np.asarray(self.energy, float)
        self.flux = np.asarray(self.flux, float)
        order = np.argsort(self.energy)
        self.energy, self.flux = self.energy[order], self.flux[order]
        if self.flux_err is not None:
            self.flux_err = np.asarray(self.flux_err, float)[order]
        self._lx = np.log(self.energy)
        self._ly = np.log(self.flux)

    @classmethod
    def from_csv(cls, path, label=None, energy_unit="GeV",
                 flux_unit: Literal["dPhi_dE", "E2Phi"] = "E2Phi", **kw):
        """Load a digitised curve.

        Columns: E, flux[, flux_err].  `energy_unit='log10GeV'` and
        `flux_unit='E2Phi'` are the natural choices for points read off the
        standard E^2 Phi vs log10(E/GeV) atmospheric-flux figure.
        """
        arr = np.genfromtxt(path, delimiter=",", comments="#")
        if arr.ndim == 1:
            arr = arr[None, :]
        e = arr[:, 0]
        if energy_unit == "log10GeV":
            e = 10.0**e
        f = arr[:, 1]
        err = arr[:, 2] if arr.shape[1] > 2 else None
        if flux_unit == "E2Phi":
            f = f / e**2
            if err is not None:
                err = err / e**2
        return cls(e, f, flux_err=err, label=label or str(Path(path).name), **kw)

    def _fallback_model(self):
        """Model used outside the tabulated range under outside='model'."""
        return self.zenith_shape_from

    def _model_total(self, e):
        """Zenith-averaged nu+nubar flux of the fallback model."""
        m = self._fallback_model()
        e1 = np.atleast_1d(np.ravel(e))
        tot = sum(m.zenith_averaged(e1, sp, self.shape_hemisphere)
                  for sp in ("nu", "nubar"))
        return tot.reshape(np.shape(e)) if np.shape(e) else tot[0]

    def _table(self, e_nu):
        """Tabulated total flux, with controlled behaviour outside the range.

        WHY THIS MATTERS.  A curve digitised off the IceCube unfolding figure
        typically starts around a few hundred GeV, where the spectrum has
        already steepened past the pion critical energy to roughly E^-3.7.
        Blindly log-log extrapolating that slope down to 10 GeV overshoots the
        true ~E^-2.7 spectrum by more than an order of magnitude, and since the
        rate integrand is dominated by the low-energy end, the answer is then
        wrong by nearly as much.  The default `outside='model'` instead hands
        off to the analytic model, rescaled to join continuously at the table
        edge, so the measured normalisation is preserved where it exists and
        the model shape takes over where it does not.
        """
        x = np.log(np.asarray(e_nu, float))
        y = np.interp(x, self._lx, self._ly)
        lo, hi = x < self._lx[0], x > self._lx[-1]
        if not np.any(lo | hi):
            return np.exp(y)

        if self.outside == "error":
            raise ValueError(
                f"energies outside the tabulated range "
                f"[{self.energy[0]:.3g}, {self.energy[-1]:.3g}] GeV; "
                "set outside='model' or extend the table")
        if self.outside == "zero":
            return np.where(lo | hi, 0.0, np.exp(y))
        if self.outside == "extrapolate" or self._fallback_model() is None:
            if np.any(lo):
                m = (self._ly[1] - self._ly[0]) / (self._lx[1] - self._lx[0])
                y = np.where(lo, self._ly[0] + m * (x - self._lx[0]), y)
            if np.any(hi):
                m = (self._ly[-1] - self._ly[-2]) / (self._lx[-1] - self._lx[-2])
                y = np.where(hi, self._ly[-1] + m * (x - self._lx[-1]), y)
            return np.exp(y)

        out = np.exp(y)
        e = np.exp(x)
        for mask, edge in ((lo, self.energy[0]), (hi, self.energy[-1])):
            if not np.any(mask):
                continue
            scale = np.interp(np.log(edge), self._lx, self._ly)
            scale = np.exp(scale) / self._model_total(edge)
            out = np.where(mask, self._model_total(e) * scale, out)
        return out

    def __call__(self, e_nu, cos_zenith, species: Species = "nu") -> np.ndarray:
        e = np.asarray(e_nu, float)
        cz = np.asarray(cos_zenith, float)
        total = self._table(e) * np.ones_like(e * cz * 1.0)
        if self.zenith_shape_from is not None:
            m = self.zenith_shape_from
            # Normalise by the SPECIES-SUMMED hemisphere average, not the
            # per-species one: `total` is already the nu+nubar table value, so
            # dividing each species by its own average would make the two
            # species each reproduce the full table and double the flux.
            # The model's nu/nubar split is inherited through the numerator.
            denom = self._model_total(e)
            return total * m(e, cz, species) / denom
        frac = self.nu_fraction if species == "nu" else 1.0 - self.nu_fraction
        return frac * total

    def fractional_error(self, e_nu) -> np.ndarray:
        """Interpolated fractional 1-sigma from the tabulated error bars.
        Falls back to 0.30 (a typical unfolding uncertainty) if none given."""
        e = np.asarray(e_nu, float)
        if self.flux_err is None:
            return np.full_like(e, 0.30)
        rel = self.flux_err / self.flux
        return np.interp(np.log(e), self._lx, rel)

    def as_metadata(self):
        return {"flux_model": self.label, "flux_n_points": len(self.energy),
                "flux_e_min_GeV": float(self.energy[0]),
                "flux_e_max_GeV": float(self.energy[-1]),
                "flux_zenith_shape": (None if self.zenith_shape_from is None
                                      else getattr(self.zenith_shape_from,
                                                   "label", "model"))}


# ---------------------------------------------------------------------------
# Astrophysical
# ---------------------------------------------------------------------------

@dataclass
class AstrophysicalPowerLaw:
    """Isotropic diffuse astrophysical nu_mu + nubar_mu, single power law.

    Default: IceCube 9.5-year northern tracks (arXiv:2111.10299),
        Phi = 1.44e-18 (E / 100 TeV)^-2.37   GeV^-1 cm^-2 s^-1 sr^-1
    https://arxiv.org/abs/2111.10299
    """

    norm: float = 1.44e-18
    e0: float = 1.0e5
    gamma: float = 2.37
    label: str = "IceCube 9.5y northern tracks SPL (2111.10299)"

    def __call__(self, e_nu, cos_zenith, species: Species = "nu") -> np.ndarray:
        e = np.asarray(e_nu, float)
        cz = np.asarray(cos_zenith, float)
        return 0.5 * self.norm * (e / self.e0) ** (-self.gamma) * np.ones_like(e * cz * 1.0)

    def as_metadata(self):
        return {"flux_astro": self.label, "flux_astro_norm": self.norm,
                "flux_astro_gamma": self.gamma}


@dataclass
class SumFlux:
    components: Sequence = field(default_factory=list)
    label: str = "sum"

    def __call__(self, e_nu, cos_zenith, species: Species = "nu"):
        out = 0.0
        for c in self.components:
            out = out + c(e_nu, cos_zenith, species)
        return out

    def fractional_error(self, e_nu):
        for c in self.components:
            if hasattr(c, "fractional_error"):
                return c.fractional_error(e_nu)
        return np.full_like(np.asarray(e_nu, float), 0.30)

    def as_metadata(self):
        md = {}
        for i, c in enumerate(self.components):
            md.update({f"c{i}_{k}": v for k, v in c.as_metadata().items()})
        return md


# ---------------------------------------------------------------------------
# Earth absorption
# ---------------------------------------------------------------------------

# Layered approximation to PREM (Dziewonski & Anderson 1981): (r_outer/km, rho)
_EARTH_SHELLS = ((1221.5, 12.95), (3480.0, 11.05), (5701.0, 4.90),
                 (6346.6, 3.55), (6371.0, 2.60))
_COLUMN_CACHE = None


def earth_column_depth(cos_zenith, n_steps: int = 2000) -> np.ndarray:
    """Column density (g/cm^2) traversed by a neutrino reaching a detector near
    the surface.  Zero for down-going."""
    cz = np.atleast_1d(np.asarray(cos_zenith, float))
    out = np.zeros_like(cz)
    up = cz < 0
    if np.any(up):
        r = _EARTH_SHELLS[-1][0]
        chord = 2.0 * r * np.abs(cz[up])
        s = np.linspace(0.0, 1.0, n_steps)[None, :]
        ell = chord[:, None] * (s - 0.5)
        b = r * np.sqrt(np.maximum(1.0 - cz[up][:, None] ** 2, 0.0))
        radius = np.sqrt(b**2 + ell**2)
        rho = np.zeros_like(radius)
        prev = 0.0
        for r_out, d in _EARTH_SHELLS:
            rho = np.where((radius > prev) & (radius <= r_out), d, rho)
            prev = r_out
        out[up] = np.trapezoid(rho, dx=1.0 / (n_steps - 1), axis=1) * chord * 1e5
    shape = np.shape(cos_zenith)
    return out.reshape(shape) if shape else out[0]


def _column_cached(cos_zenith):
    global _COLUMN_CACHE
    if _COLUMN_CACHE is None:
        g = np.linspace(-1.0, 0.0, 2001)
        _COLUMN_CACHE = (g, earth_column_depth(g))
    gx, gy = _COLUMN_CACHE
    cz = np.asarray(cos_zenith, float)
    return np.where(cz < 0, np.interp(np.clip(cz, -1.0, 0.0), gx, gy), 0.0)


def earth_transmission(e_nu, cos_zenith, cross_section, nc_ratio: float = 0.42):
    """exp(-X / lambda_int), lambda_int = m_N / sigma_total.

    Absorption only -- NC regeneration (which degrades rather than removes the
    neutrino) is neglected, so this slightly over-absorbs above ~100 TeV.
    Irrelevant below ~10 TeV.
    """
    x = _column_cached(cos_zenith)
    sigma_tot = cross_section(e_nu) * (1.0 + nc_ratio)
    return np.exp(-x * sigma_tot * N_AVOGADRO)
