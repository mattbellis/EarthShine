"""
nu_background.py -- neutrino-induced muon background at an underground detector.

Everything here is deterministic quadrature.  No muons are thrown.

--------------------------------------------------------------------------
The geometry result that makes this tractable
--------------------------------------------------------------------------
Write the rate as an integral over the rock volume, with each element weighted
by the solid angle the detector subtends from it:

    dN/dt = Int_V dV n_N Int dE_nu sigma phi_nu dOmega_det(r) P_reach(r)

With dV = r^2 dr dOmega and dOmega_det ~ A_proj / r^2, the r^2 cancels
*exactly*.  The shrinking solid angle and the growing shell volume are the same
factor with opposite signs, so there is no volume to choose: the r integral
runs from 0 to the muon range and converges by itself.

    dN/dt = Int dOmega A_proj(Omega) Phi_mu(Omega)

--------------------------------------------------------------------------
The arriving muon spectrum, in closed form
--------------------------------------------------------------------------
A muon produced at column depth X with energy E_prod arrives with
E_mu = energy_after(E_prod, X).  Converting the integral over X into an
integral over E_mu introduces the Jacobian |dE/dX|^-1 = 1/(a + b E_mu), and the
integral over the inelasticity collapses onto the survival function
S(z) = P(E_mu/E_nu > z).  The result is

    dPhi_mu/dE_mu (Omega) = N_A / (a + b E_mu)
          x Int dE_nu phi_nu(E_nu, Omega) sigma_CC(E_nu)
                       [ S(E_mu/E_nu) - S(E_cap/E_nu) ]

where E_cap = energy_before(E_mu, X_max) caps the production energy when the
available rock path X_max is finite (down-going muons under a shallow
overburden); for X_max -> inf the second term vanishes.

Two consequences worth noting:

  * the rock *density cancels* (n_N ~ rho, range ~ 1/rho), so molasse and
    standard rock give the same answer to <1%.  What matters is a and b.
  * integrating the expression above over E_mu from E_thr to infinity
    reproduces the familiar n_N sigma <R_mu> form exactly.  That identity is
    a unit test (`test_spectrum_integral_matches_range_formula`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from .muon_transport import EnergyLoss, ROCK_LOSS
from .nu_flux import ChirkinAtmospheric, SumFlux, earth_transmission
from .nu_xsec import CrossSection, Inelasticity, default_cc, sigma_uncertainty

N_AVOGADRO = 6.02214076e23
SEC_PER_YEAR = 365.25 * 24 * 3600.0
SPECIES = ("nu", "nubar")

__all__ = [
    "Rock", "STANDARD_ROCK", "MOLASSE", "Cylinder", "CMS",
    "muon_production_spectrum", "arriving_muon_spectrum",
    "integrated_muon_flux", "detector_spectrum", "background_muons",
]


# ---------------------------------------------------------------------------
# Materials and geometry
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Rock:
    """Rock target.  `density` only enters conversions to physical distance --
    the rate itself is density-independent (see module docstring)."""

    name: str = "standard rock"
    density: float = 2.65                     # g/cm^3
    loss: EnergyLoss = field(default_factory=lambda: ROCK_LOSS)

    @property
    def n_nucleon_per_cm3(self) -> float:
        return self.density * N_AVOGADRO

    def as_metadata(self):
        md = {"rock_name": self.name, "rock_density_g_cm3": self.density}
        md.update(self.loss.as_metadata())
        return md


STANDARD_ROCK = Rock()
MOLASSE = Rock("LHC molasse", 2.40)


@dataclass(frozen=True)
class Cylinder:
    """Cylindrical detector with a horizontal axis (CMS-like: axis = beam line).

    radius, half_length, depth in metres.  `depth` is the rock overburden,
    which limits the path available to down-going muons.
    """

    radius: float = 7.5
    half_length: float = 10.5
    depth: float = 100.0
    axis_azimuth: float = 0.0
    name: str = "CMS-like"

    @property
    def volume_m3(self) -> float:
        return np.pi * self.radius**2 * 2.0 * self.half_length

    def projected_area_m2(self, cos_zenith, azimuth) -> np.ndarray:
        """Shadow area presented to a track: 2 R L sin(alpha) + pi R^2 |cos(alpha)|,
        alpha being the angle between the track and the cylinder axis."""
        cz = np.asarray(cos_zenith, float)
        sz = np.sqrt(np.maximum(1.0 - cz**2, 0.0))
        ca = sz * np.cos(np.asarray(azimuth, float) - self.axis_azimuth)
        sa = np.sqrt(np.maximum(1.0 - ca**2, 0.0))
        return (2.0 * self.radius * 2.0 * self.half_length * sa
                + np.pi * self.radius**2 * np.abs(ca))

    def mean_projected_area_m2(self, cos_zenith, n_azi: int = 64) -> np.ndarray:
        az = np.linspace(0.0, 2 * np.pi, n_azi, endpoint=False)
        cz = np.atleast_1d(np.asarray(cos_zenith, float))
        return self.projected_area_m2(cz[:, None], az[None, :]).mean(axis=1)

    def max_rock_path_gcm2(self, cos_zenith, rock: Rock) -> np.ndarray:
        """Rock available along a given direction.  Unlimited from below;
        limited by the overburden from above."""
        cz = np.asarray(cos_zenith, float)
        with np.errstate(divide="ignore"):
            return np.where(cz < 0, np.inf,
                            self.depth * 100.0 * rock.density / np.maximum(cz, 1e-12))

    def as_metadata(self):
        return {"detector_name": self.name, "detector_radius_m": self.radius,
                "detector_half_length_m": self.half_length,
                "detector_depth_m": self.depth,
                "detector_volume_m3": self.volume_m3}


CMS = Cylinder()


# ---------------------------------------------------------------------------
# Defaults holder
# ---------------------------------------------------------------------------

def _defaults(flux, xsec, inel, rock):
    flux = flux if flux is not None else ChirkinAtmospheric()
    xsec = xsec if xsec is not None else {s: default_cc(s) for s in SPECIES}
    if isinstance(xsec, CrossSection):
        xsec = {s: xsec for s in SPECIES}
    inel = inel if inel is not None else Inelasticity()
    rock = rock if rock is not None else STANDARD_ROCK
    return flux, xsec, inel, rock


def _energy_grid(lo, hi, n):
    return np.geomspace(lo, hi, n)


# ---------------------------------------------------------------------------
# 1. Muon yield from a cubic metre of rock  (diagnostic)
# ---------------------------------------------------------------------------

def muon_production_spectrum(
    e_mu,
    *,
    flux=None, xsec=None, inelasticity=None, rock=None,
    e_nu_range=(10.0, 1.0e6), n_energy: int = 400,
    hemisphere: Literal["up", "down", "all"] = "up",
    n_cos: int = 40, absorb: bool = True,
) -> dict:
    """Rate of muons *produced* per cubic metre of rock, differential in the
    muon energy at the interaction point.

        d2N / (dV dE_mu dt)   [ m^-3 s^-1 GeV^-1 ]

    integrated over the chosen hemisphere of neutrino directions.  This is the
    "muons coming out of a cubic metre of rock" diagnostic: it says nothing
    about whether they reach the detector, only about production.

    Uses  dsigma/dE_mu = sigma(E_nu) f_u(E_mu/E_nu) / E_nu.
    """
    flux, xsec, inel, rock = _defaults(flux, xsec, inelasticity, rock)
    e_mu = np.atleast_1d(np.asarray(e_mu, float))
    e_nu = _energy_grid(*e_nu_range, n_energy)
    lo, hi = {"up": (-1.0, 0.0), "down": (0.0, 1.0), "all": (-1.0, 1.0)}[hemisphere]
    cz = np.linspace(lo, hi, n_cos)

    total = np.zeros_like(e_mu)
    per_species = {}
    for sp in SPECIES:
        sig = xsec[sp](e_nu)                                     # (nE,)
        phi = flux(e_nu[None, :], cz[:, None], sp)               # (nZ, nE)
        if absorb:
            phi = phi * earth_transmission(e_nu[None, :], cz[:, None], xsec[sp])
        # integrate flux over solid angle: 2 pi Int dcos
        phi_omega = 2.0 * np.pi * np.trapezoid(phi, cz, axis=0)  # (nE,)

        u = e_mu[:, None] / e_nu[None, :]
        fu = inel.pdf_u(u, sp) if inel.model != "none" else None
        if fu is None:
            # delta at u = 1: handle by direct interpolation of phi*sigma
            contrib = np.interp(e_mu, e_nu, phi_omega * sig, left=0.0, right=0.0)
        else:
            integrand = phi_omega[None, :] * sig[None, :] * fu / e_nu[None, :]
            contrib = np.trapezoid(integrand, e_nu, axis=1)
        per_species[sp] = contrib
        total = total + contrib

    # per cm^3 -> per m^3
    n_per_m3 = rock.n_nucleon_per_cm3 * 1.0e6
    md = {"hemisphere": hemisphere, "e_nu_range_GeV": tuple(e_nu_range)}
    md.update(rock.as_metadata()); md.update(inel.as_metadata())
    md.update(flux.as_metadata()); md.update(xsec["nu"].as_metadata())
    return {"e_mu": e_mu,
            "dN_dEmu_per_m3_per_s": total * n_per_m3,
            "per_species": {k: v * n_per_m3 for k, v in per_species.items()},
            "metadata": md}


# ---------------------------------------------------------------------------
# 2. Arriving muon spectrum
# ---------------------------------------------------------------------------

def arriving_muon_spectrum(
    e_mu, cos_zenith,
    *,
    flux=None, xsec=None, inelasticity=None, rock=None,
    max_path_gcm2=np.inf,
    e_nu_range=(10.0, 1.0e6), n_energy: int = 400,
    absorb: bool = True,
    flux_scale=None, xsec_scale=None,
) -> np.ndarray:
    """dPhi_mu/dE_mu at the detector, GeV^-1 cm^-2 s^-1 sr^-1.

    Shape: (len(cos_zenith), len(e_mu)).

    flux_scale / xsec_scale : optional callables of E_nu returning a
    multiplicative factor.  Used by the uncertainty budget to apply correlated
    +-1 sigma shifts without duplicating the integration code.
    """
    flux, xsec, inel, rock = _defaults(flux, xsec, inelasticity, rock)
    loss = rock.loss
    e_mu = np.atleast_1d(np.asarray(e_mu, float))
    cz = np.atleast_1d(np.asarray(cos_zenith, float))
    xmax = np.broadcast_to(np.atleast_1d(np.asarray(max_path_gcm2, float)), cz.shape)
    e_nu = _energy_grid(*e_nu_range, n_energy)

    # E_cap[z, m]: production energy that exactly exhausts the available rock
    with np.errstate(over="ignore"):
        e_cap = loss.energy_before(e_mu[None, :], xmax[:, None])   # (nZ, nMu)
    e_cap = np.where(np.isfinite(e_cap), e_cap, np.inf)

    out = np.zeros((cz.size, e_mu.size))
    for sp in SPECIES:
        sig = xsec[sp](e_nu)
        if xsec_scale is not None:
            sig = sig * xsec_scale(e_nu)
        phi = flux(e_nu[None, :], cz[:, None], sp)                 # (nZ, nE)
        if absorb:
            phi = phi * earth_transmission(e_nu[None, :], cz[:, None], xsec[sp])
        if flux_scale is not None:
            phi = phi * flux_scale(e_nu)[None, :]
        w = phi * sig[None, :]                                     # (nZ, nE)

        z_lo = e_mu[None, :, None] / e_nu[None, None, :]           # (1, nMu, nE)
        z_hi = e_cap[:, :, None] / e_nu[None, None, :]             # (nZ, nMu, nE)
        band = inel.survival(z_lo, sp) - inel.survival(z_hi, sp)
        band = np.clip(band, 0.0, None)

        out += np.trapezoid(band * w[:, None, :], e_nu, axis=2)

    # N_A per gram / (a + b E_mu) in g cm^-2 GeV^-1
    out *= N_AVOGADRO / loss.stopping_power(e_mu)[None, :]
    return out


def integrated_muon_flux(
    e_thr: float, cos_zenith,
    *, flux=None, xsec=None, inelasticity=None, rock=None,
    max_path_gcm2=np.inf, e_nu_range=(10.0, 1.0e6), n_energy: int = 400,
    absorb: bool = True,
) -> np.ndarray:
    """Phi_mu(> e_thr) in cm^-2 s^-1 sr^-1, computed directly from the range
    formula (not by integrating the spectrum).  The two must agree -- that is
    the main internal consistency check of this package.
    """
    flux, xsec, inel, rock = _defaults(flux, xsec, inelasticity, rock)
    loss = rock.loss
    cz = np.atleast_1d(np.asarray(cos_zenith, float))
    xmax = np.broadcast_to(np.atleast_1d(np.asarray(max_path_gcm2, float)), cz.shape)
    e_nu = _energy_grid(*e_nu_range, n_energy)

    # <R> over the inelasticity distribution, by quadrature in u
    u = np.linspace(1e-6, 1.0, 512)
    out = np.zeros(cz.size)
    for sp in SPECIES:
        fu = (np.ones_like(u) if inel.model == "none" else inel.pdf_u(u, sp))
        e_prod = u[:, None] * e_nu[None, :]
        r = loss.range_gcm2(e_prod, e_thr)                          # (nU, nE)
        if inel.model == "none":
            r_mean = loss.range_gcm2(e_nu, e_thr)
        else:
            r_mean = np.trapezoid(fu[:, None] * r, u, axis=0)        # (nE,)
        r_clipped = np.minimum(r_mean[None, :], xmax[:, None])

        sig = xsec[sp](e_nu)
        phi = flux(e_nu[None, :], cz[:, None], sp)
        if absorb:
            phi = phi * earth_transmission(e_nu[None, :], cz[:, None], xsec[sp])
        out += np.trapezoid(phi * sig[None, :] * r_clipped, e_nu, axis=1)

    return out * N_AVOGADRO


# ---------------------------------------------------------------------------
# 3. Detector-level spectrum and counts
# ---------------------------------------------------------------------------

def detector_spectrum(
    e_mu, detector: Cylinder = CMS,
    *, rock=None, hemisphere: Literal["up", "down", "all"] = "up",
    n_cos: int = 40, **kw,
) -> dict:
    """dN/dE_mu in s^-1 GeV^-1 for muons crossing the detector."""
    rock = rock if rock is not None else STANDARD_ROCK
    lo, hi = {"up": (-1.0, 0.0), "down": (0.0, 1.0), "all": (-1.0, 1.0)}[hemisphere]
    cz = np.linspace(lo, hi, n_cos)
    if hemisphere in ("up", "all"):
        cz = np.where(np.abs(cz) < 1e-6, -1e-6, cz)

    xmax = detector.max_rock_path_gcm2(cz, rock)
    dphi = arriving_muon_spectrum(e_mu, cz, rock=rock, max_path_gcm2=xmax, **kw)
    area_cm2 = detector.mean_projected_area_m2(cz) * 1.0e4         # (nZ,)
    integrand = dphi * area_cm2[:, None] * 2.0 * np.pi
    return {"e_mu": np.atleast_1d(np.asarray(e_mu, float)),
            "dN_dEmu_per_s": np.trapezoid(integrand, cz, axis=0),
            "cos_zenith": cz,
            "dPhi_dEmu": dphi,
            "mean_projected_area_m2": area_cm2 / 1e4}


def background_muons(
    detector: Cylinder = CMS,
    *,
    e_mu_min: float = 10.0,
    e_nu_range: tuple = (10.0, 1.0e6),
    livetime_s: float = SEC_PER_YEAR,
    rock=None, flux=None, xsec=None, inelasticity=None,
    hemisphere: Literal["up", "down", "all"] = "up",
    n_energy: int = 400, n_cos: int = 40, n_mu: int = 160,
    absorb: bool = True, uncertainty: bool = True,
) -> dict:
    """TOP-LEVEL API.

    Parameters
    ----------
    detector    : Cylinder (dimensions and overburden)
    e_mu_min    : minimum muon energy *at the detector*, GeV
    e_nu_range  : (E_min, E_max) of the neutrinos considered, GeV
    livetime_s  : exposure in seconds
    hemisphere  : "up" is the EarthShine signal region

    Returns
    -------
    dict with
        e_mu                : muon energy grid at the detector, GeV
        dN_dEmu_per_s       : arriving spectrum, s^-1 GeV^-1
        dN_dEmu_livetime    : same x livetime, GeV^-1
        n_muons             : integrated count above e_mu_min in the livetime
        n_muons_err         : total 1-sigma
        error_budget        : dict of the individual contributions
        rate_per_year       : convenience
        metadata            : flat dict, ready for a parquet footer
    """
    flux, xsec, inel, rock = _defaults(flux, xsec, inelasticity, rock)
    e_hi = min(e_nu_range[1], 1.0e9)
    e_mu = np.geomspace(e_mu_min, e_hi, n_mu)

    common = dict(flux=flux, xsec=xsec, inelasticity=inel,
                  e_nu_range=e_nu_range, n_energy=n_energy, absorb=absorb)

    spec = detector_spectrum(e_mu, detector, rock=rock, hemisphere=hemisphere,
                             n_cos=n_cos, **common)
    dn = spec["dN_dEmu_per_s"]
    n_expected = float(np.trapezoid(dn, e_mu) * livetime_s)

    budget = {}
    if uncertainty:
        # --- flux normalisation (correlated shift, energy dependent)
        if hasattr(flux, "fractional_error"):
            ferr = flux.fractional_error
        else:
            ferr = lambda e: np.full_like(np.asarray(e, float), 0.30)
        up = detector_spectrum(e_mu, detector, rock=rock, hemisphere=hemisphere,
                               n_cos=n_cos, flux_scale=lambda e: 1 + ferr(e),
                               **common)["dN_dEmu_per_s"]
        budget["flux"] = abs(float(np.trapezoid(up, e_mu)) * livetime_s - n_expected)

        # --- cross section
        up = detector_spectrum(e_mu, detector, rock=rock, hemisphere=hemisphere,
                               n_cos=n_cos,
                               xsec_scale=lambda e: 1 + sigma_uncertainty(e),
                               **common)["dN_dEmu_per_s"]
        budget["cross_section"] = abs(float(np.trapezoid(up, e_mu)) * livetime_s
                                      - n_expected)

        # --- muon energy loss (constant-b bracket)
        soft, _, hard = rock.loss.bracket()
        vals = []
        for lo_ in (soft, hard):
            r2 = Rock(rock.name, rock.density, lo_)
            s2 = detector_spectrum(e_mu, detector, rock=r2, hemisphere=hemisphere,
                                   n_cos=n_cos, **common)["dN_dEmu_per_s"]
            vals.append(float(np.trapezoid(s2, e_mu)) * livetime_s)
        budget["energy_loss"] = 0.5 * abs(vals[0] - vals[1])

        # --- inelasticity model (qpm vs flat)
        alt = Inelasticity("flat" if inel.model != "flat" else "qpm")
        c2 = dict(common); c2["inelasticity"] = alt
        s2 = detector_spectrum(e_mu, detector, rock=rock, hemisphere=hemisphere,
                               n_cos=n_cos, **c2)["dN_dEmu_per_s"]
        budget["inelasticity"] = abs(float(np.trapezoid(s2, e_mu)) * livetime_s
                                     - n_expected)

        budget["poisson"] = np.sqrt(max(n_expected, 0.0))

    total_err = float(np.sqrt(sum(v**2 for v in budget.values()))) if budget else 0.0

    md = {"stage": "nu_background", "e_mu_min_GeV": e_mu_min,
          "e_nu_min_GeV": e_nu_range[0], "e_nu_max_GeV": e_nu_range[1],
          "livetime_s": livetime_s, "hemisphere": hemisphere,
          "earth_absorption": absorb, "n_muons": n_expected,
          "n_muons_err": total_err}
    md.update(rock.as_metadata()); md.update(detector.as_metadata())
    md.update(inel.as_metadata()); md.update(flux.as_metadata())
    md.update(xsec["nu"].as_metadata())

    return {"e_mu": e_mu,
            "dN_dEmu_per_s": dn,
            "dN_dEmu_livetime": dn * livetime_s,
            "n_muons": n_expected,
            "n_muons_err": total_err,
            "error_budget": budget,
            "rate_per_year": float(np.trapezoid(dn, e_mu) * SEC_PER_YEAR),
            "cos_zenith": spec["cos_zenith"],
            "dPhi_dEmu": spec["dPhi_dEmu"],
            "metadata": md}
