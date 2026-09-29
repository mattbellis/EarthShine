"""
cosmic_muons.py -- cosmic-ray muon flux at depth.

This is NOT part of the neutrino background calculation.  It exists so that an
all-sky background estimate can show the thing that actually dominates the
down-going hemisphere: ordinary cosmic-ray muons punching through the
overburden, which outnumber neutrino-induced muons by many orders of magnitude
at shallow depth.

Method
------
Take the Gaisser parameterisation of the muon flux at the surface,

    dN/dE dOmega = 0.14 E^-2.7 [ 1/(1 + 1.1 E cos(th*)/115 GeV)
                               + 0.054/(1 + 1.1 E cos(th*)/850 GeV) ]

in cm^-2 s^-1 sr^-1 GeV^-1 (Gaisser, *Cosmic Rays and Particle Physics*, 1990),
and require that a muon arriving with E_thr had at least

    E_surface = energy_before(E_thr, X(theta))

at the surface, where X(theta) = rho * depth / cos(theta) is the slant column
density.  The same `EnergyLoss` object used for the neutrino-induced muons is
used here, so the two estimates are on a common energy-loss footing -- which is
the point: any comparison between them is then apples to apples.

Accuracy
--------
This is an ORDER-OF-MAGNITUDE tool, good to maybe a factor of 2, and that is
all it is used for.  Known limitations:

  * The Gaisser formula is stated to be valid for E_mu > 100/cos(theta) GeV and
    for theta < 70 deg.  At a 100 m overburden the surviving muons have surface
    energies of only ~60-100 GeV, right at the edge of validity.
  * Muon decay is neglected (fine above ~10 GeV).
  * No stochastic energy-loss fluctuations, so the survival threshold is sharp
    when it should be smeared.  This under-counts the high-energy tail.
  * A flat, uniform overburden is assumed.  Real sites have topography, and
    near the horizon sec(theta) diverges while the true rock path does not.
    `max_sec_theta` caps this.
  * For depths of 1-10 km.w.e. the empirical Mei & Hime relation
    (arXiv:astro-ph/0512125) is more accurate; it is NOT used here because CMS
    at ~240 m.w.e. lies below its stated validity range.

For a real number at a real site, use MUSIC/MUSUN or a full CORSIKA+propagation
chain.  This module is for establishing which background dominates where.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .muon_transport import EnergyLoss, ROCK_LOSS
from .nu_flux import cos_theta_star

__all__ = ["GaisserSurfaceMuons", "muon_flux_at_depth", "cosmic_muon_rate"]


@dataclass
class GaisserSurfaceMuons:
    """Differential muon flux at the surface, cm^-2 s^-1 sr^-1 GeV^-1."""

    norm: float = 0.14
    gamma: float = 2.7
    e_pi: float = 115.0
    e_k: float = 850.0
    k_frac: float = 0.054
    label: str = "Gaisser 1990 surface muon flux"

    def __call__(self, e_mu, cos_zenith) -> np.ndarray:
        e = np.asarray(e_mu, float)
        cs = cos_theta_star(cos_zenith)
        pion = 1.0 / (1.0 + 1.1 * e * cs / self.e_pi)
        kaon = self.k_frac / (1.0 + 1.1 * e * cs / self.e_k)
        return self.norm * e ** (-self.gamma) * (pion + kaon)

    def integral_above(self, e_min, cos_zenith, n: int = 400) -> np.ndarray:
        """Flux above e_min, cm^-2 s^-1 sr^-1."""
        e_min = np.atleast_1d(np.asarray(e_min, float))
        cz = np.broadcast_to(np.atleast_1d(np.asarray(cos_zenith, float)), e_min.shape)
        out = np.zeros_like(e_min)
        for i, (em, c) in enumerate(zip(e_min.ravel(), cz.ravel())):
            # non-finite or absurd e_min means no muon survives that much rock
            if not np.isfinite(em) or em <= 0 or em > 1e14:
                continue
            e = np.geomspace(em, min(max(em * 1e5, 1e7), 1e19), n)
            out.ravel()[i] = np.trapezoid(self(e, c), e)
        return out

    def as_metadata(self):
        return {"cr_surface_model": self.label}


def muon_flux_at_depth(
    e_thr: float, cos_zenith, depth_m: float, density: float = 2.65,
    *, loss: EnergyLoss = ROCK_LOSS, surface=None, max_sec_theta: float = 10.0,
) -> np.ndarray:
    """Cosmic-ray muon flux above `e_thr` at depth, cm^-2 s^-1 sr^-1.

    Zero for up-going directions (cos_zenith < 0): cosmic-ray muons cannot come
    through the Earth.  That is the whole basis of the up-going signal region.
    """
    surface = surface or GaisserSurfaceMuons()
    cz = np.atleast_1d(np.asarray(cos_zenith, float))
    out = np.zeros_like(cz)
    down = cz > 0
    if not np.any(down):
        return out.reshape(np.shape(cos_zenith)) if np.shape(cos_zenith) else out[0]

    sec = np.minimum(1.0 / np.maximum(cz[down], 1e-12), max_sec_theta)
    x = depth_m * 100.0 * density * sec                 # g/cm^2 of slant rock
    e_surface_min = loss.energy_before(e_thr, x)
    out[down] = surface.integral_above(e_surface_min, cz[down])
    return out.reshape(np.shape(cos_zenith)) if np.shape(cos_zenith) else out[0]


def cosmic_muon_rate(detector, e_thr: float, *, rock=None, n_cos: int = 200,
                     loss: EnergyLoss = None, surface=None) -> dict:
    """Rate of cosmic-ray muons through the detector, s^-1.

    Returns the same shape of dict as the neutrino calculation so the two can
    be plotted against each other directly.
    """
    from .nu_background import STANDARD_ROCK, SEC_PER_YEAR
    rock = rock if rock is not None else STANDARD_ROCK
    loss = loss or rock.loss
    cz = np.linspace(1e-3, 1.0, n_cos)
    phi = muon_flux_at_depth(e_thr, cz, detector.depth, rock.density,
                             loss=loss, surface=surface)
    area_cm2 = detector.mean_projected_area_m2(cz) * 1.0e4
    d_rate = phi * area_cm2 * 2.0 * np.pi
    rate = float(np.trapezoid(d_rate, cz))
    return {"cos_zenith": cz, "flux": phi, "dRate_dcos": d_rate,
            "rate_per_s": rate, "rate_per_year": rate * SEC_PER_YEAR,
            "rate_per_day": rate * 86400.0}
