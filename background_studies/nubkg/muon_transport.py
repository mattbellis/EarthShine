"""
muon_transport.py -- parameterised muon energy loss in rock.

Continuous-slowing-down approximation with

    -dE/dX = a + b E        [X in g/cm^2, E in GeV]

which integrates analytically to

    E(X)          = (E0 + a/b) exp(-b X) - a/b
    R(E0 -> Eth)  = (1/b) ln[(a + b E0) / (a + b Eth)]        [g/cm^2]

Everything downstream in this package uses only these three closed forms plus
the stopping power itself, which is why no muon-by-muon Monte Carlo is needed.

Parameter values
----------------
The relevant reference for consistency with the signal side is Feng, Smolinsky
& Tanedo, arXiv:1509.07525, Eq. (32), which uses
    a = 2.0e-3 GeV cm^2/g,  b = 4.2e-6 cm^2/g   (water/ice, rho = 1)
Standard rock is slightly harder.  Groom, Mokhov & Striganov (Atomic Data and
Nuclear Data Tables 78, 183 (2001)) tabulate muon losses in standard rock; the
usual working values are a ~ 2.2e-3 and b ~ 4.4e-6 in the 0.1-10 TeV decade.
    https://pdg.lbl.gov/2024/AtomicNuclearProperties/  (muon dE/dx tables)

CAVEAT ON CONSTANT b: b actually rises with energy (roughly 3.5e-6 at 10 GeV to
~6e-6 at 100 TeV in rock) because pair production and photonuclear losses grow.
A single constant b is therefore a compromise.  `EnergyLoss.bracket()` returns
soft/hard variants so this can be carried as an explicit systematic rather than
being buried.  Stochastic fluctuations are also ignored: the CSDA range is a
mild overestimate of the effective range above the critical energy a/b ~ 500
GeV.  Both are the sort of thing the eventual GENIE -> GEANT4 chain fixes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["EnergyLoss", "ROCK_LOSS", "WATER_LOSS"]


@dataclass(frozen=True)
class EnergyLoss:
    """Muon energy-loss parameters.  a in GeV cm^2/g, b in cm^2/g."""

    a: float = 2.2e-3
    b: float = 4.4e-6
    label: str = "standard rock (Groom et al., ~TeV decade)"

    @property
    def critical_energy(self) -> float:
        """a/b in GeV: above this, radiative losses dominate and the range
        grows only logarithmically with energy."""
        return self.a / self.b

    def stopping_power(self, e_mu) -> np.ndarray:
        """-dE/dX in GeV per g/cm^2."""
        return self.a + self.b * np.asarray(e_mu, float)

    def range_gcm2(self, e0, e_thr) -> np.ndarray:
        """Column density (g/cm^2) over which a muon stays above e_thr."""
        e0 = np.asarray(e0, float)
        eth = np.asarray(e_thr, float)
        ratio = (self.a + self.b * e0) / (self.a + self.b * eth)
        return np.where(e0 > eth, np.log(np.maximum(ratio, 1.0)) / self.b, 0.0)

    def range_cm(self, e0, e_thr, density) -> np.ndarray:
        return self.range_gcm2(e0, e_thr) / density

    def energy_after(self, e0, x_gcm2) -> np.ndarray:
        """Energy after traversing x_gcm2, clipped at zero."""
        e0 = np.asarray(e0, float)
        ab = self.a / self.b
        return np.maximum((e0 + ab) * np.exp(-self.b * np.asarray(x_gcm2, float)) - ab, 0.0)

    def energy_before(self, e_final, x_gcm2) -> np.ndarray:
        """Energy needed at production to arrive with e_final after x_gcm2.
        Inverse of `energy_after`; used to cap the production energy when the
        available rock path is limited (down-going muons, shallow detector)."""
        ef = np.asarray(e_final, float)
        ab = self.a / self.b
        return (ef + ab) * np.exp(self.b * np.asarray(x_gcm2, float)) - ab

    def bracket(self):
        """(soft, nominal, hard) variants spanning the plausible constant-b
        range for rock.  Use to quote an energy-loss systematic."""
        return (EnergyLoss(self.a, 0.8 * self.b, self.label + " [soft b]"),
                self,
                EnergyLoss(self.a, 1.3 * self.b, self.label + " [hard b]"))

    def as_metadata(self, prefix="dedx"):
        return {f"{prefix}_a_GeV_cm2_g": self.a, f"{prefix}_b_cm2_g": self.b,
                f"{prefix}_label": self.label,
                f"{prefix}_E_crit_GeV": self.critical_energy}


ROCK_LOSS = EnergyLoss()
WATER_LOSS = EnergyLoss(2.0e-3, 4.2e-6, "water/ice (1509.07525 Eq. 32)")
