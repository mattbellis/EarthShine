"""EarthShine neutrino-induced muon background package."""
from pathlib import Path as _Path

PACKAGE_DIR = _Path(__file__).resolve().parent
DATA_DIR = PACKAGE_DIR / "data"


def data_path(name: str) -> str:
    """Absolute path to a file in nubkg/data/.

    Use this instead of a relative path so notebooks and scripts work from any
    working directory:

        flux = nb.TabulatedFlux.from_csv(nb.data_path("ic59_atmospheric_numu_APPROX.csv"))
    """
    p = DATA_DIR / name
    if not p.exists():
        available = sorted(f.name for f in DATA_DIR.glob("*"))
        raise FileNotFoundError(
            f"{name!r} not found in {DATA_DIR}. Available: {available}")
    return str(p)

from .muon_transport import EnergyLoss, ROCK_LOSS, WATER_LOSS
from .nu_xsec import CrossSection, Inelasticity, default_cc, sigma_uncertainty
from .nu_flux import (ChirkinAtmospheric, TabulatedFlux, AstrophysicalPowerLaw,
                      SumFlux, cos_theta_star, earth_column_depth,
                      earth_transmission)
from .nu_background import (Rock, STANDARD_ROCK, MOLASSE, Cylinder, CMS,
                            muon_production_spectrum, arriving_muon_spectrum,
                            integrated_muon_flux, detector_spectrum,
                            background_muons, SEC_PER_YEAR)
__all__ = [n for n in dir() if not n.startswith("_")]

