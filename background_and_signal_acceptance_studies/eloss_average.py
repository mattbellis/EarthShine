"""
Muon Energy Loss Through Rock

Calculates the energy/momentum of muons after traversing rock using the
continuous energy loss approximation (average energy loss).

Physics:
--------
The average energy loss for high-energy muons is parameterized as:
    -dE/dx = a(E) + b(E) * E

where:
    a(E) = ionization energy loss (Bethe-Bloch, slowly varying with E)
    b(E) = fractional energy loss from radiative processes:
           - Bremsstrahlung
           - Pair production (e+e-)
           - Photonuclear interactions

For "standard rock" (Z=11, A=22, rho=2.65 g/cm³), using PDG values:
    a ≈ 2.0 MeV/(g/cm²)  [minimum ionizing]
    b ≈ 4.0 × 10⁻⁶ (g/cm²)⁻¹

The critical energy E_c = a/b ≈ 500 GeV marks where radiative losses
equal ionization losses.

References:
-----------
- PDG Review: "Passage of Particles Through Matter"
- Groom, Mokhov, Striganov, "Muon Stopping Power and Range Tables"
  Atomic Data and Nuclear Data Tables 78 (2001) 183-356
"""

import numpy as np
from typing import Tuple, Optional


# Physical constants
MUON_MASS_GEV = 0.1056583755  # Muon mass in GeV/c²

# Standard rock parameters (Z=11, A=22)
ROCK_DENSITY = 2.65  # g/cm³ (standard rock)
ROCK_Z = 11
ROCK_A = 22

# Energy loss parameters for standard rock (PDG values)
# These are averages; in reality they vary slowly with energy
A_IONIZATION = 2.0e-3  # GeV/(g/cm²) - ionization loss
B_RADIATIVE = 4.0e-6   # (g/cm²)⁻¹ - radiative loss coefficient

# Critical energy where ionization = radiative losses
E_CRITICAL = A_IONIZATION / B_RADIATIVE  # ~500 GeV


def momentum_to_energy(p):
    """
    Convert momentum to total energy for a muon.
    
    Parameters
    ----------
    p : float or array-like
        Momentum in GeV/c
        
    Returns
    -------
    float or ndarray
        Total energy in GeV
    """
    p = np.asarray(p)
    return np.sqrt(p**2 + MUON_MASS_GEV**2)


def energy_to_momentum(E):
    """
    Convert total energy to momentum for a muon.
    
    Parameters
    ----------
    E : float or array-like
        Total energy in GeV
        
    Returns
    -------
    float or ndarray
        Momentum in GeV/c (returns 0 if E < muon mass)
    """
    E = np.asarray(E)
    p_squared = np.maximum(E**2 - MUON_MASS_GEV**2, 0.0)
    return np.sqrt(p_squared)


def muon_energy_loss(E0, depth_gcm2, 
                     a: float = A_IONIZATION, 
                     b: float = B_RADIATIVE):
    """
    Calculate muon energy after traversing a given depth in rock.
    
    Uses the analytical solution to -dE/dx = a + b*E:
        E(x) = (E0 + a/b) * exp(-b*x) - a/b
    
    Parameters
    ----------
    E0 : float or array-like
        Initial muon energy in GeV
    depth_gcm2 : float or array-like
        Depth traversed in g/cm²
    a : float, optional
        Ionization loss parameter in GeV/(g/cm²)
    b : float, optional
        Radiative loss parameter in (g/cm²)⁻¹
        
    Returns
    -------
    float or ndarray
        Final energy in GeV (minimum is muon rest mass)
    """
    E0 = np.asarray(E0)
    depth_gcm2 = np.asarray(depth_gcm2)
    
    # Analytical solution
    E_final = (E0 + a/b) * np.exp(-b * depth_gcm2) - a/b
    
    # Muon cannot have less than rest mass energy
    return np.maximum(E_final, MUON_MASS_GEV)


def meters_to_gcm2(distance_m: float, density: float = ROCK_DENSITY) -> float:
    """
    Convert distance in meters to depth in g/cm².
    
    Parameters
    ----------
    distance_m : float
        Distance in meters
    density : float, optional
        Material density in g/cm³
        
    Returns
    -------
    float
        Depth in g/cm²
    """
    # 1 m = 100 cm, so depth = distance_m * 100 cm * density g/cm³
    return distance_m * 100.0 * density


def propagate_muon(p_initial, distance_m,
                   density: float = ROCK_DENSITY,
                   a: float = A_IONIZATION,
                   b: float = B_RADIATIVE):
    """
    Calculate muon momentum after traversing rock.
    
    This is the main function for typical use.
    
    Parameters
    ----------
    p_initial : float or array-like
        Initial muon momentum in GeV/c
    distance_m : float or array-like
        Distance through rock in meters
    density : float, optional
        Rock density in g/cm³ (default: 2.65, standard rock)
    a : float, optional
        Ionization loss parameter in GeV/(g/cm²)
    b : float, optional
        Radiative loss parameter in (g/cm²)⁻¹
        
    Returns
    -------
    float or ndarray
        Final muon momentum in GeV/c (0 if muon stopped)
    """

    # Convert to energy
    E_initial = momentum_to_energy(p_initial)
    
    # Convert distance to g/cm²
    depth = meters_to_gcm2(distance_m, density)
    
    # Calculate energy loss
    E_final = muon_energy_loss(E_initial, depth, a, b)
    
    # Convert back to momentum
    return energy_to_momentum(E_final)

##########################################################################
# Derived from CMSSW code
# 
# https://github.com/cms-sw/cmssw/blob/master/GeneratorInterface/CosmicMuonGenerator/src/SingleParticleEvent.cc#L232C1-L244C2
#
# Look at
#
# SingleParticleEvent
#
##########################################################################
def propagate_muon_CMSSW(p_initial, distance_m,
                   density: float = ROCK_DENSITY,
                   ):

    #double L10E = log10(E);
    # parameters for standard rock (PDG 2004, page 230)
    #double A = (1.91514 + 0.254957 * L10E) / 1000.;                              // a [GeV g^-1 cm^2]
    #double B = (0.379763 + 1.69516 * L10E - 0.175026 * L10E * L10E) / 1000000.;  // b [g^-1 cm^2]
    #double EPS = A / B;                                                          // epsilon [GeV]
    #E = (E + EPS) * exp(-B * waterEquivalents) - EPS;

    E = p_initial
    #
    # Water equivalent = distance (cm) * density (gm/c^3)
    #
    # distance comes in as meters so we need to convert to cm
    #
    waterEquivalents = (distance_m * 100) * ROCK_DENSITY

    L10E = np.log10(E)
    A = (1.91514 + 0.254957 * L10E) / 1000.;                              # a [GeV g^-1 cm^2]
    B = (0.379763 + 1.69516 * L10E - 0.175026 * L10E * L10E) / 1000000.;  # b [g^-1 cm^2]
    EPS = A / B;                                                          # epsilon [GeV]
    E = (E + EPS) * np.exp(-B * waterEquivalents) - EPS;

    return E

##########################################################################

def muon_range(p_initial, density: float = ROCK_DENSITY,
               a: float = A_IONIZATION, b: float = B_RADIATIVE):
    """
    Calculate the range (stopping distance) of a muon in rock.
    
    The range is found by solving E(x) = m_mu:
        x = (1/b) * ln[(E0 + a/b) / (m_mu + a/b)]
    
    Parameters
    ----------
    p_initial : float or array-like
        Initial muon momentum in GeV/c
    density : float, optional
        Rock density in g/cm³
    a : float, optional
        Ionization loss parameter
    b : float, optional
        Radiative loss parameter
        
    Returns
    -------
    float or ndarray
        Range in meters
    """
    E0 = momentum_to_energy(p_initial)
    E_ratio = a / b
    
    # Depth in g/cm² where muon stops
    range_gcm2 = (1.0 / b) * np.log((E0 + E_ratio) / (MUON_MASS_GEV + E_ratio))
    
    # Convert to meters
    return range_gcm2 / (100.0 * density)


def print_energy_loss_table(p_initial: float, distances: list,
                            density: float = ROCK_DENSITY) -> None:
    """
    Print a table of muon momentum vs distance.
    
    Parameters
    ----------
    p_initial : float
        Initial momentum in GeV/c
    distances : list
        List of distances in meters
    density : float, optional
        Rock density in g/cm³
    """
    E_initial = momentum_to_energy(p_initial)
    max_range = muon_range(p_initial, density)
    
    print(f"\nMuon Energy Loss in Standard Rock (ρ = {density:.2f} g/cm³)")
    print(f"Initial momentum: {p_initial:.2f} GeV/c")
    print(f"Initial energy:   {E_initial:.2f} GeV")
    print(f"Critical energy:  {E_CRITICAL:.0f} GeV")
    print(f"Maximum range:    {max_range:.1f} m")
    print("-" * 55)
    print(f"{'Distance (m)':>12} {'Depth (g/cm²)':>14} {'p (GeV/c)':>12} {'E (GeV)':>12}")
    print("-" * 55)
    
    for d in distances:
        p_final = propagate_muon(p_initial, d, density)
        E_final = momentum_to_energy(p_final)
        depth = meters_to_gcm2(d, density)
        
        if p_final > 0:
            print(f"{d:12.1f} {depth:14.0f} {p_final:12.2f} {E_final:12.2f}")
        else:
            print(f"{d:12.1f} {depth:14.0f} {'STOPPED':>12} {'STOPPED':>12}")


