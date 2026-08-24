"""
Unit tests for the neutrino background package.   Run with:  pytest -v

Three kinds of test:
  * closed-form identities the code must satisfy exactly (analytic checks)
  * internal consistency between two independently-coded routes to the same
    number (the spectrum vs. the range formula)
  * comparisons against published measurements (marked `literature`)
"""

import numpy as np
import pytest

from . import (CMS, ChirkinAtmospheric, CrossSection, Cylinder, Inelasticity,
               MOLASSE, ROCK_LOSS, Rock, STANDARD_ROCK, TabulatedFlux,
               arriving_muon_spectrum, background_muons, cos_theta_star,
               default_cc, detector_spectrum, earth_column_depth,
               integrated_muon_flux, muon_production_spectrum)
from .muon_transport import EnergyLoss

E_NU_RANGE = (10.0, 1.0e6)


# ---------------------------------------------------------------------------
# Muon transport: analytic identities
# ---------------------------------------------------------------------------

def test_range_and_energy_after_are_inverses():
    loss = ROCK_LOSS
    for e0 in (20.0, 100.0, 1e3, 1e5):
        for eth in (1.0, 10.0, 0.5 * e0):
            x = loss.range_gcm2(e0, eth)
            assert loss.energy_after(e0, x) == pytest.approx(eth, rel=1e-9)


def test_energy_before_inverts_energy_after():
    loss = ROCK_LOSS
    x = np.array([1e2, 1e4, 1e5])
    e_final = 50.0
    e0 = loss.energy_before(e_final, x)
    assert loss.energy_after(e0, x) == pytest.approx(e_final, rel=1e-9)


def test_range_is_zero_below_threshold():
    assert ROCK_LOSS.range_gcm2(5.0, 10.0) == 0.0


def test_low_energy_range_is_ionisation_only():
    """Well below the critical energy, R -> (E0 - Eth)/a."""
    loss = ROCK_LOSS
    e0, eth = 20.0, 1.0
    assert loss.range_gcm2(e0, eth) == pytest.approx((e0 - eth) / loss.a, rel=0.03)


def test_critical_energy():
    assert ROCK_LOSS.critical_energy == pytest.approx(500.0, rel=0.02)


# ---------------------------------------------------------------------------
# Inelasticity: normalisation and survival functions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model", ["qpm", "flat"])
@pytest.mark.parametrize("species", ["nu", "nubar"])
def test_inelasticity_pdf_normalised(model, species):
    inel = Inelasticity(model)
    u = np.linspace(0, 1, 20001)
    assert np.trapezoid(inel.pdf_u(u, species), u) == pytest.approx(1.0, rel=1e-6)


@pytest.mark.parametrize("species", ["nu", "nubar"])
def test_survival_is_the_pdf_integral(species):
    """S(z) = int_z^1 f_u(u) du -- the two must agree, since the spectrum uses
    S and the production spectrum uses f_u."""
    inel = Inelasticity("qpm")
    for z in (0.0, 0.2, 0.5, 0.9):
        u = np.linspace(z, 1.0, 20001)
        assert inel.survival(z, species) == pytest.approx(
            np.trapezoid(inel.pdf_u(u, species), u), rel=1e-5, abs=1e-9)


def test_qpm_mean_inelasticity():
    inel = Inelasticity("qpm")
    assert inel.mean_u("nu") == pytest.approx(0.5)      # <y> = 0.5
    assert inel.mean_u("nubar") == pytest.approx(0.75)  # <y> = 0.25


# ---------------------------------------------------------------------------
# Cross sections
# ---------------------------------------------------------------------------

def test_xsec_linear_in_measured_regime():
    """sigma/E must be flat at the PDG world-average value for 30-350 GeV."""
    s = default_cc("nu")
    e = np.array([30.0, 100.0, 350.0])
    ratio = s(e) / e
    assert ratio == pytest.approx(0.677e-38, rel=0.03)
    sb = default_cc("nubar")
    assert sb(e) / e == pytest.approx(0.334e-38, rel=0.03)


def test_xsec_high_energy_anchor():
    """sigma_CC(nu) at 10 TeV should sit at the NLO QCD value used by
    GQRS / CSMS and plotted in IceCube's Nature 551, 596 measurement."""
    assert default_cc("nu")(1.0e4) == pytest.approx(4.6e-35, rel=0.05)


def test_xsec_breaks_linearity_at_high_energy():
    """The W propagator must suppress sigma relative to linear scaling."""
    s = default_cc("nu")
    assert s(1e6) / 1e6 < 0.3 * s(100.0) / 100.0


def test_xsec_monotonic_and_positive():
    s = default_cc("nu")
    e = np.geomspace(10, 1e8, 500)
    v = s(e)
    assert np.all(v > 0)
    assert np.all(np.diff(v) > 0)


def test_xsec_roundtrip_csv(tmp_path):
    s = default_cc("nu")
    p = tmp_path / "x.csv"
    s.to_csv(p)
    s2 = CrossSection.from_csv(p)
    e = np.geomspace(20, 1e5, 50)
    assert s2(e) == pytest.approx(s(e), rel=1e-6)


def test_xsec_from_sigma_over_E_csv(tmp_path):
    """Digitising Formaggio & Zeller Fig 9 gives sigma/E, not sigma."""
    p = tmp_path / "y.csv"
    p.write_text("# E,sigma/E\n30,0.677e-38\n100,0.677e-38\n300,0.677e-38\n")
    s = CrossSection.from_csv(p, sigma_unit="cm2_per_GeV")
    assert s(100.0) == pytest.approx(0.677e-36, rel=1e-6)


# ---------------------------------------------------------------------------
# Flux
# ---------------------------------------------------------------------------

def test_cos_theta_star_limits():
    assert cos_theta_star(1.0) == pytest.approx(1.0, rel=1e-6)
    assert 0.0 < cos_theta_star(0.0) < 0.2
    assert cos_theta_star(-0.5) == pytest.approx(cos_theta_star(0.5))


def test_flux_horizontal_exceeds_vertical_at_high_energy():
    """Above the pion critical energy the flux is enhanced near the horizon."""
    f = ChirkinAtmospheric()
    assert f(1.0e4, 0.0) > 5.0 * f(1.0e4, 1.0)


def test_flux_falls_like_E_minus_2p7_at_low_energy():
    """Well below the pion critical energy the neutrino spectrum tracks the
    primary cosmic-ray spectrum, ~E^-2.7.  Tested at the horizon, where
    cos(theta*) ~ 0.1 keeps the meson-decay denominators near unity; at
    cos(zenith) = 0.5 the steepening has already begun by 20 GeV."""
    f = ChirkinAtmospheric()
    e = np.array([5.0, 10.0])
    slope = np.diff(np.log(f(e, 0.0))) / np.diff(np.log(e))
    assert slope == pytest.approx(-2.7, abs=0.1)


def test_tabulated_flux_interpolates_its_own_points():
    e = np.geomspace(100, 1e5, 10)
    phi = 1e-6 * e**-3
    t = TabulatedFlux(e, phi, nu_fraction=1.0)
    assert t(e[3], -0.5) == pytest.approx(phi[3], rel=1e-6)


def test_tabulated_zenith_shape_preserves_normalisation():
    """With zenith_shape_from, the hemisphere average must reproduce the table."""
    e = np.geomspace(100, 1e5, 12)
    phi = 1e-6 * e**-3
    model = ChirkinAtmospheric()
    t = TabulatedFlux(e, phi, zenith_shape_from=model, shape_hemisphere="up")
    cz = np.linspace(-1, -1e-6, 400)
    tot = sum(np.trapezoid(t(np.full_like(cz, 1000.0), cz, sp), cz) / (cz[-1] - cz[0])
              for sp in ("nu", "nubar"))
    assert tot == pytest.approx(np.interp(1000.0, e, phi), rel=0.02)


def test_earth_column_depth_scale():
    """Straight up through the Earth: ~1.1e10 g/cm^2; horizontal: ~0."""
    assert earth_column_depth(-1.0) == pytest.approx(1.1e10, rel=0.15)
    assert earth_column_depth(0.0) == 0.0
    assert earth_column_depth(0.5) == 0.0


def test_earth_absorbs_pev_neutrinos_from_below():
    """A vertical chord is ~1.1e10 g/cm^2.  At 1 PeV the interaction length is
    ~5e9 g/cm^2, so transmission is ~2 interaction lengths, i.e. ~10% -- the
    Earth is translucent, not black, at PeV.  Below ~1 TeV it is transparent."""
    from .nu_flux import earth_transmission
    s = default_cc("nu")
    assert 0.02 < earth_transmission(1e6, -1.0, s) < 0.25
    # even at 100 GeV a full Earth chord removes ~0.6%; at 10 GeV, ~0.06%
    assert earth_transmission(1e2, -1.0, s) == pytest.approx(0.994, abs=0.005)
    assert earth_transmission(1e1, -1.0, s) == pytest.approx(1.0, abs=1e-3)
    assert earth_transmission(1e6, 0.5, s) == pytest.approx(1.0)  # down-going


# ---------------------------------------------------------------------------
# THE key consistency test
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("e_thr", [10.0, 100.0, 1000.0])
@pytest.mark.parametrize("model", ["qpm", "flat"])
def test_spectrum_integral_matches_range_formula(e_thr, model):
    """Integrating dPhi/dE_mu from E_thr upward must reproduce the
    n_N sigma <R_mu> form.  These are two independently coded routes: one
    integrates the Jacobian 1/(a+bE) with the survival function, the other
    integrates the range directly.  Agreement is the strongest check that the
    arriving-spectrum derivation is right."""
    inel = Inelasticity(model)
    cz = np.array([-0.8, -0.3])
    e_mu = np.geomspace(e_thr, 1.0e6, 900)

    dphi = arriving_muon_spectrum(e_mu, cz, inelasticity=inel,
                                  e_nu_range=E_NU_RANGE, n_energy=800)
    integrated = np.trapezoid(dphi, e_mu, axis=1)

    direct = integrated_muon_flux(e_thr, cz, inelasticity=inel,
                                  e_nu_range=E_NU_RANGE, n_energy=800)
    assert integrated == pytest.approx(direct, rel=0.02)


def test_spectrum_integral_matches_with_finite_rock_path():
    """Same identity, but with the overburden cap active (down-going)."""
    cz = np.array([0.5])
    xmax = np.array([2.0e4])          # g/cm^2, ~75 m of rock
    e_mu = np.geomspace(10.0, 1.0e5, 900)
    dphi = arriving_muon_spectrum(e_mu, cz, max_path_gcm2=xmax,
                                  e_nu_range=E_NU_RANGE, n_energy=800)
    direct = integrated_muon_flux(10.0, cz, max_path_gcm2=xmax,
                                  e_nu_range=E_NU_RANGE, n_energy=800)
    assert np.trapezoid(dphi, e_mu, axis=1) == pytest.approx(direct, rel=0.03)


# ---------------------------------------------------------------------------
# Physical behaviour
# ---------------------------------------------------------------------------

def test_rate_is_independent_of_rock_density():
    """n_N ~ rho and R_mu ~ 1/rho, so the rate must not depend on density."""
    a = background_muons(rock=STANDARD_ROCK, uncertainty=False)["n_muons"]
    b = background_muons(rock=MOLASSE, uncertainty=False)["n_muons"]
    assert a == pytest.approx(b, rel=0.01)


def test_rate_decreases_with_threshold():
    prev = np.inf
    for e_thr in (10.0, 100.0, 1000.0):
        n = background_muons(e_mu_min=e_thr, uncertainty=False)["n_muons"]
        assert n < prev
        prev = n


def test_ignoring_inelasticity_overestimates():
    """E_mu = E_nu should inflate the rate by roughly a factor 2."""
    base = background_muons(uncertainty=False)["n_muons"]
    naive = background_muons(inelasticity=Inelasticity("none"),
                             uncertainty=False)["n_muons"]
    assert 1.5 < naive / base < 2.5


def test_rate_scales_with_livetime():
    a = background_muons(livetime_s=1e7, uncertainty=False)["n_muons"]
    b = background_muons(livetime_s=2e7, uncertainty=False)["n_muons"]
    assert b == pytest.approx(2 * a, rel=1e-6)


def test_rate_scales_with_detector_area():
    """Doubling both transverse dimensions should roughly double the rate."""
    small = Cylinder(radius=5.0, half_length=10.0)
    big = Cylinder(radius=10.0, half_length=10.0)
    a = background_muons(small, uncertainty=False)["n_muons"]
    b = background_muons(big, uncertainty=False)["n_muons"]
    # slightly super-linear: A_proj = 2RL sin(a) + pi R^2 |cos(a)|, so the
    # end-cap term grows as R^2 while the barrel term grows as R
    assert 2.0 < b / a < 2.6


def test_overburden_limits_vertical_downgoing_muons():
    """Straight down, only production energies within one overburden's worth of
    range above threshold contribute, so the flux is suppressed relative to the
    unlimited (up-going) case.  Give the detector a 1 km overburden instead and
    the suppression goes away.

    Note this is a *vertical* statement.  Integrated over the whole down-going
    hemisphere the suppression is mild (ratio ~0.7), because the 100 m depth is
    a vertical figure: near-horizontal down-going directions see effectively
    unlimited rock and also present the largest projected area.  Down-going is
    swamped by cosmic-ray muons regardless; this test is about the geometry.
    """
    e_thr = 500.0
    up = integrated_muon_flux(e_thr, np.array([-1.0]))[0]
    shallow = integrated_muon_flux(
        e_thr, np.array([1.0]),
        max_path_gcm2=CMS.max_rock_path_gcm2(np.array([1.0]), STANDARD_ROCK))[0]
    deep = integrated_muon_flux(
        e_thr, np.array([1.0]),
        max_path_gcm2=Cylinder(depth=1000.0).max_rock_path_gcm2(
            np.array([1.0]), STANDARD_ROCK))[0]
    assert shallow < 0.5 * up
    assert deep > 0.9 * up


def test_downgoing_over_upgoing_ratio_falls_with_threshold():
    """The overburden cap bites harder as the required muon range grows, so the
    hemisphere-integrated down/up ratio falls with threshold.

    KNOWN ODDITY, not fully explained: the ratio stops falling and turns back
    up somewhere between 1 and 5 TeV.  Three candidate causes have been tested
    and excluded -- it persists with Earth absorption disabled, with a
    perfectly isotropic flux (so it is not the horizon enhancement), and with
    the E_nu integral extended to 1e8 GeV (so it is not truncation).  The
    remaining suspect is the saturation of R_mu: above the critical energy
    R -> (1/b) ln(E/E_thr), which depends only on the *ratio* of energies, so
    <R> stops growing with E_thr while the overburden cap X stays fixed and
    min(R, X)/R tends to a constant.  Chase this before the numbers go in a
    paper.  Down-going is swamped by cosmic-ray muons anyway, so it does not
    affect the EarthShine result; see notebook 01.
    """
    def ratio(e_thr):
        a = background_muons(hemisphere="down", e_mu_min=e_thr,
                             uncertainty=False)["n_muons"]
        b = background_muons(hemisphere="up", e_mu_min=e_thr,
                             uncertainty=False)["n_muons"]
        return a / b
    assert ratio(1000.0) < ratio(10.0)


def test_production_spectrum_positive_and_falling():
    d = muon_production_spectrum(np.geomspace(10, 1e4, 40))
    y = d["dN_dEmu_per_m3_per_s"]
    assert np.all(y > 0)
    assert y[0] > 100 * y[-1]


# ---------------------------------------------------------------------------
# Literature comparisons
# ---------------------------------------------------------------------------

@pytest.mark.literature
def test_upgoing_muon_flux_matches_superk_macro():
    """Up-going through-going muon flux above ~1.6 GeV.

    Super-Kamiokande and MACRO both measure ~1.5-2 x 10^-13 cm^-2 s^-1 sr^-1
    averaged over the up-going hemisphere.  This is the single most useful
    external check available, because it convolves flux x sigma x range into
    one measured number.  Note the measurements include nu_mu -> nu_tau
    oscillation (a ~20% suppression for through-going muons) which this
    calculation does not apply, so agreement on the low side is expected.
    """
    cz = np.linspace(-1.0, -0.02, 40)
    phi = integrated_muon_flux(1.6, cz, e_nu_range=(1.0, 1e6), n_energy=600)
    mean = np.trapezoid(phi, cz) / (cz[-1] - cz[0])
    assert 1.0e-13 < mean < 2.5e-13


@pytest.mark.literature
def test_probability_of_detectable_muon_at_1_TeV():
    """P_mu = n_N sigma <R_mu>, the probability a TeV neutrino makes a muon
    that reaches a detector.  Standard treatments put this at ~1e-6."""
    inel = Inelasticity("qpm")
    u = np.linspace(1e-6, 1, 4001)
    r = np.trapezoid(inel.pdf_u(u, "nu") * ROCK_LOSS.range_gcm2(u * 1e3, 1.0), u)
    from .nu_background import N_AVOGADRO
    p = N_AVOGADRO * default_cc("nu")(1e3) * r
    assert 2e-7 < p < 3e-6


@pytest.mark.literature
def test_muon_range_matches_feng_smolinsky_tanedo():
    """1509.07525 Eq. (32): a 1 TeV muon with E_th = 50 GeV travels 2.5 km in
    a rho = 1 g/cm^3 medium."""
    from .muon_transport import WATER_LOSS
    r_cm = WATER_LOSS.range_cm(1.0e3, 50.0, density=1.0)
    assert r_cm / 1e5 == pytest.approx(2.5, rel=0.05)
