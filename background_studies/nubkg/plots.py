"""
plots.py -- diagnostic and paper figures.

Every function takes an optional `ax` and returns it, so figures can be
composed.  Nothing here computes physics; it all calls the other modules.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .nu_background import (CMS, STANDARD_ROCK, SEC_PER_YEAR,
                            arriving_muon_spectrum, detector_spectrum,
                            integrated_muon_flux, muon_production_spectrum)
from .nu_flux import ChirkinAtmospheric
from .nu_xsec import default_cc, sigma_uncertainty

__all__ = [
    "plot_flux", "plot_flux_e2", "plot_cross_sections", "plot_inelasticity",
    "plot_muon_range", "plot_rock_yield", "plot_arriving_spectrum",
    "plot_zenith_dependence", "plot_integrand_decomposition",
    "plot_threshold_scan", "plot_error_budget",
]

_KW = dict(lw=1.6)


def _finish(ax, xlabel, ylabel, title=None, logx=True, logy=True):
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.25, which="both")
    return ax


# ---------------------------------------------------------------------------
# 1. Fluxes
# ---------------------------------------------------------------------------

def plot_flux(fluxes, e_nu=None, cos_zenith=-0.5, ax=None, species=("nu", "nubar")):
    """dPhi/dE vs E_nu for one or more flux models."""
    e_nu = np.geomspace(10, 1e6, 300) if e_nu is None else e_nu
    ax = ax or plt.gca()
    for f in np.atleast_1d(fluxes):
        tot = sum(f(e_nu, cos_zenith, s) for s in species)
        ax.plot(e_nu, tot, label=getattr(f, "label", str(f)), **_KW)
    _finish(ax, r"$E_\nu$ [GeV]",
            r"$d\Phi/dE_\nu$ [GeV$^{-1}$cm$^{-2}$s$^{-1}$sr$^{-1}$]",
            f"muon neutrino flux, cos(zenith) = {cos_zenith}")
    ax.legend(fontsize=8)
    return ax


def plot_flux_e2(fluxes, e_nu=None, cos_zenith=-0.5, ax=None,
                 species=("nu", "nubar"), tabulated=None):
    """E^2 dPhi/dE -- the presentation used in the IceCube unfolding figure,
    which is the one to compare digitised points against."""
    e_nu = np.geomspace(1, 1e7, 400) if e_nu is None else e_nu
    ax = ax or plt.gca()
    for f in np.atleast_1d(fluxes):
        tot = sum(f(e_nu, cos_zenith, s) for s in species)
        ax.plot(e_nu, e_nu**2 * tot, label=getattr(f, "label", str(f)), **_KW)
    if tabulated is not None:
        t = tabulated
        yerr = None if t.flux_err is None else t.energy**2 * t.flux_err
        ax.errorbar(t.energy, t.energy**2 * t.flux, yerr=yerr, fmt="o", ms=4,
                    color="crimson", capsize=2, label=t.label, zorder=5)
    _finish(ax, r"$E_\nu$ [GeV]",
            r"$E_\nu^2\,\Phi_\nu$ [GeV cm$^{-2}$s$^{-1}$sr$^{-1}$]",
            f"cos(zenith) = {cos_zenith}")
    ax.legend(fontsize=8)
    return ax


def plot_zenith_dependence(flux=None, energies=(30, 300, 3000, 3e4), ax=None):
    flux = flux or ChirkinAtmospheric()
    ax = ax or plt.gca()
    cz = np.linspace(-1, 1, 300)
    for e in energies:
        tot = sum(flux(np.full_like(cz, e), cz, s) for s in ("nu", "nubar"))
        ax.plot(cz, tot / tot[len(cz) // 2], label=f"{e:g} GeV", **_KW)
    _finish(ax, r"$\cos\theta_{\rm zenith}$", "flux / flux(horizontal)",
            "horizon enhancement", logx=False)
    ax.legend(fontsize=8)
    return ax


# ---------------------------------------------------------------------------
# 2. Cross sections and inelasticity
# ---------------------------------------------------------------------------

def plot_cross_sections(xsecs=None, e_nu=None, ax=None, per_energy=False,
                        band=True):
    """sigma_CC (or sigma/E) vs E_nu, with the uncertainty band."""
    e_nu = np.geomspace(10, 1e7, 400) if e_nu is None else e_nu
    ax = ax or plt.gca()
    xsecs = xsecs or {s: default_cc(s) for s in ("nu", "nubar")}
    if not isinstance(xsecs, dict):
        xsecs = {"": xsecs}
    for name, s in xsecs.items():
        y = s(e_nu) / (e_nu if per_energy else 1.0)
        line, = ax.plot(e_nu, y, label=f"{name} {s.label}".strip(), **_KW)
        if band:
            f = sigma_uncertainty(e_nu)
            ax.fill_between(e_nu, y * (1 - f), y * (1 + f), alpha=0.18,
                            color=line.get_color(), lw=0)
    if per_energy:
        ax.axhline(0.677e-38, ls=":", c="k", lw=1)
        ax.axhline(0.334e-38, ls=":", c="gray", lw=1)
    _finish(ax, r"$E_\nu$ [GeV]",
            r"$\sigma_{\rm CC}/E_\nu$ [cm$^2$/GeV]" if per_energy
            else r"$\sigma_{\rm CC}$ [cm$^2$ / nucleon]",
            "CC cross section per nucleon (isoscalar)")
    ax.legend(fontsize=8)
    return ax


def plot_inelasticity(inel, ax=None):
    ax = ax or plt.gca()
    u = np.linspace(0, 1, 400)
    for sp, ls in (("nu", "-"), ("nubar", "--")):
        ax.plot(u, inel.pdf_u(u, sp), ls,
                label=rf"{sp}:  $\langle u\rangle$ = {inel.mean_u(sp):.2f}", **_KW)
    _finish(ax, r"$u = E_\mu/E_\nu = 1-y$", r"$f(u)$",
            f"inelasticity, model = '{inel.model}'", logx=False, logy=False)
    ax.legend(fontsize=8)
    return ax


def plot_muon_range(losses, e_thr=(1.0, 10.0, 100.0), density=2.65, ax=None):
    ax = ax or plt.gca()
    e = np.geomspace(10, 1e6, 300)
    for loss in np.atleast_1d(losses):
        for th in e_thr:
            ax.plot(e, loss.range_cm(e, th, density) / 1e5,
                    label=f"{loss.label.split('(')[0].strip()}, "
                          rf"$E_{{\rm th}}$={th:g} GeV", **_KW)
    _finish(ax, r"$E_\mu$ at production [GeV]", "range [km]",
            f"CSDA muon range, rho = {density} g/cm$^3$")
    ax.legend(fontsize=7)
    return ax


# ---------------------------------------------------------------------------
# 3. Rock yield and arriving spectrum
# ---------------------------------------------------------------------------

def plot_rock_yield(e_mu=None, ax=None, per_year=True, **kw):
    """Muons produced per cubic metre of rock, differential in muon energy."""
    e_mu = np.geomspace(10, 1e5, 120) if e_mu is None else e_mu
    d = muon_production_spectrum(e_mu, **kw)
    scale = SEC_PER_YEAR if per_year else 1.0
    ax = ax or plt.gca()
    ax.plot(d["e_mu"], d["e_mu"] * d["dN_dEmu_per_m3_per_s"] * scale,
            color="k", label="total", **_KW)
    for sp, ls in (("nu", "--"), ("nubar", ":")):
        ax.plot(d["e_mu"], d["e_mu"] * d["per_species"][sp] * scale, ls,
                label=sp, lw=1.2)
    unit = "yr" if per_year else "s"
    _finish(ax, r"$E_\mu$ at production [GeV]",
            rf"$E_\mu\,dN/dE_\mu$  [m$^{{-3}}$ {unit}$^{{-1}}$]",
            "muon production in 1 m$^3$ of rock")
    ax.legend(fontsize=8)
    return ax, d


def plot_arriving_spectrum(result, ax=None, per_year=True, cumulative=False):
    """Arriving muon spectrum at the detector, from `background_muons`."""
    ax = ax or plt.gca()
    e, dn = result["e_mu"], result["dN_dEmu_per_s"]
    scale = result["metadata"]["livetime_s"] if per_year else 1.0
    if cumulative:
        y = np.concatenate([[np.trapezoid(dn, e)],
                            np.trapezoid(dn, e) - np.cumsum(
                                0.5 * (dn[1:] + dn[:-1]) * np.diff(e))])
        ax.plot(e, y * scale, color="crimson", **_KW)
        ylab = "N($>E_\\mu$)"
    else:
        ax.plot(e, e * dn * scale, color="crimson", **_KW)
        ylab = r"$E_\mu\,dN/dE_\mu$"
    _finish(ax, r"$E_\mu$ at the detector [GeV]", ylab,
            f"background muons, {result['metadata']['detector_name']}, "
            f"{scale/SEC_PER_YEAR:.2f} yr")
    return ax


def plot_integrand_decomposition(e_thr=10.0, cos_zenith=-0.5, ax=None,
                                 e_nu_range=(10.0, 1e6), **kw):
    """Where in E_nu does the background actually come from?

    Plots the integrand of Phi_mu(>E_thr) per decade of E_nu, plus the
    cumulative fraction.  This is the plot that tells you which part of the
    flux and cross-section tables the answer is actually sensitive to.
    """
    ax = ax or plt.gca()
    e_nu = np.geomspace(*e_nu_range, 300)
    # finite-difference the integral w.r.t. the lower E_nu limit
    vals = np.array([integrated_muon_flux(e_thr, np.array([cos_zenith]),
                                          e_nu_range=(x, e_nu_range[1]),
                                          n_energy=250, **kw)[0] for x in e_nu])
    total = vals[0]
    frac = 1.0 - vals / total
    dens = -np.gradient(vals, np.log10(e_nu)) / total
    ax.plot(e_nu, dens, color="navy", label="fraction per decade", **_KW)
    ax2 = ax.twinx()
    ax2.plot(e_nu, frac, color="darkorange", ls="--", lw=1.3,
             label="cumulative")
    ax2.set_ylabel("cumulative fraction")
    ax2.set_ylim(0, 1.02)
    _finish(ax, r"$E_\nu$ [GeV]", r"$d(\Phi_\mu)/d\log_{10}E_\nu\ /\ \Phi_\mu$",
            rf"origin of the background, $E_\mu > {e_thr:g}$ GeV,"
            rf" $\cos\theta$ = {cos_zenith}", logy=False)
    for q in (0.5, 0.9):
        ax.axvline(np.interp(q, frac, e_nu), color="gray", ls=":", lw=1)
    return ax, {"e_nu": e_nu, "fraction_per_decade": dens, "cumulative": frac,
                "median_E_nu": float(np.interp(0.5, frac, e_nu)),
                "e_nu_90pct": float(np.interp(0.9, frac, e_nu))}


def plot_threshold_scan(thresholds=None, detector=CMS, ax=None,
                        livetime_s=SEC_PER_YEAR, **kw):
    """Integrated background count vs muon energy threshold, with errors."""
    from .nu_background import background_muons
    thresholds = thresholds if thresholds is not None else np.geomspace(10, 5e3, 12)
    ax = ax or plt.gca()
    n, err = [], []
    for t in thresholds:
        r = background_muons(detector, e_mu_min=float(t),
                             livetime_s=livetime_s, **kw)
        n.append(r["n_muons"]); err.append(r["n_muons_err"])
    n, err = np.array(n), np.array(err)
    ax.plot(thresholds, n, "o-", color="crimson", ms=4, **_KW)
    ax.fill_between(thresholds, np.maximum(n - err, 1e-6), n + err,
                    alpha=0.2, color="crimson", lw=0)
    _finish(ax, r"$E_\mu$ threshold [GeV]",
            f"background muons in {livetime_s/SEC_PER_YEAR:.2f} yr",
            f"{detector.name}: R={detector.radius} m, "
            f"L={2*detector.half_length} m, depth={detector.depth} m")
    return ax, {"thresholds": thresholds, "n": n, "err": err}


def plot_error_budget(result, ax=None):
    ax = ax or plt.gca()
    b = result["error_budget"]
    keys = sorted(b, key=lambda k: -b[k])
    ax.barh(keys, [b[k] for k in keys], color="steelblue")
    ax.set_xlabel(f"contribution to sigma(N)  [N = {result['n_muons']:.1f}]")
    ax.set_title("error budget", fontsize=10)
    ax.grid(alpha=0.25, axis="x")
    return ax
