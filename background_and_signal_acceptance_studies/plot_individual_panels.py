"""
plot_individual_panels.py
=========================
Emit each panel of the Category-A paper figures as its own single-plot PDF
(instead of the merged multi-panel figures).  Panels reproduce exactly the
curves/ranges of the merged versions in diagnostics_v2.py.

Outputs (into --plotdir), for the paper's Figures/:

  ePtAng{Core,Floating}_{theta,eta,phi,Elin,Elog,pTlin,pTlog}.pdf   (7 each)
  openAngleSep{Core,Floating}_{angle,sep}.pdf                       (2 each)
  density_profile_{1,10,100}TeV_{lin,log}.pdf                       (2 each)
  emu_at_det_{E,pT}.pdf                                             (2)

    venv/bin/python plot_individual_panels.py --data-root data --plotdir figures_pdf/individual
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import numpy as np

import earthshine_io as eio
from plot_category_A_figures import find_dataset

# consolidated-overlay colours: one colour per DM model
MODEL_COLOR = {"core": "#1f77b4", "floating": "#d62728"}

# ----------------------------------------------------------------------------
# Shared "pretty" style: large fonts (these panels are shown small, 0.3-0.48
# column width, so text must be big to stay legible), light grid, thick lines.
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "figure.figsize": (5.4, 5.4),
    "figure.dpi": 150,
    "font.size": 17,
    "axes.labelsize": 20,
    "axes.titlesize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 15,
    "legend.frameon": True,
    "legend.framealpha": 0.85,
    "lines.linewidth": 3.6,
    "patch.linewidth": 3.6,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
    "savefig.bbox": "tight",
})

YLABEL = "Normalized entries"


def _panel():
    fig = plt.figure()
    return fig, plt.gca()


def _finish(ax, xlabel, ylabel=YLABEL, logy=False, logx=False, legend=False,
            xsize=None):
    ax.set_xlabel(xlabel, fontsize=xsize)
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    if logx:
        ax.set_xscale("log")
    if legend:
        ax.legend()


def _save(fig, outdir, name):
    out = Path(outdir) / f"{name}.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print(f"  -> {out}")


def _select(df, mass):
    df = eio.derive_columns(df)
    mask, _ = eio.apply_selection(df, mass=mass, min_efinal=1.0)
    return df[mask]


# ---------------------------------------------------------------- ePtAng
def eptang_panels(df, mass, prefix, outdir):
    d = _select(df, mass)
    print(f"{prefix}: {len(d)} events")
    xr = (0, 1.1 * mass / 2)

    fig, ax = _panel()
    d['theta1'].hist(bins=64, range=(0, 3.2), histtype="step", density=True)
    _finish(ax, r'Muon polar angle $\theta_{\mu}$ [rad]')
    _save(fig, outdir, f"{prefix}_theta")

    fig, ax = _panel()
    d['eta1'].hist(bins=64, range=(0, 3), density=True, histtype="step")
    _finish(ax, r'Muon pseudorapidity $\eta_{\mu}$', logy=True)
    _save(fig, outdir, f"{prefix}_eta")

    fig, ax = _panel()
    d['phi1'].hist(bins=64, range=(0, 3.2), density=True, histtype="step")
    _finish(ax, r'Muon azimuthal angle $\phi_{\mu}$ [rad]', logy=True)
    _save(fig, outdir, f"{prefix}_phi")

    for logy, tag in ((False, "Elin"), (True, "Elog")):
        fig, ax = _panel()
        d['e1'].hist(bins=50, range=xr, histtype="step", density=True,
                     label='At generation')
        d['efinal_mu1'].hist(bins=50, range=xr, histtype="step", density=True,
                             label='At detector')
        _finish(ax, r'Muon energy $E_{\mu}$ [GeV]', logy=logy, legend=True)
        _save(fig, outdir, f"{prefix}_{tag}")

    for logy, tag in ((False, "pTlin"), (True, "pTlog")):
        fig, ax = _panel()
        d['pt1_detector_acceptance'].hist(bins=50, range=xr, density=True,
                                          histtype="step", label='Ignoring eloss')
        d['pt1_detector_acceptance_eloss'].hist(bins=50, range=xr, density=True,
                                                histtype="step", label='With eloss')
        _finish(ax, r'Muon $p_{T}$ [GeV]', logy=logy, legend=True)
        _save(fig, outdir, f"{prefix}_{tag}")


# -------------------------------------------------- ePtAng overlay (models)
def eptang_overlay_panels(cat, mass, outdir):
    """Consolidated ePtAng figure: both DM models overlaid on the same axes.
    Line style encodes the production model (core = solid, floating = dashed),
    consistent with the other signal figures. For the energy and pT panels
    colour encodes the stage: at generation vs at the detector (after energy
    loss in the rock)."""

    # For displaying in a Jupyter notebook
    figs,axes=[],[]

    model_ls = {"core": "-", "floating": "--"}
    # this figure is a single mass -> reuse that mass's colour (green for
    # 10 TeV) for the nominal/at-generation curves so it matches the mass-coded
    # figures; the at-detector stage uses gold, a colour no mass uses.
    mass_c = MASS_COLOR.get(mass, "#2ca02c")
    gold = "#DAA520"
    ang_color = mass_c                    # angular panels: colour = the mass
    stage_color = {"gen": mass_c, "det": gold}
    data = {}
    for model in ("core", "floating"):
        ds = find_dataset(cat, model, mass)
        if ds:
            data[model] = _select(ds[0], mass)
            print(f"eptang_overlay {model}: {len(data[model])} events")
    if not data:
        return
    xr = (0, 1.1 * mass / 2)

    # --- angular panels: one curve per model, model = line style -------------
    ang = [("theta1", (0, 3.2), r'Muon polar angle $\theta_{\mu}$ [rad]',
            True, "theta"),
           ("eta1", (0, 3), r'Muon pseudorapidity $\eta_{\mu}$', True, "eta"),
           ("phi1", (0, 3.2), r'Muon azimuthal angle $\phi_{\mu}$ [rad]',
            True, "phi")]
    for col, rng, xlabel, logy, tag in ang:
        fig, ax = _panel()
        for model, d in data.items():
            d[col].hist(bins=64, range=rng, histtype="step", density=True,
                        color=ang_color, linestyle=model_ls[model],
                        label=model.capitalize())
        _finish(ax, xlabel, logy=logy, legend=True)
        _save(fig, outdir, f"eptang_{tag}")

        figs.append(fig), axes.append(ax)

    # --- energy / pT panels: model = line style, stage = colour --------------
    twolayer = [(("e1", "efinal_mu1"), r'Muon energy $E_{\mu}$ [GeV]', "E"),
                (("pt1_detector_acceptance", "pt1_detector_acceptance_eloss"),
                 r'Muon $p_{T}$ [GeV]', "pT")]
    for (gen_col, det_col), xlabel, tag in twolayer:
        fig, ax = _panel()
        for model, d in data.items():
            ls = model_ls[model]
            d[gen_col].hist(bins=50, range=xr, histtype="step", density=True,
                            color=stage_color["gen"], linestyle=ls)
            d[det_col].hist(bins=50, range=xr, histtype="step", density=True,
                            color=stage_color["det"], linestyle=ls)
        _finish(ax, xlabel, logy=True)
        model_h = [Line2D([0], [0], color="0.25", lw=3.6, ls="-", label="Core"),
                   Line2D([0], [0], color="0.25", lw=3.6, ls="--",
                          label="Floating")]
        stage_h = [Line2D([0], [0], color=stage_color["gen"], lw=3.6, ls="-",
                          label="At generation"),
                   Line2D([0], [0], color=stage_color["det"], lw=3.6, ls="-",
                          label="At detector")]
        ax.legend(handles=model_h + stage_h)
        _save(fig, outdir, f"eptang_{tag}")
        figs.append(fig), axes.append(ax)

    return figs,axes

# ------------------------------------------------------------ openAngleSep
def openanglesep_panels(df, mass, prefix, outdir,
                        opening_angle_max=1e-2, sep_max=30.0):
    df = eio.derive_columns(df).copy()
    # transverse muon-PAIR separation at the detector = opening angle x flight
    # distance (small-angle). NOT the ip_*0/ip_*1 distance, which is a single
    # muon's entry->exit chord through the detector cylinder.
    df['sep'] = (np.deg2rad(df['opening angle'])
                 * df['distance_to_detector'] * 1000.0)
    mask, _ = eio.apply_selection(df, mass=mass, min_efinal=1.0)
    d = df[mask]
    print(f"{prefix}: {len(d)} events")

    fig, ax = _panel()
    d['opening angle'].hist(bins=50, range=(0, opening_angle_max),
                            histtype="step", density=True)
    _finish(ax, r'Opening angle [deg]')
    _save(fig, outdir, f"{prefix}_angle")

    fig, ax = _panel()
    d['sep'].hist(bins=50, range=(0, sep_max), density=True, histtype="step")
    _finish(ax, r'Muon separation at detector [mm]', logy=True)
    _save(fig, outdir, f"{prefix}_sep")


# -------------------------------------------------- openAngleSep overlay
def openanglesep_overlay_panels(cat, masses, outdir):
    """Consolidated opening-angle / separation figure overlaying several DM
    masses.  Colour encodes the mass; line style encodes the production model
    (core = solid, floating = dashed).  The opening angle scales as ~1/E, so
    the distributions span decades across masses -> use log x-axes."""
    figs,axes = [],[]
    model_ls = {"core": "-", "floating": "--"}
    # gather selected data for every (mass, model) present in the catalog
    data = {}
    for m in masses:
        for model in ("core", "floating"):
            ds = find_dataset(cat, model, m)
            if not ds:
                continue
            df = eio.derive_columns(ds[0]).copy()
            # transverse muon-PAIR separation = opening angle x flight distance
            df['sep'] = (np.deg2rad(df['opening angle'])
                         * df['distance_to_detector'] * 1000.0)
            mask, _ = eio.apply_selection(df, mass=m, min_efinal=1.0)
            data[(m, model)] = df[mask]
            print(f"openanglesep_overlay {MASS_TEV[m]} {model}: "
                  f"{len(data[(m, model)])} events")
    if not data:
        return
    present = [m for m in masses if any((m, mo) in data
                                        for mo in ("core", "floating"))]

    def _panel_for(col, xlabel, name):
        vals = np.concatenate([d[col].to_numpy() for d in data.values()])
        vals = vals[vals > 0]
        lo, hi = np.quantile(vals, 0.001), vals.max()
        bins = np.logspace(np.log10(lo), np.log10(hi), 50)
        fig, ax = _panel()
        for (m, model), d in data.items():
            d[col].hist(bins=bins, histtype="step", density=True,
                        color=MASS_COLOR[m], linestyle=model_ls[model])
        mass_h = [Line2D([0], [0], color=MASS_COLOR[m], lw=3.6,
                         label=MASS_TEV[m]) for m in present]
        style_h = [Line2D([0], [0], color="0.25", lw=3.6, ls="-", label="Core"),
                   Line2D([0], [0], color="0.25", lw=3.6, ls="--",
                          label="Floating")]
        ax.legend(handles=mass_h + style_h, ncol=2, fontsize=13)
        _finish(ax, xlabel, logy=True, logx=True)
        _save(fig, outdir, name)

        figs.append(fig), axes.append(ax)

    _panel_for('opening angle', r'Opening angle [deg]', "openAngleSep_angle")
    _panel_for('sep', r'Muon separation at detector [mm]', "openAngleSep_sep")

    return figs,axes

# ------------------------------------------------------------ density depth
def density_panels(df, mass, prefix, outdir):
    d = _select(df, mass)
    print(f"{prefix}: {len(d)} events")
    # y0 is the decay's vertical coordinate measured from the detector at y=0;
    # decays happen below the detector (y0 < 0), so the depth below the
    # detector is -y0 (a positive distance). Use it in both panels so the two
    # axes share the same convention.
    depth = -d['y0']
    yrange = 4000.0

    fig, ax = _panel()
    ax.hist(depth, bins=100, range=(0, yrange), density=True,
            histtype="step")
    _finish(ax, r'Depth below detector [m]')
    _save(fig, outdir, f"{prefix}_lin")

    fig, ax = _panel()
    bins = np.logspace(np.log10(1), np.log10(4000), 34)
    ax.hist(depth, bins=bins, density=True, histtype="step")
    _finish(ax, r'Depth below detector [m] (log)', logx=True, xsize=16)
    _save(fig, outdir, f"{prefix}_log")


# ------------------------------------------------ density depth overlay (mass)
MASS_TEV = {100.0: "0.1 TeV", 1000.0: "1 TeV", 10000.0: "10 TeV",
            100000.0: "100 TeV"}
MASS_COLOR = {100.0: "#9467bd", 1000.0: "#1f77b4", 10000.0: "#2ca02c",
              100000.0: "#d62728"}


def density_overlay_panels(cat, masses, outdir, dm_model='floating'):
    """Consolidated decay-depth figure: all DM masses overlaid on a single
    linear panel and a single log panel (colour encodes the DM mass)."""

    figs,axes = [],[]

    data = {}
    for m in masses:
        #ds = find_dataset(cat, "floating", m)
        ds = find_dataset(cat, dm_model, m)
        if ds:
            d = _select(ds[0], m)
            data[m] = -d['y0']   # depth below detector (positive distance)
            print(f"density_overlay {MASS_TEV[m]}: {len(d)} events")
    if not data:
        return
    yrange = 4000.0

    fig, ax = _panel()
    for m, depth in data.items():
        depth.hist(bins=100, range=(0, yrange), density=True, histtype="step",
                   color=MASS_COLOR[m], label=MASS_TEV[m])
    _finish(ax, r'Depth below detector [m]', legend=True)
    _save(fig, outdir, f"density_profile_lin_{dm_model}")
    figs.append(fig), axes.append(ax)

    bins = np.logspace(np.log10(1), np.log10(4000), 34)
    fig, ax = _panel()
    for m, depth in data.items():
        depth.hist(bins=bins, density=True, histtype="step",
                   color=MASS_COLOR[m], label=MASS_TEV[m])
    _finish(ax, r'Depth below detector [m] (log)', logx=True,
            legend=True, xsize=16)
    _save(fig, outdir, f"density_profile_log_{dm_model}")
    figs.append(fig), axes.append(ax)

    return figs,axes

# ------------------------------------------------------------ emu_at_det
def emu_panels(df, mass, outdir):
    d = _select(df, mass)
    print(f"emu_at_det: {len(d)} events")
    xr = (0, 1.5 * mass / 2)

    # stage colours consistent with the ePtAng figure: 10 TeV green at
    # generation, gold at the detector (after energy loss)
    gen_c, det_c = MASS_COLOR.get(mass, "#2ca02c"), "#DAA520"

    figs,axes = [],[]

    fig, ax = _panel()
    d['e1'].hist(bins=50, range=xr, density=True, histtype="step",
                 color=gen_c, label='At generation')
    d['efinal_mu1'].hist(bins=50, range=xr, density=True, histtype="step",
                         color=det_c, label='At detector')
    _finish(ax, r'Muon energy $E_{\mu}$ [GeV]', logy=True, legend=True)
    _save(fig, outdir, "emu_at_det_E")
    figs.append(fig)
    axes.append(ax)

    fig, ax = _panel()
    d['pt1_detector_acceptance'].hist(bins=50, range=xr, density=True,
                                      histtype="step", color=gen_c,
                                      label='Ignoring eloss')
    d['pt1_detector_acceptance_eloss'].hist(bins=50, range=xr, density=True,
                                            histtype="step", color=det_c,
                                            label='With eloss')
    _finish(ax, r'Muon $p_{T}$ [GeV]', logy=True, legend=True)
    _save(fig, outdir, "emu_at_det_pT")
    figs.append(fig)
    axes.append(ax)

    return figs,axes


# ------------------------------------------------------ pT overlap (core vs floating)
TEV = {1000: "1TeV", 10000: "10TeV", 100000: "100TeV"}


def pt_overlap_panel(cat, mass, outdir):
    """Overlay the muon pT-at-detector (with eloss) for the core and floating
    models at one mass -- the momentum-resolution comparison figure."""
    figs,axes = [],[]
    xr = (0, 1.1 * mass / 2)
    fig, ax = _panel()
    n = 0
    for model in ("core", "floating"):
        ds = find_dataset(cat, model, mass)
        if ds is None:
            continue
        d = _select(ds[0], mass)
        d['pt1_detector_acceptance_eloss'].hist(bins=50, range=xr, density=True,
                                                histtype="step", label=model)
        n += 1
        print(f"  pt_overlap {TEV[int(mass)]} {model}: {len(d)} events")
    if n:
        ax.set_title(rf'$m_\chi$ = {TEV[int(mass)]}')
        _finish(ax, r'Muon $p_{T}$ at detector [GeV]', logy=True, legend=True)
        _save(fig, outdir, f"pt_overlap_{TEV[int(mass)]}")
        figs.append(fig), axes.append(ax)

    return figs,axes

def pt_overlay_all(cat, masses, outdir):
    """Single consolidated momentum-resolution figure: the muon pT-at-detector
    (with eloss) for every DM mass, colour-coded by mass, with the core model
    as a solid line and the floating model as a dashed line. Log-x so the very
    different pT scales of the light and heavy masses all fit on one axis."""
    figs,axes = [],[]
    masses = [m for m in masses if find_dataset(cat, "core", m)
              or find_dataset(cat, "floating", m)]
    if not masses:
        return
    xr = (1.0, 1.1 * max(masses) / 2)
    bins = np.logspace(np.log10(xr[0]), np.log10(xr[1]), 50)
    style = {"core": "-", "floating": "--"}
    fig, ax = _panel()
    present = []
    for m in masses:
        for model in ("core", "floating"):
            ds = find_dataset(cat, model, m)
            if ds is None:
                continue
            d = _select(ds[0], m)
            d['pt1_detector_acceptance_eloss'].hist(
                bins=bins, density=True, histtype="step",
                color=MASS_COLOR[m], linestyle=style[model])
            print(f"  pt_overlay {MASS_TEV[m]} {model}: {len(d)} events")
        present.append(m)
    mass_h = [Line2D([0], [0], color=MASS_COLOR[m], lw=3.6, label=MASS_TEV[m])
              for m in present]
    style_h = [Line2D([0], [0], color="0.25", lw=3.6, ls="-", label="Core"),
               Line2D([0], [0], color="0.25", lw=3.6, ls="--", label="Floating")]
    ax.legend(handles=mass_h + style_h, ncol=2, fontsize=13)
    _finish(ax, r'Muon $p_{T}$ at detector [GeV]', logy=True, logx=True)
    _save(fig, outdir, "pt_overlay_all")
    figs.append(fig), axes.append(ax)

    return figs,axes


# ------------------------------------------------------ muon transit time
C_M_PER_NS = 0.299792458   # speed of light, m/ns
M_MU = 0.1056583755        # muon mass, GeV


def _transit_time_ns(d):
    """Entry->exit crossing time: chord between the muon's detector entry and
    exit points (ip_*0 -> ip_*1) divided by beta c, with beta from the muon
    energy at the detector."""
    chord = np.sqrt((d.ip_x0 - d.ip_x1) ** 2
                    + (d.ip_y0 - d.ip_y1) ** 2
                    + (d.ip_z0 - d.ip_z1) ** 2).to_numpy()      # m
    E = d.efinal_mu1.to_numpy()
    beta = np.sqrt(np.clip(E ** 2 - M_MU ** 2, 0.0, None)) / E
    return chord / (beta * C_M_PER_NS)


def transit_time_overlay_panel(cat, masses, outdir, tmax=110.0):
    """Consolidated transit-time figure in the style of pt_overlay_all: the time
    a signal muon takes to cross the detector (entry->exit chord over beta c) for
    every DM mass, colour-coded by mass, core = solid, floating = dashed. The
    muons are ultra-relativistic (beta ~ 1) so the curves are mass-independent;
    the shape is set by the crossing geometry -- vertical (core) tracks cut off
    at the detector diameter, inclined (floating) tracks tail out toward the
    full length."""
    figs,axes = [],[]
    style = {"core": "-", "floating": "--"}
    bins = np.linspace(0, tmax, 60)
    fig, ax = _panel()
    present = []
    for m in masses:
        drawn = False
        for model in ("core", "floating"):
            ds = find_dataset(cat, model, m)
            if ds is None:
                continue
            t = _transit_time_ns(_select(ds[0], m))
            ax.hist(t, bins=bins, density=True, histtype="step",
                    color=MASS_COLOR[m], linestyle=style[model])
            print(f"  transit {MASS_TEV[m]} {model}: {len(t)} events, "
                  f"median={np.median(t):.1f} ns")
            drawn = True
        if drawn:
            present.append(m)
    mass_h = [Line2D([0], [0], color=MASS_COLOR[m], lw=3.6, label=MASS_TEV[m])
              for m in present]
    style_h = [Line2D([0], [0], color="0.25", lw=3.6, ls="-", label="Core"),
               Line2D([0], [0], color="0.25", lw=3.6, ls="--", label="Floating")]
    # give ~50% headroom above the tallest bar so the legend sits at the top
    # without overlapping the histogram
    ax.set_ylim(top=ax.get_ylim()[1] * 1.5)
    ax.legend(handles=mass_h + style_h, ncol=2, fontsize=13, loc="upper center")
    _finish(ax, r'Muon transit time [ns]')
    _save(fig, outdir, "muon_transit_time_vs_mass")
    figs.append(fig), axes.append(ax)

    return figs,axes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data")
    ap.add_argument("--plotdir", default="figures_pdf/individual")
    ap.add_argument("--ptang-mass", type=float, default=10000.0)
    ap.add_argument("--density-masses", type=float, nargs="+",
                    default=[1000.0, 10000.0, 100000.0])
    ap.add_argument("--emu-mass", type=float, default=10000.0)
    args = ap.parse_args()

    cat = eio.read_catalog(args.data_root)

    print("ePtAng overlay (core vs floating):")
    eptang_overlay_panels(cat, args.ptang_mass, args.plotdir)

    print("openAngleSep overlay (all masses, core vs floating):")
    openanglesep_overlay_panels(cat, [100.0, 1000.0, 10000.0, 100000.0],
                                args.plotdir)

    print("density depth overlay (all masses):")
    density_overlay_panels(cat, args.density_masses, args.plotdir)

    ds = find_dataset(cat, "floating", args.emu_mass)
    if ds:
        emu_panels(ds[0], args.emu_mass, args.plotdir)

    print("pt overlay (all masses, core vs floating):")
    pt_overlay_all(cat, [100.0, 1000.0, 10000.0, 100000.0], args.plotdir)

    print("transit time overlay (all masses, core vs floating):")
    transit_time_overlay_panel(cat, args.density_masses, args.plotdir)

    print(f"\nDone.  Individual panels in {args.plotdir}/")


if __name__ == "__main__":
    main()
