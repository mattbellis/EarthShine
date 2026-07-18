"""
plot_transit_time.py
====================
Time for a signal muon to cross the cylindrical detector, entry to exit.
The crossing path is the chord between the muon's detector entry and exit
points (ip_*0 -> ip_*1), and the time is chord / (beta c) with beta taken
from the muon energy at the detector.

Because the muons are ultra-relativistic (beta ~ 1) at every DM mass, the
transit time is essentially independent of m_chi; the shape is instead set
by the incidence angle, so it is compared here between the (near-vertical)
core model and the (inclined) floating model.

    venv/bin/python plot_transit_time.py --data-root data --plotdir figures_pdf/individual
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt
import numpy as np

import earthshine_io as eio
from plot_category_A_figures import find_dataset

plt.rcParams.update({
    "figure.figsize": (6.4, 4.8), "figure.dpi": 150,
    "font.size": 17, "axes.labelsize": 20, "axes.titlesize": 18,
    "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 16,
    "legend.frameon": True, "legend.framealpha": 0.85,
    "lines.linewidth": 3.6, "patch.linewidth": 3.6,
    "axes.grid": True, "axes.axisbelow": True,
    "grid.alpha": 0.3, "grid.linestyle": ":", "savefig.bbox": "tight",
})

C_M_PER_NS = 0.299792458   # speed of light, m/ns
M_MU = 0.1056583755        # muon mass, GeV
COL = {"core": "C0", "floating": "C1"}


def transit_time_ns(d):
    chord = np.sqrt((d.ip_x0 - d.ip_x1) ** 2
                    + (d.ip_y0 - d.ip_y1) ** 2
                    + (d.ip_z0 - d.ip_z1) ** 2).to_numpy()      # m
    E = d.efinal_mu1.to_numpy()
    beta = np.sqrt(np.clip(E ** 2 - M_MU ** 2, 0.0, None)) / E
    return chord / (beta * C_M_PER_NS)


def plot_by_mass(cat, masses, tmax, plotdir):
    """Overlay the transit-time distribution for several DM masses, one panel
    per model.  Demonstrates that the transit time is set by the crossing
    geometry (chord length), not by m_chi: because the muons are
    ultra-relativistic (beta ~ 1) at every mass, the curves nearly coincide,
    and the only visible shift comes from the incidence angle (core = vertical,
    floating = inclined)."""
    cmap = plt.get_cmap("viridis")
    colors = {m: cmap(i / max(len(masses) - 1, 1)) for i, m in enumerate(masses)}

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), sharey=True)
    for ax, model in zip(axes, ("core", "floating")):
        for mass in masses:
            ds = find_dataset(cat, model, mass)
            if ds is None:
                print(f"[skip] no {model} dataset at {mass}")
                continue
            df = eio.derive_columns(ds[0])
            m, _ = eio.apply_selection(df, mass=mass, min_efinal=1.0)
            if m.sum() == 0:
                print(f"[skip] {model} {mass}: no events after selection")
                continue
            t = transit_time_ns(df[m])
            print(f"{model:9s} m={mass/1000:>7.1f} TeV  N={m.sum():>6d}  "
                  f"median={np.median(t):5.1f} ns  p95={np.quantile(t,.95):5.1f} ns")
            label = (f"{mass/1e6:g} PeV" if mass >= 1e6 else
                     f"{mass/1e3:g} TeV" if mass >= 1e3 else f"{mass:g} GeV")
            ax.hist(t, bins=60, range=(0, tmax), density=True,
                    histtype="step", color=colors[mass], label=label)
        ax.set_title(f"{model} model")
        ax.set_xlabel("Muon transit time across detector [ns]")
        ax.legend(title=r"$m_\chi$", ncol=1)
    axes[0].set_ylabel("Normalized entries")
    fig.tight_layout()
    out = Path(plotdir) / "muon_transit_time_vs_mass.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"  -> {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data")
    ap.add_argument("--plotdir", default="figures_pdf/individual")
    ap.add_argument("--mass", type=float, default=10000.0,
                    help="representative DM mass (GeV); transit time is nearly "
                         "mass-independent")
    ap.add_argument("--by-mass", action="store_true",
                    help="overlay several DM masses (one panel per model) "
                         "instead of comparing the two models at one mass")
    ap.add_argument("--masses", type=float, nargs="+",
                    default=[1000.0, 10000.0, 100000.0, 1000000.0],
                    help="DM masses (GeV) to overlay with --by-mass")
    ap.add_argument("--tmax", type=float, default=120.0)
    args = ap.parse_args()

    cat = eio.read_catalog(args.data_root)

    if args.by_mass:
        plot_by_mass(cat, args.masses, args.tmax, args.plotdir)
        return

    fig, ax = plt.subplots()
    for model in ("core", "floating"):
        ds = find_dataset(cat, model, args.mass)
        if ds is None:
            print(f"[skip] no {model} dataset at {args.mass}")
            continue
        df = eio.derive_columns(ds[0])
        m, _ = eio.apply_selection(df, mass=args.mass, min_efinal=1.0)
        t = transit_time_ns(df[m])
        print(f"{model}: N={len(t)} median={np.median(t):.1f} ns "
              f"p95={np.quantile(t,.95):.1f} ns max={t.max():.1f} ns")
        ax.hist(t, bins=60, range=(0, args.tmax), density=True,
                histtype="step", color=COL[model], label=model)

    ax.set_xlabel("Muon transit time across detector [ns]")
    ax.set_ylabel("Normalized entries")
    ax.legend()
    fig.tight_layout()
    out = Path(args.plotdir) / "muon_transit_time.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
