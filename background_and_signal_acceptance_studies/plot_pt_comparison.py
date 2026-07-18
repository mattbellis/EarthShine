"""
plot_pt_comparison.py
=====================
Muon pT-at-detector (with energy loss) comparison, CORE vs FLOATING,
one PDF per DM mass: 1, 10, 100, 1000 TeV.  Answers the paper's
"compare pT for both core and floating" request in the Momentum-resolution
subsection.

    venv/bin/python plot_pt_comparison.py --data-root data --plotdir figures_pdf/individual
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
    "figure.figsize": (6.0, 4.6),
    "figure.dpi": 150,
    "font.size": 17,
    "axes.labelsize": 20,
    "axes.titlesize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
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

# filename tag (no dot, to keep LaTeX graphicx happy) and display label
TEV = {100.0: "0p1TeV", 1000.0: "1TeV", 10000.0: "10TeV",
       100000.0: "100TeV", 1000000.0: "1000TeV"}
LABEL = {100.0: "0.1 TeV", 1000.0: "1 TeV", 10000.0: "10 TeV",
         100000.0: "100 TeV", 1000000.0: "1000 TeV"}
COL = {"core": "C0", "floating": "C1"}


def _sel(df, mass):
    df = eio.derive_columns(df)
    m, _ = eio.apply_selection(df, mass=mass, min_efinal=1.0)
    return df[m]


def compare_pt(cat, mass, outdir):
    series = {}
    for model in ("core", "floating"):
        ds = find_dataset(cat, model, mass)
        if ds is None:
            print(f"  [skip {mass}] no {model} dataset")
            return
        d = _sel(ds[0], mass)
        series[model] = d["pt1_detector_acceptance_eloss"].to_numpy()
        print(f"  {model}: {len(d)} events, median pT={np.median(series[model]):.0f} GeV")

    # common x-range from the 99th percentile of both models
    xmax = 1.05 * max(np.quantile(v, 0.99) for v in series.values())

    fig, ax = plt.subplots()
    for model in ("core", "floating"):
        ax.hist(series[model], bins=50, range=(0, xmax), density=True,
                histtype="step", color=COL[model], label=model)
    ax.set_xlabel(r"Muon $p_{T}$ at detector [GeV]")
    ax.set_ylabel("Normalized entries")
    ax.set_yscale("log")
    ax.legend(title=f"$m_\\chi$ = {LABEL[mass]}")
    fig.tight_layout()
    out = Path(outdir) / f"pTcompare_{TEV[mass]}.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    print(f"  -> {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data")
    ap.add_argument("--plotdir", default="figures_pdf/individual")
    ap.add_argument("--masses", type=float, nargs="+",
                    default=[100.0, 1000.0, 10000.0, 100000.0, 1000000.0])
    args = ap.parse_args()

    cat = eio.read_catalog(args.data_root)
    for mass in args.masses:
        print(f"[pT compare] m_chi = {TEV.get(mass, mass)}")
        compare_pt(cat, mass, args.plotdir)


if __name__ == "__main__":
    main()
