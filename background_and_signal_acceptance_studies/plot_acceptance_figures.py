"""
plot_acceptance_figures.py
==========================
Geometric-acceptance-vs-mass figure (paper label detectorAcceptanceCoreFloating),
for core and floating, as vector PDFs.

Reuses the acceptance_v2 logic from PROTOTYPE_calculate_efficiencies.ipynb
(counts and generated totals summed across all files at a given (volume, mass)
before any fraction is formed; volume read from the parquet footer) and the
scatter-plot styling from that notebook's plotting cell.

    venv/bin/python plot_acceptance_figures.py --data-root data --plotdir figures_pdf
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt

import earthshine_io as eio

plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 15,
    "axes.labelsize": 18,
    "axes.titlesize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 12,
    "legend.frameon": True,
    "legend.framealpha": 0.85,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
    "savefig.bbox": "tight",
})

FOOTER_KEY = b"earthshine"
_PATH_RE = re.compile(r"depth=(m?[0-9p]+)-(m?[0-9p]+)_diskR=(m?[0-9p]+)")


def _tok_to_float(tok):
    sign = -1.0 if tok.startswith("m") else 1.0
    return sign * float(tok.lstrip("m").replace("p", "."))


def get_volume_params(path):
    path = Path(path)
    try:
        meta = pq.read_schema(path).metadata or {}
        if FOOTER_KEY in meta:
            p = json.loads(meta[FOOTER_KEY])
            if all(p.get(k) is not None
                   for k in ("depth_min", "depth_max", "disk_radius")):
                return {"depth_min": float(p["depth_min"]),
                        "depth_max": float(p["depth_max"]),
                        "disk_radius": float(p["disk_radius"])}
    except Exception:
        pass
    m = _PATH_RE.search(str(path))
    if m is None:
        raise ValueError(f"no volume params for {path}")
    d1, d2 = _tok_to_float(m.group(1)), _tok_to_float(m.group(2))
    return {"depth_min": min(d1, d2), "depth_max": max(d1, d2),
            "disk_radius": _tok_to_float(m.group(3))}


def rock_volume_m3(vp):
    return (vp["depth_max"] - vp["depth_min"]) * np.pi * vp["disk_radius"] ** 2


def get_df_of_acceptances(files, energycut=(10, 100, 1000), verbose=True):
    if isinstance(files, (str, Path)):
        files = [files]
    energycut = list(energycut)
    records = []
    for f in files:
        vp = get_volume_params(f)
        vol = rock_volume_m3(vp)
        df = pd.read_parquet(f)
        if verbose:
            print(f"{f}\n  volume={vol/1e9:.4g} km^3")
        for mass, sub in df.groupby("M_DM"):
            org = sub["total_org_nevents"].iloc[0]
            rec = {"volume_m3": vol, "M_DM": mass, "n_generated": org}
            hit_id = sub["hit_inner_detector"] == True   # noqa: E712
            for ecut in energycut:
                passed = sub["efinal_mu1"] > ecut
                rec[f"count_ecut{ecut}"] = int(passed.sum())
                rec[f"count_hit_id_ecut{ecut}"] = int((passed & hit_id).sum())
            records.append(rec)
    per_file = pd.DataFrame(records)
    keys = ["volume_m3", "M_DM"]
    sums = ["n_generated"] + [c for c in per_file.columns
                              if c.startswith("count_")]
    dfacc = per_file.groupby(keys, as_index=False)[sums].sum()
    for ecut in energycut:
        dfacc[f"frac_ecut{ecut}"] = dfacc[f"count_ecut{ecut}"] / dfacc["n_generated"]
        dfacc[f"frac_hit_id_ecut{ecut}"] = (
            dfacc[f"count_hit_id_ecut{ecut}"] / dfacc["n_generated"])
    return dfacc.sort_values(["volume_m3", "M_DM"]).reset_index(drop=True)


def plot_model(dfacc, model, plotdir):
    # data-driven y-range with ~2 decades of headroom above the highest
    # point so the (log-y) legend sits in clear space instead of overlapping
    frac_cols = [c for c in dfacc.columns if c.startswith("frac_")]
    vals = dfacc[frac_cols].to_numpy().ravel()
    vals = vals[vals > 0]
    ylim = (vals.min() / 3.0, vals.max() * 1e3)
    fig = plt.figure(figsize=(6.5, 6.5))
    ax = plt.gca()
    styles = [("frac_ecut10", "b", "o", r"$E_{\mu}$ > 10 GeV"),
              ("frac_ecut100", "r", "^", r"$E_{\mu}$ > 100 GeV"),
              ("frac_ecut1000", "g", "v", r"$E_{\mu}$ > 1000 GeV"),
              ("frac_hit_id_ecut10", "b", "s", r"$E_{\mu}$ > 10 GeV (ID)"),
              ("frac_hit_id_ecut100", "r", "P", r"$E_{\mu}$ > 100 GeV (ID)"),
              ("frac_hit_id_ecut1000", "g", ">", r"$E_{\mu}$ > 1000 GeV (ID)")]
    for col, c, mk, lab in styles:
        d = dfacc[dfacc[col] > 0]
        ax.scatter(d["M_DM"], d[col], color=c, marker=mk, s=55, label=lab)
    ax.set_title(f"DM model: {model}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"DM mass $M_{DM}$ [GeV/$c^2$]")
    ax.set_ylabel("Geometric acceptance [fraction]")
    ax.set_ylim(*ylim)
    ax.set_xlim(5e0, 2e5)
    ax.legend(loc="upper center", ncol=2, columnspacing=1.0,
              handletextpad=0.3, framealpha=0.9, fontsize=11)
    plt.tight_layout()
    out = Path(plotdir) / f"acc_dm_model_{model}.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"  -> {out}")

    return fig,ax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data")
    ap.add_argument("--plotdir", default="figures_pdf")
    args = ap.parse_args()

    for model in ("core", "floating"):
        print(f"[{model}] acceptance")
        sm = eio.summary(args.data_root, stage="combined", dm_model=model,
                         full_paths=True)
        if len(sm) == 0:
            print(f"  [skip] no {model} datasets in catalog")
            continue
        files = list(sm["path"])
        dfacc = get_df_of_acceptances(files, energycut=[10, 100, 1000])
        plot_model(dfacc, model, args.plotdir)
        dfacc.to_parquet(Path(args.plotdir) / f"acc_dm_model_{model}.parquet")


if __name__ == "__main__":
    main()
