"""
plot_category_A_figures.py
==========================
Driver for the Category-A (signal kinematics / geometry) paper figures,
per PLOTTING_PLAN.md.  Produces, at production statistics, as vector PDFs:

    ePtAngCore / ePtAngFloating          kinematic_diagnostic @ 10 TeV
    openAngleSepCore / openAngleSepFloating  opening_angle_plots @ 10 TeV
    density_profiles                     origin_plots @ 1 / 10 / 100 TeV (floating)
    emu_at_det                           e_and_pt_plots (energy at detector)

Reads the partitioned tree written by generate_signal_events_CL_v2.py
(--data-root data) via earthshine_io.  Run after the SLURM generation jobs
finish:

    venv/bin/python plot_category_A_figures.py --data-root data --plotdir figures_pdf

Missing (dm_model, mass) points are skipped with a warning instead of
crashing, so a partial data tree still yields whatever is ready.
"""

import argparse

import matplotlib
matplotlib.use("Agg")

import earthshine_io as eio
import diagnostics_v2 as diag


def find_dataset(cat, dm_model, mass, stage="combined"):
    """Return (df, params) for the unique combined dataset of this model
    whose [mDM_min, mDM_max] contains `mass`.  Returns None if absent."""
    rows = cat[
        (cat["dm_model"] == dm_model)
        & (cat["stage"] == stage)
        & (cat["mDM_min"] <= mass)
        & (cat["mDM_max"] >= mass)
    ]
    if len(rows) == 0:
        print(f"  [skip] no {dm_model} dataset contains mDM={mass}")
        return None
    if len(rows) > 1:
        # prefer the tightest mass range (most specific file), then the most
        # statistics (largest n_generated) -- so when several runs at the same
        # physics point coexist (e.g. a 5-trial and a 500-trial rerun), the
        # higher-stats file wins
        ng = rows["n_generated"] if "n_generated" in rows.columns else 0
        rows = rows.assign(_span=rows["mDM_max"] - rows["mDM_min"], _ng=ng) \
                   .sort_values(["_span", "_ng"], ascending=[True, False])
        print(f"  [note] {len(rows)} {dm_model} datasets contain mDM={mass}; "
              f"using range "
              f"[{rows.iloc[0]['mDM_min']}, {rows.iloc[0]['mDM_max']}] "
              f"n_generated={rows.iloc[0].get('n_generated')}")
    path = rows.iloc[0]["path"]
    import pandas as pd
    return pd.read_parquet(path), eio.read_params(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data")
    ap.add_argument("--plotdir", default="figures_pdf")
    ap.add_argument("--ptang-mass", type=float, default=10000.0,
                    help="mass (GeV) for ePtAng + openAngleSep (paper: 10 TeV)")
    ap.add_argument("--density-masses", type=float, nargs="+",
                    default=[1000.0, 10000.0, 100000.0],
                    help="masses (GeV) for density_profiles (paper: 1/10/100 TeV)")
    ap.add_argument("--emu-mass", type=float, default=10000.0,
                    help="mass (GeV) for the emu_at_det energy panel")
    args = ap.parse_args()

    eio.DEFAULT_FIG_EXT = "pdf"          # vector output for the paper
    cat = eio.read_catalog(args.data_root)
    print(f"catalog: {len(cat)} datasets under {args.data_root!r}\n")

    # ---- ePtAng + openAngleSep, core & floating, @ ptang-mass -------------
    for model in ("core", "floating"):
        print(f"[{model}] ePtAng / openAngleSep @ mDM={args.ptang_mass}")
        ds = find_dataset(cat, model, args.ptang_mass)
        if ds is None:
            continue
        df, params = ds
        diag.kinematic_diagnostic(df, params, masses=[args.ptang_mass],
                                  plotdir=args.plotdir)
        diag.opening_angle_plots(df, params, masses=[args.ptang_mass],
                                 plotdir=args.plotdir)

    # ---- density_profiles (floating), 1/10/100 TeV -----------------------
    print("\n[floating] density_profiles")
    for mass in args.density_masses:
        ds = find_dataset(cat, "floating", mass)
        if ds is None:
            continue
        df, params = ds
        diag.origin_plots(df, params, masses=[mass], plotdir=args.plotdir,
                          save_histograms=False)

    # ---- emu_at_det (floating) -------------------------------------------
    print("\n[floating] emu_at_det")
    ds = find_dataset(cat, "floating", args.emu_mass)
    if ds is not None:
        df, params = ds
        diag.e_and_pt_plots(df, params, masses=[args.emu_mass],
                            plotdir=args.plotdir)

    print(f"\nDone.  PDFs in {args.plotdir}/")


if __name__ == "__main__":
    main()
