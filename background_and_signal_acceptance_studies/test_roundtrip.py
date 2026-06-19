"""End-to-end test of earthshine_io + migration on synthetic data."""
import sys
sys.path.insert(0, ".")

import numpy as np
import pandas as pd

import earthshine_io as eio
from migrate_existing_files import parse_legacy_name, LEGACY_RE

# ----------------------------------------------------------------------
# 1. legacy filename parser against real-world patterns
# ----------------------------------------------------------------------
names = [
    # the exact example from the conversation
    "generated_data_depth_-8.0--4000.0_diskR_4000.0_mDM_10000.0-90000.0_mA_0.22_dmModel_momentum_constrained_HIT_DETECTOR_ave_eloss__COMBINED.parquet",
    # a shard with trial number (double underscore from the tuple bug era)
    "generated_data_depth_-8.0--4000.0_diskR_4000.0_mDM_10000.0-90000.0_mA_0.22_dmModel_core_HIT_DETECTOR_ave_eloss__000042.parquet",
    # GEN_EVENTS.sh style: depth -8 2000 -> '-8.0-2000.0'
    "generated_data_depth_-8.0-2000.0_diskR_2000.0_mDM_10000.0-90000.0_mA_0.22_dmModel_floating_HIT_DETECTOR_ave_eloss__COMBINED.parquet",
    # diagnostics-era default tag style: dm_model, no flags, int-ish depth
    "generated_data_depth_-8_diskR_500_mDM_200-10000_mA_0.22_dm_model_floating.parquet",
    # single mass, single depth, no eloss
    "generated_data_depth_-10.0_diskR_40.0_mDM_20000.0_mA_0.22_dmModel_core_HIT_DETECTOR_COMBINED.parquet",
]

print("=" * 70)
print("PARSER TESTS")
print("=" * 70)
ok = True
for n in names:
    p = parse_legacy_name(n)
    if p is None:
        print(f"FAIL: {n}")
        ok = False
        continue
    print(f"\n{n[:74]}...")
    print(f"  model={p['dm_model']!r} stage={p['stage']!r} trial={p['trial']}")
    print(f"  depth=[{p['depth_min']}, {p['depth_max']}]  "
          f"diskR={p['disk_radius']}")
    print(f"  mDM=[{p['mDM_min']}, {p['mDM_max']}]  "
          f"mA=[{p['mA_min']}, {p['mA_max']}]")
    print(f"  selection={p['selection']!r} eloss={p['eloss']!r}")
    print(f"  -> new path: {eio.rel_path(p)}")

# spot checks
p = parse_legacy_name(names[0])
assert p["dm_model"] == "momentum_constrained", p["dm_model"]
assert p["depth_min"] == -4000.0 and p["depth_max"] == -8.0
assert p["mDM_min"] == 10000.0 and p["mDM_max"] == 90000.0
assert p["mA_min"] == 0.22 and p["stage"] == "combined"
assert p["selection"] == "hit_detector" and p["eloss"] == "ave"

p = parse_legacy_name(names[1])
assert p["stage"] == "shard" and p["trial"] == 42 and p["dm_model"] == "core"

p = parse_legacy_name(names[3])
assert p["dm_model"] == "floating" and p["selection"] == "all"
assert p["stage"] == "combined"  # no trial number -> standalone file
assert p["eloss"] == "none" and p["mDM_min"] == 200.0

p = parse_legacy_name(names[4])
assert p["mDM_min"] == p["mDM_max"] == 20000.0 and p["eloss"] == "none"

print("\nAll parser assertions passed." if ok else "PARSER FAILURES ABOVE")

# ----------------------------------------------------------------------
# 2. synthetic legacy files -> migrate -> catalog -> load -> figure
# ----------------------------------------------------------------------
print("\n" + "=" * 70)
print("ROUND-TRIP TEST")
print("=" * 70)

import shutil, os
from pathlib import Path

shutil.rmtree("test_legacy", ignore_errors=True)
shutil.rmtree("test_data", ignore_errors=True)
os.makedirs("test_legacy")

rng = np.random.default_rng(0)


def fake_df(n, masses):
    m = rng.choice(masses, n)
    df = pd.DataFrame({
        "M_DM": m, "M_A": 0.22,
        "px1": rng.normal(0, 10, n), "py1": rng.normal(0, 10, n),
        "pz1": rng.normal(0, 10, n),
        "e1": m / 2 * (1 + 0.001 * rng.normal(size=n)),
        "x0": rng.uniform(-100, 100, n), "y0": rng.uniform(-2000, -8, n),
        "z0": rng.uniform(-100, 100, n),
        "ip_x0": rng.normal(0, 3, n), "ip_y0": rng.uniform(-5, 5, n),
        "ip_z0": rng.normal(0, 6, n),
        "ip_inner_x0": rng.normal(0, 0.7, n),
        "ip_inner_y0": rng.normal(0, 0.7, n),
        "ip_inner_z0": rng.normal(0, 1.5, n),
        "theta1": rng.uniform(0, np.pi, n),
        "eta1": rng.exponential(0.8, n),
        "phi1": rng.uniform(0, np.pi, n),
        "hit_inner_detector": rng.random(n) < 0.3,
        "efinal_mu1": m / 2 * rng.random(n),
        "pt1_detector_acceptance": rng.exponential(50, n),
        "pt1_detector_acceptance_eloss": rng.exponential(45, n),
        "org_nevents": 1_000_000,
        "total_org_nevents": 2_000_000_000,
    })
    return df


legacy_names = [
    ("generated_data_depth_-8.0--4000.0_diskR_4000.0_mDM_10000.0-90000.0"
     "_mA_0.22_dmModel_momentum_constrained_HIT_DETECTOR_ave_eloss"
     "__COMBINED.parquet", [10000., 50000., 90000.]),
    ("generated_data_depth_-8.0--4000.0_diskR_4000.0_mDM_10000.0-90000.0"
     "_mA_0.22_dmModel_core_HIT_DETECTOR_ave_eloss__COMBINED.parquet",
     [10000., 50000., 90000.]),
    ("generated_data_depth_-8.0-2000.0_diskR_2000.0_mDM_10000.0-90000.0"
     "_mA_0.22_dmModel_floating_HIT_DETECTOR_ave_eloss__COMBINED.parquet",
     [10000., 90000.]),
]
for name, masses in legacy_names:
    fake_df(5000, masses).to_parquet(f"test_legacy/{name}")

from migrate_existing_files import migrate
migrate("test_legacy", "test_data", dry_run=False)

# catalog
cat = eio.read_catalog("test_data")
print("\ncatalog columns:", sorted(cat.columns)[:12], "...")
print(cat[["dm_model", "depth_min", "depth_max", "disk_radius",
           "n_generated", "n_rows"]].to_string())

# load exactly one
df, params = eio.load(cat, dm_model="core")
assert params["n_generated"] == 2_000_000_000
print(f"\nload(dm_model='core'): {len(df)} rows, "
      f"n_generated={params['n_generated']}")

# ambiguous load should raise helpfully
try:
    eio.load(cat, mA_min=0.22)
    raise SystemExit("should have raised")
except ValueError as e:
    print(f"\nambiguous load correctly raised:\n{str(e)[:200]}...")

# float-tolerant filtering
assert len(eio.select(cat, mA_min=0.22)) == 3

# selection + figure naming
mask, sel = eio.apply_selection(df, mass=50000.0, inner_detector=True)
print(f"\nselection: {mask.sum()} / {len(df)} events pass")
fname = eio.figure_name("e_and_pt", params, sel)
print(f"figure name: {fname}")
assert "." not in fname.replace(".png", ""), "dots in figure name!"
assert " " not in fname

# save_figure writes png + sidecar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig = plt.figure()
plt.hist(df[mask]["efinal_mu1"], bins=30)
out = eio.save_figure(fig, "e_and_pt", params, sel, plotdir="test_plots")
import json
sidecar = json.load(open(str(out).replace(".png", ".json")))
assert sidecar["params"]["dm_model"] == "core"
assert sidecar["selection"]["mass"] == 50000.0
print(f"figure + sidecar written: {out}")

# diagnostics_v2 smoke tests
import diagnostics_v2 as diag
_, figs = diag.e_and_pt_plots(df, params, masses=[50000.0],
                              inner_detector=True, plotdir="test_plots")
print(f"\ndiagnostics_v2.e_and_pt_plots produced {len(figs)} figure(s)")

_, figs = diag.kinematic_diagnostic(df, params, masses=[50000.0],
                                    plotdir="test_plots")
assert len(figs) == 2
print(f"diagnostics_v2.kinematic_diagnostic produced {len(figs)} figure(s)")

_, figs = diag.origin_plots(df, params, masses=[10000.0, 50000.0],
                            plotdir="test_plots", out_root="test_data")
assert len(figs) == 4   # 2 figures x 2 masses
print(f"diagnostics_v2.origin_plots produced {len(figs)} figure(s)")

_, figs = diag.detector_hits_plots(df, params, masses=[50000.0],
                                   plotdir="test_plots")
assert len(figs) == 1
print(f"diagnostics_v2.detector_hits_plots produced {len(figs)} figure(s)")

# compare two datasets (core vs momentum_constrained from this test tree)
df_b, params_b = eio.load(cat, dm_model="momentum_constrained")
_, figs = diag.compare_pt_plots([(df, params), (df_b, params_b)],
                                masses=[50000.0], plotdir="test_plots")
assert len(figs) == 1
import glob as _glob
cmp_pngs = _glob.glob("test_plots/compare_pt_*core-vs-momentum_constrained*.png")
assert cmp_pngs, "comparison figure not named with both labels"
cmp_sidecar = json.load(open(cmp_pngs[0].replace(".png", ".json")))
assert len(cmp_sidecar["params"]["compare"]) == 2
print(f"diagnostics_v2.compare_pt_plots produced {len(figs)} figure(s) "
      f"with multi-dataset sidecar")

# crash-fix regressions: masses=None and max_xval both work now
_, figs = diag.compare_pt_plots([(df, params), (df_b, params_b)],
                                max_xval=30000.0, plotdir="test_plots")
assert len(figs) >= 1
print("compare_pt_plots: masses=None and max_xval= both OK")

# the derived depth histograms must now be discoverable through the catalog
cat2 = eio.build_catalog("test_data", verbose=False)
derived = cat2[cat2["stage"] == "derived"]
assert len(derived) == 1, derived
assert derived.iloc[0]["kind"] == "depth_histogram"
# mass-subset identity: histogram entry reflects requested masses, not file range
assert derived.iloc[0]["mDM_min"] == 10000.0
assert derived.iloc[0]["mDM_max"] == 50000.0
dh = pd.read_parquet(derived.iloc[0]["path"])
assert "10000 count" in dh.columns and "50000 count" in dh.columns
print(f"derived depth histograms found in catalog: {derived.iloc[0]['path']}")

# different selection cuts must produce a DIFFERENT derived file, not overwrite
diag.save_depth_histograms(df, params, masses=[10000.0, 50000.0],
                           out_root="test_data", inner_detector=True)
cat3 = eio.build_catalog("test_data", verbose=False)
assert (cat3["stage"] == "derived").sum() == 2, "selection cuts collided!"
print("distinct cuts -> distinct derived files: OK")

print("\n" + "=" * 70)
print("ALL TESTS PASSED")
print("=" * 70)
