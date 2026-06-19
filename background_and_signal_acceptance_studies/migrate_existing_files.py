"""
migrate_existing_files.py
=========================
One-time migration of legacy EarthShine MC files into the new layout.

For each generated_data_*.parquet in --input-dir:
  1. parse the legacy filename into a params dict   (the LAST time these
     filenames are ever parsed)
  2. read the file; pull n_generated from org_nevents / total_org_nevents
  3. write a copy into the partitioned tree under --data-root with the
     params embedded in the parquet footer
  4. rebuild catalog.parquet

Originals are left in place unless --delete-originals is given.
Run with --dry-run first to see the planned mapping and any parse failures.

Usage:
    python migrate_existing_files.py --input-dir OUTPUT_FILES --data-root data --dry-run
    python migrate_existing_files.py --input-dir OUTPUT_FILES --data-root data
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

import earthshine_io as eio

# ----------------------------------------------------------------------------
# legacy filename parser
#
# Handles names like:
#   generated_data_depth_-8.0--4000.0_diskR_4000.0_mDM_10000.0-90000.0_mA_0.22
#       _dmModel_momentum_constrained_HIT_DETECTOR_ave_eloss__COMBINED.parquet
#   generated_data_depth_-8_diskR_2000_mDM_10000-90000_mA_0.22_dmModel_core
#       _HIT_DETECTOR_ave_eloss__000042.parquet                       (shards)
#   generated_data_depth_-8_diskR_500_mDM_200-10000_mA_0.22_dm_model_floating
#       .parquet                                       (older diagnostics era)
#
# Number grammar matches return_tag(): scalar or 'first-last', where each can
# be negative, so depth ranges look like '-8.0--4000.0'.
# ----------------------------------------------------------------------------
_NUM = r"-?\d+(?:\.\d+)?"

LEGACY_RE = re.compile(
    rf"^generated_data_"
    rf"depth_(?P<d1>{_NUM})(?:-(?P<d2>{_NUM}))?"
    rf"_diskR_(?P<r1>{_NUM})(?:-(?P<r2>{_NUM}))?"
    rf"_mDM_(?P<m1>{_NUM})(?:-(?P<m2>{_NUM}))?"
    rf"_mA_(?P<a1>{_NUM})(?:-(?P<a2>{_NUM}))?"
    rf"_(?:dmModel|dm_model)_(?P<model>.+?)"
    rf"(?P<hit>_HIT_DETECTOR)?"
    rf"(?P<eloss>_ave_eloss)?"
    rf"(?:_+(?P<trial>\d+))?"          # trial numbers: padded (000042) or not (997)
    rf"(?:_*(?P<combined>COMBINED))?"
    rf"\.parquet$"
)

# strings that must never appear inside a parsed model name; if they do, the
# lazy model group has swallowed flags it shouldn't have, and the file must
# go to the unparseable list (loud failure) rather than be silently mangled
_MODEL_FORBIDDEN = ("HIT_DETECTOR", "ave_eloss", "COMBINED")


def parse_legacy_name(name):
    """Parse a legacy filename -> params dict, or None if it doesn't match
    or if the parse looks corrupted (model containing flag strings or a
    trailing run number)."""
    m = LEGACY_RE.match(name)
    if m is None:
        return None
    g = m.groupdict()

    # validation net: a model name containing flag markers or ending in a
    # bare number means the catch-all group ate something it shouldn't have
    model = g["model"]
    if any(tok in model for tok in _MODEL_FORBIDDEN):
        return None
    if re.search(r"_\d+$", model):
        return None

    depth = [float(g["d1"])] + ([float(g["d2"])] if g["d2"] else [])
    diskR = float(g["r1"])  # ranges in diskR never occurred, but r2 tolerated
    mDM = [float(g["m1"])] + ([float(g["m2"])] if g["m2"] else [])
    mA = [float(g["a1"])] + ([float(g["a2"])] if g["a2"] else [])

    # shards always carry a _NNNNNN trial number; a file with neither a
    # trial number nor a COMBINED marker is a standalone single-run file
    # and belongs in the combined bucket
    trial = int(g["trial"]) if g["trial"] else None
    stage = "shard" if trial is not None and not g["combined"] else "combined"

    return eio.make_params(
        dm_model=g["model"],
        masses_dm=mDM,           # legacy names only kept first-last; the
        masses_a=mA,             # footer min/max is what matters downstream
        depth=depth,
        disk_radius=diskR,
        selection="hit_detector" if g["hit"] else "all",
        eloss="ave" if g["eloss"] else "none",
        stage=stage,
        trial=trial,
        extra={"legacy_filename": name},
    )


# ----------------------------------------------------------------------------
def migrate(input_dir, data_root, dry_run=True, delete_originals=False,
            overwrite=False):
    input_dir, data_root = Path(input_dir), Path(data_root)
    files = sorted(input_dir.glob("generated_data_*.parquet"))
    if not files:
        print(f"no generated_data_*.parquet found in {input_dir}")
        return

    plan, failures, existing = [], [], []
    for f in files:
        params = parse_legacy_name(f.name)
        if params is None:
            failures.append(f)
            continue
        dst = data_root / eio.rel_path(params)
        if dst.exists() and not overwrite:
            existing.append((f, dst))
            continue
        plan.append((f, dst, params))

    # ---- report ----
    print(f"{len(files)} files found: {len(plan)} to migrate, "
          f"{len(existing)} already migrated (skipped), "
          f"{len(failures)} NOT parseable\n")
    if existing and len(existing) <= 5:
        for src, dst in existing:
            print(f"  exists, skipping: {dst}")
    if failures:
        print("UNPARSEABLE (handle manually or extend LEGACY_RE):")
        for f in failures:
            print(f"  {f.name}")
        print()

    # collision check: distinct sources mapping to the same destination
    dests = {}
    for src, dst, _ in plan:
        dests.setdefault(str(dst), []).append(src.name)
    collisions = {d: s for d, s in dests.items() if len(s) > 1}
    if collisions:
        print("COLLISIONS -- multiple legacy files map to one destination.")
        print("(Identical params from separate runs; merge or tag manually.)")
        for d, srcs in collisions.items():
            print(f"  {d}")
            for s in srcs:
                print(f"      <- {s}")
        print("Aborting. Resolve collisions and rerun.")
        sys.exit(1)

    for src, dst, _ in plan[:8]:
        print(f"  {src.name}\n    -> {dst}")
    if len(plan) > 8:
        print(f"  ... and {len(plan) - 8} more")

    if dry_run:
        print("\n--dry-run: nothing written.")
        return

    # ---- execute ----
    for src, dst, params in plan:
        df = pd.read_parquet(src)
        # promote legacy bookkeeping columns into footer metadata
        if "total_org_nevents" in df.columns and len(df):
            params["n_generated"] = int(df["total_org_nevents"].iloc[0])
        elif "org_nevents" in df.columns and len(df):
            params["n_generated"] = int(df["org_nevents"].iloc[0])
        eio.write_parquet(df, dst, params)
        print(f"wrote {dst}  ({len(df)} rows)")
        if delete_originals:
            src.unlink()

    eio.build_catalog(data_root)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--data-root", default="data")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--delete-originals", action="store_true")
    ap.add_argument("--overwrite", action="store_true",
                    help="rewrite destinations that already exist")
    a = ap.parse_args()
    migrate(a.input_dir, a.data_root, dry_run=a.dry_run,
            delete_originals=a.delete_originals, overwrite=a.overwrite)
