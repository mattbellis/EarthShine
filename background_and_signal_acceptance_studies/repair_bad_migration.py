"""
repair_bad_migration.py
=======================
Removes files that were mis-migrated by the original (too-permissive)
legacy filename parser -- identifiable by dm_model values containing
flag strings (e.g. 'momentum_constrained_HIT_DETECTOR_ave_eloss_997')
-- then rebuilds the catalog.

After running this, re-run the migration with the FIXED
migrate_existing_files.py; it now skips destinations that already exist,
so only the previously-corrupted files get re-migrated (this time with
correct params, or into the unparseable list for manual review).

Usage:
    python repair_bad_migration.py --data-root data --dry-run
    python repair_bad_migration.py --data-root data
    python migrate_existing_files.py --input-dir OUTPUT_FILES --data-root data
"""

import argparse
import re
from pathlib import Path

import earthshine_io as eio

BAD_TOKENS = ("HIT_DETECTOR", "ave_eloss", "COMBINED")


def is_corrupted(params):
    model = str(params.get("dm_model", ""))
    if any(tok in model for tok in BAD_TOKENS):
        return True
    if re.search(r"_\d+$", model):
        return True
    return False


def repair(data_root, dry_run=True):
    data_root = Path(data_root)
    bad = []
    for p in sorted(data_root.rglob("*.parquet")):
        if p.name == "catalog.parquet":
            continue
        params = eio.read_params(p)
        if params and is_corrupted(params):
            bad.append((p, params))

    if not bad:
        print("No corrupted entries found -- nothing to do.")
        return

    print(f"{len(bad)} corrupted file(s) found:\n")
    for p, params in bad:
        legacy = params.get("legacy_filename", "<unknown>")
        print(f"  {p}")
        print(f"      dm_model = {params['dm_model']!r}")
        print(f"      from legacy file: {legacy}")

    if dry_run:
        print("\n--dry-run: nothing deleted.")
        print("Rerun without --dry-run, then re-run migrate_existing_files.py")
        return

    for p, _ in bad:
        p.unlink()
        print(f"deleted {p}")
        # remove now-empty partition dirs up to data_root
        d = p.parent
        while d != data_root and not any(d.iterdir()):
            d.rmdir()
            d = d.parent

    eio.build_catalog(data_root)
    print("\nDone. Now re-run the (fixed) migration to re-import the "
          "affected legacy files:\n"
          "  python migrate_existing_files.py --input-dir <legacy dir> "
          "--data-root", data_root)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default="data")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    repair(a.data_root, dry_run=a.dry_run)
