"""
acceptance_v2 -- multi-file, multi-run muon detection efficiencies for EarthShine.

Key change from the original get_df_of_acceptances: counts (numerators) and
generated totals (denominators) are SUMMED across all input files sharing the
same (rock volume, M_DM) before any fraction is computed. This makes it safe
to feed in multiple runs at the same physics point (same masses, different
run_ids) as well as files covering different masses.

Files generated over different rock volumes are kept as separate rows, keyed
by volume_m3 (plus depth_min/depth_max/disk_radius for readability), so you
can cross-check that efficiencies agree between volumes.

Volume parameters are read from the parquet footer metadata (key b"earthshine")
when present, with a fallback that parses the hive-style path, e.g.
    .../depth=m100-m8_diskR=40/combined-add8072c-....parquet
where 'm' denotes a minus sign and 'p' a decimal point.
"""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

FOOTER_KEY = b"earthshine"

# e.g. depth=m100-m8_diskR=40  or  depth=m4000-m8_diskR=4000p5
_PATH_RE = re.compile(r"depth=(m?[0-9p]+)-(m?[0-9p]+)_diskR=(m?[0-9p]+)")


def _tok_to_float(tok):
    """'m100' -> -100.0, '4000p5' -> 4000.5"""
    sign = -1.0 if tok.startswith("m") else 1.0
    return sign * float(tok.lstrip("m").replace("p", "."))


def get_volume_params(path):
    """
    Return dict(depth_min, depth_max, disk_radius) for a file.

    Prefers the embedded footer metadata; falls back to parsing the path.
    Convention: depth_max is the numerically larger (shallower) boundary.
    """
    path = Path(path)

    # 1) Footer metadata (source of truth)
    try:
        meta = pq.read_schema(path).metadata or {}
        if FOOTER_KEY in meta:
            p = json.loads(meta[FOOTER_KEY])
            if all(p.get(k) is not None for k in ("depth_min", "depth_max", "disk_radius")):
                return {
                    "depth_min": float(p["depth_min"]),
                    "depth_max": float(p["depth_max"]),
                    "disk_radius": float(p["disk_radius"]),
                }
    except Exception:
        pass  # fall through to path parsing

    # 2) Path fallback
    m = _PATH_RE.search(str(path))
    if m is None:
        raise ValueError(
            f"Could not determine volume parameters for {path}: "
            "no 'earthshine' footer metadata and no 'depth=..._diskR=...' in path."
        )
    d1 = _tok_to_float(m.group(1))
    d2 = _tok_to_float(m.group(2))
    return {
        "depth_min": min(d1, d2),
        "depth_max": max(d1, d2),
        "disk_radius": _tok_to_float(m.group(3)),
    }


def rock_volume_m3(vp):
    """Cylinder volume in m^3 from depth_min/depth_max/disk_radius (all in m)."""
    span = vp["depth_max"] - vp["depth_min"]
    return span * np.pi * vp["disk_radius"] ** 2


def get_df_of_acceptances(files, energycut=(10, 100, 1000), verbose=True):
    """
    Compute per-mass detection efficiencies from one or many parquet files.

    Parameters
    ----------
    files : str, Path, or list of str/Path
        Input parquet file(s). May freely mix files with the same masses
        (multiple runs at one physics point) and different masses. Files
        with different rock volumes are aggregated separately.
    energycut : sequence of floats
        Energy cuts (GeV) applied to efinal_mu1.
    verbose : bool
        Print per-file summary lines.

    Returns
    -------
    DataFrame with one row per (volume, M_DM):
        depth_min, depth_max, disk_radius, volume_m3, M_DM, n_generated,
        count_ecut{E}, frac_ecut{E},
        count_hit_id_ecut{E}, frac_hit_id_ecut{E}   for each E in energycut.

    Notes
    -----
    - Fractions are computed only AFTER summing counts and generated totals
      across files, which is the statistically correct combination for
      multiple runs at the same mass point.
    - n_generated for a mass is the sum of total_org_nevents over the files
      in which that mass appears. A mass can only be discovered from the
      rows present in a file; a run in which a mass had zero surviving
      events would be invisible here (astronomically unlikely for your
      statistics, but worth knowing).
    """
    if isinstance(files, (str, Path)):
        files = [files]
    energycut = list(energycut)

    records = []
    for f in files:
        vp = get_volume_params(f)
        vol = rock_volume_m3(vp)

        df = pd.read_parquet(f)

        if verbose:
            print(f"{f}")
            print(
                f"  depth=[{vp['depth_min']}, {vp['depth_max']}] m, "
                f"diskR={vp['disk_radius']} m, volume={vol:.4g} m^3 "
                f"({vol / 1e9:.4g} km^3)"
            )

        for mass, sub in df.groupby("M_DM"):
            org_nevents = sub["total_org_nevents"].iloc[0]
            rec = {
                "depth_min": vp["depth_min"],
                "depth_max": vp["depth_max"],
                "disk_radius": vp["disk_radius"],
                "volume_m3": vol,
                "M_DM": mass,
                "n_generated": org_nevents,
            }
            hit_id = sub["hit_inner_detector"] == True  # noqa: E712 (handles NaN/object)
            for ecut in energycut:
                passed = sub["efinal_mu1"] > ecut
                rec[f"count_ecut{ecut}"] = int(passed.sum())
                rec[f"count_hit_id_ecut{ecut}"] = int((passed & hit_id).sum())
            records.append(rec)

            if verbose:
                print(f"    M_DM={mass}: n_generated={org_nevents}")

    per_file = pd.DataFrame(records)

    group_keys = ["depth_min", "depth_max", "disk_radius", "volume_m3", "M_DM"]
    sum_cols = ["n_generated"] + [c for c in per_file.columns if c.startswith("count_")]
    dfacc = per_file.groupby(group_keys, as_index=False)[sum_cols].sum()

    # Fractions computed only after summation
    for ecut in energycut:
        dfacc[f"frac_ecut{ecut}"] = dfacc[f"count_ecut{ecut}"] / dfacc["n_generated"]
        dfacc[f"frac_hit_id_ecut{ecut}"] = (
            dfacc[f"count_hit_id_ecut{ecut}"] / dfacc["n_generated"]
        )

    return dfacc.sort_values(["volume_m3", "M_DM"]).reset_index(drop=True)


def compare_across_volumes(dfacc, frac_col):
    """
    Cross-check helper: pivot one fraction column so each rock volume is a
    column, indexed by M_DM. Consistent physics should give consistent
    efficiencies across volumes (within Poisson errors).
    """
    return dfacc.pivot_table(index="M_DM", columns="volume_m3", values=frac_col)
