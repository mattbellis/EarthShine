"""
earthshine_io.py
================
I/O, metadata, naming, and catalog utilities for EarthShine MC parquet files.

Design (three layers):
  1. Every parquet file carries its full generation parameters in its footer
     metadata under the key b"earthshine".  Files are standalone.
  2. Filenames are short and dumb (run/content-hash ids); directories give a
     human-browsable partition tree.  Directory names are REDUNDANT with the
     footer metadata -- the footer always wins.
  3. catalog.parquet is a GENERATED index (one row per file), rebuilt by
     scanning footers.  Never edited by hand, always regenerable.

Typical notebook use:

    import earthshine_io as eio
    cat = eio.read_catalog("data")
    df, params = eio.load(cat, dm_model="core", mA=0.22, mDM_min=10000)
    mask, sel = eio.apply_selection(df, inner_detector=True)
    ...
    eio.save_figure("e_and_pt", params, sel, plotdir="plots")
"""

from __future__ import annotations

import datetime
import getpass
import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

META_KEY = b"earthshine"
SCHEMA_VERSION = 1

# Parameters that define the *physics identity* of a dataset.  Used for the
# content hash (file naming) and the default plot tag.  Bookkeeping fields
# (timestamps, git hash, n_generated, ...) are stored too but excluded here
# so that regenerating with more statistics maps to the same identity.
IDENTITY_KEYS = (
    "dm_model",
    "mA_min", "mA_max",
    "mDM_min", "mDM_max",
    "depth_min", "depth_max",
    "disk_radius",
    "selection", "eloss", "stage",
)


# ----------------------------------------------------------------------------
# params construction / normalization
# ----------------------------------------------------------------------------
def _minmax(v):
    """Accept scalar, or list/tuple/array; return (min, max) as floats."""
    if v is None:
        return None, None
    if np.iterable(v) and not isinstance(v, str):
        v = [float(x) for x in v]
        return min(v), max(v)
    return float(v), float(v)


def make_params(
    *,
    dm_model,
    masses_dm,                 # scalar or list (GeV)
    masses_a,                  # scalar or list (GeV)
    depth,                     # scalar or [min, max] (m; negative = below ground)
    disk_radius,               # m (None -> store as None, see note)
    selection="hit_detector",  # "hit_detector" | "all"
    eloss="ave",               # "ave" | "geant_spline" | "none"
    stage="combined",          # "combined" | "shard"
    trial=None,                # shard index, if stage == "shard"
    detector_radius=None, detector_half_len=None,
    inner_detector_radius=None, inner_detector_half_len=None,
    n_generated=None,          # total thrown events this file represents
    n_events_per_trial=None, n_trials=None,
    seed=None,
    extra=None,                # dict of anything else worth keeping
):
    """Build the canonical metadata dict stored in every file footer."""
    dmin, dmax = _minmax(depth)
    mmin, mmax = _minmax(masses_dm)
    amin, amax = _minmax(masses_a)

    p = {
        "schema_version": SCHEMA_VERSION,
        "dm_model": str(dm_model),
        "mDM_min": mmin, "mDM_max": mmax,
        "mDM_values": ([float(x) for x in np.atleast_1d(masses_dm)]
                       if masses_dm is not None else None),
        "mA_min": amin, "mA_max": amax,
        "mA_values": ([float(x) for x in np.atleast_1d(masses_a)]
                      if masses_a is not None else None),
        "depth_min": dmin, "depth_max": dmax,
        "disk_radius": (float(disk_radius) if disk_radius is not None else None),
        "selection": selection,
        "eloss": eloss,
        "stage": stage,
        "trial": trial,
        "detector_radius": detector_radius,
        "detector_half_len": detector_half_len,
        "inner_detector_radius": inner_detector_radius,
        "inner_detector_half_len": inner_detector_half_len,
        "n_generated": n_generated,
        "n_events_per_trial": n_events_per_trial,
        "n_trials": n_trials,
        "seed": seed,
        "code_git_hash": git_hash(),
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "created_by": getpass.getuser(),
    }
    if extra:
        p.update(extra)
    return p


def git_hash(short=True):
    """Git hash of the repo we are running from, or None outside a repo."""
    try:
        cmd = ["git", "rev-parse"] + (["--short"] if short else []) + ["HEAD"]
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        return out.stdout.strip() or None
    except Exception:
        return None


def params_hash(params, n=8):
    """Deterministic short hash of the physics-identity subset of params."""
    ident = {k: params.get(k) for k in IDENTITY_KEYS}
    blob = json.dumps(ident, sort_keys=True, default=str).encode()
    return hashlib.sha1(blob).hexdigest()[:n]


# ----------------------------------------------------------------------------
# naming:  LaTeX/shell-safe tokens, tags, partitioned paths
# ----------------------------------------------------------------------------
def fmt_token(v):
    """Format a value for use in file/dir names: no dots, no leading '-'.
    0.22 -> 0p22 ; -8.0 -> m8 ; 10000.0 -> 10000 ; 'momentum_constrained' kept.
    Only [A-Za-z0-9_pm-] survive."""
    if v is None:
        return "none"
    if isinstance(v, bool):
        return "T" if v else "F"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        if float(v).is_integer():
            s = str(int(v))
        else:
            s = repr(float(v))
        s = s.replace("-", "m").replace(".", "p")
        return s
    s = str(v)
    return "".join(c if (c.isalnum() or c in "_-") else "_" for c in s)


def _range_token(lo, hi):
    if lo is None:
        return "none"
    if lo == hi:
        return fmt_token(lo)
    return f"{fmt_token(lo)}-{fmt_token(hi)}"


def tag(params, keys=None):
    """Human-readable, LaTeX-safe tag for plot filenames.
    Built ONLY from params (and optional selection params merged in by the
    caller) -- never assembled by ad hoc string surgery at call sites."""
    parts = [
        f"model-{fmt_token(params['dm_model'])}",
        f"mA-{_range_token(params['mA_min'], params['mA_max'])}",
        f"mDM-{_range_token(params['mDM_min'], params['mDM_max'])}",
        f"depth-{_range_token(params['depth_min'], params['depth_max'])}",
        f"diskR-{fmt_token(params['disk_radius'])}",
        f"eloss-{fmt_token(params.get('eloss'))}",
    ]
    if keys:
        for k in keys:
            parts.append(f"{k}-{fmt_token(params.get(k))}")
    return "_".join(parts)


def rel_path(params):
    """Partitioned relative path for a data file (filename = content hash).

    data/<stage>/dm_model=<m>/mA=<t>/mDM=<t>/depth=<t>_diskR=<t>/
        <stage>-<hash8>[-<trial>].parquet
    """
    d = Path(params["stage"])
    d = d / f"dm_model={fmt_token(params['dm_model'])}"
    d = d / f"mA={_range_token(params['mA_min'], params['mA_max'])}"
    d = d / f"mDM={_range_token(params['mDM_min'], params['mDM_max'])}"
    d = d / (f"depth={_range_token(params['depth_min'], params['depth_max'])}"
             f"_diskR={fmt_token(params['disk_radius'])}")
    h = params_hash(params)
    name = f"{params['stage']}-{h}"
    if params.get("run_id"):
        # distinct generation runs at the same physics point coexist as
        # separate files; load them together with load_many()
        name += f"-{params['run_id']}"
    if params.get("trial") is not None:
        name += f"-{int(params['trial']):06d}"
    return d / f"{name}.parquet"


# ----------------------------------------------------------------------------
# parquet read/write with embedded metadata
# ----------------------------------------------------------------------------
def write_parquet(df, path, params):
    """Write df to path with params embedded in the footer (key 'earthshine')."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(df)
    meta = dict(table.schema.metadata or {})
    meta[META_KEY] = json.dumps(params, default=str).encode()
    table = table.replace_schema_metadata(meta)
    pq.write_table(table, path)
    return path


def read_params(path):
    """Read embedded params from a parquet footer (no data pages touched).
    Returns None if the file has no earthshine metadata."""
    md = pq.read_schema(path).metadata or {}
    raw = md.get(META_KEY)
    return json.loads(raw) if raw else None


# ----------------------------------------------------------------------------
# catalog
# ----------------------------------------------------------------------------
# list-valued params are kept in footers but flattened out of the catalog table
_CATALOG_SKIP = ("mDM_values", "mA_values")


def _catalog_safe(v):
    """Footers may hold lists/dicts (merged seeds, merged_from provenance,
    selection_cuts...); catalog COLUMNS must be scalar or Arrow refuses to
    mix list and non-list values in one column.  JSON-stringify structures
    so they stay visible and grep-able in the catalog."""
    if isinstance(v, (list, tuple, dict)):
        return json.dumps(v, default=str)
    return v


def build_catalog(root, out=None, verbose=True):
    """Scan root for *.parquet, read footers, write/return catalog dataframe.

    Concurrency-tolerant: files that can't be read (e.g. mid-write by a
    parallel generation job) are skipped with a note, and the catalog is
    written atomically (temp file + rename) so a reader never sees a
    half-written catalog.  If parallel jobs race on the rebuild, the
    catalog may briefly miss the slower job's entry -- rerun
    `python earthshine_io.py <root> --rebuild` after all jobs finish.
    """
    root = Path(root)
    rows, skipped = [], []
    for p in sorted(root.rglob("*.parquet")):
        if p.name == "catalog.parquet":
            continue
        try:
            params = read_params(p)
        except Exception as e:
            skipped.append((p, f"unreadable ({type(e).__name__}): "
                               f"possibly being written by another job"))
            continue
        if params is None:
            skipped.append((p, "no earthshine metadata"))
            continue
        row = {k: _catalog_safe(v) for k, v in params.items()
               if k not in _CATALOG_SKIP}
        row["path"] = str(p)
        row["n_rows"] = pq.read_metadata(p).num_rows
        rows.append(row)
    cat = pd.DataFrame(rows)
    # safety net: if any column ended up with mixed types across files
    # (e.g. an old footer written before a schema refinement), coerce the
    # whole column to string rather than letting Arrow refuse the catalog
    for c in cat.columns[cat.dtypes == object]:
        nonnull = cat[c].dropna()
        if len({type(v) for v in nonnull}) > 1:
            cat[c] = cat[c].map(lambda v: v if pd.isna(v) else str(v))
            if verbose:
                print(f"note: catalog column {c!r} had mixed types across "
                      f"files; stored as strings")
    if out is None:
        out = root / "catalog.parquet"
    if len(cat):
        # atomic write: never leave a half-written catalog for a
        # concurrent reader
        tmp = Path(str(out) + f".tmp{os.getpid()}")
        cat.to_parquet(tmp)
        os.replace(tmp, out)
    if verbose:
        print(f"catalog: {len(cat)} files indexed -> {out}")
        if skipped:
            print(f"WARNING: {len(skipped)} parquet files skipped:")
            for s, reason in skipped[:10]:
                print(f"  {s}  [{reason}]")
            if len(skipped) > 10:
                print(f"  ... and {len(skipped) - 10} more")
    return cat


def read_catalog(root):
    """Read data/catalog.parquet (rebuild with build_catalog if missing)."""
    root = Path(root)
    f = root / "catalog.parquet"
    if not f.exists():
        return build_catalog(root)
    return pd.read_parquet(f)


def _match(col, val):
    """Equality that is float-tolerant (so mA=0.22 matches 0.22000000001)."""
    if pd.api.types.is_float_dtype(col) or isinstance(val, float):
        return np.isclose(col.astype(float), float(val), rtol=1e-9, atol=0)
    return col == val


def select(catalog, **filters):
    """Filter the catalog by params equality.  Returns matching rows."""
    sel = catalog
    for k, v in filters.items():
        if k not in sel.columns:
            raise KeyError(f"'{k}' is not a catalog column. "
                           f"Columns: {sorted(sel.columns)}")
        sel = sel[_match(sel[k], v)]
    return sel


##########################################################################

def _read_parquet_filtered(path, columns=None, masses=None):
    """pd.read_parquet with optional row filtering on M_DM pushed down to
    the parquet reader -- rows for other masses are never materialized,
    which is the difference between loading 60M rows and 7M.

    Requested columns absent from the file's schema (e.g. older files
    predating a new column) are skipped with a warning rather than
    failing the read."""
    kw = {}
    if columns is not None:
        have = set(pq.read_schema(path).names)
        missing = [c for c in columns if c not in have]
        if missing:
            print(f"note: {Path(path).name} lacks requested column(s) "
                  f"{missing}; loading the rest")
        kw["columns"] = [c for c in columns if c in have]
    if masses is not None:
        kw["filters"] = [("M_DM", "in", [float(m) for m in masses])]
    return pd.read_parquet(path, **kw)

##########################################################################

def load(catalog, columns=None, masses=None, **filters):
    """Load exactly one dataset matching filters.

    Returns (df, params) where params is the full footer dict of that file.
    Raises with a helpful listing if 0 or >1 files match.

    Memory control for large files:
      columns=[...]  read only these columns (parquet is columnar; unread
                     columns cost nothing).  See diagnostics_v2.columns_for()
                     for per-plot-function column sets.
      masses=[...]   read only rows with M_DM in this list, filtered inside
                     the parquet reader before rows are materialized.
    """
    sel = select(catalog, **filters)
    if len(sel) != 1:
        show_cols = [c for c in
                     ("path", "stage", "dm_model", "mA_min", "mDM_min",
                      "mDM_max", "depth_min", "depth_max", "disk_radius",
                      "eloss", "selection") if c in sel.columns]
        raise ValueError(
            f"{len(sel)} catalog entries match {filters}; need exactly 1.\n"
            f"{sel[show_cols].to_string() if len(sel) else '(no matches)'}"
        )
    path = sel.iloc[0]["path"]
    df = _read_parquet_filtered(path, columns=columns, masses=masses)
    return df, read_params(path)

##########################################################################
'''
def load(catalog, columns=None, **filters):
    """Load exactly one dataset matching filters.

    Returns (df, params) where params is the full footer dict of that file.
    Raises with a helpful listing if 0 or >1 files match.
    """
    sel = select(catalog, **filters)
    if len(sel) != 1:
        show_cols = [c for c in
                     ("path", "stage", "dm_model", "mA_min", "mDM_min",
                      "mDM_max", "depth_min", "depth_max", "disk_radius",
                      "eloss", "selection") if c in sel.columns]
        raise ValueError(
            f"{len(sel)} catalog entries match {filters}; need exactly 1.\n"
            f"{sel[show_cols].to_string() if len(sel) else '(no matches)'}"
        )
    path = sel.iloc[0]["path"]
    df = pd.read_parquet(path, columns=columns)
    return df, read_params(path)
'''
##########################################################################
##########################################################################

def load_many(catalog, columns=None, masses=None, **filters):
    """Load and concat all datasets matching filters.
    Returns (df, list_of_params).  Adds a 'source_path' column.

    columns= and masses= behave as in load() (read-time pruning).

    For plotting/tagging a merged multi-run dataset, collapse the params
    list with merge_params():

        df, plist = load_many(cat, dm_model="core", stage="combined")
        params = merge_params(plist)
    """
    sel = select(catalog, **filters)
    if len(sel) == 0:
        raise ValueError(f"no catalog entries match {filters}")
    dfs, plist = [], []
    for path in sel["path"]:
        d = _read_parquet_filtered(path, columns=columns, masses=masses)
        d["source_path"] = path
        dfs.append(d)
        plist.append(read_params(path))
    return pd.concat(dfs, ignore_index=True), plist

##########################################################################
##########################################################################

'''
def load_many(catalog, columns=None, **filters):
    """Load and concat all datasets matching filters.
    Returns (df, list_of_params).  Adds a 'source_path' column.

    For plotting/tagging a merged multi-run dataset, collapse the params
    list with merge_params():

        df, plist = load_many(cat, dm_model="core", stage="combined")
        params = merge_params(plist)
    """
    sel = select(catalog, **filters)
    if len(sel) == 0:
        raise ValueError(f"no catalog entries match {filters}")
    dfs, plist = [], []
    for path in sel["path"]:
        d = pd.read_parquet(path, columns=columns)
        d["source_path"] = path
        dfs.append(d)
        plist.append(read_params(path))
    return pd.concat(dfs, ignore_index=True), plist
'''


def merge_params(plist):
    """Collapse a list of params dicts (e.g. from load_many) into one dict
    describing the combined dataset.

    All physics-identity fields (IDENTITY_KEYS) must agree across the list
    -- merging different physics points is almost certainly an analysis
    error, so it raises.  Bookkeeping fields are combined sensibly:
    n_generated and n_trials are SUMMED (this is the exposure denominator
    for acceptances), run_id becomes 'merged-<N>runs', and the individual
    runs' provenance is preserved under 'merged_from'.
    """
    if isinstance(plist, dict):           # already a single params dict
        return plist
    plist = list(plist)
    if len(plist) == 0:
        raise ValueError("merge_params: empty params list")
    if len(plist) == 1:
        return plist[0]

    first = plist[0]
    for k in IDENTITY_KEYS:
        vals = [p.get(k) for p in plist]
        if any(v != vals[0] for v in vals[1:]):
            raise ValueError(
                f"merge_params: datasets disagree on identity field "
                f"{k!r}: {vals}.  Refusing to merge different physics "
                f"points; filter more tightly or merge intentionally by "
                f"hand.")

    def _sum(key):
        vals = [p.get(key) for p in plist]
        return sum(v for v in vals if v is not None) if any(
            v is not None for v in vals) else None

    merged = dict(first)
    merged["n_generated"] = _sum("n_generated")
    merged["n_trials"] = _sum("n_trials")
    merged["run_id"] = f"merged-{len(plist)}runs"
    # do NOT overload the scalar 'seed' field with a list -- every column in
    # the catalog must hold one type, and other rows have integer seeds
    merged["seed"] = None
    merged["seeds"] = [p.get("seed") for p in plist]
    merged["created"] = max(str(p.get("created")) for p in plist)
    merged["merged_from"] = [
        {"run_id": p.get("run_id"), "n_generated": p.get("n_generated"),
         "created": p.get("created")} for p in plist
    ]
    return merged


def as_params(params):
    """Accept a params dict, or a list of them (as load_many returns),
    and hand back a single dict -- merging via merge_params if needed."""
    if isinstance(params, (list, tuple)):
        return merge_params(params)
    return params


def restamp(path, data_root=None, **updates):
    """Correct a file's footer metadata: read it, apply the keyword
    updates, rewrite the file, and (if data_root is given) move it to the
    partition path implied by the corrected params, pruning empty dirs.

        eio.restamp("data/combined/.../combined-xxxx.parquet",
                    data_root="data", disk_radius=6000.0,
                    disk_radius_auto=True)

    Returns the (possibly new) path.  Rebuild the catalog afterwards.
    """
    path = Path(path)
    params = read_params(path)
    if params is None:
        raise ValueError(f"{path} has no earthshine metadata to update")
    params.update(updates)
    df = pd.read_parquet(path)

    if data_root is not None:
        new_path = Path(data_root) / rel_path(params)
    else:
        new_path = path

    write_parquet(df, new_path, params)
    if new_path != path:
        path.unlink()
        d = path.parent
        root = Path(data_root)
        while d != root and d.exists() and not any(d.iterdir()):
            d.rmdir()
            d = d.parent
        print(f"restamped and moved: {path} -> {new_path}")
    else:
        print(f"restamped in place: {path}")
    return new_path


# ----------------------------------------------------------------------------
# event selection (one place; the mask and its name cannot diverge)
# ----------------------------------------------------------------------------
def apply_selection(
    df,
    mass=None,                 # select one M_DM value
    min_efinal=1.0,            # baseline: muon reaches detector with E > this (GeV)
    inner_detector=False,      # require hit_inner_detector
    detector_y_max=None,       # require ip_y0 < this  (m)
    constrain_energy=False,    # |e1 - M_DM/2| / M_DM < 1%
):
    """Build the standard selection mask.  Returns (mask, sel_params).

    sel_params is a small dict describing the cuts, suitable for merging
    into a plot tag via selection_tag()."""
    mask = pd.Series(True, index=df.index)
    if min_efinal is not None:
        mask &= df["efinal_mu1"] > min_efinal
    if mass is not None:
        mask &= _match(df["M_DM"], mass)
    if constrain_energy:
        if mass is None:
            raise ValueError("constrain_energy requires mass=")
        mask &= (np.abs(df["e1"] - mass / 2) / mass) < 0.01
    if detector_y_max is not None:
        mask &= df["ip_y0"] < detector_y_max
    if inner_detector:
        mask &= df["hit_inner_detector"]

    sel_params = {
        "mass": mass,
        "min_efinal": min_efinal,
        "inner_detector": inner_detector,
        "detector_y_max": detector_y_max,
        "constrain_energy": constrain_energy,
    }
    return mask, sel_params


def selection_tag(sel_params):
    """Short tag fragment for the cuts; only non-default cuts appear."""
    parts = []
    if sel_params.get("mass") is not None:
        parts.append(f"mass-{fmt_token(sel_params['mass'])}")
    if sel_params.get("min_efinal") not in (None, 1.0):
        parts.append(f"efin-{fmt_token(sel_params['min_efinal'])}")
    if sel_params.get("inner_detector"):
        parts.append("innertk")
    if sel_params.get("detector_y_max") is not None:
        parts.append(f"ymax-{fmt_token(sel_params['detector_y_max'])}")
    if sel_params.get("constrain_energy"):
        parts.append("econ")
    return "_".join(parts)


# ----------------------------------------------------------------------------
# figures:  name from the same params that selected the data + JSON sidecar
# ----------------------------------------------------------------------------
# default file format for saved figures; override globally with e.g.
#   eio.DEFAULT_FIG_EXT = "pdf"
# so the diagnostics functions (which call save_figure without an ext) emit
# vector PDFs instead of PNGs without touching their signatures.
DEFAULT_FIG_EXT = "png"


def figure_name(stem, params, sel_params=None, ext=None):
    ext = ext or DEFAULT_FIG_EXT
    t = tag(params)
    st = selection_tag(sel_params) if sel_params else ""
    name = f"{stem}_{t}" + (f"_{st}" if st else "")
    return f"{name}.{ext}"


def save_figure(fig, stem, params, sel_params=None, plotdir="plots", dpi=150,
                ext=None):
    """Save fig as <plotdir>/<stem>_<tag>[_<seltag>].<ext> plus a .json sidecar
    recording the full provenance (params, cuts, source file, time).  ext
    defaults to DEFAULT_FIG_EXT ('png' unless overridden globally)."""
    plotdir = Path(plotdir)
    plotdir.mkdir(parents=True, exist_ok=True)
    name = figure_name(stem, params, sel_params, ext=ext)
    out = plotdir / name
    fig.savefig(out, dpi=dpi)

    sidecar = {
        "figure": str(out),
        "params": params,
        "selection": sel_params,
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    with open(out.with_suffix(".json"), "w") as f:
        json.dump(sidecar, f, indent=2, default=str)
    return out


# ----------------------------------------------------------------------------
# unique parameter combinations
# ----------------------------------------------------------------------------
# shorthand names that expand to their min/max column pairs in combos()
_COMBO_GROUPS = {
    "depth": ("depth_min", "depth_max"),
    "mDM": ("mDM_min", "mDM_max"),
    "mA": ("mA_min", "mA_max"),
}

# a sensible default "run matrix" view
_COMBO_DEFAULT_FIELDS = ("dm_model", "mA", "mDM", "depth", "disk_radius")


def combos(catalog_or_root="data", *fields, **filters):
    """Unique combinations of the named parameter fields, with per-combo
    aggregates: n_files, total n_generated (exposure), total n_rows.

    fields may be raw catalog columns ("dm_model", "disk_radius", "run_id",
    ...) or the shorthands "depth", "mDM", "mA", which expand to their
    min/max pairs and also get a collapsed display column ('-1000 (disc)'
    vs '[-4000, -8]').  With no fields, shows the full run matrix
    (dm_model x mA x mDM x depth x disk_radius).

    Keyword args are equality filters, same as select()/summary().

    Examples:
        eio.combos(cat, "mDM", dm_model="core")
        eio.combos(cat, "depth", "mDM", dm_model="core", stage="combined")
        eio.combos("data")                      # the whole run matrix
        eio.combos(cat, "depth", "run_id", dm_model="core")

    Returns a plain dataframe -- chain further with .query()/.sort_values()
    as needed.
    """
    if isinstance(catalog_or_root, (str, Path)):
        cat = read_catalog(catalog_or_root)
    else:
        cat = catalog_or_root
    if len(cat) == 0:
        print("catalog is empty")
        return cat

    sel = select(cat, **filters) if filters else cat
    if len(sel) == 0:
        print(f"no catalog entries match {filters}")
        return sel

    if not fields:
        fields = _COMBO_DEFAULT_FIELDS

    group_cols, display_pairs = [], []
    for f in fields:
        if f in _COMBO_GROUPS:
            lo, hi = _COMBO_GROUPS[f]
            for c in (lo, hi):
                if c not in sel.columns:
                    raise KeyError(f"catalog has no column {c!r}")
                if c not in group_cols:
                    group_cols.append(c)
            display_pairs.append((f, lo, hi))
        elif f in sel.columns:
            if f not in group_cols:
                group_cols.append(f)
        else:
            raise KeyError(
                f"{f!r} is neither a combo shorthand "
                f"({sorted(_COMBO_GROUPS)}) nor a catalog column.  "
                f"Columns: {sorted(sel.columns)}")

    agg = {"n_files": ("path", "count")}
    if "n_generated" in sel.columns:
        agg["n_generated"] = ("n_generated", "sum")
    if "n_rows" in sel.columns:
        agg["n_rows"] = ("n_rows", "sum")

    out = (sel.groupby(group_cols, dropna=False)
              .agg(**agg)
              .reset_index()
              .sort_values(group_cols)
              .reset_index(drop=True))

    # collapsed display columns for the min/max shorthands, shown first
    front = []
    for name, lo, hi in display_pairs:
        disc = " (disc)" if name == "depth" else ""
        out[name] = [_fmt_range(a, b, disc_label=disc)
                     for a, b in zip(out[lo], out[hi])]
        front.append(name)
    rest = [c for c in out.columns if c not in front]
    return out[front + rest]



# ----------------------------------------------------------------------------
# dataset summary view (notebook + command line)
# ----------------------------------------------------------------------------
# default column order for the human-facing summary; anything missing from
# the catalog is silently skipped
SUMMARY_COLS = (
    "stage", "kind", "dm_model",
    "mA_min", "mA_max", "mDM_min", "mDM_max",
    "depth_min", "depth_max", "disk_radius",
    "selection", "eloss",
    "n_generated", "n_rows", "created", "path",
)


def _fmt_range(lo, hi, disc_label=""):
    """'-3500 (disc)' for lo==hi, '[-4000, -8]' for a volume/range."""
    if lo is None or (isinstance(lo, float) and np.isnan(lo)):
        return ""
    def f(v):
        v = float(v)
        return str(int(v)) if v.is_integer() else str(v)
    if lo == hi:
        return f(lo) + disc_label
    return f"[{f(lo)}, {f(hi)}]"


def summary(catalog_or_root="data", stage=None, full_paths=False, **filters):
    """One row per dataset, interesting columns only, sorted for browsing.

    Accepts either a catalog dataframe or a data-root path (str/Path).
    Extra keyword args are equality filters, e.g. summary("data",
    dm_model="core").  Returns the dataframe (also nicely displayable in a
    notebook cell).

    Range-valued parameters are collapsed into single readable columns:
    a single value means the quantity was fixed (e.g. a flat disc at that
    depth); '[a, b]' means events were generated over that interval.
    NOTE on sign conventions: depths are negative (below ground), so the
    underlying catalog columns depth_min/depth_max are the NUMERIC min and
    max -- depth_max is the SHALLOW boundary.  Filter on those columns
    accordingly; this summary view shows the interval explicitly to avoid
    the ambiguity.
    """
    if isinstance(catalog_or_root, (str, Path)):
        cat = read_catalog(catalog_or_root)
    else:
        cat = catalog_or_root
    if len(cat) == 0:
        print("catalog is empty")
        return cat

    if stage is not None:
        filters["stage"] = stage
    sel = select(cat, **filters) if filters else cat

    out = sel.copy()
    # collapse min/max pairs into single readable columns
    for name, lo, hi, disc in (
        ("depth", "depth_min", "depth_max", " (disc)"),
        ("mDM", "mDM_min", "mDM_max", ""),
        ("mA", "mA_min", "mA_max", ""),
    ):
        if lo in out.columns and hi in out.columns:
            out[name] = [
                _fmt_range(a, b, disc_label=disc)
                for a, b in zip(out[lo], out[hi])
            ]

    cols = [c for c in
            ("stage", "kind", "dm_model", "mA", "mDM", "depth",
             "disk_radius", "selection", "eloss",
             "n_generated", "n_rows", "run_id", "created", "path")
            if c in out.columns]
    sort_by = [c for c in ("stage", "dm_model", "mDM_min", "depth_min")
               if c in out.columns]
    if sort_by:
        out = out.sort_values(sort_by)
    out = out[cols]
    if not full_paths and "path" in out.columns:
        out["path"] = out["path"].map(lambda p: Path(p).name)
    return out.reset_index(drop=True)


def _cli():
    """Command line:  python earthshine_io.py [data_root] [options]

      python earthshine_io.py data                      # summary table
      python earthshine_io.py data --stage combined
      python earthshine_io.py data --filter dm_model=core --filter mA_min=0.22
      python earthshine_io.py data --rebuild             # rescan footers
      python earthshine_io.py data --full-paths
    """
    import argparse
    ap = argparse.ArgumentParser(
        description="List EarthShine MC datasets and their parameters.",
        epilog=_cli.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_root", nargs="?", default="data")
    ap.add_argument("--stage", default=None,
                    help="combined | shard | derived")
    ap.add_argument("--filter", action="append", default=[],
                    metavar="KEY=VALUE",
                    help="params equality filter; repeatable")
    ap.add_argument("--rebuild", action="store_true",
                    help="rescan footers and rewrite catalog.parquet first")
    ap.add_argument("--full-paths", action="store_true")
    args = ap.parse_args()

    if args.rebuild:
        cat = build_catalog(args.data_root)
    else:
        cat = read_catalog(args.data_root)

    filters = {}
    for f in args.filter:
        k, _, v = f.partition("=")
        if not _:
            raise SystemExit(f"--filter needs KEY=VALUE, got {f!r}")
        # best-effort numeric coercion so mA_min=0.22 compares as a float
        try:
            v = float(v) if "." in v or "e" in v.lower() else int(v)
        except ValueError:
            pass
        filters[k] = v

    out = summary(cat, stage=args.stage, full_paths=args.full_paths, **filters)
    with pd.option_context("display.max_rows", None,
                           "display.max_columns", None,
                           "display.width", None,
                           "display.max_colwidth", 60):
        print(out.to_string())
    print(f"\n{len(out)} dataset(s)")


# ----------------------------------------------------------------------------
# derived kinematic columns (computed in exactly one place)
# ----------------------------------------------------------------------------
def derive_columns(df):
    """Add the derived columns the diagnostics use.  Idempotent.

    eta is computed as arctanh(pz/|p|) with |p| recomputed from the
    components (the old p_eta_CMS expression was unused test code).
    """
    out = df
    out["rho0_origin"] = np.sqrt(out["x0"] ** 2 + out["y0"] ** 2)
    out["phi0_origin"] = np.arctan2(out["y0"], out["x0"])
    for i in (1, 2):
        if f"px{i}" not in out.columns:
            continue
        pmag = np.sqrt(out[f"px{i}"] ** 2 + out[f"py{i}"] ** 2
                       + out[f"pz{i}"] ** 2)
        out[f"pt{i}_xy"] = np.sqrt(out[f"px{i}"] ** 2 + out[f"py{i}"] ** 2)
        out[f"phi{i}_xy"] = np.arctan2(out[f"py{i}"], out[f"px{i}"])
        with np.errstate(divide="ignore", invalid="ignore"):
            out[f"eta{i}_derived"] = np.arctanh(out[f"pz{i}"] / pmag)
    return out


if __name__ == "__main__":
    _cli()
