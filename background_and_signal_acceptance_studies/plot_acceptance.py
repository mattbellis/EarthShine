"""
plot_acceptance -- acceptance vs M_DM plots for EarthShine, old style, new files.

Workflow:
    1. Select files (via eio.summary or an explicit list).
    2. acceptance_v2.get_df_of_acceptances pools counts across runs, keyed
       by rock volume.
    3. plot_acceptances makes one figure per rock volume in the legacy style
       (blue/red/green = 10/100/1000 GeV cuts; circles/triangles = all muons,
       squares/plus/right-triangles = inner-detector hits), and saves both a
       .png and the underlying .parquet next to it.

Typical use:

    import earthshine_io as eio
    from acceptance_v2 import get_df_of_acceptances
    from plot_acceptance import plot_acceptances, make_acceptance_plots

    # One-liner driver (uses eio.summary under the hood):
    make_acceptance_plots("data", dm_model="core",
                          depth_min=-100, depth_max=-8, disk_radius=40,
                          masses=[20, 50, 90], energycut=[10, 100, 1000])

    # Or fully manual:
    files = eio.summary("data", stage="combined", dm_model="core",
                        disk_radius=40, full_paths=True)["path"]
    dfacc = get_df_of_acceptances(files, energycut=[10, 100, 1000])
    plot_acceptances(dfacc, dm_model="core", masses=[20, 50, 90])
"""

import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from acceptance_v2 import FOOTER_KEY, get_df_of_acceptances


def _write_parquet_with_footer(df, path, params):
    """Write df to parquet with EarthShine metadata embedded in the footer."""
    table = pa.Table.from_pandas(df, preserve_index=False)
    meta = dict(table.schema.metadata or {})
    meta[FOOTER_KEY] = json.dumps(params).encode()
    pq.write_table(table.replace_schema_metadata(meta), path)


def read_footer(path):
    """
    Return the decoded EarthShine footer dict from any parquet file
    (generation output or acceptance summary), or None if absent.
    Reads only the schema -- cheap even for large files.
    """
    meta = pq.read_schema(path).metadata or {}
    if FOOTER_KEY not in meta:
        return None
    return json.loads(meta[FOOTER_KEY])

# ---------------------------------------------------------------------------
# Legacy styling
# ---------------------------------------------------------------------------

# ecut -> (color, marker) for the "all muons above cut" series
_STYLES = {10: ("b", "o"), 100: ("r", "^"), 1000: ("g", "v")}
# ecut -> (color, marker) for the "hit inner detector" series
_STYLES_ID = {10: ("b", "s"), 100: ("r", "P"), 1000: ("g", ">")}

# fallbacks for non-standard ecut lists
_FALLBACK_COLORS = ["b", "r", "g", "m", "c", "y", "k"]
_FALLBACK_MARKERS = ["o", "^", "v", "D", "*", "X", "d"]
_FALLBACK_MARKERS_ID = ["s", "P", ">", "<", "p", "h", "8"]

# dm_model -> (label used in title/filename, default ylim)
MODEL_INFO = {
    "core": ("core", (1e-4, 0.2)),
    "floating": ("floating", (2e-9, 2e-5)),
    "momentum": ("mono-energetic", (2e-9, 2e-5)),
    "momentum_constrained": ("mono-energetic", (2e-9, 2e-5)),
}


def _style(idx, ecut, id_series):
    table = _STYLES_ID if id_series else _STYLES
    if ecut in table:
        return table[ecut]
    color = _FALLBACK_COLORS[idx % len(_FALLBACK_COLORS)]
    markers = _FALLBACK_MARKERS_ID if id_series else _FALLBACK_MARKERS
    return color, markers[idx % len(markers)]


def _fmt_num(x):
    """4000.0 -> '4000', 0.22 -> '0.22'"""
    return f"{x:g}"


def _path_token(x):
    """-4000.0 -> 'm4000', 0.22 -> '0p22'  (matches the data path convention)"""
    return _fmt_num(x).replace("-", "m").replace(".", "p")


def _depth_title(depth_min, depth_max):
    if depth_min == depth_max:
        return f"depth={_fmt_num(depth_min)} m"
    return f"depth=[{_fmt_num(depth_min)}, {_fmt_num(depth_max)}] m"


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_acceptances(
    dfacc,
    dm_model,
    energycut=(10, 100, 1000),
    masses=None,
    ylim=None,
    xlim=None,
    outdir=".",
    save=True,
    user_file_tag="",
    mass_rtol=1e-3,
):
    """
    Make one legacy-style acceptance figure per rock volume in dfacc.

    Parameters
    ----------
    dfacc : DataFrame
        Output of acceptance_v2.get_df_of_acceptances. May contain several
        rock volumes; each gets its own figure.
    dm_model : str
        'core', 'floating', or 'momentum_constrained'/'momentum'. Sets the
        title tag and default ylim.
    energycut : sequence
        Energy cuts to plot; must match columns present in dfacc.
    masses : sequence or None
        If given, restrict the plot to these M_DM points (matched with
        relative tolerance mass_rtol, so 20 matches 20.0).
    ylim, xlim : tuple or None
        Overrides. ylim=None uses the per-model legacy default; xlim=None
        autoscales.
    outdir : str or Path
        Where to write the .png and .parquet outputs.
    save : bool
        Write outputs to disk (in addition to returning the figures).

    Returns
    -------
    list of dicts, one per rock volume, with keys:
        'figure'  : the matplotlib Figure
        'png'     : Path to the saved image (None if save=False)
        'parquet' : Path to the saved summary parquet (None if save=False)
        'footer'  : the metadata dict (also embedded in the parquet footer)
        'data'    : the plotted DataFrame slice for that volume
    """
    tag, default_ylim = MODEL_INFO.get(dm_model, (dm_model, None))
    if ylim is None:
        ylim = default_ylim

    df = dfacc
    if masses is not None:
        masses = np.asarray(masses, dtype=float)
        keep = df["M_DM"].apply(
            lambda m: bool(np.any(np.isclose(m, masses, rtol=mass_rtol)))
        )
        df = df[keep]
        if df.empty:
            raise ValueError(
                f"No rows left after mass selection {list(masses)}; "
                f"available masses: {sorted(dfacc['M_DM'].unique())}"
            )

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    results = []
    vol_keys = ["depth_min", "depth_max", "disk_radius", "volume_m3"]

    n_volumes = df.groupby(vol_keys).ngroups
    if n_volumes > 1:
        warnings.warn(
            f"Selection spans {n_volumes} rock volumes; producing one figure/parquet "
            "per volume. Downstream rate scaling (e.g. DarkCapPy) expects a single "
            "volume -- tighten depth/diskR filters if this was unintentional.",
            stacklevel=2,
        )

    for (dmin, dmax, radius, volume), sub in df.groupby(vol_keys):
        sub = sub.sort_values("M_DM")

        fig = plt.figure(figsize=(8, 4))
        ax = plt.gca()

        for idx, ecut in enumerate(energycut):
            c, m = _style(idx, ecut, id_series=False)
            ax.scatter(sub["M_DM"], sub[f"frac_ecut{ecut}"], color=c, marker=m,
                       s=30, label=rf"$E_{{\mu}}$ > {ecut} GeV")
        for idx, ecut in enumerate(energycut):
            c, m = _style(idx, ecut, id_series=True)
            ax.scatter(sub["M_DM"], sub[f"frac_hit_id_ecut{ecut}"], color=c,
                       marker=m, s=30, label=rf"$E_{{\mu}}$ > {ecut} GeV (ID)")

        ax.set_title(
            f"DM model: {tag}     {_depth_title(dmin, dmax)}, "
            f"diskR={_fmt_num(radius)} m, volume={volume / 1e9:.3f} km$^3$"
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$M_{DM}$ GeV/$c^2$", fontsize=16)
        ax.set_ylabel("acceptance (frac)", fontsize=16)
        if ylim is not None:
            ax.set_ylim(*ylim)
        if xlim is not None:
            ax.set_xlim(*xlim)
        ax.legend(fontsize=10)
        plt.tight_layout()

        footer = {
            "stage": "acceptance_summary",
            "dm_model": dm_model,
            "depth_min": float(dmin),
            "depth_max": float(dmax),
            "disk_radius": float(radius),
            "volume_m3": float(volume),
            "energycut": [float(e) for e in energycut],
            "masses": [float(m) for m in sorted(sub["M_DM"].unique())],
            "created": datetime.now(timezone.utc).isoformat(),
        }

        png_path = pq_path = None
        if save:
            stem = (
                f"acc_dm_model_{tag}_depth={_path_token(dmin)}-{_path_token(dmax)}"
                f"_diskR={_path_token(radius)}{user_file_tag}"
            )
            png_path = outdir / f"{stem}.png"
            fig.savefig(png_path)

            pq_path = outdir / f"{stem}.parquet"
            _write_parquet_with_footer(sub, pq_path, footer)
            print(f"Wrote {png_path} and {pq_path}")

        results.append(
            {"figure": fig, "png": png_path, "parquet": pq_path,
             "footer": footer, "data": sub}
        )

    return results


# ---------------------------------------------------------------------------
# Driver: selection -> efficiencies -> plots in one call
# ---------------------------------------------------------------------------

def make_acceptance_plots(
    data_dir="data",
    dm_model="core",
    depth_min=None,
    depth_max=None,
    disk_radius=None,
    inner_detector_radius=None,
    inner_detector_half_len=None,
    masses=None,
    energycut=(10, 100, 1000),
    files=None,
    user_file_tag="",
    outdir=".",
    verbose=True,
    **plot_kwargs,
):
    """
    Select combined files, compute pooled efficiencies, and plot.

    Selection is done with eio.summary(data_dir, stage='combined', ...),
    passing dm_model / depth_min / depth_max / disk_radius filters when
    given. Alternatively, pass an explicit `files` list to skip eio.

    Note the depth convention: depth_max is the numerically larger
    (shallower) boundary, e.g. depth_min=-4000, depth_max=-8.

    Returns (dfacc, results) where results is the list of per-volume dicts
    from plot_acceptances (figure, png/parquet paths, footer, data).
    """
    if files is None:
        import earthshine_io as eio

        filters = {"dm_model": dm_model}
        if depth_min is not None:
            filters["depth_min"] = depth_min
        if depth_max is not None:
            filters["depth_max"] = depth_max
        if disk_radius is not None:
            filters["disk_radius"] = disk_radius
        if inner_detector_radius is not None:
            filters["inner_detector_radius"] = inner_detector_radius
        if inner_detector_half_len is not None:
            filters["inner_detector_half_len"] = inner_detector_half_len


        files = list(
            eio.summary(data_dir, stage="combined", full_paths=True, **filters)["path"]
        )
        if verbose:
            print(f"Selected {len(files)} file(s):")
            for f in files:
                print(f"  {f}")
        if not files:
            raise FileNotFoundError(f"No combined files matched {filters} in {data_dir}")

    dfacc = get_df_of_acceptances(files, energycut=energycut, verbose=verbose)
    figs = plot_acceptances(
        dfacc, dm_model, energycut=energycut, masses=masses, outdir=outdir, user_file_tag=user_file_tag,
        **plot_kwargs,
    )
    return dfacc, figs
