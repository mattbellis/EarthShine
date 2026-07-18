"""
diagnostics_v2.py
=================
Refactored diagnostics following the new pattern.  The key change relative
to the old diagnostics.py:

    OLD:  function(depth=-8, diskR=500, tag='mDM_..._dm_model_floating', ...)
          -> rebuilds a filename fragment, reads the file itself, appends
             selection suffixes with per-function string surgery

    NEW:  function(df, params, <selection kwargs>, plotdir=...)
          -> I/O lives in earthshine_io.load(); selection lives in
             earthshine_io.apply_selection(); plot names come from
             earthshine_io.save_figure() using the SAME params that
             selected the data.

Only e_and_pt_plots was the original template; kinematic_diagnostic and
origin_plots are now ported as well.  All three share the same signature
pattern:  (df, params, <selection kwargs>, plotdir=...).

Notebook usage:

    import earthshine_io as eio
    import diagnostics_v2 as diag

    cat = eio.read_catalog("data")
    cat[["dm_model", "mDM_min", "mDM_max", "depth_min", "depth_max",
         "disk_radius", "n_rows", "path"]]          # <- "what do I have?"

    df, params = eio.load(cat, dm_model="floating", mA_min=0.22,
                          stage="combined")
    diag.e_and_pt_plots(df, params, plotdir="plots")
"""

import matplotlib.pylab as plt
import numpy as np

import earthshine_io as eio


def _disk_radius(params, df):
    """disk_radius from metadata, falling back to a data-driven estimate
    for files written before auto-computed radii were recorded (the
    generator computes one from angular considerations when --disk-radius
    is omitted; older v2 files stored None).  The estimate -- the largest
    origin radius among DETECTED events -- is a lower bound on the true
    generation radius, fine for plot ranges.  Fix the file properly with
    eio.restamp(path, disk_radius=...)."""
    r = params.get("disk_radius")
    if r is None or (isinstance(r, float) and np.isnan(r)):
        r = float(np.sqrt(df["x0"] ** 2 + df["z0"] ** 2).max())
        print(f"disk_radius missing from metadata; using data-driven "
              f"~{r:.0f} m for plot ranges (a lower bound; restamp the "
              f"file to fix permanently)")
    return float(r)


def _figure_title(params, mass):
    """Standard figure title.  Single-depth (flat disc) runs show one
    value; volume runs show the interval.  This is the one place titles
    are formatted -- keep the convention here and every plot follows."""
    depth_tag = rf'depth=[{params["depth_min"]}, {params["depth_max"]}] m  '
    if params["depth_min"] == params["depth_max"]:
        depth_tag = rf'depth={params["depth_min"]} m  '
    title = (rf'$M_{{DM}}$={int(mass)} GeV/$c^2$ '
             rf'{depth_tag}'
             rf'diskR={params["disk_radius"]}')
    return title


##########################################################################
def e_and_pt_plots(df, params, masses=None,
                   inner_detector=False, detector_y_max=None,
                   constrain_energy=False, min_efinal=1.0,
                   plotdir='plots', max_xval=None):
    """Energy and pT comparison plots, one figure per DM mass.

    df, params come from earthshine_io.load().  Selection kwargs mirror the
    old INNER_DETECTOR_FILTER / DETECTOR_Y_CONSTRAINED / CONSTRAIN_ENERGY
    flags but now flow through apply_selection so the cut and the plot
    name cannot diverge.
    """
    params = eio.as_params(params)   # accept dict or load_many list
    df = eio.derive_columns(df)

    print(f"dataset: {eio.tag(params)}")
    print(f"shape:   {df.shape}")
    print(f"masses in file: {sorted(df['M_DM'].unique())}")
    if masses is None:
        masses = sorted(df['M_DM'].unique())

    figs = []
    for mass in masses:
        mask, sel = eio.apply_selection(
            df, mass=mass, min_efinal=min_efinal,
            inner_detector=inner_detector,
            detector_y_max=detector_y_max,
            constrain_energy=constrain_energy,
        )
        d = df[mask]
        print(f"M_DM = {mass}:  {len(d)} events after selection")

        ######################################################################
        fig = plt.figure(figsize=(12, 4))

        nbins = 50
        xranges = (0, max_xval if max_xval is not None else 1.5 * mass / 2)

        plt.subplot(1, 2, 1)
        d['e1'].hist(bins=nbins, range=xranges, density=True,
                     histtype="step", linewidth=2.5, label='Orig. energy')
        d['efinal_mu1'].hist(bins=nbins, range=xranges, density=True,
                             histtype="step", linewidth=2.5,
                             label='Energy at detector')
        plt.xlabel(r'$E_{\mu}$ (GeV)', fontsize=14)
        plt.legend()
        plt.yscale('log')

        plt.subplot(1, 2, 2)
        d['pt1_detector_acceptance'].hist(
            bins=nbins, range=xranges, density=True,
            histtype="step", linewidth=2.5, label='Ignoring eloss')
        d['pt1_detector_acceptance_eloss'].hist(
            bins=nbins, range=xranges, density=True,
            histtype="step", linewidth=2.5, label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=14)
        plt.legend()
        plt.yscale('log')

        title = _figure_title(params, mass)
        plt.suptitle(title, fontsize=14)
        plt.tight_layout()

        out = eio.save_figure(fig, "e_and_pt", params, sel, plotdir=plotdir)
        print(f"  -> {out}")
        figs.append(fig)

    return df, figs


##########################################################################
def kinematic_diagnostic(df, params, masses=None,
                         inner_detector=False, detector_y_max=None,
                         constrain_energy=False, min_efinal=1.0,
                         plotdir='plots'):
    """Origin/entry-point geometry figure + kinematic-distribution figure,
    one pair per DM mass.

    Port of the old kinematic_diagnostic(depth=, diskR=, tag=, ...).
    NOTE: the old version's first figure had no mass in its filename, so
    looping over masses overwrote the same PNG -- only the last mass
    survived.  Fixed here: the selection tag includes the mass.
    """
    params = eio.as_params(params)   # accept dict or load_many list
    df = eio.derive_columns(df)
    disk_radius = _disk_radius(params, df)

    print(f"dataset: {eio.tag(params)}")
    print(f"shape:   {df.shape}")
    print(f"masses in file: {sorted(df['M_DM'].unique())}")
    if masses is None:
        masses = sorted(df['M_DM'].unique())

    figs = []
    for mass in masses:
        mask, sel = eio.apply_selection(
            df, mass=mass, min_efinal=min_efinal,
            inner_detector=inner_detector,
            detector_y_max=detector_y_max,
            constrain_energy=constrain_energy,
        )
        d = df[mask]
        print(f"M_DM = {mass}:  {len(d)} events after selection")

        ######################################################################
        # Figure 1: decay origins and detector entry points
        ######################################################################
        fig = plt.figure(figsize=(12, 12))

        x, y, z = d['x0'], d['y0'], d['z0']
        R = disk_radius

        plt.subplot(3, 3, 1)
        plt.hist2d(x=x, y=y, bins=100, range=([-R, R], [-2500, 0]))
        plt.xlabel(r'Origin x (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)

        plt.subplot(3, 3, 2)
        plt.hist2d(x=z, y=y, bins=100, range=([-R, R], [-2500, 0]))
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)

        plt.subplot(3, 3, 3)
        plt.hist2d(x=z, y=x, bins=100, range=([-R, R], [-R, R]))
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin x (m)', fontsize=14)

        ######################################################################
        nbins = 50

        plt.subplot(3, 3, 4)
        d['x0'].hist(bins=nbins, range=(-50, 50), histtype="step",
                     density=False, linewidth=2.5, label='Origin')
        d['ip_x0'].hist(bins=nbins, range=(-50, 50), histtype="step",
                        density=False, linewidth=2.5, label='Detector entry')
        plt.xlabel(r'Origin x (m)', fontsize=14)
        plt.yscale('log')
        plt.legend()

        plt.subplot(3, 3, 5)
        d['y0'].hist(bins=nbins, range=(-50, 50), histtype="step",
                     density=False, linewidth=2.5, label='Origin')
        d['ip_y0'].hist(bins=nbins, range=(-50, 50), histtype="step",
                        density=False, linewidth=2.5, label='Detector entry')
        plt.xlabel(r'Origin y (m)', fontsize=14)
        plt.yscale('log')
        plt.legend()

        plt.subplot(3, 3, 6)
        d['z0'].hist(bins=200, range=(-50, 50), histtype="step",
                     density=False, linewidth=2.5, label='Origin')
        d['ip_z0'].hist(bins=200, range=(-50, 50), histtype="step",
                        density=False, linewidth=2.5, label='Detector entry')
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.yscale('log')
        plt.legend()

        ######################################################################
        inner_mask = d['hit_inner_detector']

        plt.subplot(3, 3, 7)
        plt.gca().set_aspect('equal')
        d.plot.scatter(x='ip_x0', y='ip_y0', s=1, ax=plt.gca(),
                       label='muon system')
        d.plot.scatter(x='ip_inner_x0', y='ip_inner_y0', color='r',
                       marker='s', s=1, ax=plt.gca(), label='inner tk')
        d[inner_mask].plot.scatter(x='ip_x0', y='ip_y0', color='k', s=20,
                                   ax=plt.gca(),
                                   label='muon sys. when inner tk')
        plt.legend()
        plt.xlim(-9, 9)
        plt.ylim(-9, 9)
        plt.gca().set_xlabel('x (meters)', fontsize=14)
        plt.gca().set_ylabel('y (meters)', fontsize=14)

        plt.subplot(3, 3, 8)
        plt.gca().set_aspect('equal')
        d.plot.scatter(y='ip_x0', x='ip_z0', s=1, ax=plt.gca(),
                       label='muon system')
        d.plot.scatter(y='ip_inner_x0', x='ip_inner_z0', color='r',
                       marker='s', s=1, ax=plt.gca(), label='inner tk')
        d[inner_mask].plot.scatter(y='ip_x0', x='ip_z0', color='k', s=20,
                                   ax=plt.gca(),
                                   label='muon sys. when inner tk')
        plt.legend(loc='upper left')
        plt.xlim(-20, 20)
        plt.ylim(-20, 20)
        plt.gca().set_ylabel('x (meters)', fontsize=14)
        plt.gca().set_xlabel('z (meters)', fontsize=14)

        plt.subplot(3, 3, 9)
        plt.gca().set_aspect('equal')
        if inner_mask.sum() > 0:
            plt.gca().hist2d(d[inner_mask]['ip_z0'], d[inner_mask]['ip_x0'],
                             bins=100, range=([-20, 20], [-20, 20]),
                             norm='log')
        plt.xlim(-20, 20)
        plt.ylim(-20, 20)
        plt.xlabel('z (meters)', fontsize=14)
        plt.ylabel('x (meters)', fontsize=14)

        title = _figure_title(params, mass)
        plt.suptitle(title, fontsize=14)
        
        plt.tight_layout()
        out = eio.save_figure(fig, "origin_entry", params, sel,
                              plotdir=plotdir)
        print(f"  -> {out}")
        figs.append(fig)

        ######################################################################
        # Figure 2: kinematic distributions (mosaic)
        ######################################################################
        layout = [
            ['A', 'A', 'B', 'B', 'C', 'C'],
            ['D', 'D', 'D', 'E', 'E', 'E'],
            ['F', 'F', 'F', 'G', 'G', 'G'],
        ]
        fig2, axes = plt.subplot_mosaic(layout, figsize=(12, 12))

        plt.sca(axes['A'])
        d['theta1'].hist(bins=64, range=(0, 3.2), histtype="step",
                         density=True, linewidth=2.5)
        plt.xlabel(r'$p_{\theta}$', fontsize=14)

        plt.sca(axes['B'])
        d['eta1'].hist(bins=64, range=(0, 3), density=True,
                       histtype="step", linewidth=2.5)
        plt.xlabel(r'$p_{\eta}$', fontsize=14)
        plt.yscale('log')

        plt.sca(axes['C'])
        d['phi1'].hist(bins=64, range=(0, 3.2), density=True,
                       histtype="step", linewidth=2.5)
        plt.xlabel(r'$p_{\phi}$', fontsize=14)
        plt.yscale('log')

        nbins = 50
        xranges = (0, 1.1 * mass / 2)

        plt.sca(axes['D'])
        d['e1'].hist(bins=nbins, range=xranges, histtype="step",
                     density=True, linewidth=2.5, label='Orig. energy')
        d['efinal_mu1'].hist(bins=nbins, range=xranges, histtype="step",
                             density=True, linewidth=2.5,
                             label='Energy at detector')
        plt.xlabel(r'$E_{\mu}$ (GeV)', fontsize=14)
        plt.legend()

        plt.sca(axes['E'])
        d['e1'].hist(bins=nbins, range=xranges, density=True,
                     histtype="step", linewidth=2.5, label='Orig. energy')
        d['efinal_mu1'].hist(bins=nbins, range=xranges, density=True,
                             histtype="step", linewidth=2.5,
                             label='Energy at detector')
        plt.xlabel(r'$E_{\mu}$ (GeV)', fontsize=14)
        plt.legend()
        plt.yscale('log')

        plt.sca(axes['F'])
        d['pt1_detector_acceptance'].hist(
            bins=nbins, range=xranges, density=True,
            histtype="step", linewidth=2.5, label='Ignoring eloss')
        d['pt1_detector_acceptance_eloss'].hist(
            bins=nbins, range=xranges, density=True,
            histtype="step", linewidth=2.5, label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=14)
        plt.legend()

        plt.sca(axes['G'])
        d['pt1_detector_acceptance'].hist(
            bins=nbins, range=xranges, histtype="step",
            density=True, linewidth=2.5, label='Ignoring eloss')
        d['pt1_detector_acceptance_eloss'].hist(
            bins=nbins, range=xranges, histtype="step",
            density=True, linewidth=2.5, label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=14)
        plt.legend()
        plt.yscale('log')

        plt.tight_layout()
        out = eio.save_figure(fig2, "kinematic", params, sel,
                              plotdir=plotdir)
        print(f"  -> {out}")
        figs.append(fig2)

    return df, figs


##########################################################################
def opening_angle_plots(df, params, masses=None,
                        inner_detector=False, detector_y_max=None,
                        constrain_energy=False, min_efinal=1.0,
                        plotdir='plots',
                        opening_angle_max=1e-2, sep_max=30.0):
    """Di-muon opening angle + separation-at-detector figure, one per mass.

    Port of opening_angle_plots() from
    PROTOTYPING_signal_kinematics_diagnostics_Matt.ipynb into the v2
    convention (df, params, selection kwargs, save_figure).

    Left panel:  the 'opening angle' column (degrees) between the two
    muons.  Right panel: the transverse separation of the muon PAIR at the
    detector, opening_angle x distance_to_detector, in mm.

    NOTE: do NOT use ip_*0 / ip_*1 for this -- those are the ENTRY and EXIT
    points of a SINGLE muon through the detector cylinder (see
    detector_simulation_tools.line_cylinder_intersection: both are O + t*D
    on the same ray), so their distance is that muon's path length through
    the detector (~detector size, and inclination-dependent), not the
    two-muon separation.  The physical pair separation is the small-angle
    projection opening_angle x flight distance, as the prototype used.
    opening_angle_max / sep_max set the histogram ranges (mm).
    """
    params = eio.as_params(params)   # accept dict or load_many list
    df = eio.derive_columns(df)

    # transverse muon-pair separation at the detector (deg -> rad, m -> mm)
    df = df.copy()
    df['sep'] = (np.deg2rad(df['opening angle'])
                 * df['distance_to_detector'] * 1000.0)

    print(f"dataset: {eio.tag(params)}")
    print(f"shape:   {df.shape}")
    print(f"masses in file: {sorted(df['M_DM'].unique())}")
    if masses is None:
        masses = sorted(df['M_DM'].unique())

    figs = []
    for mass in masses:
        mask, sel = eio.apply_selection(
            df, mass=mass, min_efinal=min_efinal,
            inner_detector=inner_detector,
            detector_y_max=detector_y_max,
            constrain_energy=constrain_energy,
        )
        d = df[mask]
        print(f"M_DM = {mass}:  {len(d)} events after selection")

        fig, axes = plt.subplots(1, 2, figsize=(8, 4))

        plt.sca(axes[0])
        d['opening angle'].hist(bins=50, range=(0, opening_angle_max),
                                histtype="step", density=True,
                                linewidth=2.5)
        plt.xlabel('Opening angle (degrees)', fontsize=14)

        plt.sca(axes[1])
        d['sep'].hist(bins=50, range=(0, sep_max), density=True,
                      histtype="step", linewidth=2.5)
        plt.xlabel(r'Separation (mm)', fontsize=14)
        plt.yscale('log')

        plt.suptitle(_figure_title(params, mass), fontsize=14)
        plt.tight_layout()

        out = eio.save_figure(fig, "opening_angle", params, sel,
                              plotdir=plotdir)
        print(f"  -> {out}")
        figs.append(fig)

    return df, figs


##########################################################################
def origin_plots(df, params, masses=None,
                 inner_detector=False, detector_y_max=None,
                 constrain_energy=False, min_efinal=1.0,
                 plotdir='plots', xrange=None, yrange=None, zrange=None,
                 save_histograms=True, out_root="data"):
    """2D origin maps + 1D depth projections, one set per DM mass.

    Port of the old origin_plots(depth=, diskR=, tag=, ...).  The depth
    histograms previously dumped into HISTOGRAM_FILES_FOR_RATES/ are now
    written as a stage='derived' dataset into the partitioned tree via
    save_depth_histograms() (set save_histograms=False to skip).
    """
    params = eio.as_params(params)   # accept dict or load_many list
    df = eio.derive_columns(df)

    print(f"dataset: {eio.tag(params)}")
    print(f"shape:   {df.shape}")
    print(f"masses in file: {sorted(df['M_DM'].unique())}")
    if masses is None:
        masses = sorted(df['M_DM'].unique())

    if xrange is None:
        xrange = _disk_radius(params, df)
    if yrange is None:
        yrange = 1.1 * abs(float(params["depth_min"]))
    if zrange is None:
        zrange = _disk_radius(params, df)

    figs = []
    for mass in masses:
        mask, sel = eio.apply_selection(
            df, mass=mass, min_efinal=min_efinal,
            inner_detector=inner_detector,
            detector_y_max=detector_y_max,
            constrain_energy=constrain_energy,
        )
        d = df[mask]
        print(f"M_DM = {mass}:  {len(d)} events after selection")

        x, y, z = d['x0'], d['y0'], d['z0']

        ######################################################################
        # Figure 1: 2D origin maps
        ######################################################################
        fig = plt.figure(figsize=(12, 4))
        norm = None

        plt.subplot(1, 3, 1)
        plt.hist2d(x=x, y=y, bins=50,
                   range=([-xrange, xrange], [-yrange, 0]), norm=norm)
        plt.xlabel(r'Origin x (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)
        plt.colorbar(label='Counts')

        plt.subplot(1, 3, 2)
        plt.hist2d(x=z, y=y, bins=50,
                   range=([-zrange, zrange], [-yrange, 0]), norm=norm)
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)
        plt.colorbar(label='Counts')

        plt.subplot(1, 3, 3)
        plt.hist2d(x=z, y=x, bins=50,
                   range=([-zrange, zrange], [-xrange, xrange]), norm=norm)
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin x (m)', fontsize=14)
        plt.colorbar(label='Counts')

        title = _figure_title(params, mass)
        plt.suptitle(title, fontsize=14)
        plt.tight_layout()

        out = eio.save_figure(fig, "origin", params, sel, plotdir=plotdir)
        print(f"  -> {out}")
        figs.append(fig)

        ######################################################################
        # Figure 2: 1D depth projections (linear + log-spaced zoom)
        ######################################################################
        fig2 = plt.figure(figsize=(12, 4))

        plt.subplot(1, 2, 1)
        plt.hist(y, bins=100, range=(-yrange, 0), density=True)
        plt.xlabel('Depth (m)', fontsize=18)

        plt.subplot(1, 2, 2)
        bins = np.logspace(np.log10(1), np.log10(4000), 34)
        plt.hist(-y, bins=bins, density=True)
        plt.xscale('log')
        plt.xlabel('Depth (m) [ZOOMED IN]', fontsize=18)
        plt.tight_layout()

        out = eio.save_figure(fig2, "origin_depth_1D", params, sel,
                              plotdir=plotdir)
        print(f"  -> {out}")
        figs.append(fig2)

    ##########################################################################
    # depth histograms for the rates code (replaces HISTOGRAM_FILES_FOR_RATES)
    ##########################################################################
    if save_histograms:
        save_depth_histograms(df, params, masses=masses, out_root=out_root,
                              min_efinal=min_efinal,
                              inner_detector=inner_detector,
                              detector_y_max=detector_y_max,
                              constrain_energy=constrain_energy)

    return df, figs


##########################################################################
def detector_hits_plots(df, params, masses=None,
                        inner_detector=False, detector_y_max=None,
                        constrain_energy=False, min_efinal=1.0,
                        plotdir='plots', xrange=20, zrange=20, norm=None):
    """2D map of muon entry points on the detector (z vs x), one figure
    per DM mass.

    Port of the old detector_hits_plots(depth=, diskR=, tag=, ...).
    Fixes relative to the original:
      * axis labels said 'Origin x/y' but the plotted quantities are the
        detector entry coordinates ip_z0 / ip_x0 -- labels corrected
      * the original wrote an always-empty dataframe to
        HISTOGRAM_FILES_FOR_RATES on every call (vestigial copy-paste
        from origin_plots) -- removed
      * unused position_only / yrange parameters dropped
    xrange/zrange are the half-widths (m) of the plotted window around
    the detector.
    """
    params = eio.as_params(params)   # accept dict or load_many list
    df = eio.derive_columns(df)

    print(f"dataset: {eio.tag(params)}")
    print(f"shape:   {df.shape}")
    print(f"masses in file: {sorted(df['M_DM'].unique())}")
    if masses is None:
        masses = sorted(df['M_DM'].unique())

    figs = []
    for mass in masses:
        mask, sel = eio.apply_selection(
            df, mass=mass, min_efinal=min_efinal,
            inner_detector=inner_detector,
            detector_y_max=detector_y_max,
            constrain_energy=constrain_energy,
        )
        d = df[mask]
        print(f"M_DM = {mass}:  {len(d)} events after selection")

        fig = plt.figure(figsize=(5, 4))

        plt.subplot(1, 1, 1)
        plt.hist2d(x=d['ip_z0'], y=d['ip_x0'], bins=50,
                   range=([-xrange, xrange], [-zrange, zrange]), norm=norm)
        plt.xlabel(r'Detector entry z (m)', fontsize=14)
        plt.ylabel(r'Detector entry x (m)', fontsize=14)
        plt.colorbar(label='Counts')

        title = _figure_title(params, mass)
        plt.suptitle(title, fontsize=10)
        plt.tight_layout()

        out = eio.save_figure(fig, "detector_hits", params, sel,
                              plotdir=plotdir)
        print(f"  -> {out}")
        figs.append(fig)

    return df, figs


##########################################################################
def compare_pt_plots(datasets, labels=None, masses=None,
                     inner_detector=False, detector_y_max=None,
                     constrain_energy=False, min_efinal=1.0,
                     plotdir='plots', max_xval=None):
    """Overlay the muon pT-at-detector (with eloss) distributions from two
    or more datasets, one figure per DM mass.

    datasets: list of (df, params) tuples as returned by eio.load(), e.g.

        d1 = eio.load(cat, dm_model="momentum_constrained", stage="combined")
        d2 = eio.load(cat, dm_model="core", stage="combined")
        diag.compare_pt_plots([d1, d2], masses=[20000.0])

    labels default to each dataset's dm_model (disambiguated with run_id
    if needed).

    Port of the old compare_pt_plots(infiles=[...], labels=[...]), which
    only supported exactly two files.  Fixes relative to the original:
      * masses=None crashed (referenced an undefined name); now defaults
        to the masses common to ALL datasets
      * passing max_xval crashed (plot range never assigned); fixed
      * one figure is now made PER MASS (the original piled every mass
        onto a single axes with repeated labels)
    """
    if len(datasets) < 2:
        raise ValueError("need at least two (df, params) datasets to compare")

    dfs = [eio.derive_columns(df) for df, _ in datasets]
    plist = [eio.as_params(p) for _, p in datasets]

    if labels is None:
        labels = [p["dm_model"] for p in plist]
        if len(set(labels)) < len(labels):   # same model twice -> add run id
            labels = [f"{p['dm_model']} ({p.get('run_id') or i})"
                      for i, p in enumerate(plist)]

    if masses is None:
        common = set(dfs[0]['M_DM'].unique())
        for d in dfs[1:]:
            common &= set(d['M_DM'].unique())
        masses = sorted(common)
        if not masses:
            raise ValueError("datasets share no common M_DM values; "
                             "pass masses= explicitly")

    figs = []
    for mass in masses:
        fig = plt.figure(figsize=(6, 4))
        nbins = 50
        xranges = (0, max_xval if max_xval is not None else 1.5 * mass / 2)

        sel = None
        for d, label in zip(dfs, labels):
            mask, sel = eio.apply_selection(
                d, mass=mass, min_efinal=min_efinal,
                inner_detector=inner_detector,
                detector_y_max=detector_y_max,
                constrain_energy=constrain_energy,
            )
            print(f"M_DM = {mass}, {label}: {mask.sum()} events "
                  f"after selection")
            d[mask]['pt1_detector_acceptance_eloss'].hist(
                bins=nbins, range=xranges, histtype="step",
                density=True, linewidth=2.5, label=label)

        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=14)
        plt.legend()
        plt.yscale('log')
        plt.suptitle(_figure_title(plist[0], mass), fontsize=14)
        plt.tight_layout()

        # name and sidecar reflect the comparison: the dm_model slot in the
        # tag carries all labels, and the sidecar records every dataset's
        # provenance, not just the first
        cmp_params = dict(plist[0])
        cmp_params["dm_model"] = "-vs-".join(labels)
        cmp_params["compare"] = [
            {"label": lab,
             "dm_model": p["dm_model"],
             "run_id": p.get("run_id"),
             "n_generated": p.get("n_generated"),
             "created": p.get("created")}
            for lab, p in zip(labels, plist)
        ]
        out = eio.save_figure(fig, "compare_pt", cmp_params, sel,
                              plotdir=plotdir)
        print(f"  -> {out}")
        figs.append(fig)

    return dfs, figs


##########################################################################
def save_depth_histograms(df, params, masses=None,
                          out_root="data",
                          min_efinal=1.0, **selection_kwargs):
    """Replacement for the histogram dump at the end of origin_plots().

    Writes the depth distributions as a 'derived_histogram' dataset into
    the SAME partitioned tree, with metadata, so the rates code can find
    it through the same catalog instead of a hand-named file in
    HISTOGRAM_FILES_FOR_RATES/.
    """
    import pandas as pd

    params = eio.as_params(params)   # accept dict or load_many list

    if masses is None:
        masses = sorted(df['M_DM'].unique())

    dict_hist = {}
    bins = np.logspace(np.log10(1), np.log10(4000), 34)
    for mass in masses:
        mask, sel = eio.apply_selection(df, mass=mass,
                                        min_efinal=min_efinal,
                                        **selection_kwargs)
        counts, bin_edges = np.histogram(-df[mask]['y0'], bins=bins,
                                         density=True)
        dict_hist[f'{int(mass)} bin_start'] = bin_edges[:-1]
        dict_hist[f'{int(mass)} count'] = counts

    df_hist = pd.DataFrame.from_dict(dict_hist)

    hist_params = dict(params)
    hist_params["stage"] = "derived"
    hist_params["kind"] = "depth_histogram"
    hist_params["selection_cuts"] = sel
    # the histogram contents depend on which masses were requested, not on
    # the full range in the parent file -- reflect that in the identity so
    # different subsets / cuts get distinct paths instead of overwriting
    hist_params["mDM_min"] = float(min(masses))
    hist_params["mDM_max"] = float(max(masses))

    base = eio.rel_path(hist_params)
    seltag = eio.selection_tag({**sel, "mass": None})  # cuts minus the
    if seltag:                                         # per-mass loop value
        base = base.with_name(base.stem + f"-{seltag}" + base.suffix)

    outpath = eio.write_parquet(df_hist, f'{out_root}/{base}', hist_params)
    print(f"Saved depth histograms to {outpath}")
    eio.build_catalog(out_root, verbose=False)
    return df_hist
