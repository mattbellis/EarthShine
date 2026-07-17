"""
lhe_writer.py -- Convert EarthShine muon-pair dataframes to Les Houches Event files.

Takes a dataframe with columns
    ['px1','py1','pz1','e1','px2','py2','pz2','e2','x0','y0','z0']
(one dark photon -> mu+ mu- decay per row) and writes a valid LHE file
that CMSSW's ExternalLHEProducer / Pythia8 LHE reader will accept.

Design notes
------------
* Header + <init> block: copied verbatim from a user-supplied "template" CMS
  LHE file (everything up to and including </init>). If no template is given,
  a minimal self-consistent preamble is generated instead.
* Event record layout (indices are 1-based within the event, per the LHE spec):
      [1] incoming mock particle   (status -1)   -- optional
      [2] incoming mock particle   (status -1)   -- optional
      [3] dark photon A'           (status  2)   -- optional intermediate
      [4] mu (pdg +13 by default)  (status  1)
      [5] mu (pdg -13 by default)  (status  1)
  Mother pointers are wired up appropriately for whichever of the optional
  particles are enabled.
* Mock initial state: two *massless* particles built from the light-cone
  decomposition of the total 4-momentum P:
      p_a = (E+|P|)/2 * (1,  n),   p_b = (E-|P|)/2 * (1, -n),   n = P/|P|
  so that p_a + p_b = P exactly and both are on the massless shell. This keeps
  momentum conservation exact, which Pythia's LHE reader checks. Default IDs
  are (11, -11) so no colour lines are needed; if you pass quark IDs a colour
  line 501 is added automatically.
* Vertex: the LHE standard has NO field for a production vertex. Two options
  are provided (see `vertex_mode`):
    - 'comment' (default): writes a line  `#vertex x y z`  inside each
      <event> block, after the particle lines. LHE parsers ignore '#' lines,
      so the file stays valid; a small custom CMSSW plugin (or your own
      post-GEN step) can read these back to place the primary vertex.
    - 'ctau': encodes the decay displacement as the A' proper lifetime,
      VTIMUP = |r| * mA' / |pA'|  (in mm). NOTE: this is only geometrically
      consistent if the A' is produced at the origin and (x0,y0,z0) lies
      along its momentum direction. For EarthShine events (production at the
      Earth's centre, decay near the detector) this is generally NOT the
      case, so 'comment' is the safer default. A warning is emitted per file
      if the misalignment exceeds ~1 degree on average.
    - 'none': vertex information is dropped.
* Units: momenta/energies/masses in GeV (LHE standard). VTIMUP is in mm.
  `vertex_scale` multiplies (x0,y0,z0) on the way out -- e.g. use 1000.0 if
  your dataframe stores metres and you want mm in the output.

Typical use with the EarthShine I/O layer:

    import earthshine_io as eio
    from lhe_writer import write_lhe

    df, params_list = eio.load_many(selection)
    write_lhe(df, "earthshine_events.lhe",
              template_path="cms_reference.lhe",
              vertex_scale=1000.0)          # m -> mm
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

MUON_MASS = 0.1056583745  # GeV
DEFAULT_APRIME_PDGID = 4900022  # CMS hidden-valley dark photon convention

REQUIRED_COLS = ["px1", "py1", "pz1", "e1", "px2", "py2", "pz2", "e2"]
VERTEX_COLS = ["x0", "y0", "z0"]


# ----------------------------------------------------------------------------
# Preamble handling
# ----------------------------------------------------------------------------

def extract_template_preamble(template_path: str) -> str:
    """Return everything from a reference LHE file up to and including </init>."""
    lines = []
    with open(template_path) as f:
        for line in f:
            lines.append(line)
            if line.strip().startswith("</init>"):
                return "".join(lines)
    raise ValueError(f"No </init> tag found in {template_path!r} -- "
                     "is this a complete LHE file?")


def minimal_preamble(beam_energy: float = 6500.0, comment: str = "") -> str:
    """A minimal but spec-compliant header + <init> block.

    IDWTUP = 3 (unweighted, accept all) keeps Pythia from re-sampling weights.
    The cross-section numbers are placeholders -- they only matter if you rely
    on them downstream for normalisation.
    """
    return (
        '<LesHouchesEvents version="3.0">\n'
        "<header>\n"
        f"<!-- EarthShine A' -> mu mu events. {comment} -->\n"
        "</header>\n"
        "<init>\n"
        f"    2212    2212  {beam_energy:.9e}  {beam_energy:.9e}"
        "    0    0  247000  247000    3    1\n"
        "  1.000000e+00  0.000000e+00  1.000000e+00    1\n"
        "</init>\n"
    )


# ----------------------------------------------------------------------------
# Low-level formatting
# ----------------------------------------------------------------------------

def _pline(pdg, status, m1, m2, c1, c2, px, py, pz, e, m,
           vtim=0.0, spin=9.0) -> str:
    """One particle line in the standard LHE column order."""
    return (f"    {pdg:9d} {status:5d} {m1:5d} {m2:5d} {c1:5d} {c2:5d} "
            f"{px:+.10e} {py:+.10e} {pz:+.10e} {e:+.10e} {m:+.10e} "
            f"{vtim:.4e} {spin:.4e}\n")


def _mass(px, py, pz, e) -> float:
    """Invariant mass with protection against tiny negative arguments."""
    m2 = e * e - (px * px + py * py + pz * pz)
    return float(np.sqrt(max(m2, 0.0)))


# ----------------------------------------------------------------------------
# Main writer
# ----------------------------------------------------------------------------

def write_lhe(df: pd.DataFrame,
              output_path: str,
              template_path: str | None = None,
              include_mother: bool = True,
              include_initial: bool = False,
              initial_pdgids: tuple[int, int] = (11, -11),
              muon_pdgids: tuple[int, int] = (13, -13),
              aprime_pdgid: int = DEFAULT_APRIME_PDGID,
              vertex_mode: str = "comment",
              vertex_scale: float = 1.0,
              weight_col: str | None = None,
              beam_energy: float = 6500.0) -> int:
    """Write the dataframe out as an LHE file. Returns the number of events.

    Parameters
    ----------
    df : dataframe with the EarthShine muon-pair columns.
    output_path : destination .lhe file.
    template_path : existing CMS LHE file whose header/<init> block is reused.
        If None, a minimal preamble is generated.
    include_mother : write the A' as a status-2 intermediate particle whose
        4-momentum and mass are computed from the muon pair. Recommended --
        it documents the resonance and lets Pythia/analysis code find it.
    include_initial : add two mock status -1 incoming particles. Off by
        default; turn on if your CMSSW workflow's LHE reader complains about
        missing beams.
    initial_pdgids : PDG IDs for the mock incoming pair. Leptonic IDs need no
        colour; quark IDs (|id| <= 6) get a shared colour line automatically.
    muon_pdgids : PDG IDs assigned to (muon 1, muon 2). The dataframe does not
        record charge, so by convention leg 1 -> mu- (13) and leg 2 -> mu+
        (-13); flip here if your generator uses the opposite convention.
    aprime_pdgid : PDG ID for the dark photon (default 4900022, the CMS
        hidden-valley convention).
    vertex_mode : 'comment' | 'ctau' | 'none'  (see module docstring).
    vertex_scale : multiplies (x0, y0, z0) before writing (unit conversion).
    weight_col : optional dataframe column to use as XWGTUP (default 1.0).
    beam_energy : only used when generating the minimal preamble.
    """
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Dataframe is missing required columns: {missing}")

    want_vertex = vertex_mode in ("comment", "ctau")
    if want_vertex and any(c not in df.columns for c in VERTEX_COLS):
        raise ValueError(f"vertex_mode={vertex_mode!r} needs columns {VERTEX_COLS}")
    if vertex_mode not in ("comment", "ctau", "none"):
        raise ValueError(f"Unknown vertex_mode {vertex_mode!r}")

    if template_path is not None:
        preamble = extract_template_preamble(template_path)
    else:
        preamble = minimal_preamble(beam_energy=beam_energy,
                                    comment=f"{len(df)} events")

    initial_is_coloured = all(abs(p) <= 6 for p in initial_pdgids)

    misalignments = []
    n_events = 0

    with open(output_path, "w") as out:
        out.write(preamble)

        for row in df.itertuples(index=False):
            p1 = np.array([row.px1, row.py1, row.pz1, row.e1])
            p2 = np.array([row.px2, row.py2, row.pz2, row.e2])
            P = p1 + p2  # A' four-momentum

            m1 = _mass(*p1)
            m2 = _mass(*p2)
            mA = _mass(*P)

            pmag = float(np.linalg.norm(P[:3]))
            if pmag == 0.0:
                raise ValueError("Event with zero total 3-momentum; cannot "
                                 "build light-cone initial state / ctau.")
            n_hat = P[:3] / pmag

            weight = float(getattr(row, weight_col)) if weight_col else 1.0

            # --- assemble particle lines --------------------------------
            lines = []
            idx = 0
            mother_idx = 0  # index of A' (or of nothing) for the muons

            if include_initial:
                # massless light-cone split: p_a + p_b = P exactly
                ea = 0.5 * (P[3] + pmag)
                eb = 0.5 * (P[3] - pmag)
                c1, c2 = (501, 501) if initial_is_coloured else (0, 0)
                anti = (0, 501) if initial_is_coloured else (0, 0)
                lines.append(_pline(initial_pdgids[0], -1, 0, 0, c1, 0,
                                    *(ea * n_hat), ea, 0.0))
                lines.append(_pline(initial_pdgids[1], -1, 0, 0, 0, c1,
                                    *(-eb * n_hat), eb, 0.0))
                idx += 2

            vtim = 0.0
            if vertex_mode == "ctau":
                r = np.array([row.x0, row.y0, row.z0]) * vertex_scale  # -> mm
                L = float(np.linalg.norm(r))
                vtim = L * mA / pmag if pmag > 0 else 0.0
                if L > 0:
                    cosang = float(np.dot(r, n_hat) / L)
                    misalignments.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))

            if include_mother:
                moms = (1, 2) if include_initial else (0, 0)
                lines.append(_pline(aprime_pdgid, 2, *moms, 0, 0,
                                    P[0], P[1], P[2], P[3], mA, vtim=vtim))
                idx += 1
                mother_idx = idx

            mu_moms = (mother_idx, mother_idx) if include_mother else \
                      ((1, 2) if include_initial else (0, 0))
            lines.append(_pline(muon_pdgids[0], 1, *mu_moms, 0, 0,
                                p1[0], p1[1], p1[2], p1[3], m1))
            lines.append(_pline(muon_pdgids[1], 1, *mu_moms, 0, 0,
                                p2[0], p2[1], p2[2], p2[3], m2))

            nup = len(lines)

            # --- event block --------------------------------------------
            out.write("<event>\n")
            out.write(f"    {nup:5d}     1 {weight:+.10e} {mA:.8e} "
                      f"7.2973525e-03 1.1810000e-01\n")
            out.writelines(lines)
            if vertex_mode == "comment":
                out.write(f"#vertex {row.x0 * vertex_scale:.8e} "
                          f"{row.y0 * vertex_scale:.8e} "
                          f"{row.z0 * vertex_scale:.8e}\n")
            out.write("</event>\n")
            n_events += 1

        out.write("</LesHouchesEvents>\n")

    if vertex_mode == "ctau" and misalignments:
        mean_mis = float(np.mean(misalignments))
        if mean_mis > 1.0:
            warnings.warn(
                f"vertex_mode='ctau': decay positions are misaligned with the "
                f"A' momentum by {mean_mis:.1f} deg on average. VTIMUP assumes "
                f"production at the origin with displacement along the "
                f"momentum; consider vertex_mode='comment' instead.")

    return n_events


# ----------------------------------------------------------------------------
# HepMC output (recommended for displaced vertices)
# ----------------------------------------------------------------------------
#
# Unlike LHE, HepMC vertices carry a full (x, y, z, t) position, so the muon
# production point is embedded natively and CMSSW/GEANT will start the tracks
# exactly where you say -- no comment-line hacks or custom vertex plugins.
#
# Two ASCII flavours are supported:
#   version=2 : HepMC2 "IO_GenEvent" format. This is what CMSSW's MCFileSource
#               reads, and is the safe default while the CMSSW HepMC3
#               migration is still in progress.
#   version=3 : HepMC3 "Asciiv3" format, for HepMC3-enabled workflows or
#               standalone HepMC3/Rivet-style tooling.
#
# Event topology per event: one vertex at (x0, y0, z0, t=0) with the A'
# (status 2) as its incoming particle and the two muons (status 1) outgoing.
# Units are declared as GEV / MM in the file, hence the default
# vertex_scale=1000.0 to convert metres (EarthShine convention, CMS-centred
# coordinates) to millimetres.

def _hepmc_event_lines_v2(i, ap, mus, pos, weight):
    """One event in HepMC2 IO_GenEvent format, as a list of lines."""
    nw = 0 if weight is None else 1
    wtail = "" if weight is None else f" {weight:.16e}"
    lines = [f"E {i} -1 {ap['m']:.6e} -1.0 -1.0 0 -1 1 0 0 0 {nw}{wtail}",
             "U GEV MM",
             f"V -1 0 {pos[0]:.16e} {pos[1]:.16e} {pos[2]:.16e} "
             f"{pos[3]:.16e} 1 2 0",
             f"P 10001 {ap['pdg']} {ap['px']:.16e} {ap['py']:.16e} "
             f"{ap['pz']:.16e} {ap['e']:.16e} {ap['m']:.16e} 2 0 0 -1 0"]
    for bc, p in zip((10002, 10003), mus):
        lines.append(f"P {bc} {p['pdg']} {p['px']:.16e} {p['py']:.16e} "
                     f"{p['pz']:.16e} {p['e']:.16e} {p['m']:.16e} 1 0 0 0 0")
    return lines


def _hepmc_event_lines_v3(i, ap, mus, pos, weight):
    """One event in HepMC3 Asciiv3 format, as a list of lines."""
    lines = [f"E {i} 1 3", "U GEV MM"]
    if weight is not None:
        lines.append(f"W {weight:.16e}")
    lines.append(f"P 1 0 {ap['pdg']} {ap['px']:.16e} {ap['py']:.16e} "
                 f"{ap['pz']:.16e} {ap['e']:.16e} {ap['m']:.16e} 2")
    lines.append(f"V -1 0 [1] @ {pos[0]:.16e} {pos[1]:.16e} "
                 f"{pos[2]:.16e} {pos[3]:.16e}")
    for pid, p in zip((2, 3), mus):
        lines.append(f"P {pid} -1 {p['pdg']} {p['px']:.16e} {p['py']:.16e} "
                     f"{p['pz']:.16e} {p['e']:.16e} {p['m']:.16e} 1")
    return lines


def write_hepmc(df: pd.DataFrame,
                output_path: str,
                version: int = 2,
                muon_pdgids: tuple[int, int] = (13, -13),
                aprime_pdgid: int = DEFAULT_APRIME_PDGID,
                vertex_scale: float = 1000.0,
                vertex_time: float = 0.0,
                weight_col: str | None = None) -> int:
    """Write the dataframe as a HepMC ASCII file. Returns the event count.

    Parameters
    ----------
    version : 2 for HepMC2 IO_GenEvent (read by CMSSW's MCFileSource;
        default), 3 for HepMC3 Asciiv3.
    vertex_scale : multiplies (x0, y0, z0) on output. HepMC positions are in
        mm, so the default of 1000.0 converts metre-valued EarthShine vertices
        (CMS-centred coordinates) to mm. Set to 1.0 if your dataframe is
        already in mm.
    vertex_time : the t coordinate written for every vertex, in mm (i.e.
        c*t). Default 0.0, meaning the muons are created "in time" at the
        decay point; adjust if you need a specific timing offset relative to
        the bunch crossing.
    Other parameters as in `write_lhe`. There is no `vertex_mode` here --
    the whole point of HepMC is that the vertex is a first-class citizen --
    and no `include_mother`, since the decay vertex needs the A' as its
    incoming particle to be representable in the format (see note below).
    """
    missing = [c for c in REQUIRED_COLS + VERTEX_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Dataframe is missing required columns: {missing}")
    if version not in (2, 3):
        raise ValueError(f"version must be 2 or 3, got {version!r}")
    # Note there is no include_mother option here: a vertex with zero
    # incoming particles does not survive the ASCII round trip in EITHER
    # HepMC version (verified against the reference library, whose own
    # readers reject such vertices even when written by its own writers).
    # The A' incoming particle is therefore always written.

    if version == 2:
        header = ("HepMC::Version 2.06.09\n"
                  "HepMC::IO_GenEvent-START_EVENT_LISTING\n")
        footer = "HepMC::IO_GenEvent-END_EVENT_LISTING\n"
        event_lines = _hepmc_event_lines_v2
    else:
        header = ("HepMC::Version 3.02.05\n"
                  "HepMC::Asciiv3-START_EVENT_LISTING\n")
        footer = "HepMC::Asciiv3-END_EVENT_LISTING\n"
        event_lines = _hepmc_event_lines_v3

    n_events = 0
    with open(output_path, "w") as out:
        out.write(header)
        for i, row in enumerate(df.itertuples(index=False)):
            p1 = np.array([row.px1, row.py1, row.pz1, row.e1])
            p2 = np.array([row.px2, row.py2, row.pz2, row.e2])
            P = p1 + p2

            ap = dict(pdg=aprime_pdgid, px=P[0], py=P[1], pz=P[2], e=P[3],
                      m=_mass(*P))
            mus = [dict(pdg=muon_pdgids[0], px=p1[0], py=p1[1], pz=p1[2],
                        e=p1[3], m=_mass(*p1)),
                   dict(pdg=muon_pdgids[1], px=p2[0], py=p2[1], pz=p2[2],
                        e=p2[3], m=_mass(*p2))]
            pos = (row.x0 * vertex_scale, row.y0 * vertex_scale,
                   row.z0 * vertex_scale, vertex_time)
            weight = float(getattr(row, weight_col)) if weight_col else None

            out.write("\n".join(event_lines(i, ap, mus, pos, weight)) + "\n")
            n_events += 1
        out.write(footer)

    return n_events


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="Convert an EarthShine muon-pair parquet file to LHE or HepMC.")
    ap.add_argument("input", help="Input parquet file with the muon-pair columns")
    ap.add_argument("output", help="Output file path")
    ap.add_argument("--format", choices=("lhe", "hepmc2", "hepmc3"),
                    default="hepmc2",
                    help="Output format. hepmc2 = IO_GenEvent (CMSSW "
                         "MCFileSource; default); hepmc3 = Asciiv3. HepMC "
                         "embeds the per-event vertex natively.")
    ap.add_argument("--template", default=None,
                    help="[lhe] Reference CMS LHE file to steal the header/<init> from")
    ap.add_argument("--no-mother", action="store_true",
                    help="[lhe] Do not write the A' intermediate particle "
                         "(HepMC always includes it; see write_hepmc docs)")
    ap.add_argument("--initial", action="store_true",
                    help="Add mock status -1 incoming particles")
    ap.add_argument("--initial-ids", type=int, nargs=2, default=(11, -11),
                    metavar=("ID1", "ID2"),
                    help="PDG IDs of the mock incoming pair (default: 11 -11)")
    ap.add_argument("--aprime-id", type=int, default=DEFAULT_APRIME_PDGID,
                    help=f"PDG ID for the A' (default: {DEFAULT_APRIME_PDGID})")
    ap.add_argument("--vertex-mode", choices=("comment", "ctau", "none"),
                    default="comment", help="[lhe] How to record the vertex")
    ap.add_argument("--vertex-scale", type=float, default=1000.0,
                    help="Multiply (x0,y0,z0) by this factor. Default 1000.0 "
                         "converts metres -> mm (the HepMC/LHE length unit).")
    ap.add_argument("--vertex-time", type=float, default=0.0,
                    help="[hepmc] t coordinate (c*t, mm) written for every vertex")
    ap.add_argument("--weight-col", default=None,
                    help="Dataframe column to use as the event weight")
    args = ap.parse_args()

    df = pd.read_parquet(args.input)
    if args.format == "lhe":
        n = write_lhe(df, args.output,
                      template_path=args.template,
                      include_mother=not args.no_mother,
                      include_initial=args.initial,
                      initial_pdgids=tuple(args.initial_ids),
                      aprime_pdgid=args.aprime_id,
                      vertex_mode=args.vertex_mode,
                      vertex_scale=args.vertex_scale,
                      weight_col=args.weight_col)
    else:
        n = write_hepmc(df, args.output,
                        version=2 if args.format == "hepmc2" else 3,
                        aprime_pdgid=args.aprime_id,
                        vertex_scale=args.vertex_scale,
                        vertex_time=args.vertex_time,
                        weight_col=args.weight_col)
    print(f"Wrote {n} events to {args.output} ({args.format})")


if __name__ == "__main__":
    main()
