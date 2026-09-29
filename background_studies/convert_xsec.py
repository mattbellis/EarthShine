#!/usr/bin/env python3
"""
convert_xsec.py -- turn digitised cross-section figures into one nubkg table.

Inputs (raw digitisations, kept in nubkg/data/raw/ for provenance)
------------------------------------------------------------------
pdg_cross_sections.csv
    PDG "Neutrino Cross Section Measurements" review, the sigma/E figure.
    That figure has a SPLIT x-axis: logarithmic 0.1-100 GeV, linear 100-360 GeV.
      LogUpper / LogLower       nu / nubar, log panel   (x = E in GeV)
      UpperLinear / LowerLinear nu / nubar, linear panel (x = E in GeV)
    y = sigma_CC / E in units of 1e-38 cm^2/GeV, per nucleon, isoscalar.
    The linear-panel tracks are the dashed world-average lines (0.677 / 0.334),
    i.e. constants by construction.  The log panel was traced out and back; the
    two passes agree to 0.4% rms and are simply merged.
    https://pdg.lbl.gov/2025/reviews/rpp2025-rev-nu-cross-sections.pdf

ice_cube_cross_section.csv
    IceCube, Nature 551, 596 (2017), arXiv:1711.08119, Fig. 1.
      Neutrino / Antineutrino   Cooper-Sarkar, Mertsch & Sarkar (CSMS) SM
      Weighted Combination      flux-weighted nu/nubar mix used in the analysis
      Result                    IceCube's measured value, 6.3-980 TeV
    x = log10(E/GeV), y = sigma_CC / E in 1e-38 cm^2/GeV.
    Axis interpretation is confirmed by the data themselves: Result/Weighted
    comes out at 1.298, against the published 1.30x SM.

Joining the two
---------------
At ~350 GeV the PDG world average and the CSMS curve disagree by -5.6% (nu)
and +7.9% (nubar).  Opposite signs mean this is not a digitisation offset: it is
a flat 30-350 GeV fit meeting an energy-dependent DIS calculation whose
nubar/nu ratio is already rising.  Rather than a hard switch, the two are
blended linearly in log(E) across BLEND = (340, 600) GeV; the PDG line is a
constant, so extending it to 600 GeV for the blend is harmless.  The residual
~6-8% is carried as the cross-section uncertainty in that window
(see nu_xsec.sigma_uncertainty).

Output
------
nubkg/data/xsec_cc_digitised.csv
    E_GeV, sigma_nu_cm2, sigma_nubar_cm2   on a common log grid
nubkg/data/xsec_icecube_measured.csv
    E_GeV, measured/SM   (IceCube's 1.30x, as digitised; for systematics only)

Usage:  python convert_xsec.py            (paths default to the repo layout)
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RAW = HERE / "nubkg" / "data" / "raw"
OUT = HERE / "nubkg" / "data"

UNIT = 1.0e-38            # y values are sigma/E in 1e-38 cm^2/GeV
BLEND = (340.0, 600.0)    # GeV
GRID = (0.4, 9.9e5, 20)   # E_min, E_max, points per decade


def read_multi(path: Path) -> dict[str, np.ndarray]:
    """Read a digitiser export: header row of track names over X,Y pairs."""
    rows = list(csv.reader(open(path)))
    names = [n.strip() for n in rows[0] if n.strip()]
    out = {}
    for k, name in enumerate(names):
        pts = []
        for r in rows[2:]:
            try:
                pts.append((float(r[2 * k]), float(r[2 * k + 1])))
            except (ValueError, IndexError):
                pass
        a = np.array(pts)
        out[name] = a[np.argsort(a[:, 0])]
    return out


def dedupe(x, y, tol=0.02):
    """Merge points closer than `tol` in log10(x) (out-and-back tracing),
    averaging in log(y)."""
    lx = np.log10(x)
    keep_x, keep_y = [lx[0]], [[np.log(y[0])]]
    for a, b in zip(lx[1:], y[1:]):
        if a - keep_x[-1] < tol:
            keep_y[-1].append(np.log(b))
        else:
            keep_x.append(a)
            keep_y.append([np.log(b)])
    return 10 ** np.array(keep_x), np.exp([np.mean(v) for v in keep_y])


def species_curve(pdg_log, pdg_lin, csms):
    """Return a callable sigma/E(E) for one species, in 1e-38 cm^2/GeV."""
    lo_x, lo_y = dedupe(*np.concatenate([pdg_log, pdg_lin]).T)
    ic_x, ic_y = 10 ** csms[:, 0], csms[:, 1]

    def pdg(e):
        return np.exp(np.interp(np.log(e), np.log(lo_x), np.log(lo_y)))

    def ice(e):
        return np.exp(np.interp(np.log(e), np.log(ic_x), np.log(ic_y)))

    def f(e):
        e = np.asarray(e, float)
        w = np.clip((np.log(e) - np.log(BLEND[0])) / np.log(BLEND[1] / BLEND[0]), 0, 1)
        return np.exp((1 - w) * np.log(pdg(e)) + w * np.log(ice(e)))

    return f, (lo_x, lo_y), (ic_x, ic_y)


def build(raw=RAW):
    pdg = read_multi(raw / "pdg_cross_sections.csv")
    ic = read_multi(raw / "ice_cube_cross_section.csv")
    nu, *_ = species_curve(pdg["LogUpper"], pdg["UpperLinear"], ic["Neutrino"])
    nb, *_ = species_curve(pdg["LogLower"], pdg["LowerLinear"], ic["Antineutrino"])
    n = int(np.log10(GRID[1] / GRID[0]) * GRID[2]) + 1
    e = np.geomspace(GRID[0], GRID[1], n)
    meas = ic["Result"]
    ratio = meas[:, 1] / np.interp(meas[:, 0], ic["Weighted Combination"][:, 0],
                                   ic["Weighted Combination"][:, 1])
    return {"e": e, "nu": nu(e) * e * UNIT, "nubar": nb(e) * e * UNIT,
            "meas_e": 10 ** meas[:, 0], "meas_ratio": ratio,
            "pdg": pdg, "ic": ic, "curves": (nu, nb)}


def write(res, out=OUT):
    p = out / "xsec_cc_digitised.csv"
    with open(p, "w") as fh:
        fh.write("# nu_mu / nubar_mu CC cross section per nucleon (isoscalar)\n"
                 "# digitised by MB; converted by convert_xsec.py\n"
                 "#   E < 340 GeV : PDG review sigma/E figure (world average 30-350 GeV)\n"
                 "#   E > 600 GeV : CSMS SM curves as plotted in IceCube arXiv:1711.08119 Fig 1\n"
                 f"#   {BLEND[0]:.0f}-{BLEND[1]:.0f} GeV : linear blend in log(E); "
                 "the two disagree by -5.6% (nu), +7.9% (nubar) at 350 GeV\n"
                 "# columns: E_GeV, sigma_nu_cm2, sigma_nubar_cm2\n")
        for a, b, c in zip(res["e"], res["nu"], res["nubar"]):
            fh.write(f"{a:.6e},{b:.6e},{c:.6e}\n")
    q = out / "xsec_icecube_measured.csv"
    with open(q, "w") as fh:
        fh.write("# IceCube measured / SM (flux-weighted), arXiv:1711.08119 Fig 1, as digitised\n"
                 "# published: 1.30 +0.21-0.19 (stat) +0.39-0.43 (syst), 6.3-980 TeV\n"
                 "# measurement is of CC+NC; applies as a common scale factor\n"
                 "# columns: E_GeV, ratio_to_SM\n")
        for a, b in zip(res["meas_e"], res["meas_ratio"]):
            fh.write(f"{a:.6e},{b:.6f}\n")
    return p, q


if __name__ == "__main__":
    res = build()
    p, q = write(res)
    print(f"wrote {p.relative_to(HERE)}  ({len(res['e'])} points, "
          f"{res['e'][0]:.2g}-{res['e'][-1]:.3g} GeV)")
    print(f"wrote {q.relative_to(HERE)}  (IceCube measured/SM = "
          f"{res['meas_ratio'].mean():.3f}, published 1.30)")
