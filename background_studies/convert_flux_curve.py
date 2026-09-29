#!/usr/bin/env python3
"""
convert_flux_curve.py -- convert a digitised flux curve into nubkg units.

Input  : two columns, log10(E/GeV) and a y value in unknown units.
Output : nubkg's TabulatedFlux CSV format, E[GeV], E^2*Phi, err.

The y-scale is determined by fitting a single free power of ten against an
independent model (the Chirkin CORSIKA parameterisation), rather than assumed.
A correct scale shows up as a *flat* ratio across many decades; a wrong power
of ten shows up as a ratio of 10 or 0.1.  See `identify_scale()`.

Usage:
    python convert_flux_curve.py raw.csv -o nubkg/data/converted.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nubkg as nb


def load_raw(path):
    d = np.genfromtxt(path, delimiter=",", comments="#")
    d = d[~np.isnan(d).any(axis=1)]
    d = d[np.argsort(d[:, 0])]
    return d[:, 0], d[:, 1]


def identify_scale(log_e, y, fit_range=(2.2, 4.8)):
    """Find the power of ten relating y to E^2*Phi in GeV cm^-2 s^-1 sr^-1.

    Compares against the zenith-averaged Chirkin model over a range where that
    model is well validated (600 GeV - 60 TeV, extended slightly), takes the
    median log offset, and rounds to the nearest integer power of ten.

    Returns (exponent, scatter_dex).  A small scatter means the shapes agree
    and the scale factor is real rather than a fudge absorbing a shape error.
    """
    m = (log_e >= fit_range[0]) & (log_e <= fit_range[1])
    e = 10.0 ** log_e[m]
    model = sum(nb.ChirkinAtmospheric().zenith_averaged(e, s, "up")
                for s in ("nu", "nubar")) * e**2
    offset = np.log10(y[m]) - np.log10(model)
    return int(round(np.median(offset))), float(np.std(offset))


def convert(path, frac_err=0.25, scale_exponent=None):
    log_e, y = load_raw(path)
    exp, scatter = identify_scale(log_e, y)
    if scale_exponent is not None:
        exp = scale_exponent
    e = 10.0 ** log_e
    e2phi = y * 10.0 ** (-exp)
    return e, e2phi, e2phi * frac_err, exp, scatter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--frac-err", type=float, default=0.25,
                    help="fractional uncertainty to attach (default 0.25)")
    ap.add_argument("--scale-exponent", type=int, default=None,
                    help="override the auto-detected power of ten")
    ap.add_argument("--label", default="digitised flux curve")
    a = ap.parse_args()

    e, e2phi, err, exp, scatter = convert(a.input, a.frac_err, a.scale_exponent)
    print(f"detected scale: E^2*Phi = y * 1e-{exp}   "
          f"(residual scatter {scatter:.3f} dex over the fit range)")
    print(f"energy range:   {e.min():.4g} to {e.max():.4g} GeV")

    with open(a.output, "w") as fh:
        fh.write(f"# {a.label}\n")
        fh.write(f"# converted from {Path(a.input).name} by convert_flux_curve.py\n")
        fh.write(f"# y_raw * 1e-{exp} = E^2*Phi  [GeV cm^-2 s^-1 sr^-1]\n")
        fh.write(f"# uncertainty attached as a flat {100*a.frac_err:.0f}%\n")
        fh.write("# columns: E_nu[GeV], E^2*Phi, sigma(E^2*Phi)\n")
        for a_, b_, c_ in zip(e, e2phi, err):
            fh.write(f"{a_:.6e},{b_:.6e},{c_:.6e}\n")
    print(f"wrote {a.output}  ({len(e)} points)")


if __name__ == "__main__":
    main()
