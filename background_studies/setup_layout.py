#!/usr/bin/env python3
"""
setup_layout.py -- rebuild the nubkg package layout from a flat download.

The files were shared individually, so they all landed in one directory
instead of the package structure the notebooks expect.  Run this once, from
inside that directory:

    python setup_layout.py

It creates nubkg/ and nubkg/data/, moves the module files in, writes the
placeholder flux table if it is missing, and verifies that the package
imports.  Safe to re-run: it skips anything already in place.

Result:

    .
    |- nubkg/
    |   |- __init__.py  nu_flux.py  nu_xsec.py  muon_transport.py
    |   |- nu_background.py  plots.py  test_nu_background.py  conftest.py
    |   `- data/ic59_atmospheric_numu_APPROX.csv
    |- 01_validation.ipynb
    |- 02_results.ipynb
    |- README.md
    `- pytest.ini
"""

import shutil
import subprocess
import sys
from pathlib import Path

MODULES = [
    "__init__.py",
    "conftest.py",
    "muon_transport.py",
    "nu_background.py",
    "nu_flux.py",
    "nu_xsec.py",
    "plots.py",
    "test_nu_background.py",
]

STAY_PUT = ["01_validation.ipynb", "02_results.ipynb", "README.md", "pytest.ini"]

DATA_NAME = "ic59_atmospheric_numu_APPROX.csv"

DATA_CONTENT = """\
# IC-59 unfolded atmospheric nu_mu + nubar_mu spectrum
#
# !!!! PLACEHOLDER -- REPLACE WITH YOUR OWN DIGITISATION !!!!
# These values were read BY EYE off the standard E^2 Phi figure and are
# accurate to maybe 30%.  They exist so the code runs and so the loader has
# something to be unit-tested against.  Substitute the points you scanned.
#
# Source figure / paper:
#   https://icecube.wisc.edu/news/research/2014/09/an-improved-measurement-of-atmospheric-neutrino-flux-in-icecube/
#   Aartsen et al., Eur. Phys. J. C 75, 116 (2015), arXiv:1409.4535
# Northern (up-going) sky, zenith averaged.  Includes conventional + prompt +
# whatever astrophysical component is present at the top end.
#
# NOTE: the lowest point here is ~350 GeV, but the background rate is dominated
# by E_nu ~ 100-300 GeV.  Extend this table downward (Frejus points from the
# same figure, or Honda/MCEq) before quoting an absolute number.  See README.
#
# columns: E_nu[GeV], E^2*Phi [GeV cm^-2 s^-1 sr^-1], 1-sigma on E^2*Phi
3.5e2,2.6e-4,0.9e-4
6.3e2,9.8e-5,2.6e-5
1.0e3,3.1e-5,0.8e-5
2.5e3,1.0e-5,0.2e-5
5.6e3,3.7e-6,0.7e-6
1.1e4,1.05e-6,0.20e-6
2.2e4,2.8e-7,0.6e-7
5.0e4,7.5e-8,2.0e-8
2.2e5,1.9e-8,1.0e-8
5.6e5,3.8e-9,2.5e-9
"""


def main() -> int:
    here = Path(__file__).resolve().parent
    pkg = here / "nubkg"
    data = pkg / "data"

    if pkg.exists() and (pkg / "__init__.py").exists():
        print(f"nubkg/ already exists at {pkg}")
    pkg.mkdir(exist_ok=True)
    data.mkdir(exist_ok=True)

    moved, already, missing = [], [], []
    for name in MODULES:
        src, dst = here / name, pkg / name
        if dst.exists():
            already.append(name)
            if src.exists() and src != dst:
                print(f"  NOTE: {name} exists in both places; leaving the copy "
                      f"in nubkg/ alone and not overwriting it")
        elif src.exists():
            shutil.move(str(src), str(dst))
            moved.append(name)
        else:
            missing.append(name)

    # the data file may be flat, or absent entirely (it was not shared)
    data_dst = data / DATA_NAME
    flat_data = here / DATA_NAME
    if data_dst.exists():
        data_status = "already in nubkg/data/"
    elif flat_data.exists():
        shutil.move(str(flat_data), str(data_dst))
        data_status = "moved into nubkg/data/"
    else:
        data_dst.write_text(DATA_CONTENT)
        data_status = "was missing -- placeholder written"

    print(f"\nmoved into nubkg/ : {moved or 'nothing'}")
    if already:
        print(f"already in place  : {already}")
    if missing:
        print(f"MISSING           : {missing}   <-- re-download these")
    print(f"{DATA_NAME}: {data_status}")

    for name in STAY_PUT:
        if not (here / name).exists():
            print(f"note: {name} not found here (fine if you did not download it)")

    if missing:
        print("\nCannot verify the import until the missing files are present.")
        return 1

    print("\nverifying...")
    check = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r)\n"
         "import nubkg as nb\n"
         "r = nb.background_muons(e_mu_min=10.0)\n"
         "print('  import OK')\n"
         "print('  N(E_mu>10 GeV, up-going, 1 yr) = %%.1f +/- %%.1f'\n"
         "      %% (r['n_muons'], r['n_muons_err']))\n"
         "print('  data file:', nb.data_path(%r))" % (str(here), DATA_NAME)],
        capture_output=True, text=True)
    print(check.stdout or check.stderr)
    if check.returncode:
        return check.returncode

    print("Done. Launch Jupyter from this directory and open 01_validation.ipynb.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
