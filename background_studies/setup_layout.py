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
    "cosmic_muons.py",
    "muon_transport.py",
    "nu_background.py",
    "nu_flux.py",
    "nu_xsec.py",
    "plots.py",
    "test_nu_background.py",
]

STAY_PUT = ["01_validation.ipynb", "02_results.ipynb", "03_allsky.ipynb",
            "README.md", "pytest.ini", "convert_flux_curve.py", "convert_xsec.py"]

# other data files: moved into nubkg/data/ (or data/raw/) if they arrive flat
EXTRA_DATA = ["xsec_cc_digitised.csv", "xsec_icecube_measured.csv"]
RAW_DATA = ["flux_curve.csv", "pdg_cross_sections.csv", "ice_cube_cross_section.csv"]

DATA_NAME = "atmospheric_numu_digitised.csv"

DATA_CONTENT = """\
# Frejus nu_mu model (0.12-11 GeV) + Honda H3a+ERS (125 GeV-1 PeV), digitised by MB
# converted from flux_curve.csv by convert_flux_curve.py
# y_raw * 1e-9 = E^2*Phi  [GeV cm^-2 s^-1 sr^-1]
# uncertainty attached as a flat 25%
# columns: E_nu[GeV], E^2*Phi, sigma(E^2*Phi)
1.162030e-01,2.818383e-02,7.045957e-03
1.569106e-01,3.548134e-02,8.870335e-03
2.462092e-01,4.136820e-02,1.034205e-02
3.767792e-01,4.216965e-02,1.054241e-02
5.484417e-01,4.058199e-02,1.014550e-02
7.785820e-01,3.686945e-02,9.217363e-03
1.162030e+00,3.102178e-02,7.755444e-03
1.569106e+00,2.610157e-02,6.525393e-03
2.118785e+00,2.326305e-02,5.815763e-03
2.790308e+00,1.847850e-02,4.619624e-03
3.674662e+00,1.525223e-02,3.813057e-03
4.489251e+00,1.258925e-02,3.147314e-03
5.912064e+00,1.039122e-02,2.597806e-03
7.593374e+00,8.413951e-03,2.103488e-03
9.276653e+00,6.683439e-03,1.670860e-03
1.133308e+01,5.623413e-03,1.405853e-03
1.000000e+02,5.623413e-04,1.405853e-04
1.252639e+02,4.553374e-04,1.138344e-04
1.455605e+02,3.686945e-04,9.217363e-05
1.869559e+02,3.043220e-04,7.608050e-05
2.283998e+02,2.371374e-04,5.928434e-05
2.933535e+02,1.812731e-04,4.531827e-05
3.674662e+02,1.412538e-04,3.531344e-05
4.603027e+02,1.100694e-04,2.751735e-05
5.765933e+02,8.413951e-05,2.103488e-05
7.044109e+02,6.683439e-05,1.670860e-05
8.823729e+02,5.207948e-05,1.301987e-05
1.051330e+03,3.981072e-05,9.952679e-06
1.316938e+03,3.043220e-05,7.608050e-06
1.649648e+03,2.417315e-05,6.043289e-06
1.823348e+03,2.073322e-05,5.183304e-06
2.118785e+03,1.711328e-05,4.278321e-06
2.588472e+03,1.359356e-05,3.398391e-06
3.084114e+03,1.059254e-05,2.648134e-06
3.583834e+03,8.254042e-06,2.063510e-06
4.270068e+03,6.431812e-06,1.607953e-06
4.961948e+03,5.207948e-06,1.301987e-06
6.061899e+03,3.981072e-06,9.952679e-07
7.222635e+03,3.043220e-06,7.608050e-07
8.605629e+03,2.282093e-06,5.705232e-07
1.025344e+04,1.778279e-06,4.445699e-07
1.191481e+04,1.385692e-06,3.464230e-07
1.492496e+04,1.079775e-06,2.699438e-07
1.608873e+04,8.576959e-07,2.144240e-07
1.916941e+04,6.812921e-07,1.703230e-07
2.227543e+04,5.516539e-07,1.379135e-07
2.588472e+04,4.298662e-07,1.074666e-07
3.007883e+04,3.414549e-07,8.536372e-08
3.408856e+04,2.610157e-07,6.525393e-08
4.270068e+04,1.995262e-07,4.988156e-08
4.961948e+04,1.496236e-07,3.740589e-08
6.215531e+04,1.000000e-07,2.500000e-08
7.593374e+04,6.812921e-08,1.703230e-08
9.511760e+04,5.011872e-08,1.252968e-08
1.105295e+05,3.686945e-08,9.217363e-09
1.350314e+05,2.464147e-08,6.160368e-09
1.691457e+05,1.711328e-08,4.278321e-09
1.965524e+05,1.234999e-08,3.087498e-09
2.341883e+05,9.261187e-09,2.315297e-09
2.721339e+05,7.079458e-09,1.769864e-09
3.162278e+05,5.516539e-09,1.379135e-09
3.674662e+05,4.216965e-09,1.054241e-09
4.270068e+05,3.043220e-09,7.608050e-10
5.216645e+05,2.113489e-09,5.283723e-10
6.061899e+05,1.525223e-09,3.813057e-10
7.044109e+05,1.165914e-09,2.914786e-10
7.983143e+05,8.912509e-10,2.228127e-10
9.511760e+05,6.556418e-10,1.639105e-10
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

    (data / "raw").mkdir(exist_ok=True)
    for names, dest in ((EXTRA_DATA, data), (RAW_DATA, data / "raw")):
        for name in names:
            src, dst = here / name, dest / name
            if dst.exists():
                continue
            if src.exists():
                shutil.move(str(src), str(dst))
                print(f"{name}: moved into {dst.parent.relative_to(here)}/")
            else:
                print(f"{name}: MISSING" + ("  <-- the digitised cross section; "
                      "run convert_xsec.py or re-download" if name == EXTRA_DATA[0] else ""))

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
