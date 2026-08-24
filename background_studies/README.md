# EarthShine neutrino-induced muon background

Parameterised (no per-muon Monte Carlo) estimate of muons produced by
neutrino-nucleon interactions in the rock around an underground detector.

## Layout

```
nubkg/
  nu_flux.py            atmospheric + astrophysical fluxes; tabulated & analytic
  nu_xsec.py            CC cross sections (table-driven) and inelasticity
  muon_transport.py     dE/dX = a + bE: range, energy_after, energy_before
  nu_background.py      rock yield, arriving spectrum, detector rates, top API
  plots.py              diagnostic and paper figures
  test_nu_background.py 43 unit tests (pytest --pyargs nubkg.test_nu_background)
  conftest.py           registers the `literature` marker
  data/                 flux/xsec tables -- REPLACE THE PLACEHOLDER
01_validation.ipynb     ingredient-by-ingredient checks vs literature
02_results.ipynb        rates, spectra, error budget, paper figures
```

## Getting it to import

Keep the two notebooks in the same folder as the `nubkg/` directory:

```
your_project/
  nubkg/
  01_validation.ipynb
  02_results.ipynb
```

The first cell of each notebook walks up from the working directory looking for
`nubkg/__init__.py` and prepends what it finds to `sys.path`, so the notebooks
also work from a subfolder. It prints the package root it found.

If you would rather install it once and stop thinking about paths:

```bash
pip install -e .        # with a pyproject.toml, or:
export PYTHONPATH=/path/to/your_project:$PYTHONPATH
```

Data files are located relative to the package, not the working directory:

```python
nb.data_path("ic59_atmospheric_numu_APPROX.csv")
```

so drop your own digitised CSVs into `nubkg/data/` and refer to them by name.

## Top-level call

```python
import nubkg as nb

res = nb.background_muons(
    detector   = nb.Cylinder(radius=7.5, half_length=10.5, depth=100.0),
    e_mu_min   = 10.0,            # GeV, muon energy AT the detector
    e_nu_range = (10.0, 1.0e6),   # GeV
    livetime_s = nb.SEC_PER_YEAR,
    hemisphere = "up",            # the EarthShine signal region
)
res["n_muons"], res["n_muons_err"], res["error_budget"]
res["e_mu"], res["dN_dEmu_per_s"]      # arriving spectrum
res["metadata"]                        # flat dict for a parquet footer
```

## Swapping in your own digitised curves

```python
flux = nb.TabulatedFlux.from_csv(nb.data_path("mypoints.csv"),
                                 flux_unit="E2Phi",
                                 energy_unit="log10GeV",
                                 zenith_shape_from=nb.ChirkinAtmospheric())
xsec = nb.CrossSection.from_csv(nb.data_path("myxsec.csv"),
                                sigma_unit="cm2_per_GeV")
```

`TabulatedFlux` defaults to `outside="model"`, which splices onto the analytic
model outside the tabulated range rather than log-log extrapolating. Read the
docstring before changing that -- see the caveat below.

## Known issues / to do

1. **The low-energy end of the flux table dominates and is not covered by the
   IceCube unfolding.** The rate integrand peaks at E_nu ~ 100-300 GeV for a
   10 GeV muon threshold, below where the unfolded points start. Extend the
   table downward (Frejus points from the same figure, or Honda/MCEq) before
   quoting an absolute number. Blind extrapolation of the table's steep
   high-energy slope inflates the answer by ~2x.
2. The default cross-section table is a bridge between the PDG world average
   and the 10 TeV NLO QCD anchor. Replace with a digitisation of
   arXiv:1305.7513 Fig. 9 and arXiv:1711.08119 Fig. 1.
3. Constant b in dE/dX is a compromise; the real b rises with energy. Carried
   as a systematic via `EnergyLoss.bracket()`.
4. The down-going/up-going ratio is non-monotonic in threshold above ~1 TeV
   and is not fully explained. See the docstring of
   `test_downgoing_over_upgoing_ratio_falls_with_threshold`. Does not affect
   the up-going result.
5. No oscillations, no NC regeneration, no prompt/charm flux component.

## References

* IceCube atmospheric nu_mu unfolding -- https://icecube.wisc.edu/news/research/2014/09/an-improved-measurement-of-atmospheric-neutrino-flux-in-icecube/ , arXiv:1409.4535
* Chirkin CORSIKA fit -- https://arxiv.org/abs/hep-ph/0407078
* Formaggio & Zeller cross-section review -- https://arxiv.org/abs/1305.7513
* IceCube cross-section measurement -- https://arxiv.org/abs/1711.08119
* Cooper-Sarkar, Mertsch & Sarkar (CSMS) -- https://arxiv.org/abs/1106.3723
* Gandhi, Quigg, Reno & Sarcevic -- https://arxiv.org/abs/hep-ph/9807264
* PDG neutrino cross sections -- https://pdg.lbl.gov/2025/reviews/rpp2025-rev-nu-cross-sections.pdf
* PDG muon dE/dx tables -- https://pdg.lbl.gov/2024/AtomicNuclearProperties/
* IceCube diffuse astrophysical flux -- https://arxiv.org/abs/2111.10299
* EarthShine signal paper (Feng, Smolinsky & Tanedo) -- https://arxiv.org/abs/1509.07525
