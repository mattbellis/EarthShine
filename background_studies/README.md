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
  cosmic_muons.py       CR muon flux at depth (order-of-magnitude, for the
                        all-sky case only -- NOT part of the nu calculation)
  plots.py              diagnostic and paper figures
  test_nu_background.py 61 unit tests (pytest --pyargs nubkg.test_nu_background)
  conftest.py           registers the `literature` marker
  data/
    atmospheric_numu_digitised.csv      flux (MB's digitisation)
    xsec_cc_digitised.csv               CC cross sections, nu and nubar
    xsec_icecube_measured.csv           IceCube measured/SM, for systematics
    ic59_atmospheric_numu_APPROX.csv    superseded eyeball placeholder
    raw/                                the unconverted digitiser exports
convert_flux_curve.py   digitised flux curve -> nubkg units (scale auto-detect)
convert_xsec.py         PDG + IceCube digitisations -> xsec_cc_digitised.csv
01_validation.ipynb     ingredient-by-ingredient checks vs literature
02_results.ipynb        rates, spectra, error budget, paper figures
03_allsky.ipynb         all-sky background for floating-DM models, and the
                        cosmic-ray muon comparison that dominates it
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

## Choosing the input tables

Each notebook's first cell names both inputs, and nothing else does:

```python
FLUX_FILE = "atmospheric_numu_digitised.csv"
XSEC_FILE = "xsec_cc_digitised.csv"
nb.set_default_xsec(XSEC_FILE)
```

`set_default_xsec` makes that table the default for every calculation, and
fails immediately if the file is missing. `nb.set_default_xsec("builtin")`
restores the old hand-built table for comparison -- but note it is 2x too low
at 100 TeV and 3x too low at 1 PeV.

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

## Current baseline

Using `nubkg/data/atmospheric_numu_digitised.csv` (Frejus nu_mu model below
11 GeV + Honda H3a+ERS above 125 GeV, digitised off the standard atmospheric
flux figure, with the Chirkin zenith shape grafted on):

    N(E_mu > 10 GeV, up-going, CMS-like cylinder, 1 yr) = 98 +- 29

falling to ~3/yr at a 1 TeV threshold. Error budget: flux 24, inelasticity 11,
Poisson 10, cross section 8, energy loss 4.

All-sky (floating-DM models, notebook 03): N = 175 +- 52. But see below --
that number is not the limiting background going down.

## For all-sky / floating dark matter searches

Down-going, the competition is not neutrino-induced muons but ordinary
cosmic-ray muons, which at a 100 m overburden outnumber them by ~4 x 10^7.
They cross only within ~1.1 degrees of the horizon. Practical consequence:
treat the whole down-going hemisphere as inaccessible without a cosmic-ray
veto, and use the up-going numbers for any signal region at or above the
horizon. `nubkg.cosmic_muons` is an order-of-magnitude tool for establishing
this, not for quoting a cosmic-ray rate -- see its docstring.

## Known issues / to do

1. ~~The low-energy end of the flux table dominates and is not covered.~~
   RESOLVED, and the original diagnosis was partly wrong. With the real
   digitisation in hand: neutrinos below 11 GeV contribute **0.0%** of the rate
   at a 10 GeV muon threshold (inelasticity means E_mu > 10 GeV needs
   E_nu >~ 20-30 GeV, and the muon still has to travel). What actually carries
   the rate is 11-125 GeV -- 39% of it -- which is the *gap* between the two
   digitised segments. But bridging that gap with a log-log line versus with
   the Chirkin shape changes the answer by only 2.8%, so the interpolation is
   well constrained by its endpoints. No further digitisation needed.
2. ~~The default cross-section table is a hand-built bridge.~~ RESOLVED:
   replaced by MB's digitisation (PDG below 340 GeV, CSMS via IceCube Fig. 1
   above 600 GeV, blended between). The old table was 2-3x too low above
   100 TeV; rates change by only +1% at a 10 GeV threshold but +30% at 10 TeV.
   The PDG and CSMS sources disagree by -5.6% (nu) / +7.9% (nubar) at 350 GeV,
   carried as an 8% uncertainty across the blend.
3. Constant b in dE/dX is a compromise; the real b rises with energy. Carried
   as a systematic via `EnergyLoss.bracket()`.
4. The down-going/up-going ratio is non-monotonic in threshold above ~1 TeV
   and is not fully explained. See the docstring of
   `test_downgoing_over_upgoing_ratio_falls_with_threshold`. Does not affect
   the up-going result.
5. No oscillations, no NC regeneration. (The digitised flux DOES include the
   prompt/charm component, via the ERS term in Honda H3a+ERS.)
6. The zenith shape barely affects the *integrated* rate -- the graft preserves
   the hemisphere average, and the cylinder's projected area varies only 1.12x
   with zenith, so the horizon enhancement largely cancels. It matters for the
   *differential* dN/dcos(theta), which is what a directional search needs.
7. `ChirkinAtmospheric` was fitted over 600 GeV - 60 TeV
   (`ChirkinAtmospheric.VALID_RANGE`), and the plots shade that band. It is a
   *fit to CORSIKA output*, not something CORSIKA uses. Extrapolating down to
   10 GeV costs ~20%; extrapolating up to 1 PeV costs ~60%, because the form
   contains no cosmic-ray knee. This does not propagate into the rates -- the
   model is only the zenith-shape donor there, worth ~1%.
8. `cosmic_muons.py` assumes a flat overburden and is used below the validity
   floor of both Gaisser (energy) and Mei & Hime (depth); the two bracket each
   other to within an order of magnitude at 0.265 km.w.e.

## References

* IceCube atmospheric nu_mu unfolding (IC-59; "59" = the 59-string detector
  configuration, 2009-10) -- https://icecube.wisc.edu/news/research/2014/09/an-improved-measurement-of-atmospheric-neutrino-flux-in-icecube/ , arXiv:1409.4535
* Honda atmospheric flux calculation -- https://arxiv.org/abs/astro-ph/0611418
* Enberg, Reno & Sarcevic (ERS, the prompt/charm component) -- https://arxiv.org/abs/0806.0418
* Chirkin CORSIKA fit -- https://arxiv.org/abs/hep-ph/0407078
* Formaggio & Zeller cross-section review -- https://arxiv.org/abs/1305.7513
* IceCube cross-section measurement -- https://arxiv.org/abs/1711.08119
* Cooper-Sarkar, Mertsch & Sarkar (CSMS) -- https://arxiv.org/abs/1106.3723
* Gandhi, Quigg, Reno & Sarcevic -- https://arxiv.org/abs/hep-ph/9807264
* PDG neutrino cross sections -- https://pdg.lbl.gov/2025/reviews/rpp2025-rev-nu-cross-sections.pdf
* PDG muon dE/dx tables -- https://pdg.lbl.gov/2024/AtomicNuclearProperties/
* IceCube diffuse astrophysical flux -- https://arxiv.org/abs/2111.10299
* EarthShine signal paper (Feng, Smolinsky & Tanedo) -- https://arxiv.org/abs/1509.07525
* Super-Kamiokande up-going through-going muons (the absolute validation anchor)
  -- Fukuda et al., PRL 82, 2644 (1999), https://arxiv.org/abs/hep-ex/9812014
* MACRO up-going muons (ratio only, not used as an absolute anchor)
  -- Ambrosio et al., PLB 434, 451 (1998), https://arxiv.org/abs/hep-ex/9809003
* Mei & Hime, muon depth-intensity relation -- https://arxiv.org/abs/astro-ph/0512125
* Woodley et al., cosmic ray muons deep underground -- https://arxiv.org/abs/2406.10339
