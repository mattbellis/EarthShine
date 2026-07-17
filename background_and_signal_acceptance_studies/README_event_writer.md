[Link to Claude/Fable chat](https://claude.ai/share/6f462d4b-217d-40d2-924b-54e9666270c6)

# EarthShine → LHE / HepMC converter

`lhe_writer.py` converts EarthShine muon-pair dataframes (one `A' → μ⁺μ⁻`
decay per row) into event files that can be fed into CMS' standard
simulation chain. It has no dependencies beyond `numpy` and `pandas`.

## Input

A dataframe (or parquet file) with the columns

```
px1 py1 pz1 e1   px2 py2 pz2 e2   x0 y0 z0
```

Momenta and energies are in **GeV**. `(x0, y0, z0)` is the dark photon
decay point (= muon production point) in **metres**, relative to the
centre of CMS.

## Output formats

| Format   | Function / flag       | Vertex handling | When to use |
|----------|-----------------------|-----------------|-------------|
| HepMC2 (`IO_GenEvent`) | `write_hepmc(..., version=2)` / `--format hepmc2` | **Native.** Each event's vertex carries the full (x, y, z, t) position. | **Default.** Read directly by CMSSW's `MCFileSource`. |
| HepMC3 (`Asciiv3`)     | `write_hepmc(..., version=3)` / `--format hepmc3` | Native, same as above. | HepMC3-enabled CMSSW workflows or standalone HepMC3/Rivet tooling. The CMSSW HepMC3 migration is still in progress, so check your release before relying on this. |
| LHE      | `write_lhe(...)` / `--format lhe` | **Not native.** LHE has no vertex field; the position is written as a `#vertex x y z` comment line per event (parsers ignore it; a small custom CMSSW plugin is needed to apply it). A `ctau` mode also exists but is geometrically wrong for EarthShine events — see docstring. | Only if a downstream tool insists on LHE. |

**Short version: use HepMC2.** The whole reason this module grew a HepMC
writer is that LHE cannot express per-event production vertices, and the
displaced vertex *is* the physics here. With HepMC, GEANT starts the muons
exactly where you say, with no plugins or hacks.

Common conventions across all writers:

- The A' is written as an intermediate particle (PDG **4900022**, the CMS
  hidden-valley convention; configurable via `aprime_pdgid`), with its
  4-momentum and mass computed per event from the muon pair. In HepMC it is
  the incoming particle of the decay vertex and is **mandatory** — a vertex
  with no incoming particle does not survive the ASCII round trip in either
  HepMC version (verified against the reference library).
- Leg 1 → μ⁻ (13), leg 2 → μ⁺ (−13). Flip `muon_pdgids` if dgt uses the
  opposite convention.
- `vertex_scale` defaults to **1000.0** (metres → mm, the HepMC/LHE length
  unit). Set to 1.0 if your vertices are already in mm.
- HepMC vertex time is written as t = 0 ("in time" at the decay point);
  override with `vertex_time` (units of c·t, mm) if you need a specific
  offset relative to the bunch crossing.
- Optional per-event weights via `weight_col`.

## Usage from Python

```python
import earthshine_io as eio
from lhe_writer import write_hepmc, write_lhe

df, params_list = eio.load_many(selection)

# Recommended: HepMC2 for CMSSW's MCFileSource
write_hepmc(df, "events.hepmc")                    # version=2 is the default

# HepMC3 variant
write_hepmc(df, "events_v3.hepmc", version=3)

# LHE, reusing the header/<init> block from a known-good CMS file
write_lhe(df, "events.lhe", template_path="cms_reference.lhe")
```

## Usage from the command line

The module has a `main()` entry point guarded by
`if __name__ == "__main__":`, so there are three equivalent ways to run it.

**1. Run the file directly** (works from anywhere, given the path):

```bash
python lhe_writer.py events.parquet events.hepmc                 # HepMC2 (default)
python lhe_writer.py events.parquet events_v3.hepmc --format hepmc3
python lhe_writer.py events.parquet events.lhe --format lhe --template cms_ref.lhe
python lhe_writer.py --help                                      # all options
```

**2. Run it as a module** with `-m` (the thing you were half-remembering).
This works when the file is in the current directory or on `PYTHONPATH`,
and is what you want once it lives inside a package:

```bash
python -m lhe_writer events.parquet events.hepmc --format hepmc2
```

`-m` tells Python to locate the module by *import name* (no `.py`, no
path) and execute it as a script, which triggers the same
`if __name__ == "__main__":` block. If the file later moves into a
package, e.g. `earthshine/lhe_writer.py`, the invocation becomes
`python -m earthshine.lhe_writer ...`.

**3. Call any function directly with `-c`** — useful for options the CLI
doesn't expose, or for quick one-offs without writing a script:

```bash
python -c "
import pandas as pd
from lhe_writer import write_hepmc
df = pd.read_parquet('events.parquet')
write_hepmc(df, 'events.hepmc', muon_pdgids=(-13, 13), vertex_time=25.0)
"
```

### CLI options

```
positional:   input (parquet)   output
--format {hepmc2,hepmc3,lhe}   output format (default: hepmc2)
--vertex-scale FLOAT           multiply (x0,y0,z0); default 1000.0 (m -> mm)
--vertex-time FLOAT            [hepmc] c*t (mm) written for every vertex (default 0.0)
--weight-col NAME              dataframe column to use as the event weight
--aprime-id INT                PDG ID for the A' (default 4900022)
--template PATH                [lhe] CMS LHE file to reuse the header/<init> from
--vertex-mode {comment,ctau,none}   [lhe] how to record the vertex (default: comment)
--initial / --initial-ids ID1 ID2   [lhe] add mock status -1 incoming particles
--no-mother                    [lhe] omit the A' intermediate particle
```

## Feeding CMSSW

```python
process.source = cms.Source("MCFileSource",
    fileNames = cms.untracked.vstring('file:events.hepmc'))
```

Make sure the **`VtxSmeared` step is not applied** on top of the embedded
vertices — beam-spot smearing would shift the carefully placed decay
points. For `MCFileSource` input, run GEN-SIM with the smearing producer
removed or replaced by a pass-through.

## Validation

The writers were validated by round-tripping their output through the
reference HepMC3 library (`pyhepmc`) and, for LHE, through `pylhe`:

- vertex positions and muon momenta reproduced with **zero** error;
- correct event graph (A' status 2 incoming to the vertex, muons status 1
  outgoing; correct mother pointers in LHE);
- weights preserved when `weight_col` is set;
- for LHE with `--initial`, exact momentum conservation between the mock
  incoming pair (light-cone decomposition of the total 4-momentum) and the
  final-state muons.

A quick re-check after any edit:

```bash
pip install pyhepmc
python -c "
import pyhepmc
with pyhepmc.open('events.hepmc', format='HepMC2') as f:
    ev = next(iter(f))
print(ev.vertices[0].position, [p.pid for p in ev.particles])
"
```
