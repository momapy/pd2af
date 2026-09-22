# pd2af

[![PyPI](https://img.shields.io/pypi/v/pd2af)](https://pypi.org/project/pd2af/)
[![Python](https://img.shields.io/pypi/pyversions/pd2af)](https://pypi.org/project/pd2af/)
[![License](https://img.shields.io/github/license/adrienrougny/pd2af)](https://github.com/adrienrougny/pd2af/blob/main/COPYING)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen)](https://adrienrougny.github.io/pd2af/)

**pd2af** transforms a process-description (PD) map into an activity-flow (AF) map.
It reads [CellDesigner](https://www.celldesigner.org/) and [SBGN-PD](https://www.sbgn.org) maps, and writes a CellDesigner map for CellDesigner input and an SBGN-AF map for SBGN-PD input.
It is built on top of [momapy](https://github.com/adrienrougny/momapy) for map I/O and layout, and uses [clingo](https://potassco.org/clingo/) (via [clorm](https://github.com/potassco/clorm) and [aspcompose](https://github.com/adrienrougny/aspcompose)) to derive the AF model from the PD model with a set of declarative rules.

Features of pd2af include the following:

* CellDesigner and SBGN-PD input, with the output language following the input
* five transformation modes, from a true PD→AF reduction to a CellDesigner-native one that keeps the PD topology
* an option to drop the compartments, merging the species and influences that become equal
* activities discovered from the structure of the map: phenotypes, explicit active markers, modulation sources and logical-gate inputs
* influences inferred beyond the stated modulations: multi-hop chaining across reactions, catalyst-consumes-reactant and inhibitor-spares-reactant, and activation by binding (a reactant whose partner comes out of a complex-forming reaction active)
* three layout modes: Graphviz `dot` auto-layout, reuse of the original positions, and an overlay on the original layout
* RDF/MIRIAM annotations and notes carried from the input elements to the elements they produce
* a transformation defined by declarative ASP rules, listable and individually excludable
* transformation modes contributable by third-party packages through an entry point

## Installation

pd2af is available as a Python package and can be installed with pip as follows (Python >=3.12):

`pip install pd2af`

The `dot` layout mode requires Graphviz's `dot` binary on your `PATH`.

## Usage

### Python API

```python
import momapy.io.core
import pd2af

pd_map = momapy.io.core.read("my_map.xml").obj
af_map = pd2af.transform(pd_map, mode="normal", layout_mode="auto").obj
momapy.io.core.write(af_map, "my_map_af.xml", writer="celldesigner")
```

`transform` also accepts a bare model instead of a full map.

### Command-line interface

```bash
pd2af transform my_map.xml -o my_map_af.xml
pd2af transform my_map.xml --keep-species -l plain -o my_map_af.xml
pd2af list-modes
pd2af list-groups
```

## Transformation modes and options

A mode decides what counts as an activity; an option decides which activities are treated as the same thing.

### Modes

Selectable with `-m` / `mode=`:

* **`normal`** *(default)*: keep complexes as activities of their own, routing a subunit's influences to the complex it belongs to.
* **`no-complex`**: replace a complex that has an active subunit with those subunits, promoting them to top-level activities.
* **`keep-reactions`**: keep the PD topology itself: every species is an activity, and every reaction becomes one positive influence per (reactant, product) pair. CellDesigner-only, and usually wanted with `--keep-species`, which it turns on by default.

### Options

Each option is a flag with a matching `--no-` spelling, and a matching `transform` argument:

* **`--keep-species`** / `keep_species=True`: keep each species as its own activity (template + state + compartment) instead of merging the forms of the same base entity. Off by default, on for `keep-reactions`.
* **`--drop-compartments`** / `drop_compartments=True`: drop the compartments, so that species differing only by compartment become a single activity and the influences that become equal merge in turn. The output holds one compartment for CellDesigner (the `default` one every CellDesigner map declares) and none for SBGN-AF. Off by default.

```bash
pd2af transform my_map.xml --keep-species --drop-compartments -o my_map_af.xml
```

A mode names the options it starts from and never overrules an explicit flag, so `-m keep-reactions --no-keep-species` gives the merged reading of that mode.

| before | now |
| --- | --- |
| `pd2af transform map.xml` | `pd2af transform map.xml` |
| `-m normal` | `-m normal` |
| `-m normal-no-complex` | `-m no-complex` |
| `-m keep-species` | `--keep-species` |
| `-m keep-species-no-complex` | `-m no-complex --keep-species` |
| `-m keep-reactions` | `-m keep-reactions` |
| `--no-compartment` | `--drop-compartments` |

### Layout modes

`dot` is required whenever activities merge: unless `--keep-species` is set, and always under `--drop-compartments`.

## Documentation

The documentation for pd2af is available [here](https://adrienrougny.github.io/pd2af/).
