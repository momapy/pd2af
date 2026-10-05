# pd2af

**pd2af** transforms a process-description (PD) map into an activity-flow (AF) map. It reads [CellDesigner](https://www.celldesigner.org/) and [SBGN-PD](https://www.sbgn.org) maps, and writes a CellDesigner map for CellDesigner input and an SBGN-AF map for SBGN-PD input.

It is built on top of [momapy](https://github.com/momapy/momapy) for map I/O and layout, and uses [clingo](https://potassco.org/clingo/) (via [clorm](https://github.com/potassco/clorm) and [aspcompose](https://github.com/adrienrougny/aspcompose)) to derive the AF model from the PD model with a set of declarative rules.

## Installation

pd2af is a Python package (Python >=3.12). With [uv](https://docs.astral.sh/uv/):

```bash
uv pip install pd2af
```

The `dot` layout mode requires Graphviz's `dot` binary on your `PATH`.

## Usage

### Python API

```python
import momapy.io.core
import pd2af

cd_map = momapy.io.core.read("my_map.xml").obj
af_map = pd2af.transform(cd_map, mode="normal", layout_mode="auto").obj
momapy.io.core.write(af_map, "my_map_af.xml", writer="celldesigner")
```

`transform` also accepts a bare model instead of a full map: pass
`cd_map.model` and it returns the transformed model (an SBGN-AF model for
SBGN-PD input, a CellDesigner model for CellDesigner input) via `.obj`. Model
input has no geometry, so `layout_mode` is forced to `None`.

### Command-line interface

```bash
pd2af transform my_map.xml -o my_map_af.xml
pd2af transform my_map.xml --keep-species -l plain -o my_map_af.xml
pd2af list-modes
pd2af list-groups
```

See [CLI reference](cli.md) for all options.

## Transformation modes and options

A mode decides what counts as an activity; an option decides which activities are treated as the same thing.

### Modes

Selectable with `-m` / `mode=`:

- **`normal`** *(default)*: complexes are activities of their own, and a subunit's influences are routed to the complex it belongs to. An active subunit of an activity-bearing complex is subsumed into the complex and does not appear as a separate top-level activity. Influences are inferred beyond the stated modulations: multi-hop chaining across reactions, catalyst-consumes-reactant / inhibitor-spares-reactant, and activation by binding (in `L + R -> L:R` with `R` drawn active only inside the complex, `L` is an activity that positively influences `L:R`).
- **`no-complex`**: a complex with an active subunit is replaced by those subunits, promoted to top-level activities, and the influences run through them; the same inference applies, so activation by binding gives `L -> R`.
- **`keep-reactions`**: the PD topology itself is kept. Every species is an activity — no structural signal required — and every reaction becomes one positive influence per (reactant, product) pair. Modulation arcs and reaction modifiers are kept as single-hop influences with their own kinds, and none of the inference the other modes do (multi-hop chaining across reactions, catalyst-consumes-reactant / inhibitor-spares-reactant, activation by binding) is applied. CellDesigner-only, and usually wanted with `--keep-species`, which it starts from by default.

### Options

Each option is a flag with a matching `--no-` spelling, and a matching `transform` argument. A mode names the options it starts from and never overrules an explicit flag, so `-m keep-reactions --no-keep-species` gives the merged reading of that mode.

#### `--keep-species` / `keep_species=True`

Keep each species as its own activity (template + state + compartment) instead of merging the forms of the same base entity. Off by default: the forms of the same base entity within the same compartment collapse into a single activity, which is a true PD→AF transform with no PD remnants — and the only style expressible in SBGN-AF, which forbids influences between decorated entity pools.

#### `--drop-compartments` / `drop_compartments=True`

Drop the compartments, so that species differing only by compartment become a single activity and the influences that become equal merge in turn. Off by default.

- **CellDesigner output**: a single compartment, the `default` one every CellDesigner map declares, with every species pointing at it.
- **SBGN-AF output**: no compartment at all, SBGN-AF having no default compartment, so every activity ends up with none.
- **Annotations and notes**: they follow the merge, an activity merged from several compartments gathering the metadata of every species that collapsed into it. The metadata of the removed compartments is dropped with them.

```bash
pd2af transform my_map.xml --keep-species --drop-compartments -o my_map_af.xml
```

### Finding the new spelling

| before | now |
| --- | --- |
| `pd2af transform map.xml` | `pd2af transform map.xml` |
| `-m normal` | `-m normal` |
| `-m normal-no-complex` | `-m no-complex` |
| `-m keep-species` | `--keep-species` |
| `-m keep-species-no-complex` | `-m no-complex --keep-species` |
| `-m keep-reactions` | `-m keep-reactions` |
| `--no-compartment` | `--drop-compartments` |

## Layout modes

Selectable with `-l` / `layout_mode=`:

- **`auto`** (default): pick automatically from the input — a map gets `dot`, a bare model gets `None` (no layout).
- **`dot`**: Graphviz `dot` auto-layout. Required whenever activities merge — a merged activity has no single original position to reuse — so unless `--keep-species` is set, and always under `--drop-compartments`.
- **`plain`**: reuse original positions; only model elements are kept. Available exactly when `--keep-species` is set without `--drop-compartments`.
- **`overlay`**: reuse the full original layout; non-model elements are greyed out. Available exactly when `--keep-species` is set without `--drop-compartments`, and for CellDesigner output only.

## Documentation

- [CLI reference](cli.md)
- [API reference](api_reference/index.md)
