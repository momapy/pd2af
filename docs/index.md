# pd2af

**pd2af** transforms a process-description (PD) map into an activity-flow (AF) map. It reads [CellDesigner](https://www.celldesigner.org/) and [SBGN-PD](https://www.sbgn.org) maps, and writes a CellDesigner map for CellDesigner input and an SBGN-AF map for SBGN-PD input.

It is built on top of [momapy](https://github.com/adrienrougny/momapy) for map I/O and layout, and uses [clingo](https://potassco.org/clingo/) (via [clorm](https://github.com/potassco/clorm) and [aspcompose](https://github.com/adrienrougny/aspcompose)) to derive the AF model from the PD model with a set of declarative rules.

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
pd2af transform my_map.xml -m keep-species -l plain -o my_map_af.xml
pd2af list-modes
pd2af list-groups
```

See [CLI reference](cli.md) for all options.

## Transformation modes

Selectable with `-m` / `mode=`. Four of the five modes lie on two orthogonal axes — how species are mapped to activities, and how complexes are handled:

|                       | keep complexes            | drop complexes (route through subunits) |
|-----------------------|---------------------------|------------------------------------------|
| **merge proteoforms** | `normal` *(default)*      | `normal-no-complex`                             |
| **keep each species** | `keep-species`            | `keep-species-no-complex`                |

- **merge proteoforms** modes (`normal`, `normal-no-complex`) collapse all proteoforms of the same template within the same compartment into a single activity. The result is a true PD→AF transform with no PD remnants — and the only style expressible in SBGN PD, which forbids influences between EPNs. These modes require `--layout-mode dot` (or `auto`).
- **keep each species** modes (`keep-species`, `keep-species-no-complex`) emit one activity per distinct PD species (template + state + compartment), which is only meaningful for CellDesigner.
- **drop complexes** variants (`normal-no-complex`, `keep-species-no-complex`) drop a complex when one of its subunits is independently active, routing influences through the active subunits.
- **keep complexes** variants (`normal`, `keep-species`) emit complexes as their own activities, and influences involving an active complex go through the complex. Active subunits of an activity-bearing complex are subsumed into the complex and do not appear as separate top-level activities.

One further mode steps off those axes and changes where the influences come from:

- **`keep-reactions`** keeps the PD topology itself. Every species is an activity — no structural signal required — and every reaction becomes one positive influence per (reactant, product) pair. Modulation arcs and reaction modifiers are kept as single-hop influences with their own kinds, and none of the inference the other modes do (multi-hop chaining across reactions, catalyst-consumes-reactant / inhibitor-spares-reactant) is applied. Complexes and PTM decorations are kept, as in `keep-species`. CellDesigner-only.

## Layout modes

Selectable with `-l` / `layout_mode=`:

- **`auto`** (default): pick automatically from the input — a map gets `dot`, a bare model gets `None` (no layout).
- **`dot`**: Graphviz `dot` auto-layout. Required for `normal` and `normal-no-complex`.
- **`plain`**: reuse original positions; only model elements are kept. Available for `keep-species`, `keep-species-no-complex` and `keep-reactions`.
- **`overlay`**: reuse the full original layout; non-model elements are greyed out. Available for `keep-species`, `keep-species-no-complex` and `keep-reactions`.

## Documentation

- [CLI reference](cli.md)
- [API reference](api_reference/index.md)
