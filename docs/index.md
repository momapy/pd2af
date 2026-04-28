# pd2af

**pd2af** transforms a [CellDesigner](https://www.celldesigner.org/) process-description (PD) map into an activity-flow (AF) map.

It is built on top of [momapy](https://github.com/adrienrougny/momapy) for map I/O and layout, and uses [clingo](https://potassco.org/clingo/) (via [clorm](https://github.com/potassco/clorm) and [aspcompose](https://github.com/adrienrougny/aspcompose)) to derive the AF model from the PD model with a set of declarative rules.

## Installation

pd2af is a Python package (Python >=3.12). With [uv](https://docs.astral.sh/uv/):

```bash
uv pip install pd2af
```

The `auto` layout mode requires Graphviz's `dot` binary on your `PATH`.

## Usage

### Python API

```python
import momapy.io.core
import pd2af

cd_map = momapy.io.core.read("my_map.xml").obj
af_map = pd2af.transform(cd_map, mode="normal", layout_mode="plain")
momapy.io.core.write(af_map, "my_map_af.xml", writer="celldesigner")
```

### Command-line interface

```bash
pd2af my_map.xml -o my_map_af.xml
pd2af my_map.xml -m pure-af -l auto -o my_map_af.xml
```

See [CLI reference](cli.md) for all options.

## Transformation modes

pd2af supports three transformation modes, selectable with `-m` / `mode=`:

- **`normal`** (default) — emits one activity per distinct PD species (template + state + compartment). Influences are derived from PD reactions and modulations.
- **`no-complex`** — like `normal`, but does not emit complex activities; complex membership is flattened into influences between subunit activities.
- **`pure-af`** — merges all proteoforms of the same template (within a compartment) into a single activity. Requires `--layout auto` because positions from the original PD map can no longer be reused.

## Layout modes

Selectable with `-l` / `layout_mode=`:

- **`plain`** (default) — reuse original positions, only model elements are kept.
- **`overlay`** — reuse the full original layout; non-model elements are greyed out.
- **`auto`** — Graphviz `dot` auto-layout (required for `pure-af`).

## Documentation

- [CLI reference](cli.md)
- [API reference](api_reference/index.md)
