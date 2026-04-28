# CLI Reference

## Overview

The `pd2af` command-line interface transforms a CellDesigner process-description (PD) map into an activity-flow (AF) map.

The output map is written to stdout as a [momapy](https://github.com/adrienrougny/momapy) pickle (preserving layout styling) so it can be piped into `momapy visualize`. With `-o`, the writer is chosen from the output file extension.

## Synopsis

```bash
pd2af <input_file> [-m {normal,no-complex,pure-af}] [-l {plain,overlay,auto}] [-o <output_file>]
```

## Arguments

| Argument | Description |
|----------|-------------|
| `input_file` | Input CellDesigner XML file |

## Options

| Option | Short | Description |
|--------|-------|-------------|
| `--mode` | `-m` | Transformation mode (default: `normal`); see below |
| `--layout` | `-l` | Layout strategy (default: `plain`); see below |
| `--output` | `-o` | Write output to this file instead of stdout |

## Transformation modes (`-m`)

| Mode | Description |
|------|-------------|
| `normal` | One activity per distinct PD species (template + state + compartment). |
| `no-complex` | Like `normal`, but complexes are not emitted as activities; complex membership is flattened into influences between subunit activities. |
| `pure-af` | All proteoforms of the same template (within a compartment) are merged into a single activity. Requires `--layout auto`. |

## Layout modes (`-l`)

| Mode | Description |
|------|-------------|
| `plain` | Reuse original positions; only model elements are kept (default). |
| `overlay` | Reuse the full original layout; non-model elements greyed out. |
| `auto` | Graphviz `dot` auto-layout. Requires `dot` on `PATH`. Required for `pure-af`. |

## Output writers

The writer used with `-o` is selected from the file extension:

| Extension | Writer |
|-----------|--------|
| `.xml`, `.sbml` | CellDesigner XML |
| `.pickle`, `.pkl` | momapy pickle |
| (other) | momapy pickle |

When `-o` is omitted, the map is always written to stdout as a momapy pickle.

## Examples

### Basic transformation to stdout (pipe to momapy)

```bash
pd2af my_map.xml | momapy visualize -
```

### Transform to CellDesigner XML

```bash
pd2af my_map.xml -o my_map_af.xml
```

### Pure-AF mode (merge proteoforms, auto layout)

```bash
pd2af my_map.xml -m pure-af -l auto -o my_map_af.xml
```

### Drop complex activities

```bash
pd2af my_map.xml -m no-complex -o my_map_af.xml
```

### Overlay layout (keep non-model elements as greyed background)

```bash
pd2af my_map.xml -l overlay -o my_map_af.xml
```

## Getting help

```bash
pd2af --help
```
