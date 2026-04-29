# CLI Reference

## Overview

The `pd2af` command-line interface transforms a CellDesigner process-description (PD) map into an activity-flow (AF) map.

The output map is written to stdout as a [momapy](https://github.com/adrienrougny/momapy) pickle (preserving layout styling) so it can be piped into `momapy visualize`. With `-o`, the writer is chosen from the output file extension.

## Synopsis

```bash
pd2af <input_file> [-m {normal,no-complex,keep-species,keep-species-no-complex}] [-l {plain,overlay,auto}] [-o <output_file>]
```

## Arguments

| Argument | Description |
|----------|-------------|
| `input_file` | Input CellDesigner XML file |

## Options

| Option | Short | Description |
|--------|-------|-------------|
| `--mode` | `-m` | Transformation mode (default: `normal`); see below |
| `--layout` | `-l` | Layout strategy (default: `auto`); see below |
| `--output` | `-o` | Write output to this file instead of stdout |

## Transformation modes (`-m`)

The four modes lie on two orthogonal axes — species treatment and complex treatment:

|                       | keep complexes            | drop complexes (route through subunits) |
|-----------------------|---------------------------|------------------------------------------|
| **merge proteoforms** | `normal` *(default)*      | `no-complex`                             |
| **keep each species** | `keep-species`            | `keep-species-no-complex`                |

| Mode | Description |
|------|-------------|
| `normal` | Merge proteoforms of the same template (and compartment) into a single activity, but keep complexes. Influences involving an active complex route through the complex; active subunits of an activity-bearing complex are subsumed into the complex. True PD→AF transform — the only style expressible in SBGN PD. Requires `--layout auto`. |
| `no-complex` | Like `normal`, but drop any complex that has an active subunit and route influences through the subunits. Requires `--layout auto`. |
| `keep-species` | Emit one activity per distinct PD species (template + state + compartment). Keep complexes; active subunits of an activity-bearing complex are subsumed into the complex. CellDesigner-only. |
| `keep-species-no-complex` | Like `keep-species`, but drop complexes with an active subunit. CellDesigner-only. |

The merging modes (`normal`, `no-complex`) require `--layout auto` because positions from the original PD map cannot be reused for synthesized merged-proteoform activities.

## Layout modes (`-l`)

| Mode | Description |
|------|-------------|
| `auto` | Graphviz `dot` auto-layout. Requires `dot` on `PATH`. **Default**, and required for `normal` and `no-complex`. |
| `plain` | Reuse original positions; only model elements are kept. Available for `keep-species` and `keep-species-no-complex`. |
| `overlay` | Reuse the full original layout; non-model elements greyed out. Available for `keep-species` and `keep-species-no-complex`. |

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

### Transform to CellDesigner XML (default `normal` mode, auto layout)

```bash
pd2af my_map.xml -o my_map_af.xml
```

### Drop complexes via the no-complex mode

```bash
pd2af my_map.xml -m no-complex -o my_map_af.xml
```

### Keep each PD species as its own activity, reusing original layout

```bash
pd2af my_map.xml -m keep-species -l plain -o my_map_af.xml
```

### Drop complexes but keep each remaining species, with overlay layout

```bash
pd2af my_map.xml -m keep-species-no-complex -l overlay -o my_map_af.xml
```

## Getting help

```bash
pd2af --help
```
