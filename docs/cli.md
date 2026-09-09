# CLI Reference

## Overview

The `pd2af` command-line interface transforms a process-description (PD) map — CellDesigner or SBGN-PD — into an activity-flow (AF) map.

The output map is written to stdout as a [momapy](https://github.com/adrienrougny/momapy) pickle (preserving layout styling) so it can be piped into `momapy visualize`. With `-o`, the writer is chosen from the output file extension.

An unreadable input file, an unknown element id, or a mode the input or the layout does not support is reported as a short `error: ...` message on stderr, with exit status 1 and no output file written.

## Synopsis

```bash
pd2af transform <input_file> [-m {normal,normal-no-complex,keep-species,keep-species-no-complex,keep-reactions}] [-l {plain,overlay,dot,auto}] [-p {cross,nearest}] [-a <id> ...] [-i <id> ...] [-A | -I] [--exclude-group <group> ...] [--exclude-rule <rule> ...] [-o <output_file>]
pd2af list-modes [--json]
pd2af list-groups [--json]
pd2af --version
```

The CLI is organised into subcommands:

| Subcommand | Description |
|------------|-------------|
| `transform` | Transform a process-description map into an activity-flow map |
| `list-modes` | List transformation modes, layout modes, their compatibilities, and input languages |
| `list-groups` | List the rule groups each transformation mode uses, marked excludable or mandatory |

## `transform` arguments

| Argument | Description |
|----------|-------------|
| `input_file` | Input process-description map: CellDesigner XML or SBGN-ML (reads from stdin if omitted) |

## `transform` options

| Option | Short | Description |
|--------|-------|-------------|
| `--transformation-mode` | `-m` | Transformation mode (default: `normal`); see below |
| `--layout-mode` | `-l` | Layout strategy (default: `auto`); see below |
| `--influence-pairing` | `-p` | How to draw an influence whose source/target maps to several glyphs: `cross` (default) or `nearest` |
| `--set-active` | `-a` | Mark the element with this `id_` (a species or entity pool) as active, surfacing it as an activity even when the map gives it no structural activity signal. Repeatable: `-a sa1 -a sa2`. Wins over `--set-all-inactive` for these ids. Unknown ids raise an error |
| `--set-inactive` | `-i` | Mark the element with this `id_` as NOT active, suppressing any activity for it and overriding the automatic discovery. Repeatable. Wins over `--set-all-active` for these ids; passing an id to both `-a` and `-i` is an error |
| `--set-all-active` | `-A` | Mark every top-level species / entity pool as active (subunits excluded). Per-id `-i` overrides it. Mutually exclusive with `-I` |
| `--set-all-inactive` | `-I` | Suppress activity for every element, subunits included. Per-id `-a` overrides it. Mutually exclusive with `-A` |
| `--exclude-group` |  | Drop a whole rule group, e.g. `activity:phenotype` to stop treating phenotypes as activities, or `paths:chaining` to keep only single-hop influences. Repeatable. Excluding a group another included group depends on is an error; run `pd2af list-groups` for the excludable groups per mode |
| `--exclude-rule` |  | Drop a single rule by identifier, e.g. `influences:kind:celldesigner:catalysis` — the scalpel for the table groups. Repeatable. Prefer `--exclude-group` for coherent behaviors |
| `--output` | `-o` | Write output to this file instead of stdout. Input RDF annotations and notes are carried onto the corresponding output elements for file output only; the stdout pickle cannot carry them |

## Transformation modes (`-m`)

Four of the five modes lie on two orthogonal axes — species treatment and complex treatment. The fifth, `keep-reactions`, changes where the influences come from.

|                       | keep complexes            | drop complexes (route through subunits) |
|-----------------------|---------------------------|------------------------------------------|
| **merge proteoforms** | `normal` *(default)*      | `normal-no-complex`                             |
| **keep each species** | `keep-species`            | `keep-species-no-complex`                |

| Mode | Description |
|------|-------------|
| `normal` | Merge proteoforms of the same template (and compartment) into a single activity, but keep complexes. Influences involving an active complex route through the complex; active subunits of an activity-bearing complex are subsumed into the complex. True PD→AF transform — the only style expressible in SBGN PD. Requires `--layout-mode dot` (or `auto`). |
| `normal-no-complex` | Like `normal`, but drop any complex that has an active subunit and route influences through the subunits. Requires `--layout-mode dot` (or `auto`). |
| `keep-species` | Emit one activity per distinct PD species (template + state + compartment). Keep complexes; active subunits of an activity-bearing complex are subsumed into the complex. |
| `keep-species-no-complex` | Like `keep-species`, but drop complexes with an active subunit. |
| `keep-reactions` | Keep the PD topology itself: every species is an activity (no structural activity signal required), and every reaction becomes one positive influence per (reactant, product) pair. Modulation arcs and reaction modifiers are kept as single-hop influences with their own kinds; no multi-hop chaining and no consumption/sparing inference. Complexes and PTM decorations are kept, as in `keep-species`. CellDesigner-only. |

The merging modes (`normal`, `normal-no-complex`) require `--layout-mode dot` because positions from the original PD map cannot be reused for synthesized merged-proteoform activities.

## Layout modes (`-l`)

| Mode | Description |
|------|-------------|
| `auto` | Pick automatically from the input: a map gets `dot`, a bare model gets no layout. **Default**. |
| `dot` | Graphviz `dot` auto-layout. Requires `dot` on `PATH`. Required for `normal` and `normal-no-complex`. |
| `plain` | Reuse original positions; only model elements are kept. Available for `keep-species`, `keep-species-no-complex` and `keep-reactions`. |
| `overlay` | Reuse the full original layout; non-model elements greyed out. Available for `keep-species`, `keep-species-no-complex` and `keep-reactions`. |

## Output writers

The writer used with `-o` is selected from the file extension:

| Extension | Writer |
|-----------|--------|
| `.xml`, `.sbml` | CellDesigner XML |
| `.sbgn`, `.sbgnml` | SBGN-ML |
| `.pickle`, `.pkl` | momapy pickle |
| (other) | momapy pickle |

Any other extension, including `.svg`, writes a momapy pickle: an extension
is a writer choice, never a request to render an image.

When `-o` is omitted, the map is always written to stdout as a momapy pickle.

## `list-modes` subcommand

List the transformation modes, layout modes, their compatibilities, and the
supported input languages. With `--json`, the same data is emitted as
structured JSON for scripting.

```bash
pd2af list-modes
pd2af list-modes --json
```

The JSON payload mirrors the rendered tables exactly — one entry per row, one
key per column. It has two top-level keys: `transformation_modes` (each with
`transformation_mode`, its compatible `layout_modes` and `languages`, and a
`description`) and `layout_modes` (each with `layout_mode`, its compatible
`languages`, and a `description`).

## `list-groups` subcommand

List the rule groups each transformation mode uses, each marked *excludable* (a
dependency-graph leaf `--exclude-group` can drop) or *mandatory* (depended on by
another group, so not excludable). With `--json`, the same data is emitted as
structured JSON for scripting.

```bash
pd2af list-groups
pd2af list-groups --json
```

## Examples

### Basic transformation to stdout (pipe to momapy)

```bash
pd2af transform my_map.xml | momapy visualize -
```

### Transform to CellDesigner XML (default `normal` mode, auto layout)

```bash
pd2af transform my_map.xml -o my_map_af.xml
```

### Drop complexes via the normal-no-complex mode

```bash
pd2af transform my_map.xml -m normal-no-complex -o my_map_af.xml
```

### Keep each PD species as its own activity, reusing original layout

```bash
pd2af transform my_map.xml -m keep-species -l plain -o my_map_af.xml
```

### Drop complexes but keep each remaining species, with overlay layout

```bash
pd2af transform my_map.xml -m keep-species-no-complex -l overlay -o my_map_af.xml
```

### Transform an SBGN-PD map into an SBGN-AF map

```bash
pd2af transform my_map.sbgn -o my_map_af.sbgn
```

### Mark elements active by id (input parameters)

```bash
pd2af transform my_map.xml -a sa1 -a sa2 -o my_map_af.xml
```

### Start from every species active, then silence a few

```bash
pd2af transform my_map.xml -A -i sa1 -i sa2 -o my_map_af.xml
```

### Drop a rule group

```bash
pd2af transform my_map.xml --exclude-group paths:chaining -o my_map_af.xml
```

## Getting help

```bash
pd2af --help
pd2af transform --help
pd2af list-modes --help
pd2af list-groups --help
pd2af --version
```
