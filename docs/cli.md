# CLI Reference

## Overview

The `pd2af` command-line interface transforms a process-description (PD) map — CellDesigner or SBGN-PD — into an activity-flow (AF) map.

The output map is written to stdout as a [momapy](https://github.com/adrienrougny/momapy) pickle (preserving layout styling) so it can be piped into `momapy visualize`. With `-o`, the writer is chosen from the output file extension. With `-V`, the output map is opened in the momapy browser viewer, and nothing is written to stdout.

An unreadable input file, an unknown element id, or a mode the input or the layout does not support is reported as a short `error: ...` message on stderr, with exit status 1 and no output file written.

## Synopsis

```bash
pd2af transform <input_file> [-m {normal,no-complex,keep-reactions}] [-l {plain,overlay,dot,auto}] [-p {cross,nearest}] [--keep-species | --no-keep-species] [--drop-compartments | --no-drop-compartments] [-a <id> ...] [-i <id> ...] [-A | -I] [--exclude-group <group> ...] [--exclude-rule <rule> ...] [-o <output_file>] [-V]
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

`transform` is the implicit subcommand: when the first argument is neither a
subcommand nor `-h`/`--help`/`--version`, it is inserted, so `pd2af map.xml -m
no-complex` runs exactly like `pd2af transform map.xml -m no-complex`.

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
| `--keep-species` |  | Keep each species as its own activity instead of merging the forms of the same base entity. Off by default, on for `keep-reactions`; `--no-keep-species` forces it off; see below |
| `--drop-compartments` |  | Drop the compartments, so that species differing only by compartment become a single activity and the influences that become equal merge in turn. Off by default; works with every transformation mode; requires `--layout-mode dot` (or `auto`); see below |
| `--set-active` | `-a` | Mark the element with this `id_` (a species or entity pool) as active, surfacing it as an activity even when the map gives it no structural activity signal. Repeatable: `-a sa1 -a sa2`. Wins over `--set-all-inactive` for these ids. Unknown ids raise an error |
| `--set-inactive` | `-i` | Mark the element with this `id_` as NOT active, suppressing any activity for it and overriding the automatic discovery. Repeatable. Wins over `--set-all-active` for these ids; passing an id to both `-a` and `-i` is an error |
| `--set-all-active` | `-A` | Mark every top-level species / entity pool as active (subunits excluded). Per-id `-i` overrides it. Mutually exclusive with `-I` |
| `--set-all-inactive` | `-I` | Suppress activity for every element, subunits included. Per-id `-a` overrides it. Mutually exclusive with `-A` |
| `--exclude-group` |  | Drop a whole rule group, e.g. `activity:phenotype` to stop treating phenotypes as activities, or `paths:chaining` to keep only single-hop influences. Repeatable. Excluding a group another included group depends on is an error; run `pd2af list-groups` for the excludable groups per mode |
| `--exclude-rule` |  | Drop a single rule by identifier, e.g. `influences:kind:celldesigner:catalysis` — the scalpel for the table groups. Repeatable. Prefer `--exclude-group` for coherent behaviors |
| `--output` | `-o` | Write output to this file instead of stdout. Input RDF annotations and notes are carried onto the corresponding output elements for file output only; the stdout pickle cannot carry them |
| `--visualize` | `-V` | Open the output map in the momapy browser viewer. With `-o` the map is also written to the file; without `-o` nothing is written to stdout |

## Transformation modes (`-m`)

A mode decides what counts as an activity; an option decides which activities are treated as the same thing.

| Mode | Description |
|------|-------------|
| `normal` | Keep complexes as activities of their own, routing a subunit's influences to the complex it belongs to. An active subunit of an activity-bearing complex is subsumed into the complex. Infers influences beyond the stated modulations: multi-hop chaining, consumption/sparing and activation by binding (see below). **Default**. |
| `no-complex` | Replace a complex that has an active subunit with those subunits, promoted to top-level activities, and route the influences through them. Infers the same influences as `normal`. |
| `keep-reactions` | Keep the PD topology itself: every species is an activity (no structural activity signal required), and every reaction becomes one positive influence per (reactant, product) pair. Modulation arcs and reaction modifiers are kept as single-hop influences with their own kinds; no multi-hop chaining, consumption/sparing or binding-activation inference. CellDesigner-only; starts from `--keep-species`. |

## Transformation options

Each option is a flag with a matching `--no-` spelling. A mode names the options it starts from and never overrules an explicit flag, so `-m keep-reactions --no-keep-species` gives the merged reading of that mode.

### `--keep-species`

Keep each species as its own activity (template + state + compartment) instead of merging the forms of the same base entity. Off by default: the forms of the same base entity within the same compartment collapse into a single activity, which is a true PD→AF transform with no PD remnants — the only style expressible in SBGN-AF. On by default for `keep-reactions`.

### `--drop-compartments`

Drop the compartments, so that species differing only by compartment become a single activity and the influences that become equal merge in turn.

CellDesigner output keeps a single compartment, the `default` one every CellDesigner map declares, and every species points at it; SBGN-AF output carries no compartment at all, SBGN-AF having no default compartment. Annotations and notes follow the merge: an activity merged from several compartments gathers the metadata of every species that collapsed into it, and the metadata of the removed compartments is dropped with them.

```bash
pd2af transform my_map.xml --keep-species --drop-compartments -o my_map_af.xml
```

## Layout modes (`-l`)

| Mode | Description |
|------|-------------|
| `auto` | Pick automatically from the input: a map gets `dot`, a bare model gets no layout. **Default**. |
| `dot` | Graphviz `dot` auto-layout. Requires `dot` on `PATH`. Required whenever activities merge: unless `--keep-species` is set, and always under `--drop-compartments`. |
| `plain` | Reuse original positions; only model elements are kept. Available exactly when `--keep-species` is set without `--drop-compartments`. |
| `overlay` | Reuse the full original layout; non-model elements greyed out. Available exactly when `--keep-species` is set without `--drop-compartments`, and for CellDesigner output only. |

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

### Activation by binding (`influences:binding_activation`)

In `normal` and `no-complex`, a reactant that comes out of a complex-forming
reaction as an active subunit was activated by the binding: in
`L + R -> L:R` with `R` drawn active only inside `L:R`, `L` becomes an activity
and positively influences the activity carrying `R`, the complex `L:R` in
`normal` and the promoted `R` in `no-complex`. Two forms are the same entity
when they share a template, else the same class and name (CellDesigner), or the
same entity kind and label (SBGN-PD); state, multimer cardinality and
compartment are ignored. An active reactant may still be a source
(`Ras:GTP + Raf -> Ras:GTP:Raf` gives `Ras:GTP -> Ras:GTP:Raf`), two entities
activated by the same reaction draw no edge between each other, and the
binding edge is direct: it does not chain through later reactions. Drop it
with `--exclude-group influences:binding_activation`.

## Examples

### Basic transformation to stdout (pipe to momapy)

```bash
pd2af transform my_map.xml | momapy visualize -
```

### Transform and open the result in the viewer

```bash
pd2af transform my_map.xml -V
```

### Transform to CellDesigner XML (default `normal` mode, auto layout)

```bash
pd2af transform my_map.xml -o my_map_af.xml
```

### Drop complexes via the no-complex mode

```bash
pd2af transform my_map.xml -m no-complex -o my_map_af.xml
```

### Keep each PD species as its own activity, reusing original layout

```bash
pd2af transform my_map.xml --keep-species -l plain -o my_map_af.xml
```

### Drop complexes but keep each remaining species, with overlay layout

```bash
pd2af transform my_map.xml -m no-complex --keep-species -l overlay -o my_map_af.xml
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

### Keep only the influences the map states

```bash
pd2af transform my_map.xml --exclude-group paths:chaining --exclude-group influences:consumption --exclude-group influences:binding_activation -o my_map_af.xml
```

## Getting help

```bash
pd2af --help
pd2af transform --help
pd2af list-modes --help
pd2af list-groups --help
pd2af --version
```
