# API reference

The pd2af public API is organized into the following modules.

## Top-level

- [Core](core.md): the public `transform()` entry point
- [CLI](cli.md): command-line entry point
- [Languages](languages.md): input-language tokens, layout-mode vocabulary and inference from the input map type
- [Modes](modes.md): the transformation modes and the entry point that contributes new ones

## Transformation pipeline

- [Solver](solver.md): clingo solver wrapper, resolves the ASP program into activity/influence atoms
- [Rules](rules.md): the rule groups and the composition of a mode's ASP program
- [Predicates](predicates.md): clorm predicate definitions used by the ASP program
- [Ontology](ontology.md): domain ontology (species/reaction/modulation kinds)

## Build

- [Build](build.md): coordinator that drives the two-phase model/layout builder pipeline
- [CellDesigner model](celldesigner_building_model.md): build the AF model from clingo atoms (CellDesigner output)
- [CellDesigner layout](celldesigner_building_layout.md): build the AF layout (plain / overlay / dot), CellDesigner output
- [SBGN-AF model](sbgn_building_model.md): build the SBGN-AF model from clingo atoms
- [SBGN-AF layout](sbgn_building_layout.md): build the SBGN-AF layout (plain / dot)
- [SBGN-AF labels](sbgn_labels.md): build an SBGN-AF activity's label from an SBGN-PD entity pool

## Auxiliary

- [Utils](utils.md): miscellaneous utilities
