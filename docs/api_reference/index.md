# API reference

The pd2af public API is organized into the following modules.

## Top-level

- [Core](core.md) — the public `transform()` entry point
- [CLI](cli.md) — command-line entry point

## Transformation pipeline

- [Solver](solver.md) — clingo solver wrapper, builds the new CellDesigner model
- [Rules](rules.md) — Python helpers that emit ASP rules from a CellDesigner map
- [Predicates](predicates.md) — clorm predicate definitions used by the ASP program
- [Ontology](ontology.md) — domain ontology (species/reaction/modulation kinds)

## Layout

- [Layouts](layouts.md) — build the AF map's layout (plain / overlay / auto)

## Auxiliary

- [CASQ](casq.md) — CASQ-based comparison/utility helpers
- [Utils](utils.md) — miscellaneous utilities
