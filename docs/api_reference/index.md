# API reference

pd2af does one thing in three steps, and the package follows them: `core`
checks the arguments and runs the transformation, `asp` asks clingo what the
activities and influences are, and `building` turns the answer into a map.

## Top level

- [Core](core.md): the public `transform()` entry point and the `TransformerResult` it returns
- [Modes](modes.md): everything a user chooses from, the input languages, the layout modes, the influence pairings and the transformation modes
- [CLI](cli.md): command-line entry point

## `asp`: the question asked to clingo

- [Solver](asp_solver.md): build the clingo control, solve, and return the answer with the element registry
- [Rules](asp_rules.md): the rule groups and the composition of a mode's ASP program
- [Predicates](asp_predicates.md): the clorm predicates the program and the builders share

## `building`: turning the answer into a map

Shared by both output languages:

- [Context](building_context.md): the slots the model pass fills and the layout pass reads back
- [Model](building_model.md): the model pass steps that are the same for either language
- [Layout](building_layout.md): the layout helpers that are the same for either language, styling, arcs and dot
- [Provenance](building_provenance.md): where each output element comes from, and the annotations and notes it inherits

CellDesigner output:

- [CellDesigner model](building_celldesigner_model.md): build the AF model from the clingo atoms
- [CellDesigner layout](building_celldesigner_layout.md): build the AF layout, plain, overlay or dot

SBGN-AF output:

- [SBGN-AF model](building_sbgn_model.md): build the SBGN-AF model from the clingo atoms
- [SBGN-AF layout](building_sbgn_layout.md): build the SBGN-AF layout, plain or dot
- [SBGN-AF labels](building_sbgn_labels.md): build an SBGN-AF activity's label from an SBGN-PD entity pool
