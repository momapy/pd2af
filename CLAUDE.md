# pd2af — design notes for Claude

## Conceptual frame: pd2af is a PD → AF transformation

pd2af conceptually transforms an SBGN-PD map into an SBGN-AF map. SBGN-AF
has only **activities**: opaque nodes with no internal structure. An
activity may carry a *unit of information* labelling the type of entity
performing it (e.g. "complex", "macromolecule"), but it has no subunits,
no state information, no template.

CellDesigner has no AF language. Its only AF-ish feature is
*influences between species* (which roughly mimic influences between
activities performed by entity pools). So when the output format is
CellDesigner, an activity is represented by a **species** that stands
in for the entity performing it. The species type carries the
provenance of the activity.

A transformation is three choices. **A mode decides what counts as an
activity; an option decides which activities are treated as the same thing.**
Deciding what is an activity is the solver's job, so a choice that changes it
has to change the ASP program and is a mode, while a choice that only groups
the results afterwards happens at the build stage and is an option. That leaves
**whether to keep or dissolve complexes** and **where the influences come
from** as the two mode axes, and **whether to strip post-translational
decorations** (state variables in SBGN; `modifications`, `structural_states`,
template `modification_residues`/`regions` and `homomultimer` in CellDesigner)
and merge content-equal results as an option. In every mode a subunit is a
*structural component*, never an independent activity: a subunit's activity and
influences are attributed to its outermost top-level complex (subunit-level
influences are a CellDesigner artifact, not standard AF).

| mode              | complexes                                                  | influences                                            |
| ----------------- | ---------------------------------------------------------- | ----------------------------------------------------- |
| `normal`          | kept (opaque; subunits carried in the label/structure, influences routed to the complex) | inferred: modulations + multi-hop chaining + consumption + activation by binding |
| `no-complex`      | dissolved: active subunits promoted to top-level activities | inferred                                              |
| `keep-reactions`  | kept (opaque; subunits routed to the complex)              | stated only: reactant→product positive influences + direct modulations |

| option              | off (the default)                                        | on                                                    |
| ------------------- | -------------------------------------------------------- | ----------------------------------------------------- |
| `keep_species`      | PTM decorations stripped, content-equal results merged   | each species kept as its own activity                 |
| `drop_compartments` | each compartment kept                                    | every compartment merged into the default one         |

`pd2af.modes.TRANSFORMATION_OPTIONS` is the fourth vocabulary next to
`Language`, `LayoutMode` and `InfluencePairingMode`: it holds each option's
flag, default and description in one place, and the CLI builds its flags, its
`list-modes` table and its `--json` payload from it, so a third option is one
entry rather than edits in five files. A mode may name defaults, never rules:
`TransformationMode.default_options` maps option names to the value that mode
starts from (`keep-reactions` sets `{"keep_species": True}`), an explicit flag
always wins, and a key outside `TRANSFORMATION_OPTIONS` fails at mode load
alongside the name-collision and rule-group checks.

Both options are build-stage concerns and change nothing in the ASP program.
`keep_species` is read off `context.keep_species` by the two model builders and
the CellDesigner layout builder. `drop_compartments` is read off
`context.drop_compartments` by each language's "which compartment does this
element go in" helper
(`pd2af.building.celldesigner.model._compartment_for_input_species`,
`pd2af.building.sbgn.model._compartment_for_input_element`): the answer becomes
the input map's `default` compartment for CellDesigner (`find_default_compartment`,
so the output still declares the compartment its species refer to) and `None`
for SBGN-AF. Everything downstream follows from content interning: species that
become equal collapse through `register_or_reuse`, the influences built from
them collapse in turn, and provenance unions their annotations.

`pd2af.modes.get_compatible_layout_modes` is the single place the "merged
activities need `dot`" rule lives: a merged activity is synthesized from
several input species and so has no original geometry for `plain`/`overlay` to
reuse. It takes the language and the two options and no mode, because the mode
does not bear on the answer — `plain` and `overlay` are available exactly when
`keep_species` is set without `drop_compartments`.

`normal` is the canonical AF mode and the default for the SBGN-AF output.
`keep_species` and `keep-reactions` deliberately deviate: they preserve the PD
entity structure for users who want a CD-native lossy reduction rather than a
true AF view.

`keep-reactions` also deviates on activity discovery: instead of
requiring a structural signal, **every species is an activity**, and
each reaction is rendered as one positive influence per (reactant,
product) pair. It keeps the modulation arcs and reaction modifiers as
single-hop influences, and takes none of the inference layers. It is
CellDesigner-only.

Stripping is a single recursive operation over the resolved entity
(`pd2af.building.celldesigner.model.get_or_make_stripped_species`;
`pd2af.building.sbgn.labels.make_label` with `include_state_variables=False`),
applied to *every* entity when `keep_species` is off — complexes and
non-templated entities included, not just templated forms of the same base
entity. The two structural-role activity keys are `keptSpeciesKey` (a top-level
entity, or the top-level complex a subunit resolves to via the shared
`resolvesToTopLevel` ASP relation) and `promotedSubunitKey` (a subunit lifted
to top level when its complex is dissolved).

## How the package is laid out

pd2af does one thing in three steps, and the package follows them:

| place | what it holds |
| --- | --- |
| `core.py` | the transformation itself: check the arguments, solve, build the map |
| `modes.py` | everything a user chooses from |
| `asp/` | the question asked to clingo |
| `building/` | how the answer becomes a map |
| `cli.py` | the command line, on top of `core` |

Inside `building/`, a module directly in the folder is shared by both output
languages (`context`, `model`, `layout`, `provenance`); a module in
`building/celldesigner/` or `building/sbgn/` belongs to that language only, and
each language folder holds a `model` and a `layout` matching the two build
passes. Nothing imports `core`, so the import graph has no cycle.

`pd2af.modes` holds the three vocabularies a transformation mode is defined
against, next to the modes themselves. Each is a `StrEnum` whose member values
are the tokens the outside world sees: `Language`, `LayoutMode` and
`InfluencePairingMode`. Because a member is a string, a token still serialises
to JSON, resolves an aspcompose variant and compares equal to a raw string a
contributed mode writes; building one (`LayoutMode(value)`) is what validates
it, so `transform` and the CLI convert once and pass members around after that.
`LANGUAGES` carries every per-language fact (`display_name`, the
`momapy_module` whose classes seed the ontology vocabulary, and the `map_class`
/ `model_class` an input is recognised by), keyed by `Language` member in the
order the CLI and the docs list them in. Adding a language is one member plus
one `LANGUAGES` entry; `get_language_from_map_or_model` infers the member by
walking the input. `LAYOUT_MODE_DESCRIPTIONS` holds the prose the CLI lists each
layout mode with, and `LAYOUT_MODES_BY_LANGUAGE` is the join of the first two
vocabularies. `auto` is not a layout mode but a request to pick one: it stays
the `AUTO` constant, and `transform` resolves it before the build stage sees
anything.

The token is the wire format, not just a label: it is the aspcompose variant
key `pd2af.asp.rules.build_program` resolves, a segment of every
language-specific rule identifier
(`activity:core:celldesigner:from_global_activate`, public via
`--exclude-rule`), and a member of a mode's `compatible_languages`.

The per-language build dispatch lives in `core._BUILD_BEHAVIOR_BY_LANGUAGE`
(the two build-pass modules, the output map class, the `dot` auto-layout
kwargs), keyed by language token: it selects build *behavior* rather than
defines a language, and hosting builder references in `pd2af.modes` would make
it drag the whole builder tree into every importer of `pd2af.modes`. Wrapping a
bare input model in a map is pure language-table lookup, so it belongs to
`pd2af.modes.make_map_from_model` instead.

## Modes are objects, contributed through an entry point

A mode is a `pd2af.modes.TransformationMode`: its name, its `docs` (the prose
the CLI lists it with), the rule groups its program is made of
(`rule_group_references` naming registered groups, `rule_group_definitions`
carrying groups the mode brings itself), the input languages it accepts, and
its `default_options`. `pd2af.asp.rules` owns the groups
and composes the program; the mode owns the membership, so no rule group decides
which modes include it. A mode that brings its own groups names them after
itself: `keep_reactions:*`.

`pd2af.modes.get_transformation_modes()` returns the built-ins in declaration
order followed by every mode contributed through the `pd2af.modes` entry-point
group. A contributed mode may not shadow an existing name, and any failure to
load one — bad import, wrong type, name collision, a group that fails
`registry.validate()`, an unknown `default_options` key — takes down every
pd2af entry point, deliberately.
`docs/generate_rules_reference.py` passes `_BUILTIN_TRANSFORMATION_MODES` to
`build_registry` so a contributed mode never reaches the published reference.

Group lists are complete, not leaf-only. The `preparation` slot is what makes a
hand-written list safe: `preparation:complex` and `preparation:no_complex` fill
it, the three consumers of `hasActivityCarrier`/`hasActivityKey`
(`influences:core`, `influences:consumption`, `gates:core`) depend on the slot,
so a mode that omits its preparation group raises `unfilled_slot` instead of
silently emitting an activity-less program.

## Annotations and notes

momapy stores RDF/MIRIAM annotations and notes not on model elements but
in side-tables on the `ReaderResult` (`element_to_annotations`,
`element_to_notes`: `Mapping[model_element -> frozenset]`). The transform
carries them: pass the reader's side-tables to `transform(...,
element_to_annotations=..., element_to_notes=...)` and it returns
output-keyed side-tables on the `TransformerResult`, ready to hand to the
writer. The carrier is `TransformerResult.output_element_to_input_elements`, keyed to
the origin direction (`output_element -> frozenset(input_elements)`); the pure
remap lives in `pd2af.building.provenance.carry_annotations_through_provenance`,
which unions the metadata of every input that merged into a given output
(so a merged activity gathers the annotations of all its forms).

Coverage is species/activities, **complex subunits** at any depth (paired from
the species provenance by `pd2af.building.provenance.record_provenance_for_subunit_trees`,
through the builder's `input_model_element_to_canonical_model_element` map in
the merged modes and by content-equality in the kept modes), logical operators,
**compartments** (folded into `provenance` via `context.compartment_emissions`)
and the **map** itself (both the CellDesigner and SBGN writers emit map-level
annotations/notes). Deliberately dropped: modulation/influence and reaction
annotations (AF influences are synthesized and fan-out/collapse, so there is
no clean target without threading source provenance through the ASP layer);
and the stdout-pickle path (a bare-map pickle cannot hold the side-tables, so
carry needs `-o file.xml` / `.sbgn`).

## Model-element dedup invariant

momapy's dataclass-based model elements use `compare=False` on `id_`,
so `__eq__` and `__hash__` are content-based: two model elements with
identical fields and different `id_` are *equal* but have *different
identities*.

momapy's readers enforce this invariant (see
`momapy/io/utils.py::register_model_element`): when an equal element
is registered a second time, only one survives (smaller `id_` wins),
and all references to the evicted element are remapped to the
survivor (`remap_model_element`). The XML-id lookup still resolves to
the survivor, so dangling references can't happen.

**pd2af must preserve this invariant when constructing new models.**
A naive `set(...)` of model elements silently drops duplicates without
remapping references — leaving callers (e.g. `subunits`,
`modulation.source`, `proteinReference`) pointing at evicted elements,
which causes `KeyError` on read-back of the written CellDesigner XML.

The right fix when this shape of bug appears is a dedup-and-remap pass
over the constructed model (mirroring `register_model_element`), not
ad-hoc patching. See `pd2af.building.model.register_or_reuse`.

## Read-back as the integration test

The standard integration check for any pd2af change is:

```python
new_map = pd2af.transform(cd_map, mode=..., layout_mode=...).obj
momapy.io.core.write(new_map, path, writer="celldesigner")
momapy.io.core.read(path, reader="celldesigner")  # must not raise
```

CellDesigner input gives a `CellDesignerMap`, written and read with the
`celldesigner` writer/reader; SBGN-PD input gives an `SBGNAFMap`, written and
read with the `sbgnml` ones:

```python
new_map = pd2af.transform(sbgn_pd_map, mode=..., layout_mode=...).obj
momapy.io.core.write(new_map, path, writer="sbgnml")
momapy.io.core.read(path, reader="sbgnml")  # must not raise
```

The in-repo check is `pytest` plus the read-back above. A broader sweep over all
(mode, layout_mode) combinations against every map in `tests/maps/celldesigner/`
is run from a local harness that is not committed to this repo.
