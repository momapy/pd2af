# pd2af — design notes for Claude

## Conceptual frame: pd2af is a PD → AF transformation

pd2af conceptually transforms an SBGN-PD map into an SBGN-AF map. SBGN-AF
has only **activities**: opaque nodes with no internal structure. An
activity may carry a *unit of information* labelling the type of entity
performing it (e.g. "complex", "macromolecule"), but it has no subunits,
no proteoform information, no template.

CellDesigner has no AF language. Its only AF-ish feature is
*influences between species* (which roughly mimic influences between
activities performed by entity pools). So when the output format is
CellDesigner, an activity is represented by a **species** that stands
in for the entity performing it. The species type carries the
provenance of the activity.

The transformation modes encode three choices: whether to **strip
post-translational decorations** (state variables in SBGN;
`modifications`, `structural_states`, template `modification_residues`/
`regions` and `homomultimer` in CellDesigner) and merge content-equal
results, whether to **keep or dissolve complexes**, and **where the
influences come from**. The first two are orthogonal and give the four
`normal`/`keep-species` modes; the third is what sets `keep-reactions`
and `casq` apart. In every mode a subunit is a *structural component*,
never an independent activity: a subunit's activity and influences are
attributed to its outermost top-level complex (subunit-level influences
are a CellDesigner artifact, not standard AF). The first two choices are
decided by the **mode at build time**, not encoded in the activity key;
the third is decided in the ASP layer by which rule groups the mode's
profile includes.

| mode                       | PTM decorations    | complexes                                                  | influences                                            |
| -------------------------- | ------------------ | ---------------------------------------------------------- | ----------------------------------------------------- |
| `normal`                   | stripped & merged  | kept (opaque; subunits carried in the label/structure, influences routed to the complex) | inferred: modulations + multi-hop chaining + consumption |
| `normal-no-complex`               | stripped & merged  | dissolved: active subunits promoted to top-level activities | inferred                                              |
| `keep-species`             | kept               | kept (opaque; subunits routed to the complex)              | inferred                                              |
| `keep-species-no-complex`  | kept               | dissolved: active subunits promoted to top-level activities | inferred                                              |
| `keep-reactions`           | kept               | kept (opaque; subunits routed to the complex)              | stated only: reactant→product positive influences + direct modulations |
| `casq`                     | kept               | kept (CASQ-specific deletion pruning)                      | CASQ-specific direct emission with bridged rewiring    |

`normal` is the canonical AF mode and will be the default for the
future SBGN-AF output. `keep-species*`, `keep-reactions` and `casq`
deliberately deviate: they preserve PD proteoform structure for users
who want a CD-native lossy reduction rather than a true AF view.

`keep-reactions` also deviates on activity discovery: instead of
requiring a structural signal, **every species is an activity**, and
each reaction is rendered as one positive influence per (reactant,
product) pair — the same direction `casq` takes for that relationship.
It keeps the modulation arcs and reaction modifiers as single-hop
influences, and takes none of the inference layers. It is
CellDesigner-only for now.

Stripping is a single recursive operation over the resolved entity
(`pd2af.celldesigner.building_model.get_or_make_stripped_species`;
`pd2af.sbgn.labels.build_label` with `include_state_variables=False`),
applied to *every* entity in the merged modes — complexes and
non-templated entities included, not just templated proteoforms. The
two structural-role activity keys are `kept_species` (a top-level
entity, or the top-level complex a subunit resolves to via the shared
`topLevel` ASP relation) and `promoted_subunit` (a subunit lifted to top
level when its complex is dissolved). `MERGED_PROTEOFORM_MODES` in
`pd2af.languages` is the single source of truth for which modes strip.

## Annotations and notes

momapy stores RDF/MIRIAM annotations and notes not on model elements but
in side-tables on the `ReaderResult` (`element_to_annotations`,
`element_to_notes`: `Mapping[model_element -> frozenset]`). The transform
carries them: pass the reader's side-tables to `transform(...,
element_to_annotations=..., element_to_notes=...)` and it returns
output-keyed side-tables on the `TransformerResult`, ready to hand to the
writer. The carrier is `TransformerResult.provenance`, re-keyed to the
origin direction (`output_element -> frozenset(input_elements)`); the pure
remap lives in `pd2af.annotations.carry_annotations_through_provenance`,
which unions the metadata of every input that merged into a given output
(so a merged activity gathers the annotations of all its proteoforms).

Coverage is species/activities, logical operators, **compartments** (folded
into `provenance` via `context.compartment_emissions`) and the **map**
itself (both the CellDesigner and SBGN writers emit map-level
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
ad-hoc patching. See `pd2af.utils.register_or_reuse`.

## Read-back as the integration test

The standard integration check for any pd2af change is:

```python
new_map = pd2af.transform(cd_map, mode=..., layout_mode=...)
momapy.io.core.write(new_map, path, writer="celldesigner")
momapy.io.core.read(path, reader="celldesigner")  # must not raise
```

The in-repo check is `pytest` plus the read-back above. A broader sweep over all
(mode, layout_mode) combinations against every map in `tests/maps/celldesigner/`
is run from a local harness that is not committed to this repo.
