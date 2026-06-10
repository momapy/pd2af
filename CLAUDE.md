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

The transformation modes encode two orthogonal choices: whether to
**strip post-translational decorations** (state variables in SBGN;
`modifications`, `structural_states`, template `modification_residues`/
`regions` and `homomultimer` in CellDesigner) and merge content-equal
results, and whether to **keep or dissolve complexes**. In every mode a
subunit is a *structural component*, never an independent activity: a
subunit's activity and influences are attributed to its outermost
top-level complex (subunit-level influences are a CellDesigner artifact,
not standard AF). All of this is decided by the **mode at build time**,
not encoded in the activity key.

| mode                       | PTM decorations    | complexes                                                  |
| -------------------------- | ------------------ | ---------------------------------------------------------- |
| `normal`                   | stripped & merged  | kept (opaque; subunits carried in the label/structure, influences routed to the complex) |
| `no-complex`               | stripped & merged  | dissolved: active subunits promoted to top-level activities |
| `keep-species`             | kept               | kept (opaque; subunits routed to the complex)              |
| `keep-species-no-complex`  | kept               | dissolved: active subunits promoted to top-level activities |
| `casq`                     | kept               | kept (CASQ-specific deletion pruning)                      |

`normal` is the canonical AF mode and will be the default for the
future SBGN-AF output. `keep-species*` and `casq` deliberately deviate:
they preserve PD proteoform structure for users who want a CD-native
lossy reduction rather than a true AF view.

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

Sweep harness: `/tmp/pd2af_test/sweep.py` runs all
(mode, layout_mode) combos against every map in
`tests/maps/celldesigner/`.
