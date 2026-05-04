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

The transformation modes encode different choices about how PD entities
collapse into activities:

| mode                       | proteoforms                | complexes                                  |
| -------------------------- | -------------------------- | ------------------------------------------ |
| `normal`                   | merged per (template, compartment) | kept as flat complex marker        |
| `no-complex`               | merged per (template, compartment) | broken into active subunits; kept flat if no active subunit |
| `keep-species`             | each PD species → own activity | kept (PD-style, with subunits)         |
| `keep-species-no-complex`  | each PD species → own activity | broken into active subunits; otherwise kept (PD-style) |
| `casq`                     | each PD species → own activity | (CASQ-specific pruning)                |

`normal` is the canonical AF mode and will be the default for the
future SBGN-AF output. `keep-species*` and `casq` deliberately deviate:
they preserve PD structure for users who want a CD-native lossy
reduction rather than a true AF view.

## Implication: complexes are flat in `normal` and `no-complex`

In the AF view, an activity originating from a complex is just an
activity with "complex" provenance — there is no internal subunit
structure. When we represent that activity as a CD complex species,
the species must be flat: `subunits=frozenset()`.

This is enforced in `solver.make_new_cd_model` for `normal` and
`no-complex`. `keep-species*` and `casq` keep subunits intact.

This also avoids a structural bug class: when a complex retains its
subunits, those subunits drag in their original templates (`p_X`),
which can collide with synthesized merged templates
(`merged_template__id_Y`) under dataclass equality (see below).

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

The flat-complex rule in `normal`/`no-complex` removes the most common
trigger of this bug, but if other shapes of the same problem appear,
the right fix is a dedup-and-remap pass over the constructed model
(mirroring `register_model_element`), not ad-hoc patching.

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
