"""Dedup-and-remap pass mirroring momapy.io.utils.register_model_element.

momapy's dataclass-based model elements use compare=False on id_, so two
elements with identical fields and different ids are __eq__-equal but
have distinct identity. The CD reader/writer rely on the invariant that
content-equal elements share identity. When walkers yield input frozen
species objects, they can produce a model where the same content
appears under two different ids (e.g., a top-level FXYDs complex and
a content-equal FXYDs nested inside another complex's subunits — see
plans/overlay-mapping-identity-fix.md).

This module post-processes the freshly-built model to restore the
invariant: pick the smallest-id survivor for each content class
(matching the reader's tie-break), rebuild parent complexes whose
subunits referenced evicted instances, and remap modulation
source/target plus layout-model mapping values onto the survivors.
"""

import dataclasses

import momapy.celldesigner


def dedup_and_remap_model(model_builder, mapping_builder):
    canonical_by_content = {}

    def collect(species):
        existing = canonical_by_content.get(species)
        species_id = species.id_ or ""
        existing_id = existing.id_ or "" if existing is not None else None
        if existing is None or species_id < existing_id:
            canonical_by_content[species] = species
        if isinstance(species, momapy.celldesigner.Complex):
            for subunit in species.subunits:
                collect(subunit)

    for species in model_builder.species:
        collect(species)

    rebuild_cache = {}
    id_to_rebuilt = {}

    def rebuild_canonical(canonical):
        cached = rebuild_cache.get(id(canonical))
        if cached is not None:
            return cached
        # Placeholder guards against pathological cycles (model trees
        # are acyclic in practice, but cheap insurance).
        rebuild_cache[id(canonical)] = canonical
        if (
            isinstance(canonical, momapy.celldesigner.Complex)
            and canonical.subunits
        ):
            rebuilt_subunits = []
            any_changed = False
            for subunit in canonical.subunits:
                rebuilt_subunit = remap_species(subunit)
                rebuilt_subunits.append(rebuilt_subunit)
                if rebuilt_subunit is not subunit:
                    any_changed = True
            if any_changed:
                rebuilt = dataclasses.replace(
                    canonical, subunits=frozenset(rebuilt_subunits)
                )
            else:
                rebuilt = canonical
        else:
            rebuilt = canonical
        rebuild_cache[id(canonical)] = rebuilt
        return rebuilt

    def remap_species(species):
        if id(species) in id_to_rebuilt:
            return id_to_rebuilt[id(species)]
        canonical = canonical_by_content.get(species)
        if canonical is None:
            id_to_rebuilt[id(species)] = species
            return species
        rebuilt = rebuild_canonical(canonical)
        id_to_rebuilt[id(species)] = rebuilt
        return rebuilt

    new_species = set()
    for species in list(model_builder.species):
        new_species.add(remap_species(species))
    model_builder.species.clear()
    model_builder.species.update(new_species)

    modulation_remap = {}
    new_modulations = set()
    for modulation in list(model_builder.modulations):
        old_source = getattr(modulation, "source", None)
        old_target = getattr(modulation, "target", None)
        new_source = (
            remap_species(old_source) if old_source is not None else None
        )
        new_target = (
            remap_species(old_target) if old_target is not None else None
        )
        if (
            new_source is not old_source
            or new_target is not old_target
        ):
            new_modulation = dataclasses.replace(
                modulation, source=new_source, target=new_target
            )
            modulation_remap[id(modulation)] = new_modulation
        else:
            new_modulation = modulation
        # Drop modulations that collapsed onto themselves after dedup
        # (two PD species that were content-equal but distinct become the
        # same activity in AF; a modulation between them becomes a
        # meaningless self-loop and trips writer/reader geometry).
        if (
            new_source is not None
            and new_target is not None
            and new_source is new_target
        ):
            continue
        new_modulations.add(new_modulation)
    model_builder.modulations.clear()
    model_builder.modulations.update(new_modulations)

    if mapping_builder is not None:
        for key in list(mapping_builder.keys()):
            value = mapping_builder[key]
            if value is None:
                continue
            new_value = modulation_remap.get(id(value))
            if new_value is None:
                new_value = remap_species(value)
            if new_value is not value:
                mapping_builder[key] = new_value
