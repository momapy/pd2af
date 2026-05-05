import momapy.celldesigner

import pd2af.predicates
import pd2af.solver
from pd2af.walkers import (
    BuildStep,
    BuildStepKind,
    Walker,
    register,
)
from pd2af.walkers._shared import (
    collect_ancestor_compartments,
    collect_templates_from_species,
    compartments_outermost_first,
)


def _walk_keep_species(cd_map, mode):
    """Yield build steps for a keep-species-style mode.

    Calls the existing clingo solver to determine which input species
    qualify as activities and which influences exist between them. Yields
    build steps that reuse input species/templates/compartments by
    identity (the new model objects ARE the input objects). Modulations
    are fresh objects with source/target bound to the new species.
    """
    clingo_model, id_to_model_element = pd2af.solver.solve(cd_map, mode=mode)
    activity_atoms = pd2af.solver._get_activity_atoms(clingo_model)
    influence_atoms = pd2af.solver._get_influence_atoms(clingo_model)
    kept_species = []
    for atom in activity_atoms:
        key = atom.key
        if not isinstance(key, pd2af.predicates.kept_species):
            raise ValueError(
                f"keep-species walker received non-kept_species activity "
                f"key {type(key).__name__}"
            )
        kept_species.append(id_to_model_element[key.species])
    immediate_compartments = {
        species.compartment
        for species in kept_species
        if getattr(species, "compartment", None) is not None
    }
    compartments = collect_ancestor_compartments(immediate_compartments)
    templates = collect_templates_from_species(kept_species)
    input_to_new = {}
    # Content-keyed canonicalisation: dataclass model elements with
    # compare=False on id_ are __eq__-equal across distinct identities.
    # The model builder dedups by content, so the walker must too —
    # otherwise placers see duplicate build steps for one canonical
    # element. Yield the first occurrence; remap subsequent occurrences
    # to the canonical via input_to_new.
    canonical_compartments = {}
    canonical_templates = {}
    canonical_species = {}
    for compartment in compartments_outermost_first(compartments):
        existing = canonical_compartments.get(compartment)
        if existing is None:
            canonical_compartments[compartment] = compartment
            input_to_new[id(compartment)] = compartment
            yield BuildStep(
                kind=BuildStepKind.COMPARTMENT,
                new_element=compartment,
                provenances=(compartment,),
            )
        else:
            input_to_new[id(compartment)] = existing
    for template in templates:
        existing = canonical_templates.get(template)
        if existing is None:
            canonical_templates[template] = template
            input_to_new[id(template)] = template
            yield BuildStep(
                kind=BuildStepKind.TEMPLATE,
                new_element=template,
                provenances=(template,),
            )
        else:
            input_to_new[id(template)] = existing
    for species in kept_species:
        existing = canonical_species.get(species)
        if existing is not None:
            input_to_new[id(species)] = existing
            continue
        canonical_species[species] = species
        input_to_new[id(species)] = species
        top_step = BuildStep(
            kind=BuildStepKind.SPECIES,
            new_element=species,
            provenances=(species,),
        )
        yield top_step
        for subunit in getattr(species, "subunits", ()) or ():
            sub_existing = canonical_species.get(subunit)
            if sub_existing is not None:
                input_to_new[id(subunit)] = sub_existing
                continue
            canonical_species[subunit] = subunit
            input_to_new[id(subunit)] = subunit
            yield BuildStep(
                kind=BuildStepKind.SPECIES,
                new_element=subunit,
                provenances=(subunit,),
                parent=top_step,
            )
    for atom in influence_atoms:
        influence_cls = (
            pd2af.predicates.predicate_to_model_element_class[type(atom)]
        )
        if not isinstance(atom.source, pd2af.predicates.kept_species):
            raise ValueError(
                "keep-species walker received non-kept_species influence "
                f"source {type(atom.source).__name__}"
            )
        if not isinstance(atom.target, pd2af.predicates.kept_species):
            raise ValueError(
                "keep-species walker received non-kept_species influence "
                f"target {type(atom.target).__name__}"
            )
        input_source = id_to_model_element[atom.source.species]
        input_target = id_to_model_element[atom.target.species]
        new_source = input_to_new[id(input_source)]
        new_target = input_to_new[id(input_target)]
        new_modulation = influence_cls(source=new_source, target=new_target)
        yield BuildStep(
            kind=BuildStepKind.MODULATION,
            new_element=new_modulation,
            provenances=(),
        )


@register("keep-species")
class KeepSpeciesWalker(Walker):
    def walk(self, cd_map):
        yield from _walk_keep_species(cd_map, mode="keep-species")


@register("keep-species-no-complex")
class KeepSpeciesNoComplexWalker(Walker):
    def walk(self, cd_map):
        yield from _walk_keep_species(cd_map, mode="keep-species-no-complex")


@register("casq")
class CasqWalker(Walker):
    """Casq mode goes through the same clingo-driven activity-atom path
    as keep-species; the casq profile's rule set decides which species
    survive. The resulting activity atoms are all `kept_species`, so the
    same walking logic applies."""

    def walk(self, cd_map):
        yield from _walk_keep_species(cd_map, mode="casq")
