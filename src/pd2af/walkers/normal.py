import dataclasses

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
    compartments_outermost_first,
)


_FLAT_COMPLEX_MODES = frozenset({"normal", "no-complex"})


def _flatten_complex(species):
    if not isinstance(species, momapy.celldesigner.Complex):
        return species
    if not getattr(species, "subunits", None):
        return species
    return dataclasses.replace(species, subunits=frozenset())


def _resolve(
    key,
    id_to_model_element,
    key_to_species,
    key_to_provenance,
    stripped_template_cache,
    flatten_complexes,
):
    cached = key_to_species.get(key)
    if cached is not None:
        return cached, key_to_provenance[key]
    if isinstance(key, pd2af.predicates.kept_species):
        input_species = id_to_model_element[key.species]
        if flatten_complexes:
            species = _flatten_complex(input_species)
            # Flattening drops subunits, so the input layout subtree (which
            # carries subunit/state/modification glyphs) no longer matches
            # the model. Treat as having no usable input layout — the
            # placer will synthesise a stub.
            if species is input_species:
                provenance = (input_species,)
            else:
                provenance = ()
        else:
            species = input_species
            provenance = (input_species,)
    elif isinstance(key, pd2af.predicates.derived_proteoform_class):
        species = pd2af.solver._make_synthetic_species(
            key, id_to_model_element, stripped_template_cache
        )
        provenance = ()
    else:
        raise ValueError(f"unknown activity key wrapper {type(key).__name__}")
    key_to_species[key] = species
    key_to_provenance[key] = provenance
    return species, provenance


def _walk_normal(cd_map, mode):
    flatten_complexes = mode in _FLAT_COMPLEX_MODES
    clingo_model, id_to_model_element = pd2af.solver.solve(cd_map, mode=mode)
    activity_atoms = pd2af.solver._get_activity_atoms(clingo_model)
    influence_atoms = pd2af.solver._get_influence_atoms(clingo_model)
    key_to_species = {}
    key_to_provenance = {}
    stripped_template_cache = {}
    new_species_list = []
    species_provenance = []
    for atom in activity_atoms:
        species, provenance = _resolve(
            atom.key,
            id_to_model_element,
            key_to_species,
            key_to_provenance,
            stripped_template_cache,
            flatten_complexes,
        )
        new_species_list.append(species)
        species_provenance.append(provenance)
    immediate_compartments = {
        species.compartment
        for species in new_species_list
        if getattr(species, "compartment", None) is not None
    }
    compartments = collect_ancestor_compartments(immediate_compartments)
    species_templates = set()
    for species in new_species_list:
        template = getattr(species, "template", None)
        if template is not None:
            species_templates.add(template)
        for subunit in getattr(species, "subunits", ()) or ():
            subunit_template = getattr(subunit, "template", None)
            if subunit_template is not None:
                species_templates.add(subunit_template)
    canonical_compartments = {}
    canonical_templates = {}
    canonical_species = {}
    for compartment in compartments_outermost_first(compartments):
        if compartment in canonical_compartments:
            continue
        canonical_compartments[compartment] = compartment
        yield BuildStep(
            kind=BuildStepKind.COMPARTMENT,
            new_element=compartment,
            provenances=(compartment,),
        )
    for template in species_templates:
        if template in canonical_templates:
            continue
        canonical_templates[template] = template
        is_synthesized = template.id_.startswith("merged_template__")
        provenance = () if is_synthesized else (template,)
        yield BuildStep(
            kind=BuildStepKind.TEMPLATE,
            new_element=template,
            provenances=provenance,
        )
    for species, provenance in zip(new_species_list, species_provenance):
        if species in canonical_species:
            continue
        canonical_species[species] = species
        top_step = BuildStep(
            kind=BuildStepKind.SPECIES,
            new_element=species,
            provenances=provenance,
        )
        yield top_step
        # In flat-complex modes, complexes have empty subunits — no nested
        # build steps. In non-flat modes, kept complexes may carry subunits.
        if not flatten_complexes:
            for subunit in getattr(species, "subunits", ()) or ():
                if subunit in canonical_species:
                    continue
                canonical_species[subunit] = subunit
                yield BuildStep(
                    kind=BuildStepKind.SPECIES,
                    new_element=subunit,
                    provenances=(subunit,),
                    parent=top_step,
                )
    def canonicalize(species):
        return canonical_species.get(species, species)

    for atom in influence_atoms:
        cls = pd2af.predicates.predicate_to_model_element_class[type(atom)]
        source, _ = _resolve(
            atom.source,
            id_to_model_element,
            key_to_species,
            key_to_provenance,
            stripped_template_cache,
            flatten_complexes,
        )
        source = canonicalize(source)
        target, _ = _resolve(
            atom.target,
            id_to_model_element,
            key_to_species,
            key_to_provenance,
            stripped_template_cache,
            flatten_complexes,
        )
        target = canonicalize(target)
        new_modulation = cls(source=source, target=target)
        yield BuildStep(
            kind=BuildStepKind.MODULATION,
            new_element=new_modulation,
            provenances=(),
        )


@register("normal")
class NormalWalker(Walker):
    def walk(self, cd_map):
        yield from _walk_normal(cd_map, mode="normal")


@register("no-complex")
class NoComplexWalker(Walker):
    def walk(self, cd_map):
        yield from _walk_normal(cd_map, mode="no-complex")
