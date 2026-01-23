import clorm

import momapy.celldesigner.core


class species(clorm.Predicate):
    id_: str
    name: str
    active: str


class genericProtein(clorm.Predicate):
    id_: str
    name: str
    active: str


class truncatedProtein(clorm.Predicate):
    id_: str
    name: str
    active: str


class receptor(clorm.Predicate):
    id_: str
    name: str
    active: str


class ionChannel(clorm.Predicate):
    id_: str
    name: str
    active: str


class gene(clorm.Predicate):
    id_: str
    name: str
    active: str


class rna(clorm.Predicate):
    id_: str
    name: str
    active: str


class antisenseRna(clorm.Predicate):
    id_: str
    name: str
    active: str


class phenotype(clorm.Predicate):
    id_: str
    name: str
    active: str


class ion(clorm.Predicate):
    id_: str
    name: str
    active: str


class simpleMolecule(clorm.Predicate):
    id_: str
    name: str
    active: str


class drug(clorm.Predicate):
    id_: str
    name: str
    active: str


class unknown(clorm.Predicate):
    id_: str
    name: str
    active: str


class complex(clorm.Predicate):
    id_: str
    name: str
    active: str


class degraded(clorm.Predicate):
    id_: str
    name: str
    active: str


class structuralState(clorm.Predicate):
    id_: str
    value: str


class hasStructuralState(clorm.Predicate):
    species: str
    structural_state: str


class hasSubunit(clorm.Predicate):
    species: str
    subunit: str


class reaction(clorm.Predicate):
    id_: str


class stateTransition(clorm.Predicate):
    id_: str


class knownTransitionOmitted(clorm.Predicate):
    id_: str


class unknownTransition(clorm.Predicate):
    id_: str


class transcription(clorm.Predicate):
    id_: str


class translation(clorm.Predicate):
    id_: str


class transport(clorm.Predicate):
    id_: str


class heterodimerAssociation(clorm.Predicate):
    id_: str


class dissociation(clorm.Predicate):
    id_: str


class truncation(clorm.Predicate):
    id_: str


class consumes(clorm.Predicate):
    reaction: str
    reactant: str
    stoichiometry: str


class produces(clorm.Predicate):
    reaction: str
    product: str
    stoichiometry: str


class modulates(clorm.Predicate):
    source: str
    target: str


class catalyzes(clorm.Predicate):
    source: str
    target: str


class inhibits(clorm.Predicate):
    source: str
    target: str


class physicallyStimulates(clorm.Predicate):
    source: str
    target: str


class triggers(clorm.Predicate):
    source: str
    target: str


class positivelyInfluences(clorm.Predicate):
    source: str
    target: str


class negativelyInfluences(clorm.Predicate):
    source: str
    target: str


class unknownModulates(clorm.Predicate):
    source: str
    target: str


class unknownCatalyzes(clorm.Predicate):
    source: str
    target: str


class unknownInhibits(clorm.Predicate):
    source: str
    target: str


class unknownPositivelyInfluences(clorm.Predicate):
    source: str
    target: str


class unknownNegativelyInfluences(clorm.Predicate):
    source: str
    target: str


class unknownPhysicallyStimulates(clorm.Predicate):
    source: str
    target: str


class unknownTriggers(clorm.Predicate):
    source: str
    target: str


def cd_model_to_facts(cd_model):
    facts = []
    for cd_species in cd_model.species:
        facts += cd_species_to_facts(cd_species)
    for cd_reaction in cd_model.reactions:
        facts += cd_reaction_to_facts(cd_reaction)
    for cd_modulation in cd_model.modulations:
        facts += cd_modulation_to_facts(cd_modulation)
    return facts


def cd_species_to_facts(cd_species):
    facts = []
    active = "true" if cd_species.active else "false"
    species_predicate = _cd_species_type_to_predicate[type(cd_species)]
    facts.append(species_predicate(cd_species.id_, cd_species.name, active))
    if hasattr(cd_species, "structural_states"):
        for cd_structural_state in cd_species.structural_states:
            facts += cd_structural_state_to_facts(cd_structural_state)
            facts.append(hasStructuralState(cd_species.id_, cd_structural_state.id_))
    # if hasattr(cd_species, "subunits"):
    #     for cd_subunit in cd_species.subunits:
    #         facts += cd_species_to_facts(cd_subunit)
    #         facts.append(hasSubunit(cd_species.id_, cd_subunit.id_))
    return facts


def cd_structural_state_to_facts(cd_structural_state):
    facts = [structuralState(cd_structural_state.id_, cd_structural_state.value)]
    return facts


def cd_reaction_to_facts(cd_reaction):
    facts = []
    reaction_predicate = _cd_reaction_type_to_predicate[type(cd_reaction)]
    facts.append(reaction_predicate(cd_reaction.id_))
    for cd_reactant in cd_reaction.reactants:
        facts.append(
            consumes(
                cd_reaction.id_,
                cd_reactant.referred_species.id_,
                str(cd_reactant.stoichiometry),
            )
        )
    for cd_product in cd_reaction.products:
        facts.append(
            produces(
                cd_reaction.id_,
                cd_product.referred_species.id_,
                str(cd_product.stoichiometry),
            )
        )
    for cd_modifier in cd_reaction.modifiers:
        if isinstance(cd_modifier.referred_species, momapy.celldesigner.core.Species):
            modifier_predicate = _cd_modifier_type_to_predicate[type(cd_modifier)]
            facts.append(
                modifier_predicate(cd_modifier.referred_species.id_, cd_reaction.id_)
            )
    return facts


def cd_modulation_to_facts(cd_modulation):
    facts = []
    if isinstance(cd_modulation.source, momapy.celldesigner.core.Species):
        modulation_predicate = _cd_modulation_type_to_predicate[type(cd_modulation)]
        facts.append(
            modulation_predicate(cd_modulation.source.id_, cd_modulation.target.id_)
        )
    return facts


ontology_rules = [
    "species(X,Y,Z):-protein(X,Y,Z).",
    "protein(X,Y,Z):-genericProtein(X,Y,Z).",
    "protein(X,Y,Z):-truncatedProtein(X,Y,Z).",
    "protein(X,Y,Z):-receptor(X,Y,Z).",
    "protein(X,Y,Z):-ionChannel(X,Y,Z).",
    "species(X,Y,Z):-gene(X,Y,Z).",
    "species(X,Y,Z):-rna(X,Y,Z).",
    "species(X,Y,Z):-antisenseRna(X,Y,Z).",
    "species(X,Y,Z):-phenotype(X,Y,Z).",
    "species(X,Y,Z):-ion(X,Y,Z).",
    "species(X,Y,Z):-simpleMolecule(X,Y,Z).",
    "species(X,Y,Z):-drug(X,Y,Z).",
    "species(X,Y,Z):-unknown(X,Y,Z).",
    "species(X,Y,Z):-complex(X,Y,Z).",
    "species(X,Y,Z):-degraded(X,Y,Z).",
    "modulates(X,Y):-catalyzes(X,Y).",
    "modulates(X,Y):-inhibits(X,Y).",
    "modulates(X,Y):-physicallyStimulates(X,Y).",
    "modulates(X,Y):-triggers(X,Y).",
    "modulates(X,Y):-positivelyInfluences(X,Y).",
    "modulates(X,Y):-negativelyInfluences(X,Y).",
    "unknownModulates(X,Y):-unknownCatalyzes(X,Y).",
    "unknownModulates(X,Y):-unknownPositivelyInfluences(X,Y).",
    "unknownModulates(X,Y):-unknownNegativelyInfluences(X,Y).",
    "unknownModulates(X,Y):-unknownPhysicallyStimulates(X,Y).",
    "unknownModulates(X,Y):-unknownTriggers(X,Y).",
    "reaction(X):-stateTransition(X).",
    "reaction(X):-knownTransitionOmitted(X).",
    "reaction(X):-unknownTransition(X).",
    "reaction(X):-transcription(X).",
    "reaction(X):-translation(X).",
    "reaction(X):-transport(X).",
    "reaction(X):-heterodimerAssociation(X).",
    "reaction(X):-dissociation(X).",
    "reaction(X):-truncation(X).",
]


_cd_species_type_to_predicate = {
    momapy.celldesigner.core.GenericProtein: genericProtein,
    momapy.celldesigner.core.TruncatedProtein: truncatedProtein,
    momapy.celldesigner.core.Receptor: receptor,
    momapy.celldesigner.core.IonChannel: ionChannel,
    momapy.celldesigner.core.Gene: gene,
    momapy.celldesigner.core.RNA: rna,
    momapy.celldesigner.core.AntisenseRNA: antisenseRna,
    momapy.celldesigner.core.Phenotype: phenotype,
    momapy.celldesigner.core.Ion: ion,
    momapy.celldesigner.core.SimpleMolecule: simpleMolecule,
    momapy.celldesigner.core.Drug: drug,
    momapy.celldesigner.core.Unknown: unknown,
    momapy.celldesigner.core.Complex: complex,
    momapy.celldesigner.core.Degraded: degraded,
}

_cd_reaction_type_to_predicate = {
    momapy.celldesigner.core.StateTransition: stateTransition,
    momapy.celldesigner.core.KnownTransitionOmitted: knownTransitionOmitted,
    momapy.celldesigner.core.UnknownTransition: unknownTransition,
    momapy.celldesigner.core.Transcription: transcription,
    momapy.celldesigner.core.Translation: translation,
    momapy.celldesigner.core.Transport: transport,
    momapy.celldesigner.core.HeterodimerAssociation: heterodimerAssociation,
    momapy.celldesigner.core.Dissociation: dissociation,
    momapy.celldesigner.core.Truncation: truncation,
}

_cd_modifier_type_to_predicate = {
    momapy.celldesigner.core.Modulator: modulates,
    momapy.celldesigner.core.Catalyzer: catalyzes,
    momapy.celldesigner.core.Inhibitor: inhibits,
    momapy.celldesigner.core.PhysicalStimulator: physicallyStimulates,
    momapy.celldesigner.core.Trigger: triggers,
    momapy.celldesigner.core.UnknownModulator: unknownModulates,
    momapy.celldesigner.core.UnknownCatalyzer: unknownCatalyzes,
    momapy.celldesigner.core.UnknownInhibitor: unknownInhibits,
}

_cd_modulation_type_to_predicate = {
    momapy.celldesigner.core.Modulation: modulates,
    momapy.celldesigner.core.Catalysis: catalyzes,
    momapy.celldesigner.core.Inhibition: inhibits,
    momapy.celldesigner.core.PhysicalStimulation: physicallyStimulates,
    momapy.celldesigner.core.Triggering: triggers,
    momapy.celldesigner.core.PositiveInfluence: positivelyInfluences,
    momapy.celldesigner.core.NegativeInfluence: negativelyInfluences,
    momapy.celldesigner.core.UnknownModulation: unknownModulates,
    momapy.celldesigner.core.UnknownCatalysis: unknownCatalyzes,
    momapy.celldesigner.core.UnknownInhibition: unknownInhibits,
    momapy.celldesigner.core.UnknownPositiveInfluence: unknownPositivelyInfluences,
    momapy.celldesigner.core.UnknownNegativeInfluence: unknownNegativelyInfluences,
    momapy.celldesigner.core.UnknownPhysicalStimulation: unknownPhysicallyStimulates,
    momapy.celldesigner.core.UnknownTriggering: unknownTriggers,
}
