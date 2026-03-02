% active state (dashed)
hasActivity(SPECIES, active):-species(SPECIES), hasActive(SPECIES, 1).
% structural state with value "active"
hasActivity(SPECIES, structural_state_active):-species(SPECIES), hasStructuralState(SPECIES, STRUCTURAL_STATE), hasValue(STRUCTURAL_STATE, "active").
% phenotype
hasActivity(PHENOTYPE, phenotype):-phenotype(PHENOTYPE).
% source of a modulation
hasActivity(SOURCE, modulates(SOURCE, TARGET)):-species(SOURCE), modulation(MODULATION), hasSource(MODULATION, SOURCE), hasTarget(MODULATION, TARGET).
hasActivity(SOURCE, modulates(SOURCE, TARGET)):-species(SOURCE), modulator(MODULATOR), hasReferredSpecies(MODULATOR, SOURCE), hasModifier(TARGET, MODULATOR).
% inhibition by complexation
%hasActivity(Y, inhibits_by_binding(X)):-species(X, _, _), heterodimerAssociation(Z), consumes(Z, X, _), consumes(Z, Y, _), X!=Y, 1{ hasActivity(X, active); hasActivity(X, structural_state_active); hasActivity(X,  modulates(_)); hasActivity(X,  modulated_by(_)); hasActivity(X,  modulates_reaction(_)); hasActivity(X, phenotype)}, {not hasActivity(W, _):produces(Z, W, _)}.

hasActiveSubunit(SPECIES):-complex(SPECIES), hasActivity(SPECIES, _), hasSubunit(SPECIES, SUBUNIT), hasActivity(SUBUNIT, _).
new(activity(SPECIES)):-hasActivity(SPECIES, _), not hasActiveSubunit(SPECIES).

% positive modulator to product
path(START_SPECIES, END_SPECIES, positive):-reaction(REACTION), hasModifier(REACTION, MODIFIER),
                                      catalyzer(MODIFIER), hasReferredSpecies(MODIFIER, START_SPECIES),
                                      hasProduct(REACTION, PRODUCT), hasReferredSpecies(PRODUCT, END_SPECIES).
path(START_SPECIES, END_SPECIES, positive):-reaction(REACTION), hasModifier(REACTION, MODIFIER),
                                      physicalStimulator(MODIFIER), hasReferredSpecies(MODIFIER, START_SPECIES),
                                      hasProduct(REACTION, PRODUCT), hasReferredSpecies(PRODUCT, END_SPECIES).
path(START_SPECIES, END_SPECIES, positive):-reaction(REACTION), hasModifier(REACTION, MODIFIER),
                                      trigger(MODIFIER), hasReferredSpecies(MODIFIER, START_SPECIES),
                                      hasProduct(REACTION, PRODUCT), hasReferredSpecies(PRODUCT, END_SPECIES).

% positive modulator to target
path(START_SPECIES, END_SPECIES, positive):-catalyzis(MODULATION), hasSource(MODULATION, START_SPECIES), hasTarget(MODULATION, END_SPECIES).
path(START_SPECIES, END_SPECIES, positive):-positiveInfluence(MODULATION), hasSource(MODULATION, START_SPECIES), hasTarget(MODULATION, END_SPECIES).
path(START_SPECIES, END_SPECIES, positive):-physicalStimulation(MODULATION), hasSource(MODULATION, START_SPECIES), hasTarget(MODULATION, END_SPECIES).
path(START_SPECIES, END_SPECIES, positive):-triggering(MODULATION), hasSource(MODULATION, START_SPECIES), hasTarget(MODULATION, END_SPECIES).

% negative modulator to product
path(START_SPECIES, END_SPECIES, negative):-reaction(REACTION), hasModifier(REACTION, MODIFIER),
                                      inhibitor(MODIFIER), hasReferredSpecies(MODIFIER, START_SPECIES),
                                      hasProduct(REACTION, PRODUCT), hasReferredSpecies(PRODUCT, END_SPECIES).

% negative modulator to target
path(START_SPECIES, END_SPECIES, negative):-inhibition(MODULATION), hasSource(MODULATION, START_SPECIES), hasTarget(MODULATION, END_SPECIES).
path(START_SPECIES, END_SPECIES, negative):-negativeInfluence(MODULATION), hasSource(MODULATION, START_SPECIES), hasTarget(MODULATION, END_SPECIES).

% reactant to product (transitivity)
path(START_SPECIES, END_SPECIES, SIGN):-path(START_SPECIES, INTERMEDIATE_SPECIES, SIGN), reaction(REACTION),
                                      hasReactant(REACTION, REACTANT), hasReferredSpecies(REACTANT, INTERMEDIATE_SPECIES),
                                      hasProduct(REACTION, PRODUCT), hasReferredSpecies(PRODUCT, END_SPECIES).

% complex to subunit
path(START_SPECIES, SUBUNIT, SIGN):-path(START_SPECIES, END_SPECIES, SIGN), complex(END_SPECIES), hasSubunit(END_SPECIES, SUBUNIT).
path(SUBUNIT, END_SPECIES, SIGN):-path(START_SPECIES, END_SPECIES, SIGN), complex(START_SPECIES), hasSubunit(START_SPECIES, SUBUNIT).

% positive path to modulation
new(positivelyInfluences(START_SPECIES, END_SPECIES)):-path(START_SPECIES, END_SPECIES, positive), new(activity(START_SPECIES)), new(activity(END_SPECIES)).

% negative path to modulation
new(negativelyInfluences(START_SPECIES, END_SPECIES)):-path(START_SPECIES, END_SPECIES, negative), new(activity(START_SPECIES)), new(activity(END_SPECIES)).

% modulation of consumption
new(negativelyInfluences(START_SPECIES, END_SPECIES)):-reaction(REACTION),
                                                        hasModifier(REACTION, MODIFIER), catalyzer(MODIFIER), hasReferredSpecies(MODIFIER, START_SPECIES),
                                                        hasReactant(REACTION, REACTANT), hasReferredSpecies(REACTANT, END_SPECIES),
                                                        new(activity(START_SPECIES)), new(activity(END_SPECIES)).
new(negativelyInfluences(START_SPECIES, END_SPECIES)):-reaction(REACTION),
                                                        hasModifier(REACTION, MODIFIER), physicalStimulator(MODIFIER), hasReferredSpecies(MODIFIER, START_SPECIES),
                                                        hasReactant(REACTION, REACTANT), hasReferredSpecies(REACTANT, END_SPECIES),
                                                        new(activity(START_SPECIES)), new(activity(END_SPECIES)).
