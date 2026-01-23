% active state (dashed)
hasActivity(SPECIES, active):-species(SPECIES, _, "true").
% structural state with value "active"
hasActivity(SPECIES, structural_state_active):-species(SPECIES, _, _), hasStructuralState(SPECIES, STRUCTURAL_STATE), structuralState(STRUCTURAL_STATE, "active").
% phenotype
hasActivity(PHENOTYPE, phenotype):-phenotype(PHENOTYPE, _, _).
% source of a modulation
hasActivity(SOURCE, modulates(SOURCE, TARGET)):-modulates(SOURCE, TARGET), species(SOURCE, _, _).
% inhibition by complexation
%hasActivity(Y, inhibits_by_binding(X)):-species(X, _, _), heterodimerAssociation(Z), consumes(Z, X, _), consumes(Z, Y, _), X!=Y, 1{ hasActivity(X, active); hasActivity(X, structural_state_active); hasActivity(X,  modulates(_)); hasActivity(X,  modulated_by(_)); hasActivity(X,  modulates_reaction(_)); hasActivity(X, phenotype)}, {not hasActivity(W, _):produces(Z, W, _)}.

new(activity(ENTITY)):-hasActivity(ENTITY, _).

% modulator to product
posPath(START_SPECIES, END_SPECIES):-reaction(REACTION), catalyzes(START_SPECIES, REACTION), produces(REACTION, END_SPECIES, _), species(START_SPECIES, _, _), species(END_SPECIES, _, _).
posPath(START_SPECIES, END_SPECIES):-reaction(REACTION), positivelyInfluences(START_SPECIES, REACTION), produces(REACTION, END_SPECIES, _), species(START_SPECIES, _, _), species(END_SPECIES, _, _).
posPath(START_SPECIES, END_SPECIES):-reaction(REACTION), physicallyStimulates(START_SPECIES, REACTION), produces(REACTION, END_SPECIES, _), species(START_SPECIES, _, _), species(END_SPECIES, _, _).
posPath(START_SPECIES, END_SPECIES):-reaction(REACTION), triggers(START_SPECIES, REACTION), produces(REACTION, END_SPECIES, _), species(START_SPECIES, _, _), species(END_SPECIES, _, _).

% modulator to target
posPath(START_SPECIES, END_SPECIES):-catalyzes(START_SPECIES, END_SPECIES), species(START_SPECIES, _, _), species(END_SPECIES, _, _).
posPath(START_SPECIES, END_SPECIES):-positivelyInfluences(START_SPECIES, END_SPECIES), species(START_SPECIES, _, _), species(END_SPECIES, _, _).
posPath(START_SPECIES, END_SPECIES):-physicallyStimulates(START_SPECIES, END_SPECIES), species(START_SPECIES, _, _), species(END_SPECIES, _, _).
posPath(START_SPECIES, END_SPECIES):-triggers(START_SPECIES, END_SPECIES), species(START_SPECIES, _, _), species(END_SPECIES, _, _).

% reactant to product (transitivity)
posPath(START_SPECIES, END_SPECIES):-posPath(START_SPECIES, INTERMEDIATE_SPECIES), reaction(REACTION), consumes(REACTION, INTERMEDIATE_SPECIES, _), produces(REACTION, END_SPECIES, _), species(END_SPECIES, _, _).

% modulator to product
negPath(START_SPECIES, END_SPECIES):-reaction(REACTION), inhibits(START_SPECIES, REACTION), produces(REACTION, END_SPECIES, _), species(START_SPECIES, _, _), species(END_SPECIES, _, _).
negPath(START_SPECIES, END_SPECIES):-reaction(REACTION), negativelyInfluences(START_SPECIES, REACTION), produces(REACTION, END_SPECIES, _), species(START_SPECIES, _, _), species(END_SPECIES, _, _).

% modulator to target
negPath(START_SPECIES, END_SPECIES):-inhibits(START_SPECIES, END_SPECIES), species(START_SPECIES, _, _), species(END_SPECIES, _, _).
negPath(START_SPECIES, END_SPECIES):-negativelyInfluences(START_SPECIES, END_SPECIES), species(START_SPECIES, _, _), species(END_SPECIES, _, _).

% reactant to product (transitivity)
negPath(START_SPECIES, END_SPECIES):-negPath(START_SPECIES, INTERMEDIATE_SPECIES), reaction(REACTION), consumes(REACTION, INTERMEDIATE_SPECIES, _), produces(REACTION, END_SPECIES, _), species(END_SPECIES, _, _).

% positive path to modulation
new(positivelyInfluences(START_SPECIES, END_SPECIES)):-posPath(START_SPECIES, END_SPECIES), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _).

% negative path to modulation
new(negativelyInfluences(START_SPECIES, END_SPECIES)):-negPath(START_SPECIES, END_SPECIES), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _).

% modulation of consumption
new(negativelyInfluences(START_SPECIES, END_SPECIES)):-species(START_SPECIES, _, _), species(END_SPECIES, _, _), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _), reaction(REACTION), catalyzes(START_SPECIES, REACTION), consumes(REACTION, END_SPECIES, _).


% TODELETE
% new(negativelyInfluences(START_SPECIES, END_SPECIES)):-species(START_SPECIES, _, _), species(END_SPECIES, _, _), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _), negativelyInfluences(START_SPECIES, END_SPECIES).
% new(negativelyInfluences(START_SPECIES, END_SPECIES)):-species(START_SPECIES, _, _), species(END_SPECIES, _, _), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _), inhibits(START_SPECIES, END_SPECIES).
% new(positivelyInfluences(START_SPECIES, END_SPECIES)):-species(START_SPECIES, _, _), species(END_SPECIES, _, _), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _), catalyzes(START_SPECIES, END_SPECIES).
% new(positivelyInfluences(START_SPECIES, END_SPECIES)):-species(START_SPECIES, _, _), species(END_SPECIES, _, _), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _), positivelyInfluences(START_SPECIES, END_SPECIES).
% new(positivelyInfluences(START_SPECIES, END_SPECIES)):-species(START_SPECIES, _, _), species(END_SPECIES, _, _), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _), physicallyStimulates(START_SPECIES, END_SPECIES).
% new(triggers(START_SPECIES, END_SPECIES)):-species(START_SPECIES, _, _), species(END_SPECIES, _, _), hasActivity(START_SPECIES, _), hasActivity(END_SPECIES, _), triggers(START_SPECIES, END_SPECIES).
