activeParticipates(SPECIES, REACTION):-reaction(REACTION), species(SPECIES, _, _), consumes(REACTION, SPECIES, _).
activeParticipates(SPECIES, REACTION):-reaction(REACTION), species(SPECIES, _, _), modulates(SPECIES, REACTION).
participates(SPECIES, REACTION):-activeParticipates(SPECIES, REACTION).
participates(SPECIES, REACTION):-reaction(REACTION), species(SPECIES, _, _), produces(REACTION, SPECIES, _).

% rule 0
delete(DEGRADED, rule_0):-degraded(DEGRADED, _, _).

% rule 1
delete(RECEPTOR, rule_1):-
  heterodimerAssociation(REACTION), consumes(REACTION, RECEPTOR, _), consumes(REACTION, SPECIES, _),
  receptor(RECEPTOR, _, _),
  #count{ R1 : consumes(REACTION, R1, _) } = 2, % #reactants == 2
  #count{ R2 : participates(RECEPTOR, R2) } = 1, #count{ R3 : participates(SPECIES, R3) } = 1.

% rule 2
% rule is symmetric so will be triggered two times
delete(SPECIES_1, rule_2):-
  heterodimerAssociation(REACTION), consumes(REACTION, SPECIES_1, _), consumes(REACTION, SPECIES_2, _),
  SPECIES_1 != SPECIES_2, not receptor(SPECIES_1, _, _), not receptor(SPECIES_2, _, _),
  #count{ REACTANT : consumes(REACTION, REACTANT, _) } = 2, % #reactants == 2
  #count{ R1 : activeParticipates(SPECIES_1, R1) } = 1, #count{ R2 : activeParticipates(SPECIES_2, R2) } = 1.

% rule 3
delete(REACTANT, rule_3):-
  reaction(REACTION), consumes(REACTION, REACTANT, _), produces(REACTION, PRODUCT, _),
  species(REACTANT, NAME, _), species(PRODUCT, NAME, _),
  not produces(_, REACTANT, _), not modulates(REACTANT, _), not modulates(_, REACTANT),
  #count{ P : produces(REACTION, P, _) } = 1, % #products == 1
  #count{ R : consumes(R, REACTANT, _) } = 1.

% rule 4
delete(REACTANT, rule_4):-
  transport(REACTION), consumes(REACTION, REACTANT, _), produces(REACTION, PRODUCT, _),
  species(REACTANT, NAME, _), species(PRODUCT, NAME, _),
  % #count{ R : consumes(R, REACTANT, _) } = 1,
  #count{ R : activeParticipates(REACTANT, R) } = 1,
  #count{ P : produces(REACTION, P, _) } = 1. % #products == 1

% rule 2 rewiring
% we use new so to ensure only one rules pass
new(produces(REACTION_1, PRODUCT, 1)):-
  delete(SPECIES, rule_2), 
  reaction(REACTION_1), produces(REACTION_1, SPECIES, _),
  reaction(REACTION_2), consumes(REACTION_2, SPECIES, _), produces(REACTION_2, PRODUCT, _).


% rule 4 rewiring
% we use new so to ensure only one rules pass
new(produces(REACTION_1, PRODUCT, N)):-
  delete(SPECIES, rule_4),
  reaction(REACTION_1), produces(REACTION_1, SPECIES, N),
  reaction(REACTION_2), consumes(REACTION_2, SPECIES, _), produces(REACTION_2, PRODUCT, _).


new(activity(SPECIES)):-species(SPECIES, _, _), not delete(SPECIES, _).

new(positivelyInfluences(REACTANT, PRODUCT)):-reaction(REACTION), consumes(REACTION, REACTANT, _), produces(REACTION, PRODUCT, _),
  new(activity(REACTANT)), new(activity(PRODUCT)).
new(positivelyInfluences(REACTANT, PRODUCT)):-reaction(REACTION), consumes(REACTION, REACTANT, _), new(produces(REACTION, PRODUCT, _)),
  new(activity(REACTANT)), new(activity(PRODUCT)).

new(positivelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), catalyzes(MODULATOR, REACTION), produces(REACTION, PRODUCT, _),
  new(activity(MODULATOR)), new(activity(PRODUCT)).
new(positivelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), catalyzes(MODULATOR, REACTION), new(produces(REACTION, PRODUCT, _)),
  new(activity(MODULATOR)), new(activity(PRODUCT)).

new(positivelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), physicallyStimulates(MODULATOR, REACTION), produces(REACTION, PRODUCT, _),
  new(activity(MODULATOR)), new(activity(PRODUCT)).
new(positivelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), physicallyStimulates(MODULATOR, REACTION), new(produces(REACTION, PRODUCT, _)),
  new(activity(MODULATOR)), new(activity(PRODUCT)).

new(positivelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), triggers(MODULATOR, REACTION), produces(REACTION, PRODUCT, _),
  new(activity(MODULATOR)), new(activity(PRODUCT)).
new(positivelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), triggers(MODULATOR, REACTION), new(produces(REACTION, PRODUCT, _)),
  new(activity(MODULATOR)), new(activity(PRODUCT)).

new(positivelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), positivelyInfluences(MODULATOR, REACTION), produces(REACTION, PRODUCT, _),
  new(activity(MODULATOR)), new(activity(PRODUCT)).
new(positivelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), positivelyInfluences(MODULATOR, REACTION), new(produces(REACTION, PRODUCT, _)),
  new(activity(MODULATOR)), new(activity(PRODUCT)).

new(negativelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), inhibits(MODULATOR, REACTION), produces(REACTION, PRODUCT, _),
  new(activity(MODULATOR)), new(activity(PRODUCT)).
new(negativelyInfluences(MODULATOR, PRODUCT)):-reaction(REACTION), inhibits(MODULATOR, REACTION), new(produces(REACTION, PRODUCT, _)),
  new(activity(MODULATOR)), new(activity(PRODUCT)).

new(positivelyInfluences(SOURCE, TARGET)):- new(activity(SOURCE)), new(activity(TARGET)),
  catalyzes(SOURCE, TARGET).
new(positivelyInfluences(SOURCE, TARGET)):- new(activity(SOURCE)), new(activity(TARGET)),
  physicallyStimulates(SOURCE, TARGET).
new(positivelyInfluences(SOURCE, TARGET)):- new(activity(SOURCE)), new(activity(TARGET)),
  triggers(SOURCE, TARGET).
new(positivelyInfluences(SOURCE, TARGET)):- new(activity(SOURCE)), new(activity(TARGET)),
  positivelyInfluences(SOURCE, TARGET).
new(negativelyInfluences(SOURCE, TARGET)):- new(activity(SOURCE)), new(activity(TARGET)),
  negativelyInfluences(SOURCE, TARGET).
new(negativelyInfluences(SOURCE, TARGET)):- new(activity(SOURCE)), new(activity(TARGET)),
  inhibits(SOURCE, TARGET).
