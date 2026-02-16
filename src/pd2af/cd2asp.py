import momapy_kb.clingo.core


def cd_model_to_facts(cd_model):
    return momapy_kb.clingo.core.make_facts_from_object(cd_model)
