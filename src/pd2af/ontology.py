import momapy.celldesigner
import momapy.sbgn.pd
import momapy.core.elements
import momapy.core.model


_BASES = (momapy.core.elements.ModelElement, momapy.core.model.Model)

# Input language -> the momapy package whose model-element classes seed the
# momapy_kb ontology. Functor names derive from the class names, so each
# language yields its own predicate vocabulary (no overlap between the two).
_LANGUAGE_MODULES = {
    "celldesigner": momapy.celldesigner,
    "sbgn_pd": momapy.sbgn.pd,
}


def _iter_types(module):
    for attr_name in dir(module):
        if attr_name.startswith("_"):
            continue
        attr_value = getattr(module, attr_name)
        if isinstance(attr_value, type) and issubclass(attr_value, _BASES):
            yield attr_value


def make_rules(session, language):
    module = _LANGUAGE_MODULES[language]
    rules = set()
    for type_ in _iter_types(module):
        session.get_or_make_predicate_classes_from_type(
            type_, make_predicate_classes_recursively=True
        )
        rules.update(session.make_ontology_rules_from_type(type_))
    return sorted(rules)
