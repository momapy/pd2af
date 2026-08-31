"""Derive the ASP ontology rules from a language's momapy classes."""

import momapy.core.elements
import momapy.core.model

import pd2af.languages


_BASES = (momapy.core.elements.ModelElement, momapy.core.model.Model)


def _iter_types(module):
    for attr_name in dir(module):
        if attr_name.startswith("_"):
            continue
        attr_value = getattr(module, attr_name)
        if isinstance(attr_value, type) and issubclass(attr_value, _BASES):
            yield attr_value


def make_rules(session, language):
    """The sorted ontology rules for every model class of ``language``."""
    module = pd2af.languages.LANGUAGES[language]["momapy_module"]
    rules = set()
    for type_ in _iter_types(module):
        session.get_or_make_predicate_classes_from_type(
            type_, make_predicate_classes_recursively=True
        )
        rules.update(session.make_ontology_rules_from_type(type_))
    return sorted(rules)
