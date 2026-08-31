"""Input-language tokens and inference from the input map type.

:data:`LANGUAGES` is the single source of the language set and of the order the
CLI and the docs list them in: its keyset *is* the set of input languages, so
adding one is a single literal edit here.

The input language drives both the ontology vocabulary
(:mod:`pd2af.ontology`) and which rule variant is resolved
(:mod:`pd2af.rules`). The output language is *deduced* from the input:
``celldesigner`` input -> CellDesigner output, ``sbgn_pd`` input ->
SBGN-AF output.

The token is the wire format: it is the aspcompose variant key resolved by
:func:`pd2af.rules.build_program`, a segment of every language-specific rule
identifier, and a member of a mode's ``compatible_languages``.

The layout-mode vocabulary and the concrete layout modes each language supports
live in :mod:`pd2af.layout_modes`, which imports this module (never the
reverse).
"""

import typing

import momapy.celldesigner
import momapy.sbgn.pd


CELLDESIGNER = "celldesigner"
SBGN_PD = "sbgn_pd"


# Every input language, in the order the CLI and the docs list them: the token
# -> everything pd2af knows about it. `momapy_module` is the package whose
# model-element classes seed the momapy_kb ontology (functor names derive from
# the class names, so each language yields its own predicate vocabulary, with no
# overlap between the two); `map_class` and `model_class` are what an input of
# that language is recognised by.
LANGUAGES = {
    CELLDESIGNER: {
        "display_name": "CellDesigner",
        "momapy_module": momapy.celldesigner,
        "map_class": momapy.celldesigner.CellDesignerMap,
        "model_class": momapy.celldesigner.CellDesignerModel,
    },
    SBGN_PD: {
        "display_name": "SBGN PD",
        "momapy_module": momapy.sbgn.pd,
        "map_class": momapy.sbgn.pd.SBGNPDMap,
        "model_class": momapy.sbgn.pd.SBGNPDModel,
    },
}


def get_language_from_map_or_model(map_or_model: typing.Any) -> str:
    """Infer the input language token from an input map's or model's type.

    Matching is by ``isinstance`` rather than exact type: momapy's builder
    classes are subclasses of the map classes, so a builder is accepted too.
    """
    for language, properties in LANGUAGES.items():
        if isinstance(
            map_or_model, (properties["map_class"], properties["model_class"])
        ):
            return language
    raise ValueError(
        f"unsupported input type {type(map_or_model).__name__!r}; expected "
        + " or ".join(
            momapy_class.__name__
            for properties in LANGUAGES.values()
            for momapy_class in (
                properties["map_class"],
                properties["model_class"],
            )
        )
    )
