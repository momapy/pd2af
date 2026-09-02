"""Everything a user chooses from: the vocabularies and the transformation modes.

The input languages, the layout modes and the influence pairings are the three
vocabularies a transformation mode is defined against, so they live here with
it: :data:`LANGUAGES` is the single source of the input languages and of the
order the CLI and the docs list them in, :data:`LAYOUT_MODES` of the layout
modes and their CLI descriptions, and :data:`LAYOUT_MODES_BY_LANGUAGE` is the
join of the two. A language token is also the wire format: it is the aspcompose
variant key resolved by :func:`pd2af.asp.rules.build_program`, a segment of
every language-specific rule identifier, and a member of a mode's
``compatible_languages``. The output language is *deduced* from the input:
``celldesigner`` input -> CellDesigner output, ``sbgn_pd`` input -> SBGN-AF
output.

A :class:`TransformationMode` carries everything pd2af knows about a mode:
its name, the prose the CLI lists it with, the rule groups its ASP program
is made of, the input languages it accepts, and whether it strips
post-translational decorations and merges the content-equal results.

Modes are contributed from outside the package through the ``pd2af.modes``
entry-point group. Each entry point resolves to a :class:`TransformationMode`
object::

    # pyproject.toml of the contributing distribution
    [project.entry-points."pd2af.modes"]
    my-mode = "my_package:MY_MODE"

    # my_package/__init__.py
    import pd2af

    MY_MODE = pd2af.TransformationMode(
        name="my-mode",
        docs="what this mode does, in one clause",
        rule_group_references=("influences:kind", "influences:output"),
        rule_group_definitions=MY_RULE_GROUPS,
    )

``rule_group_references`` names groups pd2af already registers;
``rule_group_definitions`` holds ``RuleGroup`` objects the mode brings with it,
which :func:`pd2af.asp.rules.build_registry` registers alongside the built-in ones. A
contributed mode may not take the name of an existing mode, and a mode that
fails to load — a bad import, a wrong type, a name collision, a group that
fails registry validation — takes down every pd2af entry point rather than
just its own: a mode that half-loads is worse than one that refuses to.
"""

import dataclasses
import functools
import importlib.metadata
import types

from aspcompose import RuleGroup

import momapy.celldesigner
import momapy.core.map
import momapy.core.model
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


# The concrete layout modes -> their one-line CLI description, in display order.
LAYOUT_MODES = {
    "plain": "reuse original positions",
    "overlay": ("reuse full original layout with unmapped layout elements dimmed"),
    "dot": "graphviz `dot` auto-layout (requires `dot` on PATH)",
}

# Not a layout mode of its own: a meta value accepted by the CLI and
# `transform`, resolved from the input before the build stage ever sees it.
AUTO = "auto"
AUTO_DESCRIPTION = "pick automatically from the input (graphviz `dot` for a map)"

# The concrete layout modes each output language supports. SBGN-AF output
# supports the curated-geometry `plain` mode and the graphviz `dot` mode; the
# `overlay` dimming is CellDesigner-only.
LAYOUT_MODES_BY_LANGUAGE = {
    CELLDESIGNER: tuple(LAYOUT_MODES),
    SBGN_PD: ("plain", "dot"),
}

# How to draw an influence whose source or target maps to several glyphs.
INFLUENCE_PAIRINGS = ("cross", "nearest")


def get_language_from_map_or_model(
    map_or_model: momapy.core.map.Map | momapy.core.model.Model,
) -> str:
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


def make_map_from_model(
    model: momapy.core.model.Model, language: str
) -> momapy.core.map.Map:
    """Wrap a bare model in a layout-less map of its language.

    A ``Map`` is a keyword-only frozen dataclass whose ``layout`` and
    ``layout_model_mapping`` default to ``None``, so a model-only map is a
    direct construction of the language's ``map_class``. A language token that
    is not in :data:`LANGUAGES` raises ``KeyError``.
    """
    return LANGUAGES[language]["map_class"](model=model)


ENTRY_POINT_GROUP = "pd2af.modes"


@dataclasses.dataclass(frozen=True)
class TransformationMode:
    """One transformation mode: its metadata and the groups it is made of.

    Attributes:
        name: The mode's only identifier, as passed to ``--transformation-mode``
            and :func:`pd2af.transform`.
        docs: One clause of prose describing what the mode does, as listed by
            the CLI.
        rule_group_references: Identifiers of already-registered rule groups
            the mode's program includes. The list is *complete*, not leaf-only:
            dependencies are named explicitly rather than auto-included, so the
            tuple is the mode's full inventory.
        rule_group_definitions: Rule groups the mode defines itself, registered
            with the built-in ones. Empty for every mode pd2af ships.
        compatible_languages: The input-language tokens the mode accepts,
            every registered language by default.
        merges_proteoforms: Whether the build stage strips post-translational
            decorations and merges content-equal results.
    """

    name: str
    docs: str
    rule_group_references: tuple[str, ...] = ()
    rule_group_definitions: tuple[RuleGroup, ...] = ()
    compatible_languages: frozenset[str] = frozenset(LANGUAGES)
    merges_proteoforms: bool = False

    @property
    def rule_group_ids(self) -> tuple[str, ...]:
        """Every group this mode is made of, defined ones included."""
        return self.rule_group_references + tuple(
            rule_group.identifier for rule_group in self.rule_group_definitions
        )

    def compatible_layout_modes(self, language: str) -> tuple[str, ...]:
        """The concrete layout modes valid for this mode on the given input language.

        Merged activities are synthesized from several input species, so they
        have no original geometry for `plain`/`overlay` to reuse.
        """
        available_layout_modes = LAYOUT_MODES_BY_LANGUAGE[language]
        if not self.merges_proteoforms:
            return available_layout_modes
        return ("dot",) if "dot" in available_layout_modes else ()


_BUILTIN_TRANSFORMATION_MODES = (
    TransformationMode(
        name="normal",
        docs=(
            "merge forms of the same base species or entity pool into a single activity"
        ),
        rule_group_references=(
            "activity:core",
            "activity:phenotype",
            "activity:active_marker",
            "activity:modulation_source",
            "activity:gate_input",
            "topology:core",
            "topology:top_level",
            "preparation:complex",
            "paths:core",
            "paths:chaining",
            "influences:kind",
            "influences:core",
            "influences:consumption",
            "influences:output",
            "gates:core",
        ),
        merges_proteoforms=True,
    ),
    TransformationMode(
        name="normal-no-complex",
        docs=(
            "merge forms of the same base species or entity pool into a "
            "single activity; additionally, replace complexes with their "
            "active subunits if any, promoting them to top-level activities"
        ),
        rule_group_references=(
            "activity:core",
            "activity:phenotype",
            "activity:active_marker",
            "activity:modulation_source",
            "activity:gate_input",
            "topology:core",
            "preparation:no_complex",
            "paths:core",
            "paths:chaining",
            "paths:complex_traversal",
            "influences:kind",
            "influences:core",
            "influences:consumption",
            "influences:output",
            "gates:core",
        ),
        merges_proteoforms=True,
    ),
    TransformationMode(
        name="keep-species",
        docs="create one activity per distinct active species or entity pool",
        rule_group_references=(
            "activity:core",
            "activity:phenotype",
            "activity:active_marker",
            "activity:modulation_source",
            "activity:gate_input",
            "topology:core",
            "topology:top_level",
            "preparation:complex",
            "paths:core",
            "paths:chaining",
            "influences:kind",
            "influences:core",
            "influences:consumption",
            "influences:output",
            "gates:core",
        ),
    ),
    TransformationMode(
        name="keep-species-no-complex",
        docs=(
            "create one activity per distinct active species or entity pool; "
            "additionally, replace complexes with their active subunits if "
            "any, promoting them to top-level activities"
        ),
        rule_group_references=(
            "activity:core",
            "activity:phenotype",
            "activity:active_marker",
            "activity:modulation_source",
            "activity:gate_input",
            "topology:core",
            "preparation:no_complex",
            "paths:core",
            "paths:chaining",
            "paths:complex_traversal",
            "influences:kind",
            "influences:core",
            "influences:consumption",
            "influences:output",
            "gates:core",
        ),
    ),
    TransformationMode(
        name="keep-reactions",
        docs=(
            "create one activity per species, and turn every reaction into a "
            "positive influence from each of its reactants to each of its "
            "products"
        ),
        rule_group_references=(
            "activity:core",
            "topology:core",
            "topology:top_level",
            "preparation:complex",
            "paths:core",
            "influences:kind",
            "influences:core",
            "influences:output",
            "gates:core",
            "keep_reactions:activity",
            "keep_reactions:influences",
        ),
        compatible_languages=frozenset({CELLDESIGNER}),
    ),
)


@functools.cache
def get_transformation_modes() -> types.MappingProxyType[str, "TransformationMode"]:
    """Return ``{name: TransformationMode}``, read-only and ordered.

    Built-in modes come first in declaration order, then the modes contributed
    through the ``pd2af.modes`` entry-point group, so ``list-modes`` and the
    ``--transformation-mode`` choices read canonically.
    """
    modes = {mode.name: mode for mode in _BUILTIN_TRANSFORMATION_MODES}
    for entry_point in importlib.metadata.entry_points(group=ENTRY_POINT_GROUP):
        mode = entry_point.load()
        if not isinstance(mode, TransformationMode):
            raise TypeError(
                f"entry point {entry_point.name!r} of group "
                f"{ENTRY_POINT_GROUP!r} ({entry_point.value}) resolved to "
                f"{type(mode).__name__}, expected a TransformationMode"
            )
        if mode.name in modes:
            raise ValueError(
                f"entry point {entry_point.name!r} of group "
                f"{ENTRY_POINT_GROUP!r} contributes transformation mode "
                f"{mode.name!r}, which already exists"
            )
        modes[mode.name] = mode
    return types.MappingProxyType(modes)


def get_transformation_mode(name: str) -> "TransformationMode":
    """Return the :class:`TransformationMode` called ``name``.

    Raises ``ValueError``, listing the available names, if there is none.
    """
    modes = get_transformation_modes()
    mode = modes.get(name)
    if mode is None:
        raise ValueError(
            f"unknown transformation mode {name!r}; available modes: "
            + ", ".join(modes)
        )
    return mode
