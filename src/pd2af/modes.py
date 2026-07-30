"""Transformation modes: what a mode *is*, and how to contribute a new one.

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
which :func:`pd2af.rules.build_registry` registers alongside the built-in ones. A
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

import pd2af.languages
import pd2af.layout_modes


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
    compatible_languages: frozenset[str] = frozenset(pd2af.languages.LANGUAGES)
    merges_proteoforms: bool = False

    @property
    def rule_group_ids(self):
        """Every group this mode is made of, defined ones included."""
        return self.rule_group_references + tuple(
            rule_group.identifier for rule_group in self.rule_group_definitions
        )

    def compatible_layout_modes(self, language):
        """The concrete layout modes valid for this mode on input of the
        given language.

        Merged activities are synthesized from several input species, so they
        have no original geometry for `plain`/`overlay` to reuse.
        """
        available_layout_modes = pd2af.layout_modes.LAYOUT_MODES_BY_LANGUAGE[language]
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
        compatible_languages=frozenset({pd2af.languages.CELLDESIGNER}),
    ),
)


@functools.cache
def get_transformation_modes():
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


def get_transformation_mode(name):
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
