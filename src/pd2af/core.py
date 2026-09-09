"""The transformation itself: check the arguments, solve, build the map.

:func:`transform` validates its arguments against the transformation mode,
solves the ASP program (:mod:`pd2af.asp.solver`) and hands the answer to
:func:`_build_map`, which runs the two build passes. The model pass walks the
clingo activity and influence atoms and fills ``context.model`` with canonical,
content-deduped elements; the layout pass -- skipped when ``layout_mode is
None`` -- fills ``context.layout`` and ``context.layout_model_mapping``.
:func:`_build_map` then assembles the map from those three slots and returns it
with the provenance mapping and the carried metadata side-tables.
"""

import dataclasses
import typing

import clorm

import momapy.builder
import momapy.celldesigner
import momapy.core.map
import momapy.core.model
import momapy.sbgn.af
import momapy.sbgn.pd
import momapy.utils

import pd2af.asp.solver
import pd2af.building.celldesigner.layout
import pd2af.building.celldesigner.model
import pd2af.building.context
import pd2af.building.layout
import pd2af.building.provenance
import pd2af.building.sbgn.layout
import pd2af.building.sbgn.model
import pd2af.modes


# What building an output map of each input language takes: the two build-pass
# modules, the map class the two build slots are assembled into, and the
# keyword arguments the `dot` auto-layout needs. SBGN-AF logical operators have
# input/output connectors, so their arcs are ranked upstream
# (`reversed_arc_classes`) and re-attached to the connector tips after graphviz
# repositions (`operator_arc_resolver`); both hooks are no-ops on operator-free
# maps. Keyed by language token, so an unknown token raises `KeyError` instead
# of falling through to a language that was not asked for.
_BUILD_BEHAVIOR_BY_LANGUAGE = {
    pd2af.modes.Language.CELLDESIGNER: {
        "model_module": pd2af.building.celldesigner.model,
        "layout_module": pd2af.building.celldesigner.layout,
        "map_class": momapy.celldesigner.CellDesignerMap,
        "auto_layout_arguments": {},
    },
    pd2af.modes.Language.SBGN_PD: {
        "model_module": pd2af.building.sbgn.model,
        "layout_module": pd2af.building.sbgn.layout,
        "map_class": momapy.sbgn.af.SBGNAFMap,
        "auto_layout_arguments": {
            "compartment_layout_classes": (momapy.sbgn.af.CompartmentLayout,),
            "reversed_arc_classes": (momapy.sbgn.af.LogicArcLayout,),
            "operator_arc_resolver": (
                pd2af.building.sbgn.layout.resolve_operator_arc_segments
            ),
        },
    },
}


@dataclasses.dataclass
class TransformerResult:
    """Result of :func:`transform`: the built AF map plus provenance.

    Attributes:
        obj: The transformed map (``CellDesignerMap`` or ``SBGNAFMap``), or
            the transformed model (``CellDesignerModel`` or ``SBGNAFModel``)
            when the input was a bare model rather than a map. Named after
            momapy's ``ReaderResult.obj``, which it mirrors.
        output_element_to_input_elements: Maps each output AF model element to
            the ``frozenset`` of input PD/CD model elements it derives from --
            the "where did this come from" direction. Keys are the output
            species (or activities), their subunits at any depth, the gates (or
            logical operators) and the compartments. Several inputs map to one
            output because the merged modes content-dedup their results
            (many-to-one); the origin direction is thus the multi-valued one and
            is the forward index. Use ``.inverse`` (``id(input_element) ->
            frozenset`` of output elements) for the source-to-result direction.
            Output elements are valid keys because they are interned by
            content, so no two content-equal-but-distinct elements coexist.
        element_to_annotations: Maps each output model element to the
            ``frozenset`` of RDF annotations carried from its input source(s),
            keyed the way momapy's writer expects (``None`` unless input
            annotations were passed to :func:`transform`).
        element_to_notes: Maps each output model element to the ``frozenset``
            of notes carried from its input source(s) (``None`` unless input
            notes were passed to :func:`transform`).
    """

    obj: momapy.core.map.Map | momapy.core.model.Model | None = None
    output_element_to_input_elements: momapy.utils.FrozenIdentityMultiDict | None = None
    element_to_annotations: dict | None = None
    element_to_notes: dict | None = None


def _format_token(value: typing.Any) -> str:
    """Quote a mode or language value as the token a user passes, not as an enum."""
    return repr(str(value))


def _check_layout_mode_is_supported_by_mode(
    layout_mode: pd2af.modes.LayoutMode | None,
    mode: pd2af.modes.TransformationMode,
    language: pd2af.modes.Language,
    no_compartment: bool = False,
) -> None:
    """Raise unless a layout mode fits the transformation mode and language.

    The `None` sentinel means "build no layout at all", so it is always valid;
    every other value must be one the mode accepts on input of this language.
    `no_compartment` narrows what the mode accepts and is named in the message,
    so the error says which of the two refused the layout mode.
    """
    if layout_mode is None:
        return
    compatible_layout_modes = mode.compatible_layout_modes(language, no_compartment)
    if layout_mode in compatible_layout_modes:
        return
    description = f"transformation mode {mode.name!r}"
    if no_compartment:
        description += " with no_compartment"
    if not compatible_layout_modes:
        raise ValueError(
            f"{description} on {_format_token(language)} input "
            f"supports no layout_mode other than None, got "
            f"{_format_token(layout_mode)}"
        )
    raise ValueError(
        f"{description} on {_format_token(language)} input "
        f"supports layout_mode "
        f"{', '.join(_format_token(candidate) for candidate in compatible_layout_modes)}"
        f" or None, got {_format_token(layout_mode)}"
    )


def _make_model_and_layout_in_context(
    context: pd2af.building.context.BuilderContext,
    clingo_model: clorm.FactBase,
    language: pd2af.modes.Language,
) -> None:
    """Run the two build passes of ``language``, filling the context's slots."""
    build_behavior = _BUILD_BEHAVIOR_BY_LANGUAGE[language]
    build_behavior["model_module"].make_and_add_model(context, clingo_model)
    if context.layout_mode is not None:
        build_behavior["layout_module"].make_and_add_layout(context)


def _make_map_from_context(
    context: pd2af.building.context.BuilderContext, language: pd2af.modes.Language
) -> momapy.core.map.Map:
    """Assemble the map of ``language`` from the three filled context slots."""
    map_builder_class = momapy.builder.get_or_make_builder_cls(
        _BUILD_BEHAVIOR_BY_LANGUAGE[language]["map_class"]
    )
    map_builder = map_builder_class(
        model=context.model,
        layout=context.layout,
        layout_model_mapping=context.layout_model_mapping,
    )
    return momapy.builder.object_from_builder(map_builder)


def _build_map(
    input_map: momapy.core.map.Map,
    clingo_model: clorm.FactBase,
    clingo_id_to_model_element: dict,
    *,
    mode: pd2af.modes.TransformationMode,
    language: pd2af.modes.Language,
    layout_mode: pd2af.modes.LayoutMode | None,
    influence_pairing: pd2af.modes.InfluencePairingMode = (
        pd2af.modes.InfluencePairingMode.CROSS
    ),
    no_compartment: bool = False,
    element_to_annotations: dict | None = None,
    element_to_notes: dict | None = None,
) -> TransformerResult:
    """Run the model and layout passes and return the assembled result.

    Args:
        input_map: The input map the clingo atoms were derived from.
        clingo_model: The solved atoms, as returned by
            :func:`pd2af.asp.solver.solve`.
        clingo_id_to_model_element: The ``clingo_id -> model_element`` registry
            the build passes resolve the atoms' keys through.
        mode: The transformation mode being built.
        language: The input language token, as returned by
            :func:`pd2af.modes.get_language_from_map_or_model`.
        layout_mode: The concrete layout mode, or ``None`` to build no layout.
        influence_pairing: How to draw an influence whose source or target maps
            to several glyphs.
        no_compartment: Whether to merge every compartment into the default one.
        element_to_annotations: The reader's ``element -> annotations``
            side-table, to be carried onto the output elements.
        element_to_notes: The reader's ``element -> notes`` side-table, to be
            carried onto the output elements.

    Returns:
        A :class:`TransformerResult` holding the built map, the provenance
        mapping and the carried annotation and note side-tables.
    """
    context = pd2af.building.context.BuilderContext(
        input_map=input_map,
        layout_mode=layout_mode,
        clingo_id_to_model_element=clingo_id_to_model_element,
        influence_pairing=influence_pairing,
        no_compartment=no_compartment,
        mode=mode,
    )
    _make_model_and_layout_in_context(context, clingo_model, language)
    new_map = _make_map_from_context(context, language)
    # the "dot" mode (graphviz) repositions an already-built layout. The
    # compartment-layout classes differ per language (see
    # pd2af.building.layout.make_auto_layout).
    if layout_mode == pd2af.modes.LayoutMode.DOT:
        new_map = pd2af.building.layout.make_auto_layout(
            new_map,
            **_BUILD_BEHAVIOR_BY_LANGUAGE[language]["auto_layout_arguments"],
        )
    output_element_to_input_elements = (
        pd2af.building.provenance.make_provenance_from_context(context)
    )
    (
        output_element_to_annotations,
        output_element_to_notes,
    ) = pd2af.building.provenance.carry_annotations_through_provenance(
        output_element_to_input_elements,
        element_to_annotations,
        element_to_notes,
        input_map=input_map,
        output_map=new_map,
    )
    return TransformerResult(
        obj=new_map,
        output_element_to_input_elements=output_element_to_input_elements,
        element_to_annotations=output_element_to_annotations,
        element_to_notes=output_element_to_notes,
    )


def transform(
    map_or_model: momapy.core.map.Map | momapy.core.model.Model,
    mode: str = "normal",
    layout_mode: pd2af.modes.LayoutMode | typing.Literal["auto"] | None = "auto",
    influence_pairing: pd2af.modes.InfluencePairingMode = (
        pd2af.modes.InfluencePairingMode.CROSS
    ),
    no_compartment: bool = False,
    set_active: list[str] | None = None,
    set_inactive: list[str] | None = None,
    set_all_active: bool = False,
    set_all_inactive: bool = False,
    exclude_groups: tuple[str, ...] = (),
    exclude_rules: tuple[str, ...] = (),
    element_to_annotations: dict | None = None,
    element_to_notes: dict | None = None,
) -> TransformerResult:
    """Transform a process-description map or model into an activity-flow one.

    Args:
        map_or_model: The input CellDesigner or SBGN-PD map, or a bare model of
            either language. A bare model has no geometry, so ``layout_mode``
            must be ``"auto"`` or ``None`` and no layout is built.
        mode: Name of the transformation mode; see
            :func:`pd2af.modes.get_transformation_modes`.
        layout_mode: How to lay the output out: ``"dot"`` (graphviz),
            ``"plain"`` (reuse the original positions), ``"overlay"`` (reuse the
            full original layout, dimming what is not in the model), ``None``
            (no layout), or ``"auto"`` to pick from the input.
        influence_pairing: How to draw an influence whose source or target maps
            to several glyphs: ``"cross"`` (one arc per pair) or ``"nearest"``
            (a single arc between the closest pair).
        no_compartment: Merge every compartment into the default one, so that
            species differing only by compartment become a single activity and
            the influences that become equal merge in turn. Works with every
            mode; because activities merge, it requires ``layout_mode`` ``"dot"``
            (or ``"auto"``) or ``None``.
        set_active: Ids of elements to surface as activities whatever the map's
            structural signals say.
        set_inactive: Ids of elements to suppress, overriding the automatic
            activity discovery. Wins over ``set_all_active``; an id passed to
            both ``set_active`` and ``set_inactive`` is an error.
        set_all_active: Make every top-level species / entity pool an activity,
            subunits excluded.
        set_all_inactive: Suppress activity for every element, subunits
            included.
        exclude_groups: Identifiers of rule groups to drop from the program.
        exclude_rules: Identifiers of individual rules to drop.
        element_to_annotations: The reader's ``element -> annotations``
            side-table, to be carried onto the output elements.
        element_to_notes: The reader's ``element -> notes`` side-table, to be
            carried onto the output elements.

    Returns:
        A :class:`TransformerResult` holding the output map (or model, for
        model input), the provenance mapping, and the carried annotation and
        note side-tables.
    """
    transformation_mode = pd2af.modes.get_transformation_mode(mode)
    language = pd2af.modes.get_language_from_map_or_model(map_or_model)
    influence_pairing = pd2af.modes.InfluencePairingMode(influence_pairing)
    is_model_input = isinstance(map_or_model, momapy.core.model.Model)
    if is_model_input:
        if layout_mode not in (None, pd2af.modes.AUTO):
            raise ValueError(
                f"model input supports only layout_mode 'auto' or None, "
                f"got {_format_token(layout_mode)}"
            )
        layout_mode = None
        input_map = pd2af.modes.make_map_from_model(map_or_model, language)
    else:
        input_map = map_or_model
        if layout_mode == pd2af.modes.AUTO:
            layout_mode = pd2af.modes.LayoutMode.DOT
        elif layout_mode is not None:
            layout_mode = pd2af.modes.LayoutMode(layout_mode)
    if language not in transformation_mode.compatible_languages:
        raise ValueError(
            f"transformation mode {mode!r} does not support "
            f"{_format_token(language)} input; it supports "
            + ", ".join(
                _format_token(candidate)
                for candidate in sorted(transformation_mode.compatible_languages)
            )
        )
    _check_layout_mode_is_supported_by_mode(
        layout_mode, transformation_mode, language, no_compartment
    )
    clingo_model, clingo_id_to_model_element = pd2af.asp.solver.solve(
        input_map,
        transformation_mode,
        set_active=set_active,
        set_inactive=set_inactive,
        set_all_active=set_all_active,
        set_all_inactive=set_all_inactive,
        exclude_groups=exclude_groups,
        exclude_rules=exclude_rules,
    )
    result = _build_map(
        input_map,
        clingo_model,
        clingo_id_to_model_element,
        mode=transformation_mode,
        language=language,
        layout_mode=layout_mode,
        influence_pairing=influence_pairing,
        no_compartment=no_compartment,
        element_to_annotations=element_to_annotations,
        element_to_notes=element_to_notes,
    )
    if is_model_input:
        return dataclasses.replace(result, obj=result.obj.model)
    return result
