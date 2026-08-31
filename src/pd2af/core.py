import dataclasses
import typing

import momapy.celldesigner
import momapy.core.model
import momapy.sbgn.pd
import momapy.utils

import pd2af.annotations
import pd2af.build
import pd2af.languages
import pd2af.modes
import pd2af.solver


@dataclasses.dataclass
class TransformerResult:
    """Result of :func:`transform`: the built AF map plus provenance.

    Attributes:
        obj: The transformed map (``CellDesignerMap`` or ``SBGNAFMap``), or
            the transformed model (``CellDesignerModel`` or ``SBGNAFModel``)
            when the input was a bare model rather than a map.
        provenance: Maps each output AF model element to the ``frozenset`` of
            input PD/CD model elements it derives from -- the "where did this
            come from" direction. Keys are the output species (or activities),
            their subunits at any depth, the gates (or logical operators) and
            the compartments. Several inputs map to one output because the
            merged modes content-dedup their results (many-to-one); the origin
            direction is thus the multi-valued one and is the forward index.
            Use ``.inverse`` (``id(input_element) -> frozenset`` of output
            elements) for the source-to-result direction. Output elements are
            valid keys because they are interned by content, so no two
            content-equal-but-distinct elements coexist.
        element_to_annotations: Maps each output model element to the
            ``frozenset`` of RDF annotations carried from its input source(s),
            keyed the way momapy's writer expects (``None`` unless input
            annotations were passed to :func:`transform`).
        element_to_notes: Maps each output model element to the ``frozenset``
            of notes carried from its input source(s) (``None`` unless input
            notes were passed to :func:`transform`).
    """

    obj: typing.Any = None
    provenance: momapy.utils.FrozenIdentityMultiDict | None = None
    element_to_annotations: dict | None = None
    element_to_notes: dict | None = None


INFLUENCE_PAIRINGS = ("cross", "nearest")


def _normalize_layout_mode(layout_mode):
    if layout_mode == "none":
        return None
    return layout_mode


def _validate_layout_mode(layout_mode, mode, language):
    """Check a concrete layout mode against the transformation mode and language.

    The `None` sentinel means "build no layout at all", so it is always valid;
    every other value must be one the mode accepts on input of this language.
    """
    if layout_mode is None:
        return
    compatible_layout_modes = mode.compatible_layout_modes(language)
    if layout_mode not in compatible_layout_modes:
        raise ValueError(
            f"transformation mode {mode.name!r} on {language!r} input supports "
            f"layout_mode {', '.join(repr(candidate) for candidate in compatible_layout_modes)}"
            f" or None, got {layout_mode!r}"
        )


def _wrap_model_in_map(model, language):
    """Wrap a bare input model in a layout-less map of the matching language.

    ``Map`` is a keyword-only frozen dataclass whose ``layout`` and
    ``layout_model_mapping`` default to ``None``, so a model-only map is a
    direct construction.
    """
    if language == pd2af.languages.SBGN_PD:
        return momapy.sbgn.pd.SBGNPDMap(model=model)
    return momapy.celldesigner.CellDesignerMap(model=model)


def transform(
    map_or_model,
    mode: str = "normal",
    layout_mode: typing.Literal["auto", "dot", "plain", "overlay"] | None = "auto",
    influence_pairing: typing.Literal["cross", "nearest"] = "cross",
    set_active: list[str] | None = None,
    set_inactive: list[str] | None = None,
    set_all_active: bool = False,
    set_all_inactive: bool = False,
    exclude_groups: tuple[str, ...] = (),
    exclude_rules: tuple[str, ...] = (),
    element_to_annotations=None,
    element_to_notes=None,
):
    transformation_mode = pd2af.modes.get_transformation_mode(mode)
    layout_mode = _normalize_layout_mode(layout_mode)
    is_model_input = isinstance(map_or_model, momapy.core.model.Model)
    if is_model_input:
        if layout_mode not in (None, "auto"):
            raise ValueError(
                f"model input supports only layout_mode 'auto' or None, "
                f"got {layout_mode!r}"
            )
        layout_mode = None
        map_ = _wrap_model_in_map(
            map_or_model,
            pd2af.languages.get_language_from_map_or_model(map_or_model),
        )
    else:
        map_ = map_or_model
        if layout_mode == "auto":
            layout_mode = "dot"
    language = pd2af.languages.get_language_from_map_or_model(map_)
    if language not in transformation_mode.compatible_languages:
        raise ValueError(
            f"transformation mode {mode!r} does not support {language!r} "
            f"input; it supports "
            + ", ".join(sorted(transformation_mode.compatible_languages))
        )
    _validate_layout_mode(layout_mode, transformation_mode, language)
    if influence_pairing not in INFLUENCE_PAIRINGS:
        raise ValueError(
            f"influence_pairing must be one of {list(INFLUENCE_PAIRINGS)}, "
            f"got {influence_pairing!r}"
        )
    clingo_model, clingo_id_to_model_element = pd2af.solver.solve(
        map_,
        mode,
        set_active=set_active,
        set_inactive=set_inactive,
        set_all_active=set_all_active,
        set_all_inactive=set_all_inactive,
        exclude_groups=exclude_groups,
        exclude_rules=exclude_rules,
    )
    result = pd2af.build.build_map(
        map_,
        layout_mode,
        clingo_model,
        clingo_id_to_model_element,
        influence_pairing,
        mode=transformation_mode,
    )
    (
        output_element_to_annotations,
        output_element_to_notes,
    ) = pd2af.annotations.carry_annotations_through_provenance(
        result.provenance,
        element_to_annotations,
        element_to_notes,
        input_map=map_,
        output_map=result.obj,
    )
    result.element_to_annotations = output_element_to_annotations
    result.element_to_notes = output_element_to_notes
    if is_model_input:
        return TransformerResult(
            obj=result.obj.model,
            provenance=result.provenance,
            element_to_annotations=output_element_to_annotations,
            element_to_notes=output_element_to_notes,
        )
    return result
