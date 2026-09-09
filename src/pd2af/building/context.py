"""The shared state of the two build passes.

``BuilderContext`` is a plain dataclass of slots: :func:`pd2af.core._build_map`
creates one, the model pass fills the model and the pass-1 scratch slots, and
the layout pass reads them back to build the layout and the layout-model
mapping. It lives in its own module so both language builders can name it in
their annotations without importing the orchestrator that imports them.
"""

import dataclasses

import pd2af.modes


@dataclasses.dataclass
class BuilderContext:
    """The slots the model pass fills and the layout pass reads back."""

    # --- inputs ---
    input_map: object
    layout_mode: pd2af.modes.LayoutMode | None
    clingo_id_to_model_element: dict
    mode: pd2af.modes.TransformationMode
    influence_pairing: pd2af.modes.InfluencePairingMode = (
        pd2af.modes.InfluencePairingMode.CROSS
    )
    no_compartment: bool = False

    # --- outputs being built ---
    model: object = None
    layout: object = None
    layout_model_mapping: object = None

    # --- Pass-1 -> Pass-2 handoff ---
    activity_emissions: list = dataclasses.field(default_factory=list)
    input_model_element_to_canonical_model_element: dict = dataclasses.field(
        default_factory=dict
    )
    # (input_compartment, output_compartment) pairs, feeding the provenance
    # mapping so compartment annotations/notes carry to their AF compartment.
    compartment_emissions: list = dataclasses.field(default_factory=list)

    # --- Pass-1 scratch ---
    cache: dict = dataclasses.field(default_factory=dict)
    # The compartment every activity is put in when `no_compartment` is on: the
    # input map's default compartment for CellDesigner, `None` for SBGN-AF.
    default_compartment: object = None
    subunit_to_top_level: dict = None
    activity_atoms: list = dataclasses.field(default_factory=list)
    influence_atoms: list = dataclasses.field(default_factory=list)
    key_to_activity: dict = dataclasses.field(default_factory=dict)

    # --- Logical-operator scratch ---
    # The activity / operator key maps and emission lists hold the CellDesigner
    # species and gates or the SBGN-AF activities and logical operators,
    # depending on the language being built.
    operator_atoms: list = dataclasses.field(default_factory=list)
    operator_input_atoms: list = dataclasses.field(default_factory=list)
    key_to_operator: dict = dataclasses.field(default_factory=dict)
    operator_emissions: list = dataclasses.field(default_factory=list)

    # --- Pass-2 scratch ---
    model_element_to_layout_elements: dict = dataclasses.field(default_factory=dict)
    object_to_builder: dict = dataclasses.field(default_factory=dict)
    synthetic_index: int = 0

    # --- SBGN-AF pass scratch ---
    input_compartment_to_af_compartment: dict = dataclasses.field(default_factory=dict)
    af_compartment_to_input_compartment: dict = dataclasses.field(default_factory=dict)
