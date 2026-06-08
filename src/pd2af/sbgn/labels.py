"""Build the label of an SBGN-AF activity from an SBGN-PD entity pool.

SBGN-AF activities are opaque, so two distinct PD proteoforms (same type and
name, different state) would collapse under content-based model equality unless
their distinguishing fields are folded into the activity's content. We do that
through the activity *label*, built from the entity's units of information,
state variables and name:

    [unit|of|information][state@variables]name

Rules:

- The entity **type** and its **compartment** are not part of the label (the
  type is carried by the activity's unit-of-information glyph).
- **Units of information** are rendered ``prefix:value`` (or ``value`` with no
  prefix) and sorted (they carry no order).
- **State variables** are rendered ``value@variable`` (``@variable`` with no
  value, ``value`` with no variable, empty when neither) and listed in their
  ``order`` so multi-site proteoforms stay positionally distinct without
  showing the order number. An all-blank state-variable set renders no bracket.
- A complex's **subunits** are *not* part of its name when the complex has its
  own label; for an **unlabelled** complex the name is its subunits, built
  recursively and joined by ``:``.

There are no spaces between elements.
"""


def _state_variable_token(state_variable):
    if state_variable.variable is not None:
        if state_variable.value:
            return f"{state_variable.value}@{state_variable.variable}"
        return f"@{state_variable.variable}"
    return state_variable.value or ""


def _state_variable_sort_key(state_variable):
    return (
        state_variable.order if state_variable.order is not None else 0,
        state_variable.value or "",
        state_variable.variable or "",
    )


def _unit_of_information_token(unit_of_information):
    if unit_of_information.prefix:
        return f"{unit_of_information.prefix}:{unit_of_information.value}"
    return unit_of_information.value


def build_label(entity, include_state_variables=True):
    """Return the SBGN-AF activity label for an SBGN-PD ``entity`` (entity pool
    or subunit), per the rules in the module docstring.

    When ``include_state_variables`` is ``False`` the state-variable block is
    omitted so that distinct proteoforms (same type/name/units, different state)
    collapse into a single merged activity. Used by the merged ``normal`` /
    ``no-complex`` modes; the keep-species modes keep the default ``True``.
    """
    decorations = ""
    units_of_information = getattr(entity, "units_of_information", None)
    if units_of_information:
        decorations += (
            "["
            + "|".join(
                sorted(
                    _unit_of_information_token(unit)
                    for unit in units_of_information
                )
            )
            + "]"
        )
    state_variables = getattr(entity, "state_variables", None)
    if include_state_variables and state_variables:
        tokens = [
            _state_variable_token(state_variable)
            for state_variable in sorted(
                state_variables, key=_state_variable_sort_key
            )
        ]
        if any(tokens):  # omit an all-blank state-variable bracket
            decorations += "[" + "|".join(tokens) + "]"
    name = getattr(entity, "label", None) or ""
    if not name:
        subunits = getattr(entity, "subunits", None)
        if subunits:
            name = ":".join(
                sorted(
                    build_label(subunit, include_state_variables)
                    for subunit in subunits
                )
            )
    if decorations and name:
        return f"{decorations}{name}"
    return decorations or name
