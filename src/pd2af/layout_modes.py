"""The layout-mode vocabulary and which mode each output language supports.

:data:`LAYOUT_MODES` is the single source of the concrete layout modes and of
their CLI descriptions; :data:`LAYOUT_MODES_BY_LANGUAGE` is the join of layout
modes with the input languages of :mod:`pd2af.languages`, which this module
imports (never the reverse).
"""

import pd2af.languages


# The concrete layout modes -> their one-line CLI description, in display order.
LAYOUT_MODES = {
    "plain": "reuse original positions",
    "overlay": (
        "reuse full original layout with unmapped layout elements dimmed"
    ),
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
    pd2af.languages.CELLDESIGNER: tuple(LAYOUT_MODES),
    pd2af.languages.SBGN_PD: ("plain", "dot"),
}
