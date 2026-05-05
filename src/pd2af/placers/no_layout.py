from pd2af.placers import PlaceContext, Placer, register
from pd2af.walkers import BuildStep


@register(None)
class NoLayoutPlacer(Placer):
    """Emits no layout. Build steps still flow so the model is built."""

    def make_layout_builders(self, cd_map):
        return None, None

    def place(self, build_step: BuildStep, context: PlaceContext) -> None:
        return None
