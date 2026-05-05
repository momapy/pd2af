import abc
import dataclasses
import typing

import momapy.builder
import momapy.core.layout
import momapy.core.mapping

from pd2af.walkers import BuildStep


@dataclasses.dataclass
class PlaceContext:
    cd_map: typing.Any
    layout_mode: str | None
    new_layout_builder: typing.Any
    new_mapping_builder: typing.Any
    object_to_builder: dict = dataclasses.field(default_factory=dict)
    parent_subunit_index: dict = dataclasses.field(default_factory=dict)


class Placer(abc.ABC):
    @abc.abstractmethod
    def place(self, build_step: BuildStep, context: PlaceContext) -> None: ...

    def make_layout_builders(self, cd_map):
        layout_builder_cls = momapy.builder.get_or_make_builder_cls(
            momapy.core.layout.Layout
        )
        return (
            layout_builder_cls(),
            momapy.core.mapping.LayoutModelMappingBuilder(),
        )

    def finalize(self, context: PlaceContext, new_map_builder) -> None:
        return None

    def post_finalize(self, context: PlaceContext, new_map):
        return new_map


_REGISTRY: dict[str | None, type[Placer]] = {}


def register(layout_mode: str | None):
    def decorate(placer_cls):
        _REGISTRY[layout_mode] = placer_cls
        return placer_cls
    return decorate


def placer_for(layout_mode: str | None) -> Placer | None:
    placer_cls = _REGISTRY.get(layout_mode)
    if placer_cls is None:
        return None
    return placer_cls()
