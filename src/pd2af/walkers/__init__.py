import abc
import dataclasses
import enum
import typing


class BuildStepKind(enum.Enum):
    COMPARTMENT = "compartment"
    TEMPLATE = "template"
    SPECIES = "species"
    MODULATION = "modulation"


@dataclasses.dataclass(frozen=True)
class BuildStep:
    kind: BuildStepKind
    new_element: typing.Any
    provenances: tuple
    parent: typing.Optional["BuildStep"] = None
    extras: dict = dataclasses.field(default_factory=dict)


class Walker(abc.ABC):
    @abc.abstractmethod
    def walk(self, cd_map) -> typing.Iterator[BuildStep]: ...

    def get_new(self, input_element):
        raise NotImplementedError


_REGISTRY: dict[str, type[Walker]] = {}


def register(mode: str):
    def decorate(walker_cls):
        _REGISTRY[mode] = walker_cls
        return walker_cls
    return decorate


def walker_for(mode: str) -> Walker | None:
    walker_cls = _REGISTRY.get(mode)
    if walker_cls is None:
        return None
    return walker_cls()
