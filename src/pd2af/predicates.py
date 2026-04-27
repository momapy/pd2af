import clorm

import momapy.celldesigner


class activity(clorm.Predicate):
    name: clorm.ConstantStr


class positivelyInfluences(clorm.Predicate):
    source: clorm.ConstantStr
    target: clorm.ConstantStr


class negativelyInfluences(clorm.Predicate):
    source: clorm.ConstantStr
    target: clorm.ConstantStr


class new(clorm.Predicate):
    object_: activity | positivelyInfluences | negativelyInfluences


predicate_to_model_element_class = {
    positivelyInfluences: momapy.celldesigner.PositiveInfluence,
    negativelyInfluences: momapy.celldesigner.Inhibition,
}

model_element_class_to_layout_element_class = {
    momapy.celldesigner.PositiveInfluence: momapy.celldesigner.PositiveInfluenceLayout,
    momapy.celldesigner.Inhibition: momapy.celldesigner.InhibitionLayout,
}
