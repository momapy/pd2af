class Rule(object):
    def __init__(self, id_: str, code: str, doc: str | None = None):
        self.id_ = id_
        self.code = code
        self.doc = doc


class Control(object):
    def __init__(
        self, rules: list[Rule] | None = None, set_active: list[str] | None = None
    ):
        self.rules = rules if rules is not None else []
        self.set_active = set_active if set_active is not None else []

    def add_rule(self, rule: Rule):
        self.rules.append(rule)

    def add_active(self, id_: str):
        self.set_active.append(id_)


if __name__ == "__main__":
    rule_1 = Rule(id_="rule_1", code="a:-b.", doc="if b then a")
    rule_2 = Rule(id_="rule_2", code="c:-d.", doc="if d then c")

    DEFAULT_RULES: list[Rule] = [rule_1, rule_2]
    control = Control(DEFAULT_RULES)
