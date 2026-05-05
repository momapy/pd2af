"""Helpers shared by walkers."""


def compartments_outermost_first(compartments):
    def depth(compartment):
        result = 0
        seen = set()
        current = compartment.outside
        while current is not None and id(current) not in seen:
            seen.add(id(current))
            result += 1
            current = current.outside
        return result
    return sorted(compartments, key=depth)


def collect_ancestor_compartments(compartments):
    expanded = set(compartments)
    frontier = expanded
    while True:
        next_frontier = set(
            compartment.outside
            for compartment in frontier
            if compartment.outside is not None and compartment.outside not in expanded
        )
        if not next_frontier:
            break
        expanded |= next_frontier
        frontier = next_frontier
    return expanded


def collect_templates_from_species(species_iterable):
    collected = set()

    def visit(species):
        template = getattr(species, "template", None)
        if template is not None:
            collected.add(template)
        for subunit in getattr(species, "subunits", ()) or ():
            visit(subunit)

    for species in species_iterable:
        visit(species)
    return collected
