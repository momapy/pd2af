"""Build the Rules reference page (`rules.md`) from the live rule registry.

Run as a `mkdocs-gen-files` script at every `mkdocs build`, so the page can
never drift from `src/pd2af/rules.py`. The page has three parts that share a
single canonical definition of each rule:

1. a by-mode spine listing each mode's applicable groups and rules as tables
   (the primary, most-important view);
2. a canonical "Rule groups" reference where each rule's full ``docs`` prose is
   written exactly once, in a per-group table, with a stable anchor that the
   other parts link into;
3. an alphabetical index of every rule (identifier, group, variant, modes and a
   one-line summary).

Every one-line summary is the first sentence of a ``docs`` string; those first
sentences are written as pure natural language in ``pd2af.rules``.
"""

import mkdocs_gen_files

from pd2af.cli import _INPUT_LANGUAGE_DISPLAY_NAMES, _MODE_CHOICES
from pd2af.rules import build_registry

PAGE_PATH = "rules.md"

# Order in which base and per-language variant rules are presented within a
# group: language-agnostic base rules first, then each language variant in the
# order the CLI lists its input languages.
_VARIANT_ORDER = tuple(_INPUT_LANGUAGE_DISPLAY_NAMES)


def make_rule_anchor(rule_identifier):
    """Stable in-page anchor for a rule's canonical row."""
    return "rule-" + rule_identifier.replace(":", "-")


def make_first_sentence(docs):
    """The first sentence of a ``docs`` string, whitespace-collapsed."""
    collapsed = " ".join(docs.split())
    sentence_end = collapsed.find(". ")
    if sentence_end == -1:
        return collapsed
    return collapsed[: sentence_end + 1]


def flatten_docs(docs):
    """A ``docs`` string collapsed to a single line for a table cell."""
    return escape_table_cell(" ".join(docs.split()))


def escape_table_cell(text):
    """Escape a string so it is safe inside a Markdown table cell."""
    return text.replace("|", "\\|")


def variant_label(variant_key):
    """Human-facing label for a variant key, or empty for a base rule."""
    if variant_key is None:
        return ""
    return _INPUT_LANGUAGE_DISPLAY_NAMES[variant_key]


def collect_group_rules(group):
    """Return ``[(rule, variant_key)]`` for a group in presentation order.

    ``variant_key`` is ``None`` for language-agnostic base rules, otherwise the
    input-language token (e.g. ``"celldesigner"``) the variant rule belongs to.
    """
    rules_with_variant = [(rule, None) for rule in group.rules]
    for variant_key in _VARIANT_ORDER:
        for rule in group.variants.get(variant_key, ()):
            rules_with_variant.append((rule, variant_key))
    return rules_with_variant


def mode_to_profile(mode):
    """Map a hyphenated transform-mode name to its underscored ASP profile."""
    return mode.replace("-", "_")


def modes_for_group(group):
    """The modes a group is used in, as ``"all"`` or an enumerated cell."""
    modes = [
        mode for mode in _MODE_CHOICES if mode_to_profile(mode) in group.profiles
    ]
    if len(modes) == len(_MODE_CHOICES):
        return "all"
    return ", ".join(f"`{mode}`" for mode in modes)


def write_intro(page):
    page.write("# Rules reference\n\n")
    page.write(
        "The pd2af transformation is driven by a set of rules (encoded in "
        "[ASP](https://en.wikipedia.org/wiki/Answer_set_programming)). Every "
        "rule is language-agnostic (a *base* rule) or specific to an input "
        "language (a *variant* rule, for CellDesigner or SBGN PD). Rules are "
        "organized in groups, each representing a coherent functional unit.\n\n"
    )
    page.write(
        "This page is generated from the live registry in `pd2af.rules`. It "
        "presents a [by-mode](#by-mode) view of which rules each mode uses, a "
        "[canonical reference](#rule-groups) documenting each rule once, and an "
        "alphabetical [index](#rule-index) of every rule.\n\n"
    )


def write_by_mode(page, registry):
    page.write("## By mode {#by-mode}\n\n")
    for mode in _MODE_CHOICES:
        profile = mode_to_profile(mode)
        group_identifiers = registry.profiles.get(profile, frozenset())
        page.write(f"### `{mode}`\n\n")
        for group in registry.groups.values():
            if group.identifier not in group_identifiers:
                continue
            page.write(f"**`{group.identifier}`**\n\n")
            page.write(f"{make_first_sentence(group.docs)}\n\n")
            page.write("| Rule | Variant |\n")
            page.write("| --- | --- |\n")
            for rule, variant_key in collect_group_rules(group):
                anchor = make_rule_anchor(rule.identifier)
                identifier_cell = f"[`{rule.identifier}`](#{anchor})"
                page.write(f"| {identifier_cell} | {variant_label(variant_key)} |\n")
            page.write("\n")


def write_canonical_reference(page, registry):
    page.write("## Rule groups {#rule-groups}\n\n")
    for group in registry.groups.values():
        page.write(f"### `{group.identifier}`\n\n")
        if group.docs:
            page.write(f"{group.docs}\n\n")
        page.write("| Rule | Variant | Documentation |\n")
        page.write("| --- | --- | --- |\n")
        for rule, variant_key in collect_group_rules(group):
            anchor = make_rule_anchor(rule.identifier)
            identifier_cell = f'<a id="{anchor}"></a>`{rule.identifier}`'
            page.write(
                f"| {identifier_cell} | {variant_label(variant_key)} "
                f"| {flatten_docs(rule.docs)} |\n"
            )
        page.write("\n")


def write_rule_index(page, registry):
    page.write("## Rule index {#rule-index}\n\n")
    page.write(
        "All rules in alphabetical order (identifiers are group- and "
        "variant-prefixed).\n\n"
    )
    page.write("| Rule | Group | Variant | Modes | Summary |\n")
    page.write("| --- | --- | --- | --- | --- |\n")
    indexed_rules = []
    for group in registry.groups.values():
        for rule, variant_key in collect_group_rules(group):
            indexed_rules.append((rule, variant_key, group))
    for rule, variant_key, group in sorted(
        indexed_rules, key=lambda entry: entry[0].identifier
    ):
        anchor = make_rule_anchor(rule.identifier)
        identifier_cell = f"[`{rule.identifier}`](#{anchor})"
        summary_cell = escape_table_cell(make_first_sentence(rule.docs))
        page.write(
            f"| {identifier_cell} | `{group.identifier}` "
            f"| {variant_label(variant_key)} | {modes_for_group(group)} "
            f"| {summary_cell} |\n"
        )
    page.write("\n")


def main():
    registry = build_registry()
    with mkdocs_gen_files.open(PAGE_PATH, "w") as page:
        write_intro(page)
        write_by_mode(page, registry)
        write_canonical_reference(page, registry)
        write_rule_index(page, registry)
    mkdocs_gen_files.set_edit_path(PAGE_PATH, "generate_rules_reference.py")


main()
