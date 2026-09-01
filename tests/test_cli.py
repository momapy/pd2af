import argparse
import io
import json
import os
import unittest.mock

import pytest

import momapy.io.core

import pd2af
import pd2af.cli
import pd2af.modes

from tests._helpers import MAPS_DIR, has_dot_binary


def _parse_transform_args(argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    pd2af.cli._add_transform_parser(subparsers)
    return parser.parse_args(argv)


class TestWriterForOutput:
    @pytest.mark.parametrize(
        "filename,expected",
        [
            ("foo.xml", "celldesigner"),
            ("foo.sbml", "celldesigner"),
            ("foo.pickle", "pickle"),
            ("foo.pkl", "pickle"),
            ("foo.bin", "pickle"),
            ("noext", "pickle"),
            ("FOO.XML", "celldesigner"),
            ("FOO.PICKLE", "pickle"),
        ],
    )
    def test_extension_mapping(self, filename, expected):
        assert pd2af.cli._writer_for_output(filename) == expected


class TestModeAndLayoutChoices:
    def test_mode_choices_lists_supported_modes(self):
        assert {mode.name for mode in pd2af.modes._BUILTIN_TRANSFORMATION_MODES} == {
            "normal",
            "normal-no-complex",
            "keep-species",
            "keep-species-no-complex",
            "keep-reactions",
        }

    def test_layout_choices_includes_documented_modes(self):
        assert set(pd2af.cli._LAYOUT_CHOICES) == {"plain", "overlay", "dot", "auto"}


class TestCliMainOutputFile:
    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_writes_xml_output_when_xml_extension(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.xml"
        pd2af.cli.main(["transform", example_map_path, "-o", str(out_path)])
        assert out_path.exists()
        assert out_path.stat().st_size > 0

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_writes_pickle_output_when_pickle_extension(
        self, tmp_path, example_map_path
    ):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(["transform", example_map_path, "-o", str(out_path)])
        assert out_path.exists()
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        assert len(roundtrip.model.species) > 0

    def test_keep_species_no_complex_mode_flag(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [
                "transform",
                example_map_path,
                "-m",
                "keep-species-no-complex",
                "-l",
                "plain",
                "-o",
                str(out_path),
            ]
        )
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        names = sorted(s.name for s in roundtrip.model.species)
        # Under keep-species-no-complex, complex D drops out in favour of
        # subunit C.
        assert "D" not in names

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_normal_no_complex_mode_with_auto_layout(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [
                "transform",
                example_map_path,
                "-m",
                "normal-no-complex",
                "-l",
                "auto",
                "-o",
                str(out_path),
            ]
        )
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        names = sorted(s.name for s in roundtrip.model.species)
        # Under normal-no-complex (merge proteoforms + drop complexes), D drops out.
        assert "D" not in names

    @pytest.mark.parametrize("mode", ["normal", "normal-no-complex"])
    def test_merged_modes_with_plain_layout_raise(self, example_map_path, mode):
        with pytest.raises(ValueError):
            pd2af.cli.main(["transform", example_map_path, "-m", mode, "-l", "plain"])

    @pytest.mark.parametrize(
        "mode",
        ["keep-species", "keep-species-no-complex", "keep-reactions"],
    )
    def test_per_species_modes_accept_plain_layout(
        self, tmp_path, example_map_path, mode
    ):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [
                "transform",
                example_map_path,
                "-m",
                mode,
                "-l",
                "plain",
                "-o",
                str(out_path),
            ]
        )
        assert out_path.exists()

    def test_invalid_mode_choice_exits(self, example_map_path, capsys):
        with pytest.raises(SystemExit):
            pd2af.cli.main(["transform", example_map_path, "-m", "bogus"])

    def test_invalid_layout_choice_exits(self, example_map_path, capsys):
        with pytest.raises(SystemExit):
            pd2af.cli.main(["transform", example_map_path, "-l", "bogus"])


class TestCliMainStdout:
    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_no_output_arg_writes_pickle_to_stdout_buffer(self, example_map_path):
        captured = io.BytesIO()
        # Patch sys.stdout to expose a `.buffer` attribute, mirroring the
        # interface CLI code uses for byte writes to stdout.
        fake_stdout = unittest.mock.MagicMock()
        fake_stdout.buffer = captured
        with unittest.mock.patch("sys.stdout", new=fake_stdout):
            pd2af.cli.main(["transform", example_map_path])
        assert captured.tell() > 0


class TestActiveFlag:
    def test_active_flag_defaults_to_none(self):
        args = _parse_transform_args(["transform", "map.xml"])
        assert args.set_active is None

    def test_active_flag_is_repeatable_and_accumulates(self):
        args = _parse_transform_args(["transform", "map.xml", "-a", "s1", "-a", "s3"])
        assert args.set_active == ["s1", "s3"]

    def test_active_flag_marks_species_active(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [
                "transform",
                example_map_path,
                "-m",
                "keep-species",
                "-l",
                "plain",
                "-a",
                "s1",
                "-o",
                str(out_path),
            ]
        )
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        names = sorted(s.name for s in roundtrip.model.species)
        # Species A (id `s1`) is not an activity without the flag; -a surfaces it.
        assert "A" in names


class TestInactiveFlag:
    def test_inactive_flag_defaults_to_none(self):
        args = _parse_transform_args(["transform", "map.xml"])
        assert args.set_inactive is None

    def test_inactive_flag_is_repeatable_and_accumulates(self):
        args = _parse_transform_args(["transform", "map.xml", "-i", "s1", "-i", "s3"])
        assert args.set_inactive == ["s1", "s3"]

    def test_inactive_flag_suppresses_species(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [
                "transform",
                example_map_path,
                "-m",
                "keep-species",
                "-l",
                "plain",
                "-i",
                "s2",
                "-o",
                str(out_path),
            ]
        )
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        names = sorted(s.name for s in roundtrip.model.species)
        # Species B (id `s2`) is active by default; -i suppresses it.
        assert "B" not in names


class TestGlobalActivityFlags:
    def test_global_flags_default_to_false(self):
        args = _parse_transform_args(["transform", "map.xml"])
        assert args.set_all_active is False
        assert args.set_all_inactive is False

    def test_set_all_active_short_flag(self):
        args = _parse_transform_args(["transform", "map.xml", "-A"])
        assert args.set_all_active is True
        assert args.set_all_inactive is False

    def test_set_all_inactive_short_flag(self):
        args = _parse_transform_args(["transform", "map.xml", "-I"])
        assert args.set_all_inactive is True
        assert args.set_all_active is False

    def test_both_global_flags_is_a_parse_error(self):
        with pytest.raises(SystemExit):
            _parse_transform_args(["transform", "map.xml", "-A", "-I"])

    def test_set_all_active_marks_every_species(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [
                "transform",
                example_map_path,
                "-m",
                "keep-species",
                "-l",
                "plain",
                "-A",
                "-o",
                str(out_path),
            ]
        )
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        names = sorted(s.name for s in roundtrip.model.species)
        # Species A (id `s1`) is not an activity by default; -A surfaces it.
        assert "A" in names

    def test_set_all_inactive_drops_every_species(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [
                "transform",
                example_map_path,
                "-m",
                "keep-species",
                "-l",
                "plain",
                "-I",
                "-o",
                str(out_path),
            ]
        )
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        # Species B (id `s2`) is active by default; -I suppresses everything.
        assert not roundtrip.model.species


class TestListModes:
    def test_text_output_lists_modes_and_layouts(self, capsys):
        pd2af.cli.main(["list-modes"])
        out = capsys.readouterr().out
        for mode_name in pd2af.modes.get_transformation_modes():
            assert mode_name in out
        for layout_mode in pd2af.cli._LAYOUT_CHOICES:
            assert layout_mode in out

    def test_json_output_has_expected_structure(self, capsys):
        pd2af.cli.main(["list-modes", "--json"])
        data = json.loads(capsys.readouterr().out)
        assert set(data) == {"transformation_modes", "layout_modes"}
        # One entry per table row, one key per column.
        names = {mode["transformation_mode"] for mode in data["transformation_modes"]}
        assert names == set(pd2af.modes.get_transformation_modes())
        for mode in data["transformation_modes"]:
            assert set(mode) == {
                "transformation_mode",
                "layout_modes",
                "languages",
                "description",
            }
        layout_names = {row["layout_mode"] for row in data["layout_modes"]}
        assert layout_names == set(pd2af.cli._LAYOUT_CHOICES)

    def test_json_layout_modes_reflect_validation(self, capsys):
        pd2af.cli.main(["list-modes", "--json"])
        modes = {
            mode["transformation_mode"]: mode
            for mode in json.loads(capsys.readouterr().out)["transformation_modes"]
        }
        # Merged-proteoform modes only accept `dot`; the per-species modes
        # accept all three concrete layout modes.
        assert modes["normal"]["layout_modes"] == ["dot"]
        assert modes["normal-no-complex"]["layout_modes"] == ["dot"]
        assert set(modes["keep-species"]["layout_modes"]) == {
            "plain",
            "overlay",
            "dot",
        }

    def test_json_layout_modes_reflect_language_support(self, capsys):
        pd2af.cli.main(["list-modes", "--json"])
        rows = {
            row["layout_mode"]: row
            for row in json.loads(capsys.readouterr().out)["layout_modes"]
        }
        every_language = [
            properties["display_name"] for properties in pd2af.modes.LANGUAGES.values()
        ]
        # `overlay` dimming is CellDesigner-only; `auto` is a meta value, not a
        # concrete layout mode, so no language rejects it.
        assert rows["overlay"]["languages"] == ["CellDesigner"]
        assert rows["plain"]["languages"] == every_language
        assert rows["dot"]["languages"] == every_language
        assert rows["auto"]["languages"] == every_language

    def test_json_matches_rendered_tables(self, capsys):
        # The tables are an exact view of the JSON payload: every cell value
        # appears in the rendered text.
        pd2af.cli.main(["list-modes", "--json"])
        data = json.loads(capsys.readouterr().out)
        pd2af.cli.main(["list-modes"])
        text = capsys.readouterr().out
        for mode in data["transformation_modes"]:
            assert mode["transformation_mode"] in text
            for language in mode["languages"]:
                assert language in text


class TestAnnotationCarryCli:
    """The CLI threads input annotations end-to-end: a MIRIAM resource carried
    onto a top-level species appears in the written CellDesigner XML."""

    def test_carries_annotations_to_xml_output(self, tmp_path):
        jnk_path = os.path.join(MAPS_DIR, "JNK_pathway.xml")
        reader_result = momapy.io.core.read(jnk_path)
        # Compute a resource the carry places on a top-level output species.
        result = pd2af.transform(
            reader_result.obj,
            mode="keep-species",
            layout_mode="plain",
            element_to_annotations=reader_result.element_to_annotations,
            element_to_notes=reader_result.element_to_notes,
        )
        species = set(result.obj.model.species)
        carried_resources = [
            resource
            for output_element, annotations in result.element_to_annotations.items()
            if output_element in species
            for annotation in annotations
            for resource in annotation.resources
        ]
        assert carried_resources
        expected_resource = carried_resources[0]

        out_path = tmp_path / "out.xml"
        pd2af.cli.main(
            [
                "transform",
                jnk_path,
                "-m",
                "keep-species",
                "-l",
                "plain",
                "-o",
                str(out_path),
            ]
        )
        assert expected_resource in out_path.read_text()
