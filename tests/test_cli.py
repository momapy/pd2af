import io
import json
import os
import unittest.mock

import pytest

import momapy.io.core

import pd2af.cli

from tests._helpers import has_dot_binary


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
        assert set(pd2af.cli._MODE_CHOICES) == {
            "normal",
            "normal-no-complex",
            "keep-species",
            "keep-species-no-complex",
            "casq",
        }

    def test_layout_choices_includes_documented_modes(self):
        assert set(pd2af.cli._LAYOUT_CHOICES) == {"plain", "overlay", "auto"}


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
    def test_merged_modes_with_plain_layout_raise(
        self, example_map_path, mode
    ):
        with pytest.raises(ValueError):
            pd2af.cli.main(["transform", example_map_path, "-m", mode, "-l", "plain"])

    @pytest.mark.parametrize(
        "mode", ["keep-species", "keep-species-no-complex", "casq"]
    )
    def test_per_species_modes_accept_plain_layout(
        self, tmp_path, example_map_path, mode
    ):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            ["transform", example_map_path, "-m", mode, "-l", "plain", "-o", str(out_path)]
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


class TestListModes:
    def test_text_output_lists_modes_and_layouts(self, capsys):
        pd2af.cli.main(["list-modes"])
        out = capsys.readouterr().out
        for mode in pd2af.cli._MODE_CHOICES:
            assert mode in out
        for layout_mode in pd2af.cli._LAYOUT_CHOICES:
            assert layout_mode in out

    def test_json_output_has_expected_structure(self, capsys):
        pd2af.cli.main(["list-modes", "--json"])
        data = json.loads(capsys.readouterr().out)
        assert set(data) == {"transformation_modes", "layout_modes"}
        # One entry per table row, one key per column.
        names = {mode["transformation_mode"] for mode in data["transformation_modes"]}
        assert names == set(pd2af.cli._MODE_CHOICES)
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
        # Merged-proteoform modes only accept `auto`; the per-species/casq
        # modes accept all three layout modes.
        assert modes["normal"]["layout_modes"] == ["auto"]
        assert modes["normal-no-complex"]["layout_modes"] == ["auto"]
        assert set(modes["keep-species"]["layout_modes"]) == {
            "plain",
            "overlay",
            "auto",
        }

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
