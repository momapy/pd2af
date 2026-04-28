import io
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
    def test_mode_choices_lists_three_modes(self):
        assert set(pd2af.cli._MODE_CHOICES) == {
            "normal",
            "no-complex",
            "pure-af",
        }

    def test_layout_choices_includes_documented_modes(self):
        assert set(pd2af.cli._LAYOUT_CHOICES) == {"plain", "overlay", "auto"}


class TestCliMainOutputFile:
    def test_writes_xml_output_when_xml_extension(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.xml"
        pd2af.cli.main([example_map_path, "-o", str(out_path)])
        assert out_path.exists()
        assert out_path.stat().st_size > 0

    def test_writes_pickle_output_when_pickle_extension(
        self, tmp_path, example_map_path
    ):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main([example_map_path, "-o", str(out_path)])
        assert out_path.exists()
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        assert len(roundtrip.model.species) > 0

    def test_no_complex_mode_flag(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [example_map_path, "-m", "no-complex", "-o", str(out_path)]
        )
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        names = sorted(s.name for s in roundtrip.model.species)
        # Under no-complex, complex D drops out in favour of subunit C.
        assert "D" not in names

    @pytest.mark.skipif(
        not has_dot_binary(), reason="graphviz `dot` binary not on PATH"
    )
    def test_pure_af_mode_with_auto_layout(self, tmp_path, example_map_path):
        out_path = tmp_path / "out.pickle"
        pd2af.cli.main(
            [
                example_map_path,
                "-m",
                "pure-af",
                "-l",
                "auto",
                "-o",
                str(out_path),
            ]
        )
        roundtrip = momapy.io.core.read(str(out_path), reader="pickle").obj
        names = sorted(s.name for s in roundtrip.model.species)
        # Under pure-af (no-complex inheritance), D drops out.
        assert "D" not in names

    def test_pure_af_with_plain_layout_exits_nonzero(self, example_map_path):
        with pytest.raises(ValueError):
            pd2af.cli.main([example_map_path, "-m", "pure-af", "-l", "plain"])

    def test_invalid_mode_choice_exits(self, example_map_path, capsys):
        with pytest.raises(SystemExit):
            pd2af.cli.main([example_map_path, "-m", "bogus"])

    def test_invalid_layout_choice_exits(self, example_map_path, capsys):
        with pytest.raises(SystemExit):
            pd2af.cli.main([example_map_path, "-l", "bogus"])


class TestCliMainStdout:
    def test_no_output_arg_writes_pickle_to_stdout_buffer(self, example_map_path):
        captured = io.BytesIO()
        # Patch sys.stdout to expose a `.buffer` attribute, mirroring the
        # interface CLI code uses for byte writes to stdout.
        fake_stdout = unittest.mock.MagicMock()
        fake_stdout.buffer = captured
        with unittest.mock.patch("sys.stdout", new=fake_stdout):
            pd2af.cli.main([example_map_path])
        assert captured.tell() > 0
