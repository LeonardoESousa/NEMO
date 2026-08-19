import sys
from unittest.mock import MagicMock

import numpy as np
import pytest

import nemo.__main__ as nemo_main


def _set_argv(monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", ["nemo", *args])


def test_cli_check_option_calls_susceptibility_and_exits(monkeypatch):
    _set_argv(monkeypatch, "-c", "sample.log")

    susceptibility_mock = MagicMock()
    monkeypatch.setattr(nemo_main.nemo.tools, "susceptibility_check", susceptibility_mock)

    with pytest.raises(SystemExit) as exc:
        nemo_main.main()

    assert exc.value.code == 0
    susceptibility_mock.assert_called_once_with("sample.log")


def test_cli_geom_option_prints_geometry_and_exits(monkeypatch, capsys):
    _set_argv(monkeypatch, "-g", "geom.log")

    geometry = np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
    atoms = ["C", "H"]
    monkeypatch.setattr(nemo_main.nemo.parser, "pega_geom", lambda _: (geometry, atoms))

    with pytest.raises(SystemExit) as exc:
        nemo_main.main()

    assert exc.value.code == 0

    out = capsys.readouterr().out
    assert "2" in out
    assert "C" in out
    assert "H" in out


def test_cli_geom_option_falls_back_to_interface_on_missing_file(monkeypatch):
    _set_argv(monkeypatch, "-g", "missing.log")

    def _raise_file_not_found(_):
        raise FileNotFoundError

    interface_mock = MagicMock()
    monkeypatch.setattr(nemo_main.nemo.parser, "pega_geom", _raise_file_not_found)
    monkeypatch.setattr(nemo_main, "interface", interface_mock)

    nemo_main.main()

    interface_mock.assert_called_once()


def test_cli_ensemble_file_option_runs_gather_data_for_each_state(monkeypatch):
    _set_argv(monkeypatch, "-e", "single.log", "S1,T1")

    gather_mock = MagicMock()
    monkeypatch.setattr(nemo_main, "gather_data", gather_mock)

    with pytest.raises(SystemExit) as exc:
        nemo_main.main()

    assert exc.value.code == 0
    assert gather_mock.call_count == 2
    gather_mock.assert_any_call("S1", save=True, filename="single.log")
    gather_mock.assert_any_call("T1", save=True, filename="single.log")
