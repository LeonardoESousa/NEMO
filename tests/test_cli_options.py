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
    susceptibility_mock.assert_called_once_with("sample.log", fit=None, correction=None)


@pytest.mark.parametrize("options, correction", [((), None), (("--correction", "ss"), "ss")])
def test_cli_check_option_accepts_fit_file(monkeypatch, options, correction):
    _set_argv(monkeypatch, "-c", "sample.log", *options, "fit.npy")

    susceptibility_mock = MagicMock()
    monkeypatch.setattr(nemo_main.nemo.tools, "susceptibility_check", susceptibility_mock)

    with pytest.raises(SystemExit) as exc:
        nemo_main.main()

    assert exc.value.code == 0
    susceptibility_mock.assert_called_once_with("sample.log", fit="fit.npy", correction=correction)


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


@pytest.mark.parametrize("options, correction", [((), "both"), (("--correction", "ss"), "ss")])
def test_cli_ensemble_file_option_runs_gather_data_for_each_state(monkeypatch, options, correction):
    _set_argv(monkeypatch, "-e", "single.log", "S1,T1", *options)

    gather_mock = MagicMock()
    monkeypatch.setattr(nemo_main, "gather_data", gather_mock)

    with pytest.raises(SystemExit) as exc:
        nemo_main.main()

    assert exc.value.code == 0
    assert gather_mock.call_count == 2
    gather_mock.assert_any_call("S1", save=True, filename="single.log", correction=correction)
    gather_mock.assert_any_call("T1", save=True, filename="single.log", correction=correction)


@pytest.mark.parametrize("operation, callback", [("5", "absorption"), ("6", "rates"), ("7", "gather_data")])
def test_menu_applies_one_correction_to_all_states(monkeypatch, operation, callback):
    answers = iter([operation, "ss", *([] if operation == "7" else ["y"]), "S1, T1"])
    monkeypatch.setattr("builtins.input", lambda *args: next(answers))
    monkeypatch.setattr(nemo_main.nemo.tools, "check_for_updates", MagicMock())
    monkeypatch.setattr(nemo_main.nemo.tools, "get_nr", lambda: (3.0, 1.4))
    monkeypatch.setattr(nemo_main, "export_results", MagicMock())
    run = MagicMock(return_value=(None, None))
    monkeypatch.setattr(nemo_main, callback, run)

    nemo_main.interface()

    assert [call.args[0] for call in run.call_args_list] == ["S1", "T1"]
    assert [call.kwargs["correction"] for call in run.call_args_list] == ["ss", "ss"]
