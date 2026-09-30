import numpy as np

import nemo.tools


def test_tuning_returns_mahalanobis_distance_and_selected_root(monkeypatch, tmp_path):
    alpha_opt = 1.0 / 3.0
    parser_result = (
        np.array([2.0, 3.0]),
        np.array([1.5]),
        None,
        None,
        None,
        np.array([0.2, 0.4]) * alpha_opt,
        np.array([0.1]) * alpha_opt,
        0.05,
        np.array([0.05, 0.15]),
        np.array([0.05]),
    )
    monkeypatch.setattr(
        nemo.tools.nemo.parser, "pega_energias", lambda _: parser_result
    )
    monkeypatch.setattr(nemo.tools, "fetch_nr", lambda _: (3.0, np.sqrt(2.0)))

    fit = tmp_path / "fit.npy"
    np.save(
        fit,
        {
            "E_vac": 3.0 - 0.5 * 0.2 * 0.3243,
            "chi": 0.5,
            "covariance_matrix": np.diag([0.01**2, 0.01**2]),
        },
    )

    distance, root = nemo.tools.susceptibility_check(
        "calculation.log", fit=fit, tuning=1
    )

    assert root == 2
    assert np.isclose(distance, 0.0)


def test_plain_check_preserves_report_only_mode(monkeypatch, capsys):
    parser_result = (
        np.array([2.0]),
        np.array([1.5]),
        None,
        None,
        None,
        np.array([0.1]),
        np.array([0.1]),
        0.1,
        np.array([0.1]),
        np.array([0.1]),
    )
    monkeypatch.setattr(
        nemo.tools.nemo.parser, "pega_energias", lambda _: parser_result
    )
    monkeypatch.setattr(nemo.tools, "fetch_nr", lambda _: (3.0, np.sqrt(2.0)))

    assert nemo.tools.susceptibility_check("calculation.log") is None
    assert "E_vac(eV)" in capsys.readouterr().out
