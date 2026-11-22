import numpy as np

import nemo.tools


def test_tuning_returns_mahalanobis_distance_and_selected_root(monkeypatch, tmp_path):
    alpha_opt = 1.0 / 3.0
    parser_result = (
        np.array([2.0, 3.0]),
        np.array([1.5]),
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
        np.array([0.1]),
        np.array([0.1]),
        0.1,
        np.array([0.1]),
        np.array([0.1]),
    )
    def parse(_, return_matching=False):
        matching = {'two_calculations': True}
        for spin in ('s', 't'):
            matching[spin] = dict(pcm_roots=np.array([1]), scores=np.array([1.]),
                                  matrix=np.ones((1, 1)), retained=1,
                                  raw_shift=np.array([0.1]),assignment_scores=np.array([1.]),)
        return parser_result + (matching,) if return_matching else parser_result
    monkeypatch.setattr(nemo.tools.nemo.parser, "pega_energias", parse)
    monkeypatch.setattr(nemo.tools, "fetch_nr", lambda _: (3.0, np.sqrt(2.0)))

    assert nemo.tools.susceptibility_check("calculation.log") is None
    assert "E_vac(eV)" in capsys.readouterr().out


def test_fit_discards_dominated_root_and_reports_correction(tmp_path, monkeypatch, capsys):
    alpha_opt = 1.0 / 3.0
    parsed = (
        np.array([2.568, 3.407]), np.array([]), None,
        np.array([0.005, 0.405]) * alpha_opt, np.array([]),
        0.1, np.array([0.1, 0.1]), np.array([]),
        {"s": {"ptss": np.array([0.004, 0.400]) * alpha_opt, "assignment_scores": [1, 1]},
         "t": {"ptss": np.array([]), "assignment_scores": []}},
    )
    parse = lambda *args, **kwargs: parsed
    monkeypatch.setattr(nemo.tools.nemo.parser, "pega_energias", parse)
    monkeypatch.setattr(nemo.tools, "fetch_nr", lambda _: (3.0, np.sqrt(2.0)))
    fit = tmp_path / "fit.npy"
    np.save(fit, {"E_vac": 3.0, "chi": 0.6,
                  "covariance_matrix": [[0.04, 0.038], [0.038, 0.04]]})

    nemo.tools.susceptibility_check("calculation.log", fit=fit)
    output = capsys.readouterr().out
    assert "Numerical best match: S1 [both]" in output
    assert "Selected state: S2" in output
    assert "Selected correction: --correction both" in output

    nemo.tools.susceptibility_check("calculation.log", fit=fit, correction="ss")
    output = capsys.readouterr().out
    assert "Selected state: S2" in output
    assert "Selected correction: --correction ss" in output


def test_fit_gives_vacuum_energy_more_uncertainty_than_chi(tmp_path, monkeypatch, capsys):
    parsed = (
        np.array([2.00, 2.28]), np.array([]), None,
        np.array([0.03, 0.15]) / 3.0, np.array([]),
        0.1, np.array([0.1, 0.1]), np.array([]),
        {"s": {"assignment_scores": [1, 1]}, "t": {"assignment_scores": []}},
    )
    monkeypatch.setattr(nemo.parser, "pega_energias", lambda *args, **kwargs: parsed)
    monkeypatch.setattr(nemo.tools, "fetch_nr", lambda _: (3.0, np.sqrt(2.0)))
    fit = tmp_path / "fit.npy"
    np.save(fit, {"E_vac": 2.0, "chi": 0.15, "covariance_matrix": np.eye(2) * 0.01**2})

    nemo.tools.susceptibility_check("calculation.log", fit=fit, correction="both")
    output = capsys.readouterr().out
    assert "χ = 0.043 eV, E_vac = 0.150 eV" in output
    assert "Selected state: S2" in output
