import numpy as np

import nemo.empirical as empirical
import nemo.tools


def test_fetch_next_omega_expands_then_refines():
    assert empirical.fetch_next_omega([150], [3.0], 20) == 170
    assert empirical.fetch_next_omega([150, 170], [3.0, 2.0], 20) == 190
    assert empirical.fetch_next_omega(
        [150, 170, 190], [3.0, 1.0, 2.0], 20
    ) == 180


def test_rodar_omega_optimizes_s1_then_runs_empirical_sp(monkeypatch):
    generated = []
    checks = []

    def fake_generate(template, rem, atoms, geometry, filename, **kwargs):
        generated.append((template, filename, kwargs))
        return filename

    monkeypatch.setattr(empirical, "gera_file", fake_generate)
    monkeypatch.setattr(
        empirical, "_run_input", lambda filename, batch, nproc: filename[:-3] + "log"
    )
    monkeypatch.setattr(empirical, "_archive", lambda *files: None)
    monkeypatch.setattr(
        empirical.nemo.parser,
        "pega_geom",
        lambda _: (np.ones((1, 3)), ["H"]),
    )

    def fake_check(log, fit=None, tuning=0):
        checks.append(log)
        return 1.25, 2

    monkeypatch.setattr(empirical.nemo.tools, "susceptibility_check", fake_check)

    distance, root = empirical.rodar_omega(
        "fit.npy", ["H"], np.zeros((1, 3)), 4, 150, "batch.sh", "$rem\n$end"
    )

    assert (distance, root) == (1.25, 2)
    assert [item[0] for item in generated] == [
        "state_tracking_opt",
        "empirical",
    ]
    assert generated[0][2]["state"] == 1
    assert checks == ["td-150-sp-.log"]


def test_rodar_omega_can_skip_optimization(monkeypatch):
    generated = []

    def fake_generate(template, rem, atoms, geometry, filename, **kwargs):
        generated.append(template)
        return filename

    monkeypatch.setattr(empirical, "gera_file", fake_generate)
    monkeypatch.setattr(
        empirical, "_run_input", lambda filename, batch, nproc: filename[:-3] + "log"
    )
    monkeypatch.setattr(empirical, "_archive", lambda *files: None)
    monkeypatch.setattr(
        empirical.nemo.tools,
        "susceptibility_check",
        lambda *args, **kwargs: (1.25, 1),
    )

    empirical.rodar_omega(
        "fit.npy",
        ["H"],
        np.zeros((1, 3)),
        4,
        150,
        "batch.sh",
        "$rem\n$end",
        relax=False,
    )

    assert generated == ["empirical"]


def test_initial_geometry_uses_twocalc_qchem_input_parser(monkeypatch):
    expected = (np.zeros((1, 3)), ["H"])
    monkeypatch.setattr(empirical.nemo.parser, "pega_geom_qchem", lambda _: expected)

    assert empirical._starting_geometry("molecule.in", 150, []) == expected


def test_starting_state_defaults_to_s1(monkeypatch):
    assert empirical._starting_state("fit.npy", 150, []) == 1


def test_empirical_template_generates_vacuum_and_pcm_jobs(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    filename = empirical.gera_file(
        "empirical",
        "$rem\nmethod wb97x-d\nbasis def2-svp\n$end",
        ["H"],
        np.zeros((1, 3)),
        "test.com",
        omega="150",
        cm="0 1",
        num_ex=5,
        stat="3.0",
        optic="1.96",
    )

    contents = (tmp_path / filename).read_text(encoding="utf-8")
    assert contents.count("@@@") == 1
    assert "set_state_deriv" not in contents
    assert "StateSpecific           Perturb" in contents


def test_empirical_launcher_enables_optimization_by_default(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    answers = iter(["y", "8"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))

    def fake_fetch(description, extensions):
        return {
            "input": "molecule.in",
            "file with spec2epsilon fit data": "fit.npy",
            "batch script": "batch.sh",
        }[description]

    launched = []
    monkeypatch.setattr(nemo.tools, "fetch_file", fake_fetch)
    monkeypatch.setattr(
        nemo.tools.subprocess,
        "Popen",
        lambda command, **kwargs: launched.append((command, kwargs)),
    )

    nemo.tools.empirical_omega()

    assert launched[0][0] == [
        "empirical_tuning",
        "molecule.in",
        "8",
        "0.15",
        "0.02",
        "yes",
        "batch.sh",
        "fit.npy",
    ]
