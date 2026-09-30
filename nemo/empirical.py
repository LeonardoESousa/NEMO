#!/usr/bin/env python3
"""Empirical omega tuning against a spec2epsilon fit."""

import os
import shutil
import sys
from contextlib import redirect_stdout

import numpy as np

import nemo.parser
import nemo.tools


DEFAULT_NUM_ROOTS = 5
CONFIDENCE_LIMIT = np.sqrt(2.30)


def gera_file(template, rem, atoms, geometry, filename, **kwargs):
    """Create a Q-Chem input using one of NEMO's packaged templates."""
    template_text = nemo.tools.load_template(template)
    header = template_text.format(
        basic=nemo.tools.extract_basic_rem(rem.lower().strip()),
        **kwargs,
    )
    header, bottom = header.split("#GGG#")
    nemo.tools.write_input(atoms, geometry, header, bottom, filename)
    return filename


def _run_input(filename, batch_file, nproc):
    watcher = nemo.tools.Watcher(".", key=filename)
    watcher.run(batch_file, nproc, 1)
    watcher.hold_watch()
    log_file = filename[:-3] + "log"
    if not os.path.isfile(log_file):
        nemo.parser.fatal_error(f"Calculation failed: {log_file} was not found.")
    return log_file


def _archive(*files):
    os.makedirs("Logs", exist_ok=True)
    for filename in files:
        shutil.move(filename, os.path.join("Logs", filename))


def rodar_omega(
    fit,
    atoms,
    geometry,
    nproc,
    omega,
    batch_file,
    rem,
    state=1,
    relax=True,
):
    """Optionally optimize a state, then evaluate one empirical SP input."""
    omega_label = f"{omega:03.0f}"
    files = []

    if relax:
        optimization_input = gera_file(
            "state_tracking_opt",
            rem,
            atoms,
            geometry,
            f"OPT-{omega_label}-.com",
            omega=omega_label,
            cm="0 1",
            state=state,
            num_ex=min(DEFAULT_NUM_ROOTS, state + 2),
        )
        optimization_log = _run_input(optimization_input, batch_file, nproc)
        geometry, atoms = nemo.parser.pega_geom(optimization_log)
        files.extend([optimization_input, optimization_log])

    evaluation_input = gera_file(
        "empirical",
        rem,
        atoms,
        geometry,
        f"td-{omega_label}-sp-.com",
        omega=omega_label,
        cm="0 1",
        num_ex=max(DEFAULT_NUM_ROOTS, state + 2),
        stat="3.0",
        optic="1.96",
    )
    evaluation_log = _run_input(evaluation_input, batch_file, nproc)
    distance, final_state = nemo.tools.susceptibility_check(
        evaluation_log, fit=fit, tuning=1
    )

    files.extend([evaluation_input, evaluation_log])
    _archive(*files)
    return distance, final_state


def write_tolog(omegas, distances, message):
    """Write a restartable, human-readable tuning table."""
    ordered = sorted(zip(omegas, distances))
    with open("omega.lx", "w", encoding="utf-8") as output:
        output.write(f"{'# w(10^3 bohr^-1)':<22}{'M-distance':<12}\n")
        for omega, distance in ordered:
            output.write(f"{omega:<22.0f}{distance:<12.4f}\n")
        best_omega, _ = min(ordered, key=lambda pair: pair[1])
        output.write(f"\n{message} {best_omega:3.0f}\n")


def fetch_next_omega(
    omegas, distances, initial_step, omega_min=0, omega_max=500
):
    """Continue outward around the best point, then bisect its bracket."""
    if len(omegas) != len(distances):
        raise ValueError("omegas and distances must have the same length.")
    if not omegas:
        raise ValueError("At least one omega value is required.")
    if initial_step <= 0:
        raise ValueError("initial_step must be positive.")

    sampled = {}
    for omega, distance in zip(omegas, distances):
        omega = int(round(omega))
        distance = float(distance)
        if not np.isfinite(distance):
            raise ValueError(f"Invalid Mahalanobis distance at omega {omega}.")
        sampled[omega] = min(distance, sampled.get(omega, np.inf))

    x = np.array(sorted(sampled), dtype=int)
    values = np.array([sampled[value] for value in x], dtype=float)
    evaluated = set(x.tolist())
    best_index = int(np.argmin(values))

    def valid(candidate):
        return omega_min <= candidate <= omega_max and candidate not in evaluated

    def midpoint(left, right):
        if right - left <= 1:
            return None
        center = int(round(0.5 * (left + right)))
        for offset in range(right - left):
            for candidate in (center - offset, center + offset):
                if left < candidate < right and valid(candidate):
                    return candidate
        return None

    if len(x) == 1:
        step = max(1, int(round(initial_step)))
        for candidate in (min(omega_max, x[0] + step), max(omega_min, x[0] - step)):
            if valid(int(candidate)):
                return int(candidate)
        return None

    if best_index == 0:
        spacing = max(1, int(x[1] - x[0]))
        candidate = max(omega_min, int(x[0] - spacing))
        return candidate if valid(candidate) else midpoint(int(x[0]), int(x[1]))

    if best_index == len(x) - 1:
        spacing = max(1, int(x[-1] - x[-2]))
        candidate = min(omega_max, int(x[-1] + spacing))
        return candidate if valid(candidate) else midpoint(int(x[-2]), int(x[-1]))

    left, best, right = map(
        int, (x[best_index - 1], x[best_index], x[best_index + 1])
    )
    candidates = []
    for neighbor, candidate in (
        (best_index - 1, midpoint(left, best)),
        (best_index + 1, midpoint(best, right)),
    ):
        if candidate is not None:
            candidates.append(
                (abs(int(x[neighbor]) - best), values[neighbor], candidate)
            )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (-item[0], item[1]))
    return candidates[0][2]


def _load_restart():
    try:
        data = np.loadtxt("omega.lx", dtype=float)
    except FileNotFoundError:
        return [], []
    except ValueError as error:
        raise ValueError(
            "Could not read omega.lx; its numerical tuning table is malformed."
        ) from error
    if data.size == 0:
        return [], []
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError("omega.lx must contain omega and distance columns.")
    return data[:, 0].tolist(), data[:, 1].tolist()


def _starting_geometry(geomlog, omega, sampled_omegas, relax=True):
    """Use the closest completed optimization, as in the standard tuner."""
    if relax and sampled_omegas:
        closest = min(sampled_omegas, key=lambda value: abs(value - omega))
        closest_log = f"Logs/OPT-{closest:03.0f}-.log"
        if os.path.isfile(closest_log):
            return nemo.parser.pega_geom(closest_log)
    # twocalc's generic dispatcher identifies Q-Chem from its output banner;
    # a raw .in file therefore needs the Q-Chem parser explicitly.
    return nemo.parser.pega_geom_qchem(geomlog)


def _starting_state(fit, omega, sampled_omegas):
    """Use S1 initially, then continue from the closest evaluated omega."""
    if sampled_omegas:
        closest = min(sampled_omegas, key=lambda value: abs(value - omega))
        closest_log = f"Logs/td-{closest:03.0f}-sp-.log"
        if os.path.isfile(closest_log):
            _, state = nemo.tools.susceptibility_check(
                closest_log, fit=fit, tuning=1
            )
            return state
    return 1


def main():
    if len(sys.argv) != 8:
        nemo.parser.fatal_error(
            "Usage: empirical_tuning INPUT NPROC OMEGA STEP RELAX "
            "BATCH_SCRIPT FIT.npy"
        )
    geomlog, nproc, omega, step, relax, script, fit = sys.argv[1:]
    relax = relax.lower() == "yes"
    try:
        nproc = int(nproc)
        omega = float(omega) * 1000
        step = float(step) * 1000
    except ValueError:
        nemo.parser.fatal_error("nproc, omega, and step must be numbers. Goodbye!")
    if nproc < 1 or step <= 0:
        nemo.parser.fatal_error("nproc and step must be positive.")
    if not os.path.isfile(fit):
        nemo.parser.fatal_error(f"Fit file not found: {fit}")

    rem, _, extra = nemo.parser.busca_input(geomlog)
    rem += extra + "\n"
    try:
        omegas, distances = _load_restart()
    except ValueError as error:
        nemo.parser.fatal_error(str(error))

    for _ in range(100):
        matches = [i for i, value in enumerate(omegas) if np.isclose(value, omega)]
        if not matches:
            geometry, atoms = _starting_geometry(
                geomlog, omega, omegas, relax=relax
            )
            state = _starting_state(fit, omega, omegas)
            distance, _ = rodar_omega(
                fit,
                atoms,
                geometry,
                nproc,
                omega,
                script,
                rem,
                state=state,
                relax=relax,
            )
            omegas.append(float(omega))
            distances.append(float(distance))

        ordered = sorted(zip(omegas, distances))
        omegas = [pair[0] for pair in ordered]
        distances = [pair[1] for pair in ordered]
        write_tolog(omegas, distances, "# Best value so far:")
        if min(distances) <= CONFIDENCE_LIMIT:
            break
        next_omega = fetch_next_omega(omegas, distances, step)
        if next_omega is None:
            break
        omega = next_omega

    if not distances:
        nemo.parser.fatal_error("No successful omega calculations were obtained.")
    write_tolog(omegas, distances, "# Done! Optimized value:")
    best_index = int(np.argmin(distances))
    best_omega = omegas[best_index]
    best_log = f"Logs/td-{best_omega:03.0f}-sp-.log"
    with open("omega_final.lx", "w", encoding="utf-8") as output:
        with redirect_stdout(output):
            nemo.tools.susceptibility_check(best_log, fit=fit, tuning=2)


if __name__ == "__main__":
    sys.exit(main())
