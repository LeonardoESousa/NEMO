#!/usr/bin/env python3
import os
import re
import sys
import shutil
import numpy as np
import nemo.tools
from contextlib import redirect_stdout

###############################################################
def gera_file(template, rem, atomos, geometry, filename, **kwargs):
    rem = rem.lower().strip()
    
    # Extract relevant part of rem
    basic_rem = nemo.tools.extract_basic_rem(rem)

    template = nemo.tools.load_template(template)
    
    # Inject computed + user-provided kwargs into format
    format_dict = {
        "basic": basic_rem,
        **kwargs
    }

    header = template.format(**format_dict)
    
    header, bottom = header.split("#GGG#")
    
    nemo.tools.write_input(atomos, geometry, header, bottom, filename)

    return filename




##RUNS CALCULATIONS############################################
def rodar_omega(e_vac, chi, atomos, geom, nproc, omega, batch_file, state, rem, numjobs):
    omega = f"{omega:03.0f}"
    files = []
    file = gera_file(
        'empirical',
        rem,
        atomos,
        geom,
        f"td-{omega}-sp-.com",
        omega=omega,
        cm="0 1",
        state=state,
        num_ex="5",
        stat="3.0",
        optic="1.96"
    )
    files.append(file)
    the_watcher = nemo.tools.Watcher('.',key=f"-{omega}-")
    the_watcher.run(batch_file, nproc, min(numjobs,3))
    the_watcher.hold_watch()


    try:
        os.mkdir("Logs")
    except FileExistsError:
        pass
    

    J, new_state, sign = nemo.tools.susceptibility_check(f"td-{omega}-sp-.log",E_vac_fit=e_vac, chi_fit=chi, tuning=True)

    for file in files:
        shutil.move(file, "Logs/" + file)
        shutil.move(file[:-3] + "log", "Logs/" + file[:-3] + "log")

    return J, new_state, sign
    
    

###############################################################


##WRITES LOG WITH RESULTS######################################
def write_tolog(omegas, Js, frase):
    with open("omega.lx", "w", encoding='utf-8') as f:
        # Align headers appropriately
        f.write(f"{'#w(10^3 bohr^-1)':<22}{'J':<12}\n")
        
        # Sort the values by omega
        list1, list2 = zip(*sorted(zip(omegas, Js)))
        
        for i in range(len(list1)):
            # Format the columns with proper alignment
            f.write(f"{list1[i]:<22.0f}{list2[i]:<12.4f}\n")
        
        # Find the minimum J and its corresponding omega
        min_index = list2.index(min(list2, key=abs))
        f.write(f"\n{frase} {list1[min_index]:3.0f}\n")


###############################################################


def fetch_next_omega(
    omegas,
    Js,
    initial_step,
    initial_sign=None,
    omega_min=0,
    omega_max=500,
):
    """
    Select the next omega value without assuming that J(omega) is smooth.

    Strategy
    --------
    1. With one point, use the direction supplied by susceptibility_check.
    2. If the best point is at an edge of the sampled interval, continue
       outward using the spacing to the nearest point.
    3. If the best point is bracketed, bisect the larger neighboring
       interval.
    4. Return None when no new integer omega exists around the best point.

    Parameters
    ----------
    omegas : sequence of float
        Previously evaluated omega values in units of 10^-3 bohr^-1.
    Js : sequence of float
        Corresponding non-negative distances.
    initial_step : float
        Initial omega step in units of 10^-3 bohr^-1.
    initial_sign : int or float, optional
        Direction returned by susceptibility_check:
        +1 means increase omega and -1 means decrease omega.
        Required when only one omega has been evaluated.
    omega_min, omega_max : int
        Allowed omega range.

    Returns
    -------
    int or None
        Next unevaluated integer omega, or None when the local interval
        cannot be refined further.
    """
    if len(omegas) != len(Js):
        raise ValueError("omegas and Js must have the same length.")

    if len(omegas) == 0:
        raise ValueError("At least one omega value is required.")

    # Sort points and remove duplicate omega values.
    sampled = {}

    for omega, J in zip(omegas, Js):
        omega = int(round(omega))
        J = float(J)

        # Keep the lowest J if an omega appears more than once.
        if omega not in sampled or J < sampled[omega]:
            sampled[omega] = J

    x = np.array(sorted(sampled), dtype=int)
    f = np.array([sampled[omega] for omega in x], dtype=float)

    best_index = int(np.argmin(f))
    evaluated = set(x.tolist())

    def valid(candidate):
        return (
            candidate is not None
            and omega_min <= candidate <= omega_max
            and candidate not in evaluated
        )

    # First point: use the physical direction from susceptibility_check.
    if len(x) == 1:
        if initial_sign is None or initial_sign == 0:
            raise ValueError(
                "initial_sign must be +1 or -1 when only one omega "
                "has been evaluated."
            )

        direction = 1 if initial_sign > 0 else -1
        candidate = int(
            round(x[0] + direction * abs(initial_step))
        )
        candidate = max(omega_min, min(omega_max, candidate))

        if valid(candidate):
            return candidate

        return None

    # Best point is at the lowest sampled omega: continue toward lower omega.
    if best_index == 0:
        spacing = x[1] - x[0]
        candidate = int(x[0] - spacing)
        candidate = max(omega_min, candidate)

        if valid(candidate):
            return candidate

        # The lower boundary has been reached. Refine toward the next point.
        width = x[1] - x[0]

        if width > 1:
            candidate = int(round(0.5 * (x[0] + x[1])))

            if valid(candidate):
                return candidate

        return None

    # Best point is at the highest sampled omega: continue toward higher omega.
    if best_index == len(x) - 1:
        spacing = x[-1] - x[-2]
        candidate = int(x[-1] + spacing)
        candidate = min(omega_max, candidate)

        if valid(candidate):
            return candidate

        # The upper boundary has been reached. Refine toward the previous point.
        width = x[-1] - x[-2]

        if width > 1:
            candidate = int(round(0.5 * (x[-2] + x[-1])))

            if valid(candidate):
                return candidate

        return None

    # The best point is bracketed by evaluated points.
    left = x[best_index - 1]
    best = x[best_index]
    right = x[best_index + 1]

    left_width = best - left
    right_width = right - best

    candidates = []

    if left_width > 1:
        left_candidate = int(round(0.5 * (left + best)))
        candidates.append(
            (
                left_width,
                f[best_index - 1],
                left_candidate,
            )
        )

    if right_width > 1:
        right_candidate = int(round(0.5 * (best + right)))
        candidates.append(
            (
                right_width,
                f[best_index + 1],
                right_candidate,
            )
        )

    if not candidates:
        return None

    # Prefer the wider unexplored interval. If both intervals have the
    # same width, prefer the side whose neighboring J is lower.
    candidates.sort(
        key=lambda item: (
            -item[0],
            item[1],
        )
    )

    for _, _, candidate in candidates:
        if valid(candidate):
            return candidate

    return None



def main():
    geomlog = sys.argv[1]
    nproc = sys.argv[2]
    omega1 = sys.argv[3]
    passo = sys.argv[4]
    script = sys.argv[5]
    fit = sys.argv[6]

    rem, _, extra = nemo.parser.busca_input(geomlog)
    rem += extra + "\n"

    state = 1
    numjobs = 10
    sign = None

    try:
        nproc = int(nproc)
        passo = float(passo) * 1000
        omega1 = float(omega1) * 1000
    except ValueError:
        nemo.parser.fatal_error(
            "nproc, omega, and step must be numbers. Goodbye!"
        )

    if not os.path.isfile(fit):
        nemo.parser.fatal_error(
            f"Fit file not found: {fit}"
        )

    omegas = []
    Js = []

    try:
        data = np.loadtxt(
            "omega.lx",
            dtype=float,
        )

        if data.ndim == 1:
            data = data.reshape(1, -1)

        omegas = data[:, 0].tolist()
        Js = data[:, 1].tolist()

    except FileNotFoundError:
        pass
    except ValueError:
        nemo.parser.fatal_error(
            "Could not read omega.lx. The file must contain "
            "only the omega-tuning table."
        )

    # Recover the initial direction when restarting from one point.
    if len(omegas) == 1:
        cached_omega = int(round(omegas[0]))
        cached_log = (
            f"Logs/td-{cached_omega:03d}-sp-.log"
        )

        if not os.path.isfile(cached_log):
            nemo.parser.fatal_error(
                "Could not recover the omega-tuning direction: "
                f"{cached_log} was not found."
            )

        _, _, sign = nemo.tools.susceptibility_check(
            cached_log,
            fit=fit,
            tuning=1,
        )

    # All omega calculations use the same geometry.
    G, atomos = nemo.parser.pega_geom(geomlog)

    # Joint 68% confidence limit for two fitted parameters.
    confidence_limit = np.sqrt(2.30)

    iteration = 0

    while iteration < 100:
        if omega1 in omegas:
            index = omegas.index(omega1)
            J = Js[index]

        else:
            J, new_state, sign = rodar_omega(
                fit,
                atomos,
                G,
                nproc,
                omega1,
                script,
                state,
                rem,
                numjobs,
            )

            omegas.append(omega1)
            Js.append(J)

        omegas, Js = map(
            list,
            zip(*sorted(zip(omegas, Js))),
        )

        best_index = int(np.argmin(Js))
        omega1 = omegas[best_index]

        next_omega = fetch_next_omega(
            omegas,
            Js,
            initial_step=passo,
            initial_sign=sign,
            omega_min=0,
            omega_max=500,
        )

        reached_confidence = (
            min(Js) <= confidence_limit
        )

        if next_omega is None or reached_confidence:
            break

        omega1 = next_omega

        write_tolog(
            omegas,
            Js,
            "# Best value so far:",
        )

        iteration += 1

    write_tolog(
        omegas,
        Js,
        "# Done! Optimized value:",
    )

    best_index = int(np.argmin(Js))
    best_omega = omegas[best_index]

    best_log = (
        f"Logs/td-{best_omega:03.0f}-sp-.log"
    )

    if not os.path.isfile(best_log):
        nemo.parser.fatal_error(
            f"Best-omega log file not found: {best_log}"
        )

    # Determine the best state from the actual best-omega calculation.
    J, new_state, _ = nemo.tools.susceptibility_check(
        best_log,
        fit=fit,
        tuning=1,
    )

    # If a higher state is the final target, optimize it once using
    # state tracking at the selected omega.
    if new_state > 1:
        G, atomos = nemo.parser.pega_geom(
            best_log
        )

        basic_rem = nemo.tools.extract_basic_rem(
            rem
        )

        opt_file = gera_file(
            "state_tracking_opt",
            basic_rem,
            atomos,
            G,
            f"freqs{new_state}.com",
            omega=f"{best_omega:03.0f}",
            cm="0 1",
            state=new_state,
            num_ex=str(new_state + 2),
        )

        the_watcher = nemo.tools.Watcher(
            ".",
            key=opt_file,
        )
        the_watcher.run(
            script,
            nproc,
            1,
        )
        the_watcher.hold_watch()

        opt_log = opt_file[:-3] + "log"

        if not os.path.isfile(opt_log):
            nemo.parser.fatal_error(
                f"State-tracking optimization failed: "
                f"{opt_log} was not found."
            )

        G, atomos = nemo.parser.pega_geom(
            opt_log
        )

        sp_file = gera_file(
            "empirical",
            basic_rem,
            atomos,
            G,
            f"s{new_state}_sp.com",
            omega=f"{best_omega:03.0f}",
            cm="0 1",
            state=1,
            num_ex="5",
            soc="false",
            stat="3.0",
            optic="1.96",
        )

        the_watcher = nemo.tools.Watcher(
            ".",
            key=sp_file,
        )
        the_watcher.run(
            script,
            nproc,
            1,
        )
        the_watcher.hold_watch()

        final_file = sp_file[:-3] + "log"

        if not os.path.isfile(final_file):
            nemo.parser.fatal_error(
                f"Final single-point calculation failed: "
                f"{final_file} was not found."
            )

    else:
        final_file = best_log

    # Keep the numerical tuning table separate from the final report.
    with open(
        "omega_final.lx",
        "w",
        encoding="utf-8",
    ) as output:
        with redirect_stdout(output):
            try:
                nemo.tools.susceptibility_check(
                    final_file,
                    fit=fit,
                    tuning=2,
                )
            except Exception as error:
                print(
                    "Could not perform the final "
                    "susceptibility check."
                )
                print(f"Reason: {error}")

if __name__ == "__main__":
    sys.exit(main())
