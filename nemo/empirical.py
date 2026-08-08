#!/usr/bin/env python3
import os
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
def rodar_omega(fit, atomos, geom, nproc, omega, batch_file, state, rem, numjobs):
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
    

    J, new_state = nemo.tools.susceptibility_check(f"td-{omega}-sp-.log",fit=fit, tuning=1)

    for file in files:
        shutil.move(file, "Logs/" + file)
        shutil.move(file[:-3] + "log", "Logs/" + file[:-3] + "log")

    return J, new_state
    
    

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
    omega_min=0,
    omega_max=500,
):
    """
    Select the next omega value using only previously calculated
    Mahalanobis distances.

    Strategy
    --------
    1. With one sampled point, first try increasing omega.
    2. If the current best point is at an edge of the sampled range,
       continue sampling outward.
    3. Once the best point is bracketed, bisect the larger adjacent
       interval.
    4. Stop when no unevaluated integer omega remains around the
       current best point.

    Parameters
    ----------
    omegas : sequence of float
        Previously evaluated omega values in units of
        10^-3 bohr^-1.

    Js : sequence of float
        Mahalanobis distance associated with each omega.

    initial_step : float
        Initial step in units of 10^-3 bohr^-1.

    omega_min, omega_max : int
        Allowed omega range in units of 10^-3 bohr^-1.

    Returns
    -------
    int or None
        Next unevaluated omega value. Returns None when the search
        interval cannot be refined further.
    """
    if len(omegas) != len(Js):
        raise ValueError(
            "omegas and Js must have the same length."
        )

    if len(omegas) == 0:
        raise ValueError(
            "At least one omega value is required."
        )

    if initial_step <= 0:
        raise ValueError(
            "initial_step must be positive."
        )

    omega_min = int(round(omega_min))
    omega_max = int(round(omega_max))

    if omega_min >= omega_max:
        raise ValueError(
            "omega_min must be smaller than omega_max."
        )

    # Remove duplicate omega values. If duplicates exist, keep the
    # smallest distance.
    sampled = {}

    for omega, J in zip(omegas, Js):
        omega = int(round(omega))
        J = float(J)

        if not np.isfinite(J):
            raise ValueError(
                f"Invalid Mahalanobis distance at omega {omega}."
            )

        if omega not in sampled or J < sampled[omega]:
            sampled[omega] = J

    x = np.array(
        sorted(sampled),
        dtype=int,
    )

    f = np.array(
        [sampled[omega] for omega in x],
        dtype=float,
    )

    evaluated = set(x.tolist())

    best_index = int(np.argmin(f))

    def valid(candidate):
        return (
            candidate is not None
            and omega_min <= candidate <= omega_max
            and candidate not in evaluated
        )

    def midpoint(left, right):
        """
        Return an unevaluated integer midpoint, if one exists.
        """
        if right - left <= 1:
            return None

        candidate = int(
            round(0.5 * (left + right))
        )

        if candidate <= left:
            candidate = left + 1

        if candidate >= right:
            candidate = right - 1

        if valid(candidate):
            return candidate

        # Search outward from the midpoint if rounding or previous
        # evaluations produced an occupied point.
        for offset in range(1, right - left):
            lower = candidate - offset
            upper = candidate + offset

            if left < lower < right and valid(lower):
                return lower

            if left < upper < right and valid(upper):
                return upper

        return None

    # First point: try increasing omega. If the upper boundary is
    # unavailable, try decreasing omega.
    if len(x) == 1:
        step = max(
            1,
            int(round(abs(initial_step))),
        )

        candidate = min(
            omega_max,
            int(x[0] + step),
        )

        if valid(candidate):
            return candidate

        candidate = max(
            omega_min,
            int(x[0] - step),
        )

        if valid(candidate):
            return candidate

        return None

    # The best point is the lowest sampled omega.
    if best_index == 0:
        spacing = max(
            1,
            int(x[1] - x[0]),
        )

        # Continue toward lower omega.
        candidate = max(
            omega_min,
            int(x[0] - spacing),
        )

        if valid(candidate):
            return candidate

        # If the lower boundary has been reached, refine the
        # interval between the best point and its right neighbor.
        return midpoint(
            int(x[0]),
            int(x[1]),
        )

    # The best point is the highest sampled omega.
    if best_index == len(x) - 1:
        spacing = max(
            1,
            int(x[-1] - x[-2]),
        )

        # Continue toward higher omega.
        candidate = min(
            omega_max,
            int(x[-1] + spacing),
        )

        if valid(candidate):
            return candidate

        # If the upper boundary has been reached, refine the
        # interval between the left neighbor and the best point.
        return midpoint(
            int(x[-2]),
            int(x[-1]),
        )

    # The current best point is bracketed.
    left = int(x[best_index - 1])
    best = int(x[best_index])
    right = int(x[best_index + 1])

    left_width = best - left
    right_width = right - best

    candidates = []

    left_candidate = midpoint(
        left,
        best,
    )

    if left_candidate is not None:
        candidates.append(
            {
                "width": left_width,
                "neighbor_J": f[best_index - 1],
                "omega": left_candidate,
            }
        )

    right_candidate = midpoint(
        best,
        right,
    )

    if right_candidate is not None:
        candidates.append(
            {
                "width": right_width,
                "neighbor_J": f[best_index + 1],
                "omega": right_candidate,
            }
        )

    if not candidates:
        return None

    # Prefer the larger unexplored interval. If both intervals have
    # equal width, prefer the side with the lower neighboring J.
    candidates.sort(
        key=lambda item: (
            -item["width"],
            item["neighbor_J"],
        )
    )

    return candidates[0]["omega"]


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

        if data.shape[1] < 2:
            raise ValueError

        omegas = data[:, 0].tolist()
        Js = data[:, 1].tolist()

    except FileNotFoundError:
        pass

    except ValueError:
        nemo.parser.fatal_error(
            "Could not read omega.lx. The file must contain "
            "only the numerical omega-tuning table."
        )

    # Every omega calculation uses the original geometry.
    G, atomos = nemo.parser.pega_geom(geomlog)

    # Joint 68% confidence limit for two fitted parameters.
    confidence_limit = np.sqrt(2.30)

    iteration = 0

    while iteration < 100:
        # Reuse a previously calculated omega when available.
        matching_indices = [
            i
            for i, omega in enumerate(omegas)
            if np.isclose(omega, omega1)
        ]

        if matching_indices:
            index = matching_indices[0]
            omega1 = omegas[index]
            J = Js[index]

        else:
            J, new_state = rodar_omega(
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

            omegas.append(float(omega1))
            Js.append(float(J))

        # Keep the tuning data ordered by omega.
        ordered_data = sorted(
            zip(omegas, Js),
            key=lambda pair: pair[0],
        )

        omegas = [
            pair[0]
            for pair in ordered_data
        ]
        Js = [
            pair[1]
            for pair in ordered_data
        ]

        best_index = int(np.argmin(Js))
        best_J = Js[best_index]
        best_omega = omegas[best_index]

        write_tolog(
            omegas,
            Js,
            "# Best value so far:",
        )

        # Stop when the best result lies inside the joint 68%
        # confidence region.
        if best_J <= confidence_limit:
            break

        next_omega = fetch_next_omega(
            omegas,
            Js,
            initial_step=passo,
            omega_min=0,
            omega_max=500,
        )

        # No untested integer omega remains in the current search region.
        if next_omega is None:
            break

        omega1 = next_omega
        iteration += 1

    if not Js:
        nemo.parser.fatal_error(
            "No successful omega calculations were obtained."
        )

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

    # Reevaluate the selected state directly from the best-omega log.
    J, new_state = nemo.tools.susceptibility_check(
        best_log,
        fit=fit,
        tuning=1,
    )

    # If a higher excited state is selected, optimize that state once
    # using state tracking at the final omega.
    if new_state > 1:
        G_final, atoms_final = nemo.parser.pega_geom(
            best_log
        )

        basic_rem = nemo.tools.extract_basic_rem(
            rem
        )

        opt_file = gera_file(
            "state_tracking_opt",
            basic_rem,
            atoms_final,
            G_final,
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
                "State-tracking optimization failed: "
                f"{opt_log} was not found."
            )

        G_final, atoms_final = nemo.parser.pega_geom(
            opt_log
        )

        sp_file = gera_file(
            "empirical",
            basic_rem,
            atoms_final,
            G_final,
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
                "Final single-point calculation failed: "
                f"{final_file} was not found."
            )

    else:
        final_file = best_log

    # Keep the human-readable final diagnostic separate from the
    # numerical restart table in omega.lx.
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
