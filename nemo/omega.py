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


def pega_energia(log_file):
    """
    Returns the last total energy reported on a line containing
    'Convergence criterion met'.
    """
    energia = None

    with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "Convergence criterion met" in line:
                vals = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", line)
                if vals:
                    energia = float(vals[0])

    if energia is None:
        raise ValueError(f"Could not find total energy in {log_file}")

    return energia


def pega_homo(log_file):
    """Return the HOMO energy, in hartree, from the last Q-Chem MO blocks."""
    import re

    with open(log_file, encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    def get_spin_homo(spin):
        # Matches both "Alpha MOs" and "Alpha MOs, Unrestricted"
        header = re.compile(rf"^\s*{spin}\s+MOs\b", re.IGNORECASE)

        starts = [
            i for i, line in enumerate(lines)
            if header.match(line)
        ]

        if not starts:
            return None

        occupied = False
        energies = []

        for line in lines[starts[-1] + 1:]:
            text = line.strip()

            if re.match(r"^(Alpha|Beta)\s+MOs\b", text, re.IGNORECASE):
                break

            if "Occupied" in text:
                occupied = True
                continue

            if "Virtual" in text:
                break

            if occupied:
                # Decimal point requirement ignores symmetry labels:
                # "1 A  1 B  2 A  2 B"
                values = re.findall(
                    r"[-+]?(?:\d+\.\d*|\.\d+)(?:[DEde][-+]?\d+)?",
                    line,
                )

                energies.extend(
                    float(value.replace("D", "E").replace("d", "e"))
                    for value in values
                )

        if not energies:
            raise ValueError(
                f"Could not find occupied {spin} orbitals in {log_file}"
            )

        return energies[-1]

    homo_alpha = get_spin_homo("Alpha")
    homo_beta = get_spin_homo("Beta")

    if homo_alpha is None:
        raise ValueError(f"Could not find Alpha MOs in {log_file}")

    # Restricted calculations may print only the alpha block.
    if homo_beta is None:
        return homo_alpha

    # Unrestricted calculation: highest occupied energy across both spins.
    return max(homo_alpha, homo_beta)

##RUNS CALCULATIONS############################################
def rodar_omega(atomos, geom, nproc, omega, batch_file, relax, rem, numjobs):
    omega = f"{omega:03.0f}"
    files = []
    if relax:
        file = gera_file(
            'omega_opt',
            rem,
            atomos,
            geom,
            f"OPT-{omega}-.com",
            omega=omega,
            cm="0 1",
        )
        files.append(file)
        the_watcher = nemo.tools.Watcher('.',key="OPT-")
        the_watcher.run(batch_file, nproc, 1)
        the_watcher.hold_watch()
        geom, atomos = nemo.parser.pega_geom(file[:-3] + "log")
        file2 = gera_file(
            'omega_sp',
            rem,
            atomos,
            geom,
            f"pos-{omega}-sp-.com",
            omega=omega,
            cm="+1 2",
        )
        files.append(file2)
        files3 = gera_file(
            'omega_sp',
            rem,
            atomos,
            geom,
            f"neg-{omega}-sp-.com",
            omega=omega,
            cm="-1 2",
        )
        files.append(files3)
        the_watcher = nemo.tools.Watcher('.',key="sp-")
        the_watcher.run(batch_file, nproc, min(numjobs,2))
        the_watcher.hold_watch()
    else:
        file = gera_file(
            'omega_sp',
            rem,
            atomos,
            geom,
            f"OPT-{omega}-.com",
            omega=omega,
            cm="0 1",
        )
        files.append(file)
        file2 = gera_file(
            'omega_sp',
            rem,
            atomos,
            geom,
            f"pos-{omega}-sp-.com",
            omega=omega,
            cm="+1 2",
        )
        files.append(file2)
        files3 = gera_file(
            'omega_sp',
            rem,
            atomos,
            geom,
            f"neg-{omega}-sp-.com",
            omega=omega,
            cm="-1 2",
        )
        files.append(files3)
        the_watcher = nemo.tools.Watcher('.',key=f"-{omega}-")
        the_watcher.run(batch_file, nproc, min(numjobs,3))
        the_watcher.hold_watch()


    try:
        os.mkdir("Logs")
    except FileExistsError:
        pass
    

    for file in files:
        if "pos-" in file:
            cation = pega_energia(file[:-3] + "log")
        elif "neg-" in file:
            anion = pega_energia(file[:-3] + "log")
            homo_anion = pega_homo(file[:-3] + "log")
        elif "OPT-" in file:
            neutro = pega_energia(file[:-3] + "log")
            homo_neutro = pega_homo(file[:-3] + "log")

    for file in files:
        shutil.move(file, "Logs/" + file)
        shutil.move(file[:-3] + "log", "Logs/" + file[:-3] + "log")
    J = np.sqrt(
        ((homo_neutro + cation - neutro) ** 2 + (homo_anion + neutro - anion) ** 2)
    ) * (27.2114)
    return J
    
    

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


def fetch_grad(x, f, i):
    """
    Estimate the gradient at x[i] using finite differences and Lagrange interpolation.

    Parameters:
    x (list or array): x values.
    f (list or array): Corresponding f(x) values.
    i (int): Index of the point to estimate the gradient.

    Returns:
    float: Estimated gradient at x[i].
    """
    if len(x) == 1:
        # Only one point; gradient is undefined, return 0
        return 0.0

    if len(x) == 2:
        # Two points; compute simple finite difference
        return (f[1] - f[0]) / (x[1] - x[0])

    if i == 0:
        # Forward finite difference for the first point
        return (f[1] - f[0]) / (x[1] - x[0])

    if i == len(x) - 1:
        # Backward finite difference for the last point
        return (f[-1] - f[-2]) / (x[-1] - x[-2])

    # General case: Use Lagrange interpolation for three points
    x0, x1, x2 = x[i - 1], x[i], x[i + 1]
    f0, f1, f2 = f[i - 1], f[i], f[i + 1]

    # Compute denominators for the Lagrange basis polynomials
    denom0 = (x0 - x1) * (x0 - x2)
    denom1 = (x1 - x0) * (x1 - x2)
    denom2 = (x2 - x0) * (x2 - x1)

    # First derivative (gradient) using corrected Lagrange formula
    l0p = (x1 - x2) / denom0
    l1p = (2 * x1 - x0 - x2) / denom1
    l2p = (x1 - x0) / denom2
    grad = f0 * l0p + f1 * l1p + f2 * l2p

    return grad



def main():
    geomlog = sys.argv[1]
    nproc = sys.argv[2]
    omega1 = sys.argv[3]
    passo = sys.argv[4]
    relax = sys.argv[5].lower()
    script = sys.argv[6]
    parallel = sys.argv[7].lower()
    rem, _, extra = nemo.parser.busca_input(geomlog)
    rem += extra + "\n"

    if relax != 'yes':
        relax = False
    else:
        relax = True  

    if parallel.lower() == "y":
        numjobs = 10
    else:
        numjobs = 1

    try:
        int(nproc)
        passo = float(passo) * 1000
        omega1 = float(omega1) * 1000
    except ValueError:
        nemo.parser.fatal_error("nproc, omega and step must be numbers. Goodbye!")
    omegas, Js = [], []
    try:
        data = np.loadtxt("omega.lx", dtype=float)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        omegas = data[:, 0].tolist()
        Js = data[:, 1].tolist()
    except FileNotFoundError:
        pass    

    iteration = 0
    while iteration < 100:
        if omega1 in omegas:
            ind = omegas.index(omega1)
            J = Js[ind]
        else:
            try:
                #find existing omega closest to omega1
                d_omega = [abs(omega1 - om) for om in omegas]
                geom_log = omegas[d_omega.index(min(d_omega))]
                geom_log = f"Logs/OPT-{geom_log:03.0f}-.log"
                G, atomos = nemo.parser.pega_geom(geom_log)
            except (ValueError, FileNotFoundError):
                G, atomos = nemo.parser.pega_geom(geomlog)   
            J = rodar_omega(
                atomos, G, nproc, omega1, script, relax, rem, numjobs
                )
            omegas.append(omega1)
            Js.append(J)
        omegas, Js = map(list, zip(*sorted(zip(omegas, Js))))
        #index of min J
        ind = Js.index(min(Js))
        omega1 = omegas[ind]
        grad = fetch_grad(omegas, Js, ind)

        if grad == 0:
            delta_omega = passo
        else:
            delta_omega = -1*Js[ind] / grad

        if ind == len(Js) - 1:
            max_omega = 500
        else:
            max_omega = omegas[ind+1]
        if ind == 0:
            min_omega = 0
        else:
            min_omega = omegas[ind-1]

        omega1 += delta_omega
        omega1 = np.round(omega1, 0)
        if min_omega > omega1 or omega1 > max_omega or int(omega1) in omegas:
            left = (min_omega - omegas[ind])/2
            right = (max_omega - omegas[ind])/2
            omega1 =   omegas[ind] + max(left,right,key=abs)
        omega1 = int(omega1)
        if max_omega - min_omega <= 2:
            break
    
        write_tolog(omegas, Js, f"#Best value so far:")
        iteration += 1
        

    write_tolog(omegas, Js, "#Done! Optimized value:")
    menor = omegas[Js.index(min(Js, key=abs))] 
    log = f"Logs/OPT-{menor:.0f}-.log"
    G, atomos = nemo.parser.pega_geom(log)
    rem, _, _ = nemo.parser.busca_input(log)
    gera_file("td-dft", rem, atomos, G, "s0_sp.com", cm="0 1", stat=3.0, optic=1.96, num_ex=5, soc="false")
    
    
    gera_file(
        "opt_td",
        rem,
        atomos,
        G,
        "freqs1.com",
        omega=f"{menor:.0f}",
        cm="0 1",
    )
    the_watcher = nemo.tools.Watcher('.',key="s0_sp.com")
    the_watcher.run(script, nproc, 1)
    the_watcher.hold_watch()

    with open("state_analysis.txt", "w") as f:
        with redirect_stdout(f):
            try:
                nemo.tools.susceptibility_check("s0_sp.log")
            except:
                print("Could not perform susceptibility check. Please check s0_sp.log for details.")

if __name__ == "__main__":
    sys.exit(main())
