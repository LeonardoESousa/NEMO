#!/usr/bin/env python3
import os
import sys
import shutil
import numpy as np
import nemo.tools

###############################################################

def gera_file(atomos, G, basis, functional, omega, mem):
    header = f'''$rem
mem_total             {mem}
mem_static            500
jobtype               opt
GEOM_OPT_MAX_CYCLES   100
exchange              {functional}
omega                 {omega}
basis                 {basis}
cis_n_roots             3
cis_singlets            true
cis_triplets            false
CIS_STATE_DERIV         1
CIS_MAX_CYCLES          200
MAX_SCF_CYCLES          200
RPA                     false
$end

$molecule
0 1
'''

    bottom = f"$end\n"
    
    #generate qchem input for single point tda calculation
    nemo.tools.write_input(
            atomos,
            G,
            header,
            bottom,
            f'OPT-{omega}-.com',
        )
    return f'OPT-{omega}-.com'

def gera_sp_file(atomos, G, basis, functional, omega, mem):
    header = f'''$rem
mem_total             {mem}
mem_static            500
jobtype               sp
exchange              {functional}
omega                 {omega}
basis                 {basis}
cis_n_roots             3
cis_singlets            true
cis_triplets            false
CIS_STATE_DERIV         1
RPA                     false
CIS_RELAXED_DENSITY     TRUE
CIS_MAX_CYCLES          200
MAX_SCF_CYCLES          200
solvent_method          PCM
$end

$pcm
theory                  IEFPCM
ChargeSeparation        Marcus
StateSpecific           Perturb
$end

$solvent
Dielectric              3.00
OpticalDielectric       2.22
$end

$molecule
0 1
'''
    bottom = f"$end\n"    
    #generate qchem input for single point tda calculation
    nemo.tools.write_input(
            atomos,
            G,
            header,
            bottom,
            f'SP-{omega}-.com',
        )
    return f'SP-{omega}-.com'


##RUNS CALCULATIONS############################################
def rodar_omega(atomos, geom, nproc, basis, functional, omega, batch_file, e_exp, de_exp, chi_exp, dchi_exp, relax, mem):
    omega = f"{omega:03.0f}"
    files = []
    if relax:
        file = gera_file(atomos, geom, basis, functional, omega, mem)
        files.append(file)
        the_watcher = nemo.tools.Watcher('.',key="OPT-")
        the_watcher.run(batch_file, nproc, 1)
        the_watcher.hold_watch()
        geom, atomos = nemo.parser.pega_geom(file[:-3] + "log")
    
    file2 = gera_sp_file(atomos, geom, basis, functional, omega, mem)
    files.append(file2)   
    the_watcher = nemo.tools.Watcher('.',key="SP-")
    the_watcher.run(batch_file, nproc, 1)
    the_watcher.hold_watch()
    e_vac, chi = nemo.tools.susceptibility_check(file2[:-3] + "log", tuning=True) 
    try:
        os.mkdir("Logs")
    except FileExistsError:
        pass
    for file in files:
        shutil.move(file, "Logs/" + file)
        shutil.move(file[:-3] + "log", "Logs/" + file[:-3] + "log")
    delta_e = e_vac - e_exp
    delta_chi = chi  - chi_exp
    J = abs(delta_e/de_exp) + abs(delta_chi/dchi_exp)
    return J


###############################################################


##WRITES LOG WITH RESULTS######################################
def write_tolog(omegas, Js, frase):
    with open("omega.lx", "w", encoding='utf-8') as f:
        # Align headers appropriately
        f.write(f"{'#w(10^3 bohr^-1)':<15}{'J':<10}\n")
        
        # Sort the values by omega
        list1, list2 = zip(*sorted(zip(omegas, Js)))
        
        for i in range(len(list1)):
            # Format the columns with proper alignment
            f.write(f"{list1[i]:<15.0f}{list2[i]:<10.4f}\n")
        
        # Find the minimum J and its corresponding omega
        min_index = list2.index(min(list2, key=abs))
        f.write(f"\n{frase} {list1[min_index]:3.0f}\n")


###############################################################


def fetch_grad(x, f, i):
    """
    Estimate the first and second derivatives at x[i] using finite differences
    and Lagrange interpolation.

    Parameters:
    x (list or array): x values.
    f (list or array): Corresponding f(x) values.
    i (int): Index of the point to estimate derivatives.

    Returns:
    tuple: (first_derivative, second_derivative) at x[i].
    """
    if len(x) == 1:
        # Only one point; derivatives are undefined
        return 0.0, 0.0

    if len(x) == 2:
        # Two points; compute simple finite differences for gradient
        grad = (f[1] - f[0]) / (x[1] - x[0])
        hessian = 0.0  # Insufficient data for second derivative
        return grad, hessian

    if i == 0:
        # Forward finite difference for the first point
        grad = (f[1] - f[0]) / (x[1] - x[0])
        hessian = 0.0  # Insufficient data for second derivative
        return grad, hessian

    if i == len(x) - 1:
        # Backward finite difference for the last point
        grad = (f[-1] - f[-2]) / (x[-1] - x[-2])
        hessian = 0.0  # Insufficient data for second derivative
        return grad, hessian

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

    # Second derivative (Hessian) using corrected Lagrange formula
    l0pp = 2 / denom0
    l1pp = 2 / denom1
    l2pp = 2 / denom2
    hessian = f0 * l0pp + f1 * l1pp + f2 * l2pp

    return grad, hessian


def main():
    geomlog = sys.argv[1]
    functional = sys.argv[2]
    basis = sys.argv[3]
    nproc = sys.argv[4]
    omega1 = sys.argv[5]
    passo = sys.argv[6]
    relax = sys.argv[7].lower()
    script = sys.argv[8]
    e_exp = sys.argv[9].split()
    chi_exp = sys.argv[10].split()
    mem = sys.argv[11]
    if len(e_exp) == 2 and len(chi_exp) == 2:
        e_exp, de_exp = float(e_exp[0]), float(e_exp[1])
        chi_exp, dchi_exp = float(chi_exp[0]), float(chi_exp[1])
    else:
        e_exp = float(e_exp[0])
        chi_exp = float(chi_exp[0])
        de_exp = 1
        dchi_exp = 1
    if relax != 'yes':
        relax = False
    else:
        relax = True    
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
                geom_log = f"Logs/OPT-{geom_log:3.0f}-.log"
                G, atomos = nemo.parser.pega_geom(geom_log)
            except (ValueError, FileNotFoundError):
                G, atomos = nemo.parser.pega_geom(geomlog)   
            J = rodar_omega(
                atomos, G, nproc, basis, functional, omega1, script, e_exp, de_exp, chi_exp, dchi_exp, relax, mem)
            omegas.append(omega1)
            Js.append(J)
        omegas, Js = map(list, zip(*sorted(zip(omegas, Js))))
        #index of min J
        ind = Js.index(min(Js))
        omega1 = omegas[ind]
        grad, hessian = fetch_grad(omegas, Js, ind)

        if grad == 0:
            delta_omega = passo
        elif hessian == 0:
            delta_omega = -1*Js[ind] / grad
        else:
            delta_omega = -grad / hessian

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
    
        write_tolog(omegas, Js, f"#E_vac: {e_exp:.3f} ± {de_exp:.3f} eV\n#\u03C7_exp: {chi_exp:.3f} ± {dchi_exp:.3f} eV\n#Best value so far:")
        iteration += 1
        

    write_tolog(omegas, Js, f"#E_vac: {e_exp:.3f} ± {de_exp:.3f} eV\n#\u03C7_exp: {chi_exp:.3f} ± {dchi_exp:.3f} eV\n#Done! Optimized value:")
    menor = omegas[Js.index(min(Js, key=abs))] 
    log = f"Logs/SP-{menor:03d}-.log"
    #copy log file
    shutil.copy(log, "tuned.log")

if __name__ == "__main__":
    sys.exit(main())
