#!/usr/bin/env python3
import os
import warnings
import numpy as np
import pandas as pd
import nemo.tools


LIGHT_SPEED = nemo.tools.LIGHT_SPEED
HBAR_EV = nemo.tools.HBAR_EV
HBAR_J = nemo.tools.HBAR_J
BOLTZ_EV = nemo.tools.BOLTZ_EV
E_CHARGE = nemo.tools.E_CHARGE
MASS_E = nemo.tools.MASS_E
EPSILON_0 = nemo.tools.EPSILON_0


##RETURNS LIST OF LOG FILES WITH NORMAL TERMINATION######################################
def check_normal(files):
    normal = []
    add = False
    for file in files:
        with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
            for line in log_file:
                if (
                    "TDDFT/TDA Excitation Energies" in line
                    or "TDDFT Excitation Energies" in line
                ):
                    exc = True
                elif "Excited state" in line and exc:
                    eng = float(line.split()[7])
                    if eng < 0:
                        add = False
                    else:
                        add = True
                    exc = False
                elif "Have a nice day" in line and add:
                    normal.append(file)
    return normal


#########################################################################################


##GETS ENERGIES, OSCS, AND INDICES FOR Sn AND Tn STATES##################################
def pega_energias(file):
    ss_mark = "Excited-state properties with   relaxed density"
    with open(file, "r", encoding="utf-8") as log_file:
        exc = False
        corr = False
        correction, correction2 = [], []
        for line in log_file:
            if (
                "TDDFT/TDA Excitation Energies" in line
                or "TDDFT Excitation Energies" in line
            ):
                energies, spins, oscs, ind = [], [], [], []
                exc = True
            elif ss_mark in line:
                corr = True
            elif "Solute Internal Energy" in line:
                sol_int = float(line.split()[5])
            elif "Total Free Energy" in line:
                total_free = float(line.split()[9])
            elif "Excited state" in line and exc:
                energies.append(float(line.split()[7]))
                ind.append(int(line.split()[2].replace(":", "")))
            elif "Multiplicity" in line and exc:
                spins.append(line.split()[1])
            elif "Strength" in line and exc:
                oscs.append(float(line.split()[2]))
            elif (
                "---------------------------------------------------" in line
                and exc
                and len(energies) > 0
            ):
                exc = False
            elif "SS-PCM correction" in line and corr:
                correction.append(-1 * float(line.split()[-2]))
            elif "LR-PCM correction" in line and corr:
                correction2.append(-2 * float(line.split()[-2]))
            elif (
                "------------------------ END OF SUMMARY -----------------------"
                in line
                and corr
            ):
                corr = False
            elif "Total energy in the final basis set" in line:
                line = line.split()
                total_nopcm = float(line[8])
        if len(correction) == 0:  # When run on logs that do not employ pcm
            correction = np.zeros(len(energies))
            sol_int = total_nopcm
            total_free = total_nopcm
        singlets = np.array(
            [energies[i] for i in range(len(energies)) if spins[i] == "Singlet"]
        )
        ss_s = np.array(
            [
                correction[i] + correction2[i]
                for i in range(len(correction))
                if spins[i] == "Singlet"
            ]
        )
        ind_s = np.array([ind[i] for i in range(len(ind)) if spins[i] == "Singlet"])
        oscs = np.array(
            [oscs[i] for i in range(len(energies)) if spins[i] == "Singlet"]
        )
        triplets = np.array(
            [energies[i] for i in range(len(energies)) if spins[i] == "Triplet"]
        )
        ss_t = np.array(
            [
                correction[i] + correction2[i]
                for i in range(len(correction))
                if spins[i] == "Triplet"
            ]
        )
        ind_t = np.array([ind[i] for i in range(len(ind)) if spins[i] == "Triplet"])

        oscs = np.array([x for _, x in zip(singlets, oscs)])
        ind_s = np.array([x for _, x in zip(singlets, ind_s)])
        ind_t = np.array([x for _, x in zip(triplets, ind_t)])

        order_s = np.argsort(singlets)
        order_t = np.argsort(triplets)
        singlets = np.sort(singlets)
        triplets = np.sort(triplets)
        oscs = oscs[order_s]
        ind_s = ind_s[order_s]
        ind_t = ind_t[order_t]

        return (
            singlets,
            triplets,
            oscs,
            ind_s,
            ind_t,
            ss_s,
            ss_t,
            (sol_int - total_free) * 27.2114,
        )


#########################################################################################


##GETS SOC BETWEEN Sn STATE AND TRIPLETS#################################################
def pega_soc_singlet(file, n_state):
    socs = []
    _, _, _, ind_s, ind_t, _, _, _ = pega_energias("Geometries/" + file)
    order_s = np.argsort(ind_s)
    order_t = np.argsort(ind_t)
    n_state = order_s[n_state] + 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch = False
        for line in log_file:
            if (
                "Total SOC between the S"
                + str(n_state)
                + " state and excited triplet states:"
                in line
            ):
                catch = True
            elif catch and "T" in line and "(" not in line:
                try:
                    socs.append(float(line.split()[1]))
                except (IndexError, ValueError):
                    catch = False
    socs = np.array(socs)
    socs = socs[order_t]
    return socs[np.newaxis, :] * 0.12398 / 1000


#########################################################################################


##GETS SOC BETWEEN Tn STATE AND SINGLETS#################################################
def pega_soc_triplet(file, n_state):
    socs = []
    _, _, _, ind_s, ind_t, _, _, _ = pega_energias("Geometries/" + file)
    order_s = np.argsort(ind_s)
    order_t = np.argsort(ind_t)
    n_state = order_t[n_state] + 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch = False
        for line in log_file:
            if (
                "Total SOC between the S" in line
                and "state and excited triplet states:" in line
            ):
                catch = True
            elif catch and "T" + str(n_state) + " " in line and "(" not in line:
                try:
                    socs.append(float(line.split()[1]))
                except (IndexError, ValueError):
                    catch = False
    socs = np.array(socs)
    socs = socs[order_s]
    return socs[np.newaxis, :] * 0.12398 / 1000


#########################################################################################


##GETS SOC BETWEEN Tn STATE AND S0#######################################################
def pega_soc_ground(file, n_state):
    socs = []
    # _, _, _, ind_s, ind_t, _, _, _ = pega_energias('Geometries/'+file)
    # order_s = np.argsort(ind_s)
    # order_t = np.argsort(ind_t)
    n_state += 1  # order_t[n_state] + 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch = False
        for line in log_file:
            if (
                "Total SOC between the singlet ground state and excited triplet states:"
                in line
            ):
                catch = True
            elif catch and "T" + str(n_state) + " " in line and "(" not in line:
                try:
                    socs.append(float(line.split()[1]))
                except (IndexError, ValueError):
                    catch = False
            elif len(line.split()) < 2:
                catch = False
    socs = np.array(socs)
    # socs = socs[order_s]
    return socs[np.newaxis, :] * 0.12398 / 1000


#########################################################################################


##GETS SOC BETWEEN Tn STATE AND SINGLETS#################################################
def pega_soc_triplet_triplet(file, n_state):
    socs = []
    # _, _, _, ind_s, ind_t, _, _, _ = pega_energias('Geometries/'+file)
    # order_s = np.argsort(ind_s)
    # order_t = np.argsort(ind_t)
    n_state += 1  # order_t[n_state] + 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch, catch2 = False, False
        for line in log_file:
            if (
                "Total SOC between the T"
                + str(n_state)
                + " state and excited triplet states:"
                in line
            ):
                catch2 = True
            elif (
                "Total SOC between the T" in line
                and "state and excited triplet states:" in line
            ):
                catch = True
            elif (catch and "T" + str(n_state) + " " in line and "(" not in line) or (
                catch2 and "T" in line and "(" not in line
            ):
                try:
                    socs.append(float(line.split()[1]))
                except (IndexError, ValueError):
                    catch, catch2 = False, False
            elif len(line.split()) < 2:
                catch, catch2 = False, False
    socs = np.array(socs)
    # socs = socs[order_s]
    return socs[np.newaxis, :] * 0.12398 / 1000


#########################################################################################


##DECIDES WHICH FUNCTION TO USE IN ORDER TO GET SOCS#####################################
def avg_socs(files, tipo, n_state):
    col = None
    if tipo == "singlet":
        pega_soc = pega_soc_singlet
    elif tipo == "triplet":
        pega_soc = pega_soc_triplet
    elif tipo == "ground":
        pega_soc = pega_soc_ground
    elif tipo == "tts":
        pega_soc = pega_soc_triplet_triplet
    for file in files:
        socs = pega_soc(file, n_state)
        try:
            socs = socs[:, :col]
            total_socs = np.vstack((total_socs, socs))
        except NameError:
            col = socs.shape[1]
            total_socs = socs
    return total_socs


#########################################################################################


##GETS TRANSITION DIPOLE MOMENTS#########################################################
def pega_dipolos(file, ind, frase, state):
    dipoles = np.zeros((1, 1))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        dip = False
        for line in log_file:
            if frase in line:
                dip = True
            elif dip and "--" not in line:
                line = line.split()
                if line[0] == str(ind[state]):
                    dipole_line = [float(line[i]) for i in range(1, len(line))]
                    dipole_line = np.array([dipole_line])
                    try:
                        dipoles = np.vstack((dipoles, dipole_line))
                    except NameError:
                        dipoles = np.zeros((1, len(line) - 1))
                        dipoles = np.vstack((dipoles, dipole_line))
                elif line[1] == str(ind[state]) and int(line[0]) < int(line[1]):
                    dipole_line = [float(line[0])]
                    extra_line = [float(line[i]) for i in range(2, len(line))]
                    dipole_line.extend(extra_line)
                    dipole_line = np.array([dipole_line])
                    try:
                        dipoles = np.vstack((dipoles, dipole_line))
                    except NameError:
                        dipoles = np.zeros((1, len(line) - 1))
                        dipoles = np.vstack((dipoles, dipole_line))
            elif np.shape(dipoles)[0] > 1 and "---" in line:
                dip = False
    dipoles = dipoles[1:, :]
    muf = np.zeros((1, 3))
    ind = np.array(ind).astype(float)
    col = np.shape(dipoles)[1]
    if col > 4:
        for i in ind:
            index = np.where(dipoles[:, 0] == i)
            try:
                index = index[0][0]
                stack = dipoles[index, 1:-1]
                muf = np.vstack((muf, stack))
            except IndexError:
                pass
        muf = muf[1:, :]
    else:
        muf = dipoles
    return muf


#########################################################################################


##GETS TRANSITION DIPOLE MOMENTS#########################################################
def pega_oscs(files, indices, initial):
    spin = initial[0].upper()
    num = int(initial[1:]) - 1
    mapa = {"S": "Singlet", "T": "Triplet"}
    frase = "Transition Moments Between " + mapa[spin] + " Excited States"
    for i, file in enumerate(files):
        oscs = []
        ind = indices[i, num]
        ind_s = indices[i, :]
        location = np.where(ind_s == ind)[0][0]
        ind_s = ind_s[location + 1 :]
        ind = str(ind)
        with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
            dip = False
            for line in log_file:
                if frase in line:
                    dip = True
                elif dip and "--" not in line:
                    line = line.split()
                    if (line[0] == ind and int(line[1]) in ind_s) or (
                        line[1] == ind and int(line[0]) in ind_s
                    ):
                        oscs.append(float(line[5]))
                elif len(oscs) > 0 and "---" in line:
                    dip = False
            try:
                total_oscs = np.vstack((total_oscs, np.array(oscs)[np.newaxis, :]))
            except NameError:
                total_oscs = np.array(oscs)[np.newaxis, :]
    return total_oscs


#########################################################################################


##GETS SOCS BETWEEN S0 AND EACH TRIPLET SUBLEVEL#########################################
def soc_s0(file, mqn, ind_t):
    socs = np.zeros((1))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        read = False
        for line in log_file:
            if (
                "SOC between the singlet ground state and excited triplet states (ms="
                + mqn
                in line
            ):
                read = True
            elif read:
                if "T" in line and "Total" not in line:
                    line = line.split()
                    real_part = line[1].replace("(", "").replace(")", "").replace("i", "")
                    real_part = float(
                        real_part.replace("--", "+")
                        .replace("+-", "-")
                        .replace("-+", "-")
                        .replace("++", "+")
                    )
                    img_part = line[2] + line[3].replace("(", "").replace(")", "").replace(
                        "i", ""
                    )
                    img_part = float(
                        img_part.replace("--", "+")
                        .replace("+-", "-")
                        .replace("-+", "-")
                        .replace("++", "+")
                    )
                    complex_soc = real_part + img_part * 1j
                    # Qchem gives <S0|H|Tn>. We need to conjugate to get <Tn|H|S0>
                    complex_soc = complex_soc.conjugate()
                    socs = np.vstack((socs, np.array([complex_soc])))
                else:
                    read = False
    socs = socs[1:, :]
    indice = np.argsort(ind_t)
    socs = socs[indice, :]
    return socs * 0.12398 / 1000


#########################################################################################


##GETS SOCS BETWEEN Sm AND EACH Tn SUBLEVEL##############################################
def soc_t1(file, mqn, n_triplet, ind_s):
    socs = np.zeros((1))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        read = False
        for line in log_file:
            if "SOC between the S" in line and "(ms=" + mqn in line:
                read = True
            elif read:
                if "T" + str(n_triplet + 1) + "(ms=" + mqn in line:
                    line = line.split()
                    real_part = line[1].replace("(", "").replace(")", "").replace("i", "")
                    real_part = float(
                        real_part.replace("--", "+")
                        .replace("+-", "-")
                        .replace("-+", "-")
                        .replace("++", "+")
                    )
                    img_part = line[2] + line[3].replace("(", "").replace(")", "").replace(
                        "i", ""
                    )
                    img_part = float(
                        img_part.replace("--", "+")
                        .replace("+-", "-")
                        .replace("-+", "-")
                        .replace("++", "+")
                    )
                    # Qchem gives <Sn|H|Tn>. No need to conjugate.
                    complex_soc = real_part + img_part * 1j
                    socs = np.vstack((socs, np.array([complex_soc])))
    socs = socs[1:, :]
    indice = np.argsort(ind_s)
    socs = socs[indice, :]
    return socs * 0.12398 / 1000


##CALCULATES TRANSITION DIPOLE MOMENTS FOR Tn TO S0 TRANSITIONS##########################
def moment(file, ess, ets, dipss, dipts, n_triplet, ind_s, ind_t):
    # Conversion factor between a.u. = e*bohr to SI
    conversion = 8.4783533e-30
    fake_t = np.where(np.sort(ind_t) == ind_t[n_triplet])[0][0]
    ess = np.array(ess)
    ets = np.array(ets)
    ess = np.insert(ess, 0, 0)
    moments = []
    for mqn in ["1", "-1", "0"]:
        socst1 = soc_t1(file, mqn, fake_t, ind_s)
        socss0 = soc_s0(file, mqn, ind_t)
        socst1 = np.vstack((socss0[0, :], socst1))
        # Conjugate to get <S0|H|T1>
        socst1[0] = socst1[0].conjugate()
        # Now conjugate socst1
        socst1 = socst1.conjugate()
        ess = ess[: np.shape(socst1)[0]]
        if 0 in ets[n_triplet] - ess:
            return 0
        for i in [0, 1, 2]:
            part_1 = (socss0 / (0 - ets)) * dipts[:, i]
            part_1 = np.sum(part_1)
            part_2 = (socst1 / (ets[n_triplet] - ess)) * dipss[:, i]
            part_2 = np.sum(part_2)
            complex_dipole = part_1 + part_2
            # append magnitude squared
            moments.append((complex_dipole * complex_dipole.conjugate()).real)

    moments = np.array(moments)
    moments = np.sum(moments) * (conversion**2)
    return moments


#########################################################################################


##READS NUMBER OF EXCITED STATES FROM INPUT FILE#########################################
def read_cis(file):
    file = file[:-3] + "com"
    with open("Geometries/" + file, "r", encoding="utf-8") as com_file:
        for line in com_file:
            if "cis_n_roots" in line.lower():
                line = line.split()
                for elem in line:
                    try:
                        n_state = int(elem)
                        break
                    except ValueError:
                        pass
    return n_state


#########################################################################################


def phosph_osc(file, n_state, ind_s, ind_t, singlets, triplets):
    zero = ["0"]
    zero.extend(ind_s)
    total_moments = []
    ground_dipoles = pega_dipolos(file, zero, "Electron Dipole Moments of Ground State", 0)
    ground_singlet_dipoles = pega_dipolos(
        file, zero, "Transition Moments Between Ground and Singlet Excited States", 0
    )
    ground_dipoles = np.vstack((ground_dipoles, ground_singlet_dipoles))
    for n_triplet in range(n_state):
        triplet_dipoles = pega_dipolos(
            file, ind_t, "Electron Dipole Moments of Triplet Excited State", n_triplet
        )
        triplet_triplet_dipoles = pega_dipolos(
            file, ind_t, "Transition Moments Between Triplet Excited States", n_triplet
        )
        triplet_dipoles = np.vstack((triplet_dipoles, triplet_triplet_dipoles))
        # Fixing the order
        order = np.arange(1, n_state)
        order = np.insert(order, n_triplet, 0)
        triplet_dipoles = triplet_dipoles[order, :]
        moments = moment(file, singlets, triplets, ground_dipoles,
                          triplet_dipoles, n_triplet, ind_s, ind_t)
        total_moments.append(moments)
    total_moments = np.array(total_moments)
    term = E_CHARGE * (HBAR_J**2) / triplets
    osc_strength = (2 * MASS_E) * total_moments / (3 * term)
    return osc_strength[np.newaxis, :]


def get_osc_phosph(files, singlets, triplets, ind_s, ind_t):
    n_state = read_cis(files[0])
    #removed correction from phosph_osc calculation
    eng_singlets = singlets  # - (alphast2/alphaopt1)*Ss_s
    eng_triplets = triplets  # - (alphast2/alphaopt1)*Ss_t
    for j in range(singlets.shape[0]):
        tos = phosph_osc(
            files[j], n_state, ind_s[j, :], ind_t[j, :], eng_singlets[j, :], eng_triplets[j, :]
        )
        try:
            osc_strengths = np.vstack((osc_strengths, tos))
        except NameError:
            osc_strengths = tos
    return osc_strengths


##GETS ALL RELEVANT INFORMATION FROM LOG FILES###########################################
def analysis(files):
    n_state = read_cis(files[0])
    numbers = []
    for file in files:
        singlets, triplets, oscs, ind_s, ind_t, ss_s, ss_t, ground_pol = pega_energias(
            "Geometries/" + file
        )
        singlets = np.array([singlets[:n_state]])
        triplets = np.array([triplets[:n_state]])
        oscs = np.array([oscs[:n_state]])
        ss_s = np.array([ss_s[:n_state]])
        ss_t = np.array([ss_t[:n_state]])
        ind_s = np.array([ind_s[:n_state]])
        ind_t = np.array([ind_t[:n_state]])
        ground_pol = np.array([ground_pol])
        try:
            total_singlets = np.vstack((total_singlets, singlets))
            total_triplets = np.vstack((total_triplets, triplets))
            total_oscs = np.vstack((total_oscs, oscs))
            total_ss_s = np.vstack((total_ss_s, ss_s))
            total_ss_t = np.vstack((total_ss_t, ss_t))
            total_ind_s = np.vstack((total_ind_s, ind_s))
            total_ind_t = np.vstack((total_ind_t, ind_t))
            total_ground_pol = np.append(total_ground_pol, ground_pol)
        except NameError:
            total_singlets = singlets
            total_triplets = triplets
            total_oscs = oscs
            total_ss_s = ss_s
            total_ss_t = ss_t
            total_ind_s = ind_s
            total_ind_t = ind_t
            total_ground_pol = ground_pol
        numbers.append(int(file.split("-")[1]))
    numbers = np.array(numbers)[:, np.newaxis]
    return numbers, total_singlets, total_triplets, total_oscs, total_ss_s, total_ss_t, total_ground_pol, total_ind_s, total_ind_t


#########################################################################################


##PRINTS EMISSION SPECTRUM###############################################################
def printa_espectro_emi(initial, eps, refractive_index, tdm, energy, mean_y, error):
    mean_rate, error_rate = nemo.tools.calc_emi_rate(energy, mean_y, error)
    primeira = f"{'#Energy(ev)':4s} {'diff_rate':4s} {'error':4s} TDM={tdm:.3f} au\n"
    primeira += (
        f"# Total Rate {initial} -> S0: {mean_rate:5.2e} +/- {error_rate:5.2e} s^-1\n"
    )
    primeira += f"#Epsilon: {eps:.3f} nr: {refractive_index:.3f}\n"
    arquivo = nemo.tools.naming(f"differential_rate_{initial.upper()}.lx")
    with open(arquivo, "w", encoding="utf-8") as emi_spectrum:
        emi_spectrum.write(primeira)
        for i, eng in enumerate(energy):
            text = f"{eng:.6f} {mean_y[i]:.6e} {error[i]:.6e}\n"
            emi_spectrum.write(text)
    print(f"Spectrum printed in the {arquivo} file")


#######################################################################################


###CALCULATES WEIGHTED AVERAGES WHEN POSSIBLE##########################################
def means(variable, weight, ensemble_mean=False):
    if ensemble_mean:
        try:
            mean = np.mean(variable, axis=0)
        except IndexError:
            mean = np.mean(variable)
    else:
        try:
            mean = np.average(variable, axis=0, weights=weight)
        except IndexError:
            mean = np.average(variable, axis=0)
    return mean


########################################################################################


###FORMATS RATES AND ERRORS IN THE SAME EXPONENT########################################
def format_rate(rate, delta_rate):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exp = np.nan_to_num(np.floor(np.log10(rate)))
    try:
        exp[exp < -99] = -99
    except TypeError:
        exp = max(exp, -99)
    pre_r = rate / 10**exp
    pre_dr = delta_rate / 10**exp
    return pre_r, pre_dr, exp


#########################################################################################


###SAVES ENSEMBLE DATA#################################################################
def gather_data(initial, save=True):
    files = [i for i in os.listdir("Geometries") if ".log" in i]
    files = check_normal(files)
    files = sorted(files, key=lambda pair: float(pair.split("-")[1]))
    n_state = int(initial[1:]) - 1
    eps_i, nr_i = nemo.tools.get_nr()
    kbt = nemo.tools.detect_sigma()
    if "s" in initial.lower():
        numbers, singlets, triplets, oscs, ss_s, ss_t, ground_pol, ind_s, ind_t = analysis(
            files
        )
        if "s0" == initial.lower():
            label_oscs = [f"osc_s{i+1}" for i in range(oscs.shape[1])]
        else:
            # Oscs       = Oscs[:,n_state][:,np.newaxis]
            label_oscs = [f"osce_s{n_state+1+i}" for i in range(oscs.shape[1])]
            noscs = pega_oscs(files, ind_s, initial)
            label_oscs.extend([f"osc_s{n_state+2+i}" for i in range(noscs.shape[1])])
            oscs = np.hstack((oscs, noscs))
        try:
            header7 = []
            for i in range(singlets.shape[1]):
                socs_partial = avg_socs(files, "singlet", i)
                header7.extend(
                    [f"soc_s{i+1}_t{j}" for j in range(1, 1 + socs_partial.shape[1])]
                )
                try:
                    socs_complete = np.hstack((socs_complete, socs_partial))
                except NameError:
                    socs_complete = socs_partial
        except IndexError:
            pass
    else:
        numbers, singlets, triplets, _, ss_s, ss_t, ground_pol, ind_s, ind_t = analysis(files)
        oscs = get_osc_phosph(files, singlets, triplets, ind_s, ind_t)
        # Oscs       = Oscs[:,n_state][:,np.newaxis]
        label_oscs = [f"osce_t{n_state+1+i}" for i in range(oscs.shape[1])]
        noscs = pega_oscs(files, ind_t, initial)
        oscs = np.hstack((oscs, noscs))
        label_oscs.extend([f"osc_t{n_state+2+i}" for i in range(noscs.shape[1])])
        try:
            header7 = []
            for i in range(triplets.shape[1]):
                socs_partial = np.hstack(
                    (
                        avg_socs(files, "ground", i),
                        avg_socs(files, "triplet", i),
                        avg_socs(files, "tts", i),
                    )
                )
                indices = [
                    j + 1 for j in range(triplets.shape[1]) if j != i
                ]  # Removed Tn to Tn transfers
                header7.extend([f"soc_t{i+1}_s0"])
                header7.extend(
                    [f"soc_t{i+1}_s{j}" for j in range(1, 1 + singlets.shape[1])]
                )
                header7.extend([f"soc_t{i+1}_t{j}" for j in indices])
                try:
                    socs_complete = np.hstack((socs_complete, socs_partial))
                except NameError:
                    socs_complete = socs_partial
        except IndexError:
            pass
    header = ["geometry"]
    header.extend(["e_s" + str(i) for i in range(1, 1 + singlets.shape[1])])
    header.extend(["e_t" + str(i) for i in range(1, 1 + triplets.shape[1])])
    header.extend(["d_s" + str(i) for i in range(1, 1 + ss_s.shape[1])])
    header.extend(["d_t" + str(i) for i in range(1, 1 + ss_t.shape[1])])
    header.extend(["gp"])
    header.extend(label_oscs)
    try:
        header.extend(header7)
        data = np.hstack(
            (
                numbers,
                singlets,
                triplets,
                ss_s,
                ss_t,
                ground_pol[:, np.newaxis],
                oscs,
                socs_complete,
            )
        )
    except NameError:
        data = np.hstack(
            (numbers, singlets, triplets, ss_s, ss_t, ground_pol[:, np.newaxis], oscs)
        )
    arquivo = f"Ensemble_{initial.upper()}_.lx"
    data = pd.DataFrame(data, columns=header)
    # add 'ensemble', 'kbT', 'nr', 'eps' columns with constant values
    # values are initial.upper(), kbT, nr_i, eps_i
    data["ensemble"] = initial.upper()
    data["kbT"] = kbt
    data["nr"] = nr_i
    data["eps"] = eps_i
    # make these the first columns
    cols = data.columns.tolist()
    cols = cols[-4:] + cols[:-4]
    data = data[cols]
    if save:
        data.to_csv(arquivo, index=False)
    return data


#######################################################################################


###PRINTS RATES AND EMISSION SPECTRUM##################################################
def export_results(data, emission, dielec):
    data = data.copy()
    initial = data["Transition"][0].split(">")[0][:-1]
    printa_espectro_emi(
        initial,
        dielec[0],
        dielec[1],
        emission["TDM"][0],
        emission["Energy"].values,
        emission["Diffrate"].values,
        emission["Error"].values,
    )
    pre_r, pre_dr, exp = format_rate(data["Rate(s^-1)"], data["Error(s^-1)"])
    rate = [f"{pre_r[i]:5.2f}e{exp[i]:+03.0f}" for i in range(len(pre_r))]
    error = [f"{pre_dr[i]:5.2f}e{exp[i]:+03.0f}" for i in range(len(pre_dr))]
    headers = [
        i
        for i in data.columns.values
        if i != "Rate(s^-1)" and i != "Error(s^-1)" and i != "Transition"
    ]
    for header in headers:
        if header == "Prob(%)" or header == "AvgConc(%)":
            data[header] = data[header].map("{:.1f}".format)
        else:
            data[header] = data[header].map("{:.3f}".format)
    data["Rate(s^-1)"] = rate
    data["Error(s^-1)"] = error
    arquivo = nemo.tools.naming(f"rates_{initial}_.lx")
    solvent = f"#Epsilon: {dielec[0]:.3f} nr: {dielec[1]:.3f}\n"
    with open(arquivo, "w", encoding="utf-8") as rate_file:
        rate_file.write(solvent + data.to_string(header=True, index=False))
    print(f"Rates are written in the {arquivo} file")


#######################################################################################


def reorder(initial_state, final_state, ss_i, ss_f, socs):
    argsort = np.argsort(initial_state, axis=1)
    initial_state = np.take_along_axis(initial_state, argsort, axis=1)
    ss_i = np.take_along_axis(ss_i, argsort, axis=1)
    corredor = int(np.sqrt(socs.shape[1]))
    socs_complete = socs.reshape((socs.shape[0], corredor, corredor))
    for j in range(socs_complete.shape[1]):
        socs_complete[:, j, :] = np.take_along_axis(
            socs_complete[:, j, :], argsort, axis=1
        )
    argsort = np.argsort(final_state, axis=1)
    final_state = np.take_along_axis(final_state, argsort, axis=1)
    ss_f = np.take_along_axis(ss_f, argsort, axis=1)
    for j in range(socs_complete.shape[1]):
        socs_complete[:, :, j] = np.take_along_axis(
            socs_complete[:, :, j], argsort, axis=1
        )
    return initial_state, final_state, ss_i, ss_f, socs_complete


def fix_absent_soc(data):
    columns = data.columns.values
    # check if at least one column contains soc_
    if any("soc_" in i for i in columns):
        return data
    else:
        singlets = [i.split("_")[1] for i in columns if "e_s" in i and "osc" not in i]
        triplets = [i.split("_")[1] for i in columns if "e_t" in i and "osc" not in i]
        for singlet in singlets:
            for triplet in triplets:
                data[f"soc_{singlet}_{triplet}"] = 0
    return data

def x_values(mean,std):
    left = max(np.min(mean - 2 * std), 0.01)
    right = np.max(mean + 2 * std)
    x_axis = np.linspace(left, right, int((right - left) / 0.01))
    return x_axis

def sorting_parameters(argsort,*args):
    argsort = np.argsort(args[0], axis=1)
    args = list(args)
    for arg in args:
        arg = np.take_along_axis(arg, argsort, axis=1)
    return tuple(args)

def check_number_geoms(data):
    number_geoms = data.shape[0]
    coms = nemo.tools.start_counter()
    if number_geoms != coms:
        print(
            (f"There are {coms} inputs and just {number_geoms} log files. "
            "Something is not right! Computing the rates anyway...")
        )
    return number_geoms


def fetch(data, criteria_list):
    filtered_data = data[
        [i for i in data.columns.values if all(c in i for c in criteria_list if not c.startswith('-')) and not any(c[1:] in i for c in criteria_list if c.startswith('-'))]
    ].to_numpy()
    return filtered_data

###CALCULATES ISC AND EMISSION RATES & SPECTRA#########################################
def rates(initial, dielec, data=None, ensemble_average=False, detailed=False):
    if data is None:
        data = gather_data(initial, save=True)
        eps_i, nr_i = nemo.tools.get_nr()
        kbt = nemo.tools.detect_sigma()
    else:
        eps_i = data["eps"][0]
        nr_i = data["nr"][0]
        kbt = data["kbT"][0]
    eps, refractive_index = dielec[0], dielec[1]
    alphast1 = nemo.tools.get_alpha(eps_i)
    alphast2 = nemo.tools.get_alpha(eps)
    alphaopt1 = nemo.tools.get_alpha(nr_i**2)
    alphaopt2 = nemo.tools.get_alpha(refractive_index**2)
    n_state = int(initial[1:]) - 1
    initial = initial.lower()

    data = fix_absent_soc(data)

    # Emission Calculations
    lambda_be = (alphast2 / alphast1 - alphaopt2 / alphast1) * data["gp"].to_numpy()
    l_total = np.sqrt(2 * lambda_be * kbt + kbt**2)
    energies = data.filter(regex=f"^e_{initial[0]}").to_numpy()
    delta_emi = (
        energies
        - (alphast2 / alphaopt1) * data.filter(regex=f"^d_{initial[0]}").to_numpy()
    )
    constante = (
        (refractive_index**2)
        * (E_CHARGE**2)
        / (2 * np.pi * HBAR_EV * MASS_E * (LIGHT_SPEED**3) * EPSILON_0)
    )
    if "t" in initial:
        constante *= 1 / 3
    oscs = data.filter(regex="^osce_").to_numpy()
    # socs_s0 will be used for T1>S0 ISC calculation
    socs_s0 = data[
            [i for i in data.columns.values if "soc_t" in i and "s0" in i]
        ].to_numpy()
    delta_emi, oscs, energies, *socs_s0 = sorting_parameters(delta_emi,oscs,energies,socs_s0)
    delta_emi = delta_emi[:, n_state]
    oscs = oscs[:, n_state]
    energies = energies[:, n_state]
    socs_s0 = socs_s0[:, n_state]
    espectro = constante * ((delta_emi - lambda_be) ** 2) * oscs
    tdm = nemo.tools.calc_tdm(oscs, energies, espectro)
    x_axis = x_values(delta_emi,l_total)
    y_axis = espectro[:, np.newaxis] * nemo.tools.gauss(
        x_axis, delta_emi[:, np.newaxis], l_total[:, np.newaxis]
    )
    number_geoms = y_axis.shape[0]
    mean_y = np.sum(y_axis, axis=0) / number_geoms
    error = np.sqrt(np.sum((y_axis - mean_y) ** 2, axis=0) / (number_geoms * (number_geoms - 1)))
    emi_rate, emi_error = nemo.tools.calc_emi_rate(x_axis, mean_y, error)
    gap_emi = means(delta_emi, espectro, ensemble_average)
    mean_sigma_emi = means(l_total, espectro, ensemble_average)
    mean_part_emi = (100 / number_geoms) / means(
        espectro / np.sum(espectro), espectro, ensemble_average
    )
    emi = np.hstack((x_axis[:, np.newaxis], mean_y[:, np.newaxis], error[:, np.newaxis]))
    emi = pd.DataFrame(emi, columns=["Energy", "Diffrate", "Error"])
    emi.insert(0, "TDM", tdm)

    # Checks number of logs
    if data is None:
        check_number_geoms(data)
    # Intersystem Crossing Rates
    singlets = data[
        [i for i in data.columns.values if "e_s" in i and "osc" not in i]
    ].to_numpy()
    triplets = data[
        [i for i in data.columns.values if "e_t" in i and "osc" not in i]
    ].to_numpy()
    ss_s = data[[i for i in data.columns.values if "d_s" in i]].to_numpy()
    ss_t = data[[i for i in data.columns.values if "d_t" in i]].to_numpy()
    if "s" in initial:
        initial_state = singlets - (alphast2 / alphaopt1) * ss_s
        final_state = triplets - (alphaopt2 / alphaopt1) * ss_t
        socs_complete = data[
            [i for i in data.columns.values if "soc_s" in i]
        ].to_numpy()
        initial_state, final_state, ss_s, ss_t, socs_complete = reorder(
            initial_state, final_state, ss_s, ss_t, socs_complete
        )
        initial_state = initial_state[:, n_state]
        socs_complete = socs_complete[:, n_state, :]
        delta = final_state - np.repeat(
            initial_state[:, np.newaxis], final_state.shape[1], axis=1
        )
        lambda_b = (alphast2 / alphaopt1 - alphaopt2 / alphaopt1) * ss_t
        final = [
            i.split("_")[2].upper()
            for i in data.columns.values
            if "soc_" + initial.lower() + "_" in i
        ]
        ##FOR WHEN IC IS AVAILABLE
        # socs_complete = np.hstack((socs_complete,0.0001*np.ones((Singlets.shape[0],Singlets.shape[1]-1))))
        # delta_ss = Singlets + np.repeat((alphast2/alphaopt1)*Ss_s[:,n_state][:,np.newaxis] - Singlets[:,n_state][:,np.newaxis],Singlets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_s    #Sm (final) - Sn (initial) + lambda_b
        # indices  = [i for i in range(Singlets.shape[1]) if i != n_state] #Removed Sn to Sn transfers
        # delta    = np.hstack((delta,delta_ss[:,indices]))
        # lambda_bt= (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_s
        # lambda_b = np.hstack((lambda_b,lambda_bt[:,indices]))
    elif "t" in initial:
        # Tn to Sm ISC
        initial_state = triplets - (alphast2 / alphaopt1) * ss_t
        final_state = singlets - (alphaopt2 / alphaopt1) * ss_s
        socs_complete = data[
            [
                i
                for i in data.columns.values
                if "soc_t" in i and "s0" not in i and i.count("t") == 1
            ]
        ].to_numpy()
        initial_state, final_state, ss_t, ss_s, socs_complete = reorder(
            initial_state, final_state, ss_t, ss_s, socs_complete
        )
        initial_state = initial_state[:, n_state]
        socs_complete = socs_complete[:, n_state, :]
        delta = final_state - np.repeat(
            initial_state[:, np.newaxis], final_state.shape[1], axis=1
        )
        lambda_b = (alphast2 / alphaopt1 - alphaopt2 / alphaopt1) * ss_s
        final = [
            i.split("_")[2].upper()
            for i in data.columns.values
            if "soc_" + initial.lower() + "_" in i and i.count("t") == 1
        ]
        # Tn to S0 ISC
        socs_complete = np.hstack((socs_s0[:, np.newaxis], socs_complete))
        delta = np.hstack((delta_emi[:, np.newaxis], delta))
        lambda_b = np.hstack((lambda_be[:, np.newaxis], lambda_b))
        # Tn to Tm ISC
        # delta_tt = Triplets + np.repeat((alphast2/alphaopt1)*Ss_t[:,n_state][:,np.newaxis] - Triplets[:,n_state][:,np.newaxis],Triplets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_t    #Tm (final) - Tn (initial) + lambda_b
         #Removing Tn to Tn transfers
        # indices  = [i for i in range(Triplets.shape[1]) if i != n_state]
        # delta    = np.hstack((delta,delta_tt[:,indices]))
        # lambda_bt= (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_t
        # lambda_b = np.hstack((lambda_b,lambda_bt[:,indices]))
        # final.extend([i.upper()[4:] for i in data.columns.values if 'soc_t' in i])

    sigma = np.sqrt(2 * lambda_b * kbt + kbt**2)
    y_axis = (2 * np.pi / HBAR_EV) * (socs_complete**2) * nemo.tools.gauss(delta, 0, sigma)
    # hstack y and espectro
    individual = np.hstack((espectro[:, np.newaxis], y_axis))
    individual /= individual.shape[0]
    number_geoms = y_axis.shape[0]
    rate = np.sum(y_axis, axis=0) / number_geoms
    total = emi_rate + np.sum(rate)
    # Error estimate
    error = np.sqrt(np.sum((y_axis - rate) ** 2, axis=0) / (number_geoms * (number_geoms - 1)))

    results = np.array(
        [
            [
                emi_rate,
                emi_error,
                100 * emi_rate / total,
                gap_emi,
                np.nan,
                mean_sigma_emi,
                mean_part_emi,
            ]
        ]
    )
    mean_gap = means(delta, y_axis, ensemble_average)[:, np.newaxis]
    mean_soc = 1000 * means(socs_complete, y_axis, ensemble_average)[:, np.newaxis]
    mean_sigma = means(sigma, y_axis, ensemble_average)[:, np.newaxis]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mean_part = np.nan_to_num(100 * rate / means(y_axis, y_axis, ensemble_average))
    rate = rate[:, np.newaxis]
    error = error[:, np.newaxis]
    labels = [f"{initial.upper()}->S0"] + [f"{initial.upper()}~>{j}" for j in final]

    # make a dataframe with Ss_s and Ss_t
    breakdown = pd.DataFrame(
        np.hstack((ss_s / alphaopt1, ss_t / alphaopt1)),
        columns=[f"chi_s{i+1}" for i in range(ss_s.shape[1])]
        + [f"chi_t{i+1}" for i in range(ss_t.shape[1])],
    )
    # append a columns with energies named eng
    breakdown["eng"] = delta_emi
    breakdown["sigma"] = l_total
    # append individual to df, use labels as columns
    breakdown = pd.concat([breakdown, pd.DataFrame(individual, columns=labels)], axis=1)

    labels = np.array(labels)
    results_isc = np.hstack(
        (
            rate,
            error,
            100 * rate / total,
            mean_gap,
            mean_soc,
            mean_sigma,
            mean_part[:, np.newaxis],
        )
    )
    results = np.vstack((results, results_isc))
    results = pd.DataFrame(
        results,
        columns=[
            "Rate(s^-1)",
            "Error(s^-1)",
            "Prob(%)",
            "AvgDE+L(eV)",
            "AvgSOC(meV)",
            "AvgSigma(eV)",
            "AvgConc(%)",
        ],
    )
    results.insert(0, "Transition", labels)
    if detailed:
        return results, emi, breakdown
    else:
        return results, emi


#########################################################################################

def save_absorption_spectrum(initial,eps, refractive_index, x_axis, mean_y, sigma, labels):
    arquivo = nemo.tools.naming(f"cross_section_{initial.upper()}_.lx")
    primeira = (f"{'Energy(ev)':8s} {'cross_section(A^2)':8s} {'error(A^2)':8s}\n"
                f"Absorption from State: {initial.upper()}\n"
                f"Epsilon: {eps:.3f} nr: {refractive_index:.3f}\n")
    labels = [f"{i:14s}" for i in labels]
    primeira += " ".join(labels)
    fmt = ["%14.6e" for i in range(0, mean_y.shape[1])]
    fmt = " ".join(fmt)
    np.savetxt(
        arquivo,
        np.hstack((x_axis[:, np.newaxis], mean_y, sigma[:, np.newaxis])),
        fmt="%14.6f " + fmt + " %14.6e",
        header=primeira,
    )
    print(f"Spectrum printed in the {arquivo} file")

def make_breakdown(initial, spin, num, oscs, deltae_lambda, lambda_neq, alphaopt1, l_total):
    # concatenate oscs[:,:,0] with DE[:,:,0] and ds
    colunas = [
        f"{initial.upper()}->{spin.upper()}{i}"
        for i in range(num + 1, num + oscs.shape[1] + 1)
    ]
    colunas += [
        f"eng_{spin}{i}"
        for i in range(num + 1, num + deltae_lambda.shape[1] + 1)
    ]
    colunas += [
        f"chi_{spin}{i}"
        for i in range(num + 1, num + lambda_neq.shape[1] + 1)
    ]
    colunas += [
        f"sigma_{spin}{i}"
        for i in range(num + 1, num + l_total.shape[1] + 1)
    ]
    # concatenate oscs[:,:,0] with DE[:,:,0] and ds
    breakdown = pd.DataFrame(
        np.hstack((oscs[:, :, 0], deltae_lambda[:, :, 0], lambda_neq / alphaopt1,l_total[:, :, 0])),
        columns=colunas,
    )
    return breakdown

def another_dimension(nstates,*args):
    args = list(args)
    for arg in args:
        arg = arg[:, :nstates, np.newaxis]
    return tuple(args)


###COMPUTES ABSORPTION SPECTRA###########################################################
def absorption(initial, dielec, data=None, save=False, detailed=False, nstates=-1):
    if data is None:
        data = gather_data(initial, save=True)
        _, nr_i = nemo.tools.get_nr()
        kbt = nemo.tools.detect_sigma()
    else:
        #eps_i = data["eps"][0]
        nr_i = data["nr"][0]
        kbt = data["kbT"][0]
    eps, refractive_index = dielec[0], dielec[1]
    #alphast1 = nemo.tools.get_alpha(eps_i)
    alphast2 = nemo.tools.get_alpha(eps)
    alphaopt1 = nemo.tools.get_alpha(nr_i**2)
    alphaopt2 = nemo.tools.get_alpha(refractive_index**2)
    initial = initial.lower()
    constante = (
        (np.pi * (E_CHARGE**2) * HBAR_EV)
        / (2 * refractive_index * MASS_E * LIGHT_SPEED * EPSILON_0)
        * 1e20
    )
    spin = initial[0]
    num = int(initial[1:])
    engs = [i for i in data.columns if f"e_{spin}" in i and "osc" not in i and int(i.split("_")[1][1:]) > num]
    lambda_neq = [i for i in data.columns if f"d_{spin}" in i and int(i.split("_")[1][1:]) > num]
    oscs = [i for i in data.columns if "osc_" in i and int(i.split("_")[1][1:]) > num]
    engs = data[engs].to_numpy()
    lambda_neq = data[lambda_neq].to_numpy()
    oscs = data[oscs].to_numpy()
    lambda_b = (alphast2 / alphaopt1 - alphaopt2 / alphaopt1) * lambda_neq
    if initial == "s0":
        deltae_lambda = engs - (alphaopt2 / alphaopt1) * lambda_neq
    else:
        base = data[f"e_{initial}"].to_numpy()[:, np.newaxis]
        lambda_neq_base = data[f"d_{initial}"].to_numpy()[:, np.newaxis]
        deltae_lambda = (
            engs
            - (alphaopt2 / alphaopt1) * lambda_neq
            - np.repeat(base - (alphast2 / alphaopt1) * lambda_neq_base, engs.shape[1], axis=1)
        )

    # Sorting states by energy
    deltae_lambda, oscs, lambda_b, *lambda_neq = sorting_parameters(deltae_lambda, oscs, lambda_b, lambda_neq)
    l_total = np.sqrt(2 * lambda_b * kbt + kbt**2)
    x_axis = x_values(deltae_lambda,l_total)
    if nstates == -1:
        nstates = deltae_lambda.shape[1]
    # Add extra dimension to DE and Ltotal to match x shape
    deltae_lambda, l_total, oscs, *lambda_b = another_dimension(nstates,deltae_lambda, l_total, oscs, lambda_b)
    y_axis = constante * oscs * nemo.tools.gauss(x_axis, deltae_lambda, l_total)
    number_geoms = oscs.shape[0]
    mean_y = np.sum(y_axis, axis=0) / number_geoms
    # Error estimate
    sigma = np.sqrt(np.sum((y_axis - mean_y) ** 2, axis=0) / (number_geoms * (number_geoms - 1)))
    mean_y = mean_y.T
    sigma = sigma.T
    total = np.sum(mean_y, axis=1)
    sigma = np.sum(sigma, axis=1)
    # append total to mean_y
    mean_y = np.append(mean_y, total[:, np.newaxis], axis=1)

    # make dataframe
    labels = [
        initial[0].upper() + str(int(initial[1:]) + i + 1)
        for i in range(0, mean_y.shape[1] - 1)
    ]
    labels += ["Total", "Error"]
    labels = ["Energy"] + labels
    abs_spec = pd.DataFrame(
        np.hstack((x_axis[:, np.newaxis], mean_y, sigma[:, np.newaxis])), columns=labels
    )

    if save:
        save_absorption_spectrum(initial,eps, refractive_index, x_axis, mean_y, sigma, labels)
    if detailed:
        breakdown = make_breakdown(initial, spin, num, oscs, deltae_lambda, lambda_neq, alphaopt1, l_total)
        return abs_spec, breakdown
    return abs_spec


#########################################################################################
