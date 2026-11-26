#!/usr/bin/env python3
import os
import re
import warnings
import numpy as np
import pandas as pd
import nemo.tools
import nemo.parser
import nemo.eom

# pylint: disable=unbalanced-tuple-unpacking

LIGHT_SPEED = nemo.parser.LIGHT_SPEED
HBAR_EV = nemo.parser.HBAR_EV
HBAR_J = nemo.parser.HBAR_J
BOLTZ_EV = nemo.parser.BOLTZ_EV
E_CHARGE = nemo.parser.E_CHARGE
MASS_E = nemo.parser.MASS_E
EPSILON_0 = nemo.parser.EPSILON_0


##RETURNS LIST OF LOG FILES WITH NORMAL TERMINATION######################################
def check_normal(files):
    normal, abnormal = [], []
    add = False
    for file in files:
        with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
            for line in log_file:
                #if (
                #    "TDDFT/TDA Excitation Energies" in line
                #    or "TDDFT Excitation Energies" in line
                #):
                #    exc = True
                #elif "Excited state" in line and exc:
                #    eng = float(line.split()[7])
                #    if eng < 0:
                #        add = False
                #        abnormal.append(file)
                #    else:
                #        add = True
                #    exc = False
                if "Have a nice day" in line:# and add:
                    normal.append(file)
    if len(abnormal) > 0:
        print(f"Warning! Negative transition energies detected in {len(abnormal)} files:")
        for file in abnormal:
            print(file)
        print("They will not be considered in the analysis")
    return normal


#########################################################################################


#########################################################################################


##READS NUMBER OF EXCITED STATES FROM INPUT FILE#########################################
def read_cis(file):
    file = file[:-3] + "com"
    text = open("Geometries/" + file, "r", encoding="utf-8").read().lower()
    if "cis_n_roots" in text:
        n_state = int(text.split("cis_n_roots")[1].split()[0])
        calculation_type = "tddft"
    elif "ee_singlets" in text:
        n_state = int(text.split("ee_singlets")[1].split()[0])
        calculation_type = "eom-ccsd"
    return n_state, calculation_type


#########################################################################################


def get_osc_phosph(files, singlets, triplets, n_state, phosph_osc):
    # removed correction from phosph_osc calculation
    eng_singlets = singlets  # - (alphast2/alphaopt1)*Ss_s
    eng_triplets = triplets  # - (alphast2/alphaopt1)*Ss_t
    for j in range(singlets.shape[0]):
        tos = phosph_osc(
            files[j],
            n_state,
            eng_singlets[j, :],
            eng_triplets[j, :],
        )
        try:
            osc_strengths = np.vstack((osc_strengths, tos))
        except NameError:
            osc_strengths = tos
    return osc_strengths


##GETS ALL RELEVANT INFORMATION FROM LOG FILES###########################################
def analysis(files, n_state, get_energies):
    numbers = []
    for file in files:
        (
            singlets,
            triplets,
            oscs,
            ss_s,
            ss_t,
            ground_pol,
        ) = get_energies("Geometries/" + file)
        singlets = np.array([singlets[:n_state]])
        triplets = np.array([triplets[:n_state]])
        oscs = np.array([oscs[:n_state]])
        ss_s = np.array([ss_s[:n_state]])
        ss_t = np.array([ss_t[:n_state]])
        ground_pol = np.array([ground_pol])
        try:
            total_singlets = np.vstack((total_singlets, singlets))
            total_triplets = np.vstack((total_triplets, triplets))
            total_oscs = np.vstack((total_oscs, oscs))
            total_ss_s = np.vstack((total_ss_s, ss_s))
            total_ss_t = np.vstack((total_ss_t, ss_t))
            total_ground_pol = np.append(total_ground_pol, ground_pol)
        except NameError:
            total_singlets = singlets
            total_triplets = triplets
            total_oscs = oscs
            total_ss_s = ss_s
            total_ss_t = ss_t
            total_ground_pol = ground_pol
        numbers.append(int(file.split("-")[1]))
    numbers = np.array(numbers)[:, np.newaxis]
    return (
        numbers,
        total_singlets,
        total_triplets,
        total_oscs,
        total_ss_s,
        total_ss_t,
        total_ground_pol,
    )


#########################################################################################


##PRINTS EMISSION SPECTRUM###############################################################
def printa_espectro_emi(initial, eps, refractive_index, emission):
    mean_rate, error_rate = emission.rate, emission.error
    energy = emission["Energy"].to_numpy()
    mean_y = emission["Diffrate"].to_numpy()
    error = emission["Error"].to_numpy()
    tdm = emission.tdm
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
        return np.mean(variable, axis=0)
    weight = np.asarray(weight)
    weight = np.where(weight.sum(axis=0) == 0, 1.0, weight)
    return np.average(variable, axis=0, weights=weight)


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

get_energies = {'tddft': nemo.parser.pega_energias, 'eom-ccsd': nemo.eom.pega_energias}
get_oscs = {'tddft': nemo.parser.pega_oscs, 'eom-ccsd': nemo.eom.pega_oscs}
get_avg_socs = {'tddft': nemo.parser.avg_socs, 'eom-ccsd': nemo.eom.avg_socs}
get_phosph_osc = {'tddft': nemo.parser.phosph_osc, 'eom-ccsd': nemo.eom.phosph_osc}

###SAVES ENSEMBLE DATA#################################################################
def gather_data(initial, save=True):
    formats = {}
    files = [i for i in os.listdir("Geometries") if ".log" in i]
    files = check_normal(files)
    files = sorted(files, key=lambda pair: float(pair.split("-")[1]))
    n_state = int(initial[1:]) - 1
    eps_i, nr_i = nemo.tools.get_nr()
    kbt = nemo.tools.detect_sigma()
    total_states, calculation_type = read_cis(files[0])
    if "s" in initial.lower():
        (
            numbers,
            singlets,
            triplets,
            oscs,
            ss_s,
            ss_t,
            ground_pol,
        ) = analysis(files, total_states, get_energies[calculation_type])
        if "s0" == initial.lower():
            label_oscs = [f"osc_s{i+1}" for i in range(oscs.shape[1])]
            any({formats.update({f"osc_s{i+1}": "{:.5e}"}) for i in range(oscs.shape[1])})
        else:
            label_oscs = [f"osce_s{i+1}" for i in range(oscs.shape[1])]
            any({formats.update({f"osce_s{i+1}": "{:.5e}"}) for i in range(oscs.shape[1])})
            noscs = get_oscs[calculation_type](files, initial, singlets.shape[1])
            label_oscs.extend([f"osc_s{n_state+2+i}" for i in range(noscs.shape[1])])
            any({formats.update({f"osc_s{n_state+2+i}": "{:.5e}"}) for i in range(noscs.shape[1])})
            oscs = np.hstack((oscs, noscs))
        try:
            header7 = []
            for i in range(singlets.shape[1]):
                socs_partial = get_avg_socs[calculation_type](files, "singlet", i, triplets.shape[1])
                header7.extend(
                    [f"soc_s{i+1}_t{j}" for j in range(1, 1 + socs_partial.shape[1])]
                )
                any({formats.update({f"soc_s{i+1}_t{j}": "{:.5e}"}) for j in range(1, 1 + socs_partial.shape[1])})
                try:
                    socs_complete = np.hstack((socs_complete, socs_partial))
                except NameError:
                    socs_complete = socs_partial
        except IndexError:
            pass
    else:
        numbers, singlets, triplets, _, ss_s, ss_t, ground_pol = analysis(
            files, total_states, get_energies[calculation_type]
        )
        oscs = get_osc_phosph(files, singlets, triplets, total_states, get_phosph_osc[calculation_type])
        label_oscs = [f"osce_t{i+1}" for i in range(oscs.shape[1])]
        any({formats.update({f"osce_t{i+1}": "{:.5e}"}) for i in range(oscs.shape[1])})
        noscs =  get_oscs[calculation_type](files, initial, triplets.shape[1])
        oscs = np.hstack((oscs, noscs))
        label_oscs.extend([f"osc_t{n_state+2+i}" for i in range(noscs.shape[1])])
        any({formats.update({f"osc_t{n_state+2+i}": "{:.5e}"}) for i in range(noscs.shape[1])})
        try:
            header7 = []
            for i in range(triplets.shape[1]):
                socs_partial = np.hstack(
                    (
                        get_avg_socs[calculation_type](files, "ground", i),
                        get_avg_socs[calculation_type](files, "triplet", i, singlets.shape[1]),
                        get_avg_socs[calculation_type](files, "tts", i, triplets.shape[1]-1),
                    )
                )
                indices = [
                    j + 1 for j in range(triplets.shape[1]) if j != i
                ]  # Removed Tn to Tn transfers
                header7.extend([f"soc_t{i+1}_s0"])
                formats[f"soc_t{i+1}_s0"] = "{:.5e}"
                header7.extend(
                    [f"soc_t{i+1}_s{j}" for j in range(1, 1 + singlets.shape[1])]
                )
                any({formats.update({f"soc_t{i+1}_s{j}": "{:.5e}"}) for j in range(1, 1 + singlets.shape[1])})
                header7.extend([f"soc_t{i+1}_t{j}" for j in indices])
                any({formats.update({f"soc_t{i+1}_t{j}": "{:.5e}"}) for j in indices})
                try:
                    socs_complete = np.hstack((socs_complete, socs_partial))
                except NameError:
                    socs_complete = socs_partial
        except IndexError:
            pass
    header = ["geometry"]
    formats["geometry"] = "{:.0f}"
    header.extend([f"e_s{i}" for i in range(1, 1 + singlets.shape[1])])
    any({formats.update({f"e_s{i}": "{:.4f}"}) for i in range(1, 1 + singlets.shape[1])})
    header.extend([f"e_t{i}" for i in range(1, 1 + triplets.shape[1])])
    any({formats.update({f"e_t{i}": "{:.4f}"}) for i in range(1, 1 + triplets.shape[1])})
    header.extend([f"d_s{i}" for i in range(1, 1 + ss_s.shape[1])])
    any({formats.update({f"d_s{i}": "{:.4f}"}) for i in range(1, 1 + ss_s.shape[1])})
    header.extend([f"d_t{i}" for i in range(1, 1 + ss_t.shape[1])])
    any({formats.update({f"d_t{i}": "{:.4f}"}) for i in range(1, 1 + ss_t.shape[1])})
    header.extend(["gp"])
    formats["gp"] = "{:.4f}"
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
    formats["ensemble"] = "{:s}"
    data["kbT"] = kbt
    formats["kbT"] = "{:.4f}"
    data["nr"] = nr_i
    formats["nr"] = "{:.3f}"
    data["eps"] = eps_i
    formats["eps"] = "{:.3f}"
    # make these the first columns
    cols = data.columns.tolist()
    cols = cols[-4:] + cols[:-4]
    data = data[cols]
    if save:
        # Create a temporary copy of the DataFrame
        temp_data = data.copy()
        #Apply formats
        for column, fmt in formats.items():
            temp_data[column] = temp_data[column].map(fmt.format)
        temp_data.to_csv(arquivo, index=False)
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
        emission
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
    order_i = np.argsort(initial_state, axis=1, kind="stable")
    order_f = np.argsort(final_state, axis=1, kind="stable")
    socs_complete = socs.reshape(
        (socs.shape[0], initial_state.shape[1], final_state.shape[1])
    )
    socs_complete = np.take_along_axis(socs_complete, order_i[:, :, None], axis=1)
    socs_complete = np.take_along_axis(socs_complete, order_f[:, None, :], axis=2)
    initial_state = np.take_along_axis(initial_state, order_i, axis=1)
    final_state = np.take_along_axis(final_state, order_f, axis=1)
    ss_i = np.take_along_axis(ss_i, order_i, axis=1)
    ss_f = np.take_along_axis(ss_f, order_f, axis=1)
    return initial_state, final_state, ss_i, ss_f, socs_complete


def fix_absent_soc(data, spin="s"):
    source = fetch(data, [rf"^e_{spin}\d+$"]).shape[1]
    target_spin = "t" if spin == "s" else "s"
    target = fetch(data, [rf"^e_{target_spin}\d+$"]).shape[1]
    first_target = 0 if spin == "t" else 1
    columns = [
        f"soc_{spin}{i}_{target_spin}{j}"
        for i in range(1, source + 1)
        for j in range(first_target, target + 1)
    ]
    missing = [column for column in columns if column not in data]
    if 0 < len(missing) < len(columns):
        raise ValueError("Missing SOC columns: " + ", ".join(missing))
    return data.assign(**dict.fromkeys(missing, 0.0)) if missing else data


def x_values(mean, std):
    left = max(np.min(mean - 2 * std), 0.01)
    right = max(np.max(mean + 2 * std), left + 0.01)
    return np.linspace(left, right, max(2, int((right - left) / 0.01)))


def sorting_parameters(*args):
    args = list(args)
    argsort = np.argsort(args[0], axis=1, kind="stable")
    for i in range(len(args)):
        args[i] = np.take_along_axis(args[i], argsort, axis=1)
    return args


def check_number_geoms(data):
    number_geoms = data.shape[0]
    coms = nemo.tools.start_counter()
    if number_geoms != coms:
        print(
            (
                f"There are {coms} inputs and just {number_geoms} log files. "
                "Something is not right! Computing the rates anyway..."
            )
        )
    return number_geoms


def fetch(data, criteria_list):
    regex_list = [re.compile(c) for c in criteria_list]
    columns = [i for i in data.columns.values if all(r.search(i) for r in regex_list)]
    columns.sort(key=lambda column: tuple(int(n) for n in re.findall(r"\d+", column)))
    filtered_data = data[columns].to_numpy()
    return filtered_data


def solvent_susceptibility(corrections, reference_nr):
    """Convert the combined ptSS and ptLR stabilization into a susceptibility."""
    alpha = nemo.tools.get_alpha(reference_nr**2)
    if reference_nr == 1:
        if np.any(corrections):
            raise ValueError("Nonzero PCM corrections require a reference refractive index above 1.")
        return np.zeros_like(corrections)
    return corrections / alpha


def rescale_oscillators(strengths, reference_gap, corrected_gap):
    """Use the reference transition dipoles at the corrected photon energies."""
    return np.divide(
        strengths * np.maximum(corrected_gap, 0.0),
        reference_gap,
        out=np.zeros_like(corrected_gap),
        where=reference_gap > 0,
    )


def total_reorganization_energy(lambda_b, kbt, energies=None, weights=None, fallback=0.28):
    if energies is None:
        return np.sqrt(2 * lambda_b * kbt + kbt**2)

    solvent_width = np.sqrt(np.maximum(2 * lambda_b * kbt, 0.0))
    energies = np.asarray(energies)
    weights = np.ones_like(energies) if weights is None else np.broadcast_to(weights, energies.shape)
    column_shape = energies.shape[1:]
    bandwidth = np.full(column_shape, fallback, dtype=float)
    rows = len(energies)
    for column, (values, amplitudes) in enumerate(
        zip(energies.reshape(rows, -1).T, weights.reshape(rows, -1).T)
    ):
        active = amplitudes > 0
        if np.count_nonzero(active) < 2:
            continue
        values, amplitudes = values[active], amplitudes[active]
        order = np.argsort(values, kind="stable")
        values, amplitudes = values[order], amplitudes[order]
        p = amplitudes / amplitudes.sum()
        sum_squares = np.dot(p, p)
        if sum_squares >= 1:
            continue
        std = np.sqrt(np.dot(p, (values - np.dot(p, values))**2) / (1 - sum_squares))
        positions = np.cumsum(p) - p / 2
        positions = (positions - positions[0]) / (positions[-1] - positions[0])
        q25, q75 = np.interp([0.25, 0.75], positions, values)
        iqr_scale = (q75 - q25) / 1.349
        scale = min(std, iqr_scale) if iqr_scale > 0 else std
        bandwidth.flat[column] = 0.9 * scale * sum_squares**0.2
    width = np.maximum(solvent_width, bandwidth)
    return np.where(width > 0, width, fallback)


def rate_and_uncertainty(y_axis):
    mean_y = np.mean(y_axis, axis=0)
    if len(y_axis) == 1:
        return mean_y, np.zeros_like(mean_y)
    return mean_y, np.std(y_axis, axis=0, ddof=1) / np.sqrt(len(y_axis))


def select_columns(nstate, *args):
    args = list(args)
    modified = []
    for arg in args:
        arg = arg[:, nstate]
        modified.append(arg)
    return modified


def breakdown_emi(chi_s, chi_t, delta_emi, l_total, individual, labels):
    # make a dataframe with Ss_s and Ss_t
    breakdown = pd.DataFrame(
        np.hstack((chi_s, chi_t)),
        columns=[f"chi_s{i+1}" for i in range(chi_s.shape[1])]
        + [f"chi_t{i+1}" for i in range(chi_t.shape[1])],
    )
    # append a columns with energies named eng
    breakdown["eng"] = delta_emi
    breakdown["sigma"] = l_total
    # append individual to df, use labels as columns
    breakdown = pd.concat([breakdown, pd.DataFrame(individual, columns=labels)], axis=1)
    return breakdown


###CALCULATES ISC AND EMISSION RATES & SPECTRA#########################################
def rates(initial, dielec, data=None, ensemble_average=False, detailed=False):
    if data is None:
        data = gather_data(initial, save=True)
        check_number_geoms(data)
    nr_i, kbt = data["nr"].iloc[0], data["kbT"].iloc[0]
    eps, refractive_index = dielec
    alpha_static = nemo.tools.get_alpha(eps)
    alpha_optical = nemo.tools.get_alpha(refractive_index**2)
    alpha_slow = alpha_static - alpha_optical
    initial = initial.lower()
    spin = initial[0]
    n_state = int(initial[1:]) - 1
    data = fix_absent_soc(data, spin)

    singlets = fetch(data, [r"^e_s\d+$"])
    triplets = fetch(data, [r"^e_t\d+$"])
    chi_s = solvent_susceptibility(fetch(data, [r"^d_s\d+$"]), nr_i)
    chi_t = solvent_susceptibility(fetch(data, [r"^d_t\d+$"]), nr_i)
    if spin == "s":
        energies, chi_i, final_energies, chi_f = singlets, chi_s, triplets, chi_t
        target_spin = "t"
    else:
        energies, chi_i, final_energies, chi_f = triplets, chi_t, singlets, chi_s
        target_spin = "s"

    # The spin-free PCM energies already contain the reference static response.
    initial_eq = energies - alpha_static * chi_i
    emission_gaps = initial_eq - alpha_slow * chi_i
    emission_order = np.argsort(initial_eq, axis=1, kind="stable")
    oscs = fetch(data, [rf"^osce_{spin}\d+$"])
    oscs = rescale_oscillators(oscs, energies, emission_gaps)
    selected = emission_order[:, n_state:n_state + 1]
    delta_emi, oscs, lambda_be = [
        np.take_along_axis(values, selected, axis=1)[:, 0]
        for values in (emission_gaps, oscs, alpha_slow * chi_i)
    ]
    constante = (
        refractive_index**2 * E_CHARGE**2
        / (2 * np.pi * HBAR_EV * MASS_E * LIGHT_SPEED**3 * EPSILON_0)
    )
    if spin == "t":
        constante /= 3
    espectro = constante * delta_emi**2 * oscs
    l_total = total_reorganization_energy(lambda_be, kbt, delta_emi, oscs)
    bright = espectro > 0
    tdm = nemo.tools.calc_tdm(oscs[bright], delta_emi[bright], espectro[bright]) if np.any(bright) else 0.0
    x_axis = x_values(delta_emi, l_total)
    y_axis = espectro[:, None] * nemo.tools.gauss(x_axis, delta_emi[:, None], l_total[:, None])
    mean_y, error = rate_and_uncertainty(y_axis)
    emi_rate, emi_error = rate_and_uncertainty(espectro / HBAR_EV)
    gap_emi = means(delta_emi, espectro, ensemble_average)
    mean_sigma_emi = means(l_total, espectro, ensemble_average)
    mean_part_emi = 0.0
    if np.any(bright):
        mean_part_emi = (
            100.0 if ensemble_average
            else 100 * np.sum(espectro)**2 / (len(espectro) * np.dot(espectro, espectro))
        )
    emi = pd.DataFrame({"Energy": x_axis, "Diffrate": mean_y, "Error": error})
    emi.tdm = tdm
    emi.rate = emi_rate
    emi.error = emi_error

    final_neq = final_energies - alpha_optical * chi_f
    socs = fetch(data, [rf"^soc_{spin}\d+_{target_spin}[1-9]\d*$"])
    initial_eq, final_neq, _, chi_f, socs = reorder(initial_eq, final_neq, chi_i, chi_f, socs)
    socs_complete = socs[:, n_state, :]
    delta = final_neq - initial_eq[:, n_state, None]
    lambda_b = alpha_slow * chi_f
    final = [f"{target_spin.upper()}{j+1}" for j in range(final_neq.shape[1])]
    if spin == "t":
        # Report the downward S0 gap as a positive photon energy, as in 1.4.
        socs_s0 = np.take_along_axis(fetch(data, [r"^soc_t\d+_s0$"]), selected, axis=1)
        socs_complete = np.hstack((socs_s0, socs_complete))
        delta = np.hstack((delta_emi[:, None], delta))
        lambda_b = np.hstack((lambda_be[:, None], lambda_b))
        final.insert(0, "S0")

    participation = socs_complete**2
    sigma = total_reorganization_energy(lambda_b, kbt, delta, participation, fallback=kbt or 0.026)
    y_axis = (2 * np.pi / HBAR_EV) * participation * nemo.tools.gauss(delta, 0, sigma)


    # hstack y and espectro
    individual = np.hstack((espectro[:, np.newaxis], y_axis))
    rate, error = rate_and_uncertainty(y_axis)
    total = emi_rate + np.sum(rate)
    results = np.array(
        [
            [
                emi_rate,
                emi_error,
                100 * emi_rate / total if total else 0.0,
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
        mean_part = np.nan_to_num(100 * rate / means(y_axis, y_axis, ensemble_average), posinf=100)
    rate = rate[:, np.newaxis]
    error = error[:, np.newaxis]
    labels = [f"{initial.upper()}->S0"] + [f"{initial.upper()}~>{j}" for j in final]
    labels = np.array(labels)
    results_isc = np.hstack(
        (
            rate,
            error,
            100 * rate / total if total else np.zeros_like(rate),
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
        breakdown = breakdown_emi(
            chi_s, chi_t, delta_emi, l_total, individual, labels
        )
        return results, emi, breakdown
    else:
        return results, emi


#########################################################################################


def save_absorption_spectrum(
    initial, eps, refractive_index, x_axis, mean_y, sigma, labels
):
    arquivo = nemo.tools.naming(f"cross_section_{initial.upper()}_.lx")
    primeira = (
        f"{'Energy(ev)':8s} {'cross_section(A^2)':8s} {'error(A^2)':8s}\n"
        f"Absorption from State: {initial.upper()}\n"
        f"Epsilon: {eps:.3f} nr: {refractive_index:.3f}\n"
    )
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


def make_breakdown(
    initial, spin, num, oscs, deltae_lambda, chis, l_total
):
    # concatenate oscs[:,:,0] with DE[:,:,0] and ds
    colunas = [
        f"{initial.upper()}->{spin.upper()}{i}"
        for i in range(num + 1, num + oscs.shape[1] + 1)
    ]
    colunas += [
        f"eng_{spin}{i}" for i in range(num + 1, num + deltae_lambda.shape[1] + 1)
    ]
    colunas += [f"chi_{spin}{i}" for i in range(num + 1, num + chis.shape[1] + 1)]
    colunas += [f"sigma_{spin}{i}" for i in range(num + 1, num + l_total.shape[1] + 1)]
    # concatenate oscs[:,:,0] with DE[:,:,0] and ds
    breakdown = pd.DataFrame(
        np.hstack(
            (
                oscs[:, :, 0],
                deltae_lambda[:, :, 0],
                chis,
                l_total[:, :, 0],
            )
        ),
        columns=colunas,
    )
    return breakdown


def another_dimension(nstates, *args):
    args = list(args)
    new_args = []
    for arg in args:
        new_args.append(arg[:, :nstates, np.newaxis])
    return new_args


###COMPUTES ABSORPTION SPECTRA###########################################################
def absorption(initial, dielec, data=None, save=False, detailed=False, nstates=-1):
    if data is None:
        data = gather_data(initial, save=True)
    nr_i, kbt = data["nr"].iloc[0], data["kbT"].iloc[0]
    eps, refractive_index = dielec
    alpha_static = nemo.tools.get_alpha(eps)
    alpha_optical = nemo.tools.get_alpha(refractive_index**2)
    initial = initial.lower()
    spin, num = initial[0], int(initial[1:])
    constante = (
        np.pi * E_CHARGE**2 * HBAR_EV
        / (2 * refractive_index * MASS_E * LIGHT_SPEED * EPSILON_0) * 1e20
    )
    energies = fetch(data, [rf"^e_{spin}\d+$"])
    chi = solvent_susceptibility(fetch(data, [rf"^d_{spin}\d+$"]), nr_i)
    final_neq = energies[:, num:] - alpha_optical * chi[:, num:]
    if num:
        base = energies[:, num - 1:num]
        initial_eq = base - alpha_static * chi[:, num - 1:num]
    else:
        base = initial_eq = 0.0
    reference_gap = energies[:, num:] - base
    deltae_lambda = final_neq - initial_eq
    oscs = rescale_oscillators(fetch(data, [rf"^osc_{spin}\d+$"]), reference_gap, deltae_lambda)
    deltae_lambda, oscs, chis = sorting_parameters(deltae_lambda, oscs, chi[:, num:])
    if nstates != -1:
        deltae_lambda, oscs, chis = (values[:, :nstates] for values in (deltae_lambda, oscs, chis))
    lambda_b = (alpha_static - alpha_optical) * chis
    l_total = total_reorganization_energy(lambda_b, kbt, deltae_lambda, oscs)
    x_axis = x_values(deltae_lambda, l_total)
    deltae_lambda, l_total, oscs = another_dimension(
        deltae_lambda.shape[1], deltae_lambda, l_total, oscs
    )
    y_axis = constante * oscs * nemo.tools.gauss(x_axis, deltae_lambda, l_total)
    mean_y, _ = rate_and_uncertainty(y_axis)
    total, sigma = rate_and_uncertainty(np.sum(y_axis, axis=1))
    mean_y = mean_y.T
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
        save_absorption_spectrum(
            initial, eps, refractive_index, x_axis, mean_y, sigma, labels
        )
    if detailed:
        breakdown = make_breakdown(
            initial, spin, num, constante * oscs, deltae_lambda, chis, l_total
        )
        return abs_spec, breakdown
    return abs_spec


#########################################################################################

class Ensemble(object):
    def __init__(self, file, name=''):
        data = pd.read_csv(file)
        initial = data['ensemble'][0]
        self.eps = data['eps'][0]
        self.nr = data['nr'][0]
        self.data = data
        self.initial = initial
        self.name = name

    def rate(self, dielec, ensemble_average=False):
        results, _ = rates(
            self.initial,
            dielec,
            self.data,
            ensemble_average=ensemble_average,
            detailed=False
        )
        return results

    def emission(self, dielec, wavelength=False):
        _, emi = rates(self.initial, dielec, self.data, ensemble_average=False, detailed=False)
        if wavelength:
            emi = self.emi2wavelength(emi)
        return emi

    def complete_emi(self, dielec, ensemble_average=False, wavelength=False):
        results, emi, breakdown = rates(self.initial, dielec, self.data, ensemble_average=ensemble_average, detailed=True)
        if wavelength:
            emi = self.emi2wavelength(emi)
        breakdown.insert(0, 'Geometry', self.data['geometry'].astype(int))
        return results, emi, breakdown

    def complete_abs(self, dielec, nstates=-1, wavelength=False, extinction=False):
        abs_spec, breakdown = absorption(self.initial, dielec, self.data, nstates=nstates, save=False, detailed=True)
        if wavelength:
            abs_spec = self.abs2wavelength(abs_spec)
        if extinction:
            abs_spec = self.abs2extinction(abs_spec)
        breakdown.insert(0, 'Geometry', self.data['geometry'].astype(int))
        return abs_spec, breakdown

    def absorption(self, dielec, nstates=-1, wavelength=False, extinction=False):
        abs_spec = absorption(self.initial, dielec, data=self.data, nstates=nstates, save=False, detailed=False)
        if wavelength:
            abs_spec = self.abs2wavelength(abs_spec)
        if extinction:
            abs_spec = self.abs2extinction(abs_spec)
        return abs_spec

    def breakdown(self, dielec):
        if self.initial == 'S0':
            _, breakdown = absorption(self.initial, dielec, data=self.data, save=False, detailed=True)
        else:
            _, _, breakdown = rates(self.initial, dielec, self.data, ensemble_average=False, detailed=True)
        breakdown.insert(0, 'Geometry', self.data['geometry'].astype(int))
        return breakdown

    def save(self, dielec, mode):
        if mode == 'emi':
            results, emi = rates(self.initial, dielec, self.data, ensemble_average=False, detailed=False)
            export_results(results, emi, dielec)
        elif mode == 'abs':
            _ = absorption(self.initial, dielec, data=self.data, save=True, detailed=False)

    def emi2wavelength(self, emi):
        emi_energy = emi['Energy'].to_numpy()
        emi['Energy'] = 1239.84193 / emi_energy
        emi['Diffrate'] = emi['Diffrate'] * (1239.84193 / (emi_energy ** 2))
        emi['Error'] = emi['Error'] * (1239.84193 / (emi_energy ** 2))
        return emi
    
    def abs2wavelength(self, abs_spec):
        abs_energy = abs_spec['Energy'].to_numpy()
        abs_spec['Energy'] = 1239.84193 / abs_energy
        return abs_spec

    def abs2extinction(self, abs_spec):
        # Convert cross section (A^2) to molar extinction coefficient (M^-1 cm^-1)
        NA = 6.02214076e23  # Avogadro's number
        abs_spec_ext = abs_spec.copy()
        for col in abs_spec_ext.columns:
            if col != 'Energy':
                abs_spec_ext[col] = abs_spec_ext[col] * (1e-16) * NA / (1000 * np.log(10))
        return abs_spec_ext