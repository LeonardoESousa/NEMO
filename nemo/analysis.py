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


def get_osc_phosph(files, singlets, triplets, n_state, ind_s, ind_t, phosph_osc):
    # removed correction from phosph_osc calculation
    eng_singlets = singlets  # - (alphast2/alphaopt1)*Ss_s
    eng_triplets = triplets  # - (alphast2/alphaopt1)*Ss_t
    for j in range(singlets.shape[0]):
        tos = phosph_osc(
            files[j],
            n_state,
            ind_s[j, :],
            ind_t[j, :],
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
            ind_s,
            ind_t,
            ss_s,
            ss_t,
            e_s0,
            ground_pol,
        ) = get_energies("Geometries/" + file)
        singlets = np.array([singlets[:n_state]])
        triplets = np.array([triplets[:n_state]])
        oscs = np.array([oscs[:n_state]])
        ss_s = np.array([ss_s[:n_state]])
        ss_t = np.array([ss_t[:n_state]])
        ind_s = np.array([ind_s[:n_state]])
        ind_t = np.array([ind_t[:n_state]])
        s0 = np.array([e_s0])
        ground_pol = np.array([ground_pol])
        try:
            total_singlets = np.vstack((total_singlets, singlets))
            total_triplets = np.vstack((total_triplets, triplets))
            total_oscs = np.vstack((total_oscs, oscs))
            total_ss_s = np.vstack((total_ss_s, ss_s))
            total_ss_t = np.vstack((total_ss_t, ss_t))
            total_ind_s = np.vstack((total_ind_s, ind_s))
            total_ind_t = np.vstack((total_ind_t, ind_t))
            total_s0 = np.vstack((total_s0, s0))
            total_ground_pol = np.append(total_ground_pol, ground_pol)
        except NameError:
            total_singlets = singlets
            total_triplets = triplets
            total_oscs = oscs
            total_ss_s = ss_s
            total_ss_t = ss_t
            total_ind_s = ind_s
            total_ind_t = ind_t
            total_s0 = s0
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
        total_s0,
        total_ground_pol,
        total_ind_s,
        total_ind_t,
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
        try:
            mean = np.mean(variable, axis=0)
        except (IndexError,ZeroDivisionError):
            mean = np.mean(variable)
    else:
        try:
            mean = np.average(variable, axis=0, weights=weight)
        except ZeroDivisionError:
            column_sums = np.sum(weight, axis=0)
            # Identify columns where the sum is zero
            zero_sum_columns = column_sums == 0
            # Modify columns with zero sum by setting all their elements to 1
            weight[:, zero_sum_columns] = 1
            mean = np.average(variable, axis=0, weights=weight)
        except IndexError:
            mean = np.average(variable, weights=weight)
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
    alphaopt1 = nemo.tools.get_alpha(nr_i**2)
    alphast1 = nemo.tools.get_alpha(eps_i)
    kbt = nemo.tools.detect_sigma()
    total_states, calculation_type = read_cis(files[0])
    (
    numbers,
    singlets,
    triplets,
    oscs,
    ss_s,
    ss_t,
    e_s0,
    ground_pol,
    ind_s,
    ind_t
        ) = analysis(files, total_states, get_energies[calculation_type])
    ss_s = ss_s/alphaopt1
    ss_t = ss_t/alphaopt1
    ground_pol = ground_pol/alphast1
    #start dataframe with numbers as geometry column
    data = pd.DataFrame(numbers, columns=["geometry"])

    for i in range(singlets.shape[1]):
        data[f"e_s{i+1}"] = singlets[:, i]
        formats[f"e_s{i+1}"] = "{:.4f}"
    for i in range(triplets.shape[1]):
        data[f"e_t{i+1}"] = triplets[:, i]
        formats[f"e_t{i+1}"] = "{:.4f}"
    for i in range(ss_s.shape[1]):
        data[f"chi_s{i+1}"] = ss_s[:, i]
        formats[f"chi_s{i+1}"] = "{:.4f}"
    for i in range(ss_t.shape[1]):
        data[f"chi_t{i+1}"] = ss_t[:, i]
        formats[f"chi_t{i+1}"] = "{:.4f}"

    data["e_g"] = e_s0
    formats["e_g"] = "{:.4e}"
    data["chi_s0"] = ground_pol
    formats["chi_s0"] = "{:.4f}"
        
    if "s" in initial.lower():        
        
        if "s0" == initial.lower():
            for i in range(oscs.shape[1]):
                data[f"osc_s{i+1}"] = oscs[:, i]
                formats[f"osc_s{i+1}"] = "{:.5e}"

        else:
            for i in range(oscs.shape[1]):
                data[f"osce_s{n_state+1+i}"] = oscs[:, i]
                formats[f"osce_s{n_state+1+i}"] = "{:.5e}"

            noscs = get_oscs[calculation_type](files, ind_s, initial)
            for i in range(noscs.shape[1]):
                data[f"osc_s{n_state+2+i}"] = noscs[:, i]
                formats[f"osc_s{n_state+2+i}"] = "{:.5e}"

        try:
            for i in range(singlets.shape[1]):
                socs_partial = get_avg_socs[calculation_type](files, "singlet", i, ind_s, ind_t)
                for j in range(singlets.shape[1]):
                    data[f"soc_s{i+1}_t{j+1}"] = socs_partial[:, j]
                    formats[f"soc_s{i+1}_t{j+1}"] = "{:.5e}"
                    
        except IndexError:
            pass
    else:
        
        oscs = get_osc_phosph(files, singlets, triplets, total_states, ind_s, ind_t, get_phosph_osc[calculation_type])
        
        for i in range(oscs.shape[1]):
            data[f"osce_t{n_state+1+i}"] = oscs[:, i]
            formats[f"osce_t{n_state+1+i}"] = "{:.5e}"
        
        noscs =  get_oscs[calculation_type](files, ind_t, initial)
        
        for i in range(noscs.shape[1]):
            data[f"osc_t{n_state+2+i}"] = noscs[:, i]
            formats[f"osc_t{n_state+2+i}"] = "{:.5e}"
        
        try:
            
            for i in range(triplets.shape[1]):
                
                soc_ground = get_avg_socs[calculation_type](files, "ground", i, ind_s, ind_t)
                soc_triplet = get_avg_socs[calculation_type](files, "triplet", i, ind_s, ind_t)
                soc_tts = get_avg_socs[calculation_type](files, "tts", i, ind_s, ind_t)

                for j in range(triplets.shape[1]):
                    data[f"soc_t{i+1}_s0"] = soc_ground[:, j]
                    formats[f"soc_t{i+1}_s0"] = "{:.5e}"
                    data[f"soc_t{i+1}_s{j+1}"] = soc_triplet[:, j]
                    formats[f"soc_t{i+1}_s{j+1}"] = "{:.5e}"
                    data[f"soc_t{i+1}_t{j+1}"] = soc_tts[:, j]
                    formats[f"soc_t{i+1}_t{j+1}"] = "{:.5e}"
                
        except IndexError:
            pass
    
    
    arquivo = f"Ensemble_{initial.upper()}_.lx"
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


def x_values(mean, std):
    left = max(np.min(mean - 2 * std), 0.01)
    right = np.max(mean + 2 * std)
    x_axis = np.linspace(left, right, int((right - left) / 0.01))
    return x_axis


def sorting_parameters(*args):
    args = list(args)
    argsort = np.argsort(args[0], axis=1)
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
    filtered_data = data[
        [i for i in data.columns.values if all(r.search(i) for r in regex_list)]
    ].to_numpy()
    return filtered_data


def total_reorganization_energy(lambda_b, kbt):
    return np.sqrt(2 * lambda_b * kbt + kbt**2)


def rate_and_uncertainty(y_axis):
    number_geoms = y_axis.shape[0]
    mean_y = np.sum(y_axis, axis=0) / number_geoms
    error = np.sqrt(
        np.sum((y_axis - mean_y) ** 2, axis=0) / (number_geoms * (number_geoms - 1))
    )
    return mean_y, error


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
        kbt = nemo.tools.detect_sigma()
    else:
        kbt = data["kbT"][0]
    eps, refractive_index = dielec[0], dielec[1]
    alphast2 = nemo.tools.get_alpha(eps)
    alphaopt2 = nemo.tools.get_alpha(refractive_index**2)
    
    #excited state energies
    singlets = fetch(data, ["^e_s"])
    triplets = fetch(data, ["^e_t"])
    ground = fetch(data, ["^e_g"])

    #excited state susceptibilities
    chi_s = fetch(data, ["^chi_s(?!0)"])
    chi_t = fetch(data, ["^chi_t"])
    
    #ground state susceptibility
    chi_s0 = data['chi_s0'].to_numpy()
    #fix dimension of chi_s0
    chi_s0 = chi_s0[:, np.newaxis]

    n_state = int(initial[1:]) - 1
    initial = initial.lower()

    data = fix_absent_soc(data)

    # Emission Calculations

    energies = fetch(data, [f"^e_{initial[0]}"])
    
    
    delta_emi_unsorted = energies - chi_s * alphast2 - (ground - chi_s0 * alphaopt2) 
    constante = (
        (refractive_index**2)
        * (E_CHARGE**2)
        / (2 * np.pi * HBAR_EV * MASS_E * (LIGHT_SPEED**3) * EPSILON_0)
    )
    if "t" in initial:
        constante *= 1 / 3
        #reorganization energy
        lambda_be = (alphast2 - alphaopt2) * chi_s0
    else:
        #reorganization energy
        lambda_be = (alphast2 - alphaopt2) * chi_s0
    #make dimensions match
    lambda_be = np.repeat(lambda_be, energies.shape[1], axis=1)

    oscs = fetch(data, ["^osce_"])
    print(delta_emi_unsorted.shape, oscs.shape, lambda_be.shape)
    delta_emi, oscs, lambda_be = sorting_parameters(delta_emi_unsorted, oscs, lambda_be)
    delta_emi, oscs, lambda_be = select_columns(n_state, delta_emi, oscs, lambda_be)
    
    l_total = total_reorganization_energy(lambda_be, kbt)
    
    espectro = constante * ((delta_emi - lambda_be) ** 2) * oscs
    tdm = nemo.tools.calc_tdm(oscs, delta_emi, espectro)
    x_axis = x_values(delta_emi, l_total)
    y_axis = espectro[:, np.newaxis] * nemo.tools.gauss(
        x_axis, delta_emi[:, np.newaxis], l_total[:, np.newaxis]
    )
    number_geoms = y_axis.shape[0]
    mean_y, error = rate_and_uncertainty(y_axis)
    emi_rate = np.mean(espectro, axis=0) / HBAR_EV
    emi_error = np.sqrt(np.sum((espectro /HBAR_EV - emi_rate) ** 2, axis=0) / (number_geoms * (number_geoms - 1)))
    gap_emi = means(delta_emi, espectro, ensemble_average)
    mean_sigma_emi = means(l_total, espectro, ensemble_average)
    mean_part_emi = (100 / number_geoms) / means(
        espectro / np.sum(espectro), espectro, ensemble_average
    )
    emi = np.hstack(
        (x_axis[:, np.newaxis], mean_y[:, np.newaxis], error[:, np.newaxis])
    )
    emi = pd.DataFrame(emi, columns=["Energy", "Diffrate", "Error"])
    emi.tdm = tdm
    emi.rate = emi_rate
    emi.error = emi_error

    # Checks number of logs
    if data is None:
        check_number_geoms(data)
    # Intersystem Crossing Rates

    if "s" in initial:
        initial_state = singlets - chi_s * alphast2
        final_state = triplets - chi_t * alphaopt2
        socs_complete = fetch(data, ["^soc_s"])
        initial_state, final_state, chi_s, chi_t, socs_complete = reorder(
            initial_state, final_state, chi_s, chi_t, socs_complete
        )
        initial_state = initial_state[:, n_state]
        socs_complete = socs_complete[:, n_state, :]
        delta = final_state - initial_state[:, np.newaxis]
        lambda_b = (alphast2 - alphaopt2) * chi_t
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
        initial_state = triplets - chi_t * alphast2
        final_state = singlets - chi_s * alphaopt2
        socs_complete = fetch(data, ["^soc_t.*s[1-9]"])
        initial_state, final_state, chi_t, chi_s, socs_complete = reorder(
            initial_state, final_state, chi_t, chi_s, socs_complete
        )
        initial_state = initial_state[:, n_state]
        socs_complete = socs_complete[:, n_state, :]
        delta = final_state - initial_state[:, np.newaxis]
        lambda_b = (alphast2 - alphaopt2) * chi_s
        final = [
            i.split("_")[2].upper()
            for i in data.columns.values
            if "soc_" + initial.lower() + "_" in i and i.count("t") == 1
        ]
        # Tn to S0 ISC
        socs_s0 = fetch(data, ["^soc_t.*s0"])
        delta_emi, socs_s0 = sorting_parameters(delta_emi_unsorted, socs_s0)
        delta_emi = delta_emi[:, n_state]
        socs_s0 = socs_s0[:, n_state]
        socs_complete = np.hstack((socs_s0[:, np.newaxis], socs_complete))
        delta = np.hstack((delta_emi[:, np.newaxis], delta))
        lambda_b = np.hstack((lambda_be[:, np.newaxis], lambda_b))
        # Tn to Tm ISC
        # delta_tt = Triplets + np.repeat((alphast2/alphaopt1)*Ss_t[:,n_state][:,np.newaxis] - Triplets[:,n_state][:,np.newaxis],Triplets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_t    #Tm (final) - Tn (initial) + lambda_b
        # Removing Tn to Tn transfers
        # indices  = [i for i in range(Triplets.shape[1]) if i != n_state]
        # delta    = np.hstack((delta,delta_tt[:,indices]))
        # lambda_bt= (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_t
        # lambda_b = np.hstack((lambda_b,lambda_bt[:,indices]))
        # final.extend([i.upper()[4:] for i in data.columns.values if 'soc_t' in i])

    sigma = total_reorganization_energy(lambda_b, kbt)
    y_axis = (
        (2 * np.pi / HBAR_EV) * (socs_complete**2) * nemo.tools.gauss(delta, 0, sigma)
    )
    # hstack y and espectro
    individual = np.hstack((espectro[:, np.newaxis], y_axis))
    number_geoms = y_axis.shape[0]
    rate, error = rate_and_uncertainty(y_axis)
    total = emi_rate + np.sum(rate)
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
        mean_part = np.nan_to_num(100 * rate / means(y_axis, y_axis, ensemble_average), posinf=100)
    rate = rate[:, np.newaxis]
    error = error[:, np.newaxis]
    labels = [f"{initial.upper()}->S0"] + [f"{initial.upper()}~>{j}" for j in final]
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
                chis ,
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
        kbt = nemo.tools.detect_sigma()
    else:
        kbt = data["kbT"][0]
    eps, refractive_index = dielec[0], dielec[1]
    alphast2 = nemo.tools.get_alpha(eps)
    alphaopt2 = nemo.tools.get_alpha(refractive_index**2)
    
    initial = initial.lower()
    spin = initial[0]
    num = int(initial[1:])

    # excited state energies
    engs = fetch(data, [f"^e_{spin}"])

    #ground state susceptibility
    chi_s0 = data['chi_s0'].to_numpy()
    chi_s0 = chi_s0[:, np.newaxis]
    ground = fetch(data, ["^e_g"])
    
    #excited state susceptibilities
    chis = fetch(data, [f"^chi_{spin}(?!0)"])

    #oscillator strengths
    oscs = fetch(data, ["^osc_"])

    constante = (
        (np.pi * (E_CHARGE**2) * HBAR_EV)
        / (2 * refractive_index * MASS_E * LIGHT_SPEED * EPSILON_0)
        * 1e20
    )
    
    
    engs = engs[:, num:]
    chis = chis[:, num:]
    
    lambda_b = (alphast2  - alphaopt2) * chis

    if initial == "s0":
        deltae_lambda = engs - chis * alphaopt2 - (ground - chi_s0 * alphast2)

    else:
        base = fetch(data, [rf"\be_{initial}\b"])
        chi_i = chis[:,0]
        deltae_lambda = (engs - chis * alphaopt2) - (base - chi_i[:,np.newaxis] * alphast2 )
        

    # Sorting states by energy
    deltae_lambda, oscs, lambda_b = sorting_parameters(
        deltae_lambda, oscs, lambda_b
    )
    l_total = total_reorganization_energy(lambda_b, kbt)
    x_axis = x_values(deltae_lambda, l_total)
    if nstates == -1:
        nstates = deltae_lambda.shape[1]
    # Add extra dimension to DE and Ltotal to match x shape
    deltae_lambda, l_total, oscs, lambda_b = another_dimension(
        nstates, deltae_lambda, l_total, oscs, lambda_b
    )
    y_axis = constante * oscs * nemo.tools.gauss(x_axis, deltae_lambda, l_total)
    mean_y, sigma = rate_and_uncertainty(y_axis)
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
        save_absorption_spectrum(
            initial, eps, refractive_index, x_axis, mean_y, sigma, labels
        )
    if detailed:
        breakdown = make_breakdown(
            initial, spin, num, constante*oscs, deltae_lambda, chis[:,:nstates], l_total
        )
        return abs_spec, breakdown
    return abs_spec


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

    def emission(self, dielec):
        _, emi = rates(self.initial, dielec, self.data, ensemble_average=False, detailed=False)
        return emi

    def complete_emi(self, dielec, ensemble_average=False):
        results, emi, breakdown = rates(self.initial, dielec, self.data, ensemble_average=ensemble_average, detailed=True)
        breakdown.insert(0, 'Geometry', self.data['geometry'].astype(int))
        return results, emi, breakdown

    def complete_abs(self, dielec, nstates=-1):
        abs_spec, breakdown = absorption(self.initial, dielec, self.data, nstates=nstates, save=False, detailed=True)
        breakdown.insert(0, 'Geometry', self.data['geometry'].astype(int))
        return abs_spec, breakdown

    def absorption(self, dielec, nstates=-1):
        abs_spec = absorption(self.initial, dielec, data=self.data, nstates=nstates, save=False, detailed=False)
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
