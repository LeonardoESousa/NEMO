#!/usr/bin/env python3
import os
import re
import warnings
import numpy as np
import pandas as pd
import nemo.tools
import nemo.parser
import nemo.eom
from itertools import combinations

# pylint: disable=unbalanced-tuple-unpacking

LIGHT_SPEED = nemo.parser.LIGHT_SPEED
HBAR_EV = nemo.parser.HBAR_EV
HBAR_J = nemo.parser.HBAR_J
BOLTZ_EV = nemo.parser.BOLTZ_EV
E_CHARGE = nemo.parser.E_CHARGE
MASS_E = nemo.parser.MASS_E
EPSILON_0 = nemo.parser.EPSILON_0


##RETURNS LIST OF LOG FILES WITH NORMAL TERMINATION######################################
def check_normal(folder):
    watcher = nemo.tools.Watcher(folder)
    watcher.check()
    return [i+'.log' for i in watcher.done]

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
            ground_pol,
            y_s,
            y_t,
        ) = get_energies("Geometries/" + file)
        singlets = np.array([singlets[:n_state]])
        triplets = np.array([triplets[:n_state]])
        oscs = np.array([oscs[:n_state]])
        ss_s = np.array([ss_s[:n_state]])
        ss_t = np.array([ss_t[:n_state]])
        ind_s = np.array([ind_s[:n_state]])
        ind_t = np.array([ind_t[:n_state]])
        y_s = np.array([y_s[:n_state]])
        y_t = np.array([y_t[:n_state]])
        ground_pol = np.array([ground_pol])
        try:
            total_singlets = np.vstack((total_singlets, singlets))
            total_triplets = np.vstack((total_triplets, triplets))
            total_oscs = np.vstack((total_oscs, oscs))
            total_ss_s = np.vstack((total_ss_s, ss_s))
            total_ss_t = np.vstack((total_ss_t, ss_t))
            total_ind_s = np.vstack((total_ind_s, ind_s))
            total_ind_t = np.vstack((total_ind_t, ind_t))
            total_y_s = np.vstack((total_y_s, y_s))
            total_y_t = np.vstack((total_y_t, y_t))
            total_ground_pol = np.append(total_ground_pol, ground_pol)
        except NameError:
            total_singlets = singlets
            total_triplets = triplets
            total_oscs = oscs
            total_ss_s = ss_s
            total_ss_t = ss_t
            total_ind_s = ind_s
            total_ind_t = ind_t
            total_y_s = y_s
            total_y_t = y_t
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
        total_ind_s,
        total_ind_t,
        total_y_s,
        total_y_t,
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
    files = check_normal("Geometries")
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
    ground_pol,
    ind_s,
    ind_t,
    y_s,
    y_t,
        ) = analysis(files, total_states, get_energies[calculation_type])
    ss_s = ss_s/alphaopt1
    ss_t = ss_t/alphaopt1
    gamma_s0 = ground_pol/alphast1
    gamma_s = y_s/alphast1
    gamma_t = y_t/alphast1

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
    for i in range(gamma_s.shape[1]):
        data[f"gamma_s{i+1}"] = gamma_s[:, i]
        formats[f"gamma_s{i+1}"] = "{:.4f}"
    for i in range(gamma_t.shape[1]):
        data[f"gamma_t{i+1}"] = gamma_t[:, i]
        formats[f"gamma_t{i+1}"] = "{:.4f}"

    data["gamma_s0"] = gamma_s0
    formats["gamma_s0"] = "{:.4f}"

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
                #soc_tts = get_avg_socs[calculation_type](files, "tts", i, ind_s, ind_t)

                data[f"soc_t{i+1}_s0"] = soc_ground[:, 0]
                formats[f"soc_t{i+1}_s0"] = "{:.5e}"
                for j in range(triplets.shape[1]):
                    data[f"soc_t{i+1}_s{j+1}"] = soc_triplet[:, j]
                    formats[f"soc_t{i+1}_s{j+1}"] = "{:.5e}"
                    #data[f"soc_t{i+1}_t{j+1}"] = soc_tts[:, j]
                    #formats[f"soc_t{i+1}_t{j+1}"] = "{:.5e}"

        except IndexError:
            pass

    # Check if derivative coupling calculation was performed
    dc_computation, states = nemo.parser.check_derivative_couplings(files[0])
    if dc_computation:
        ###### OBTAINS THE B PARAMETERS ########
        formats_dc = {}
        freq_log = nemo.tools.fetch_file("frequency", [".out", ".log"])
        (
            initial_state,
            final_state,
            geometry,
            mode,
            b
        ) = nemo.parser.get_derivative_couplings(states, states, files, freq_log, debug=False)

        arquivo_dc =f"Derivative_Couplings_{initial.upper()}_.lx"
        data_dc = pd.DataFrame()
        data_dc["initial_state"] = np.array(initial_state).astype(str)
        data_dc["final_state"] = np.array(final_state).astype(str)
        data_dc["geometry"] = geometry
        data_dc["mode"] = mode
        data_dc["B"] = b

        formats_dc["initial_state"] = "{:s}"
        formats_dc["final_state"] = "{:s}"
        formats_dc["geometry"] = "{:.0f}"
        formats_dc["mode"] = "{:.0f}"
        formats_dc["B"] = "{:.5e}"

        ###### OBTAINS THE V PARAMETERS ########
        Mag_file = nemo.tools.fetch_file("Magnitudes", ['Magnitudes'])
        (
            geometry_V,
            mode_V,
            v1,
            v2
        ) = nemo.parser.get_V(Mag_file, files)

        formats_V = {}
        arquivo_V =f"V_{initial.upper()}_.lx"
        data_V = pd.DataFrame()
        data_V["geometry"] = geometry_V
        data_V["mode"] = mode_V
        data_V["v1"] = v1
        data_V["v2"] = v2

        formats_V["geometry"] = "{:.0f}"
        formats_V["mode"] = "{:.0f}"
        formats_V["v1"] = "{:.5e}"
        formats_V["v2"] = "{:.5e}"
        if save:
            # Create a temporary copy of the DataFrame
            temp_data_V = data_V.copy()
            #Apply formats
            for column, fmt in formats_V.items():
                if column in temp_data_V.columns:
                    temp_data_V[column] = temp_data_V[column].map(fmt.format)
            temp_data_V.to_csv(arquivo_V, index=False)
            # Create a temporary copy of the DataFrame
            temp_data_dc = data_dc.copy()
            #Apply formats
            for column, fmt in formats_dc.items():
                if column in temp_data_dc.columns:
                    temp_data_dc[column] = temp_data_dc[column].map(fmt.format)
            temp_data_dc.to_csv(arquivo_dc, index=False)

        for i in states:
            for j in states:
                if i== j:
                    continue
                _, _, h = gather_data_derivative_couplings("s"+str(i), "s"+str(j), data, data_dc, data_V)
                data[f"IC_{i}_{j}"] = h
                formats[f"IC_{i}_{j}"] = "{:.5e}"
        #header.extend([f"IC_{i}_{j}(eV)" for i, j in combinations(states, 2)])
        #any({formats.update({f"IC_{i}_{j}(eV)": "{:.5e}"}) for i, j in combinations(states, 2)})
        #data=np.hstack((data, h))

    arquivo = f"Ensemble_{initial.upper()}_.lx"
    data["ensemble"] = initial.upper()
    formats["ensemble"] = "{:s}"
    data["kbT"] = kbt
    formats["kbT"] = "{:.4f}"
    # make these the first columns
    cols = data.columns.tolist()
    cols = cols[-2:] + cols[:-2]
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


def gather_data_derivative_couplings(initial, final, data=None, data_dc=None, data_V=None):

    if int(initial[1])==0:
        energies = fetch(data, [f"^e_{initial.lower()[0]}"])
        delta= energies
        delta = np.hstack((np.zeros_like(energies[:, -1])[:, np.newaxis], delta)) # add a second column with the same values for the ground state
    else:
        n_state = int(initial[1:]) - 1
        energies = fetch(data, [f"^e_{initial.lower()[0]}"])
        delta = energies - energies[:, n_state ][:, np.newaxis]
        #remove column corresponding to the initial state
        #In gather_data, the states are never the same
        #delta = np.delete(delta, n_state, axis=1)
        #add -energies[:, n_state][:, np.newaxis] as the first column
        delta = np.hstack((-energies[:, n_state][:, np.newaxis], delta))
    #------------------------------------#
    # Computes the H parameters
    h = 0.0
    if data is not None:
        v1, v2 = nemo.tools.V_to_vec(data_V)

        # ----- Freq
        mag_file = nemo.tools.fetch_file("Magnitudes", ['Magnitudes'])
        data_f = pd.read_csv(mag_file)
        freq_V = data_f.filter(regex="freq").dropna().to_numpy().flatten()
        freq_row = freq_V[np.newaxis,:] #rad/s


        if int(initial[1]) < int(final[1]):
            lower=int(initial[1])
            higher=int(final[1])
        else:
            higher=int(initial[1])
            lower=int(final[1])


        # Compute the geometry dependent rate for the ith transition
        # ----- Get the transition energy for the ith transition
        #e_col = fetch(data, [f"^e_{initial.lower()[0]}"])[:,i] # eV
        #e_col = e_col[:,np.newaxis]
        #e_col *= -1.0
        e_col = delta[:, int(final[1])][:, np.newaxis] # eV
        # ----- Get oscillator strength for the ith transition
        #osc_col = fetch(data, [f"^osce_{initial.lower()[0]}"])[:,i] #
        #osc_col = osc_col[:,np.newaxis]
        #constante = E_CHARGE**2 / (2.0 * np.pi * HBAR_EV * MASS_E * (LIGHT_SPEED**3.0) * EPSILON_0)
        #espectro = constante * ((e_col) ** 2 ) * osc_col
        #gammas_lorentz = espectro / 2.0# nemo.tools.detect_sigma()# espectro / 2.0
        sigma = np.std(e_col, axis=0)



        # ----_ Get the corresponding coupling for the ith transition
        b = nemo.tools.B_to_vec(data_dc, lower, higher) # J^2
        #np.savetxt(f"debug_dc.csv", b, delimiter=",", fmt='%10.5e')

        # ----- Calculate h_IC
        # ----- Gaussian line shape
        # argument_pos =  (e_col * HBAR_EV * freq_row)/(s**2)-(HBAR_EV * freq_row)**2/(2.0*s**2)
        # argument_neg = -(e_col * HBAR_EV * freq_row)/(s**2)-(HBAR_EV * freq_row)**2/(2.0*s**2)
        # term_pos= np.exp((argument_pos))
        # term_neg= np.exp((argument_neg))

        # debugging: print the arguments before exponentiation to check for overflow issues
        # print the numpy array argument to a file
        # np.savetxt("debug_argument_pos.csv", argument_pos, delimiter=",", fmt='%10.5f')
        # np.savetxt("debug_argument_neg.csv", argument_neg, delimiter=",", fmt='%10.5f')
        # print(argument_pos)
        # np.savetxt("debug_exp_pos.csv", term_pos, delimiter=",", fmt='%10.5f')
        # np.savetxt("debug_exp_neg.csv", term_neg, delimiter=",", fmt='%10.5f')
        #print(term_pos)

        # ------ Lorentzian line shape
        #term_pos = nemo.tools.voigt(-e_col+HBAR_EV*freq_row,0.0,gammas_lorentz) / nemo.tools.voigt(-e_col,0.0,gammas_lorentz)
        #term_neg = nemo.tools.voigt(-e_col-HBAR_EV*freq_row,0.0,gammas_lorentz) / nemo.tools.voigt(-e_col,0.0,gammas_lorentz)
        term_pos = nemo.tools.gauss(e_col-HBAR_EV*freq_row,sigma) 
        term_neg = nemo.tools.gauss(e_col+HBAR_EV*freq_row,sigma) 
        #print(final)
        # np.savetxt(f"debug_exp_pos_{final}.csv", term_pos, delimiter=",", fmt='%10.5e')
        # np.savetxt(f"debug_exp_neg_{final}.csv", term_neg, delimiter=",", fmt='%10.5e')
        # np.savetxt(f"debug_exp_den_{final}.csv", nemo.tools.gauss(e_col,sigma), delimiter=",", fmt='%10.5e')
        # np.savetxt(f"debug_exp_ratio_pos_{final}.csv", term_pos/nemo.tools.gauss(e_col,sigma), delimiter=",", fmt='%10.5e')
        # np.savetxt(f"debug_exp_ratio_neg_{final}.csv", term_neg/nemo.tools.gauss(e_col,sigma), delimiter=",", fmt='%10.5e')

        h = b * (v1 * term_pos + v2 * term_neg) / (nemo.tools.gauss(e_col,sigma)) # J^2
        # ----- sum on normal modes
        h=np.sum(h, axis=1)[:,np.newaxis]

        # ----- Add H for this transition pair
        try:
            total_H = np.hstack((total_H, h))
        except NameError:
            total_H = h
    #np.savetxt("debug.csv", total_H/ (E_CHARGE**2), delimiter=",", fmt='%10.5e') #
    return data_dc, data_V, total_H / (E_CHARGE**2) # J**2 to eV**2
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


def reorder(initial_state, final_state, ss_i, ss_f, gamma_i, gamma_f, socs):
    argsort = np.argsort(initial_state, axis=1)
    initial_state = np.take_along_axis(initial_state, argsort, axis=1)
    ss_i = np.take_along_axis(ss_i, argsort, axis=1)
    gamma_i = np.take_along_axis(gamma_i, argsort, axis=1)
    corredor = int(np.sqrt(socs.shape[1]))
    socs_complete = socs.reshape((socs.shape[0], corredor, corredor))
    for j in range(socs_complete.shape[1]):
        socs_complete[:, j, :] = np.take_along_axis(
            socs_complete[:, j, :], argsort, axis=1
        )
    argsort = np.argsort(final_state, axis=1)
    final_state = np.take_along_axis(final_state, argsort, axis=1)
    ss_f = np.take_along_axis(ss_f, argsort, axis=1)
    gamma_f = np.take_along_axis(gamma_f, argsort, axis=1)
    for j in range(socs_complete.shape[1]):
        socs_complete[:, :, j] = np.take_along_axis(
            socs_complete[:, :, j], argsort, axis=1
        )
    return initial_state, final_state, ss_i, ss_f, gamma_i, gamma_f, socs_complete


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

def fix_absent_triplets(data):
    columns = data.columns.values
    # check if at least one column contains t_
    if any("e_t" in i for i in columns):
        return data
    else:
        singlets = [i.split("_")[1] for i in columns if "e_s" in i and "osc" not in i]
        for singlet in singlets:
            data[f"e_t{singlet[1:]}"] = 0
            data[f"chi_t{singlet[1:]}"] = 0
            data[f"gamma_t{singlet[1:]}"] = 0
    return data


def x_values(mean, std):
    left = max(np.min(mean - 3.5 * std), 0.01)
    right = np.max(mean + 3.5 * std)
    x_axis = np.linspace(left, right, int((right - left) / 0.01))
    return x_axis


def sorting_parameters(*args, axis=1):
    args = list(args)

    argsort = np.argsort(args[0], axis=axis)

    for i in range(len(args)):
        if args[i].ndim == 2:
            args[i] = np.take_along_axis(args[i], argsort, axis=axis)

        elif args[i].ndim == 3:
            # expand argsort to match the 3rd dimension
            indexer = np.expand_dims(argsort, axis=2)
            args[i] = np.take_along_axis(args[i], indexer, axis=axis)

        else:
            raise ValueError("Unsupported number of dimensions")

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


def total_reorganization_energy(lambda_b, kbt, std):
    return np.sqrt(2 * lambda_b * kbt + std**2)


def rate_and_uncertainty(y_axis):
    with np.errstate(invalid='ignore', divide='ignore'):
        number_geoms = y_axis.shape[0]
        mean_y = np.sum(y_axis, axis=0) / number_geoms
        error = np.sqrt(
            np.sum((y_axis - mean_y) ** 2, axis=0) / (number_geoms * (number_geoms - 1))
        )
    return mean_y, np.nan_to_num(error, nan=0.0)


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


def lambda_solvent(chi_i, chi_f, alphaopt, alphast):
    chi_t = chi_i + chi_f - 2 * np.sqrt(chi_i * chi_f)
    lambda_b = chi_t * (alphast - alphaopt)  #chi_f * (alphast - alphaopt)
    return lambda_b


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

    #data = fix_absent_triplets(data)
    #data = fix_absent_soc(data)

    n_state = int(initial[1:]) - 1
    initial = initial.lower()

    #excited state energies
    singlets = fetch(data, ["^e_s"])
    triplets = fetch(data, ["^e_t"])
    
    #excited state susceptibilities
    chi_s = fetch(data, ["^chi_s(?!0)"])
    chi_t = fetch(data, ["^chi_t"])
    gamma_s = fetch(data, ["^gamma_s(?!0)"])
    gamma_t = fetch(data, ["^gamma_t"])

    # oscillator strengths
    oscs = fetch(data, ["^osce_"])
    
    #spin orbit couplings
    if "s" in initial:
        socs_complete = fetch(data, ["^soc_"+initial[0]])
    else:
        socs_complete = fetch(data, ["^soc_t.*_s[^0].*"])

    #spin orbit couplings for Tn -> S0
    socs_s0 = fetch(data, ["^soc_t.*_s0"])

    #internal conversion
    h_ic = fetch(data, ["^IC_(?!0)"])    

    #ground state susceptibility
    gamma_s0 = data['gamma_s0'].to_numpy()
    #fix dimension of gamma_s0
    gamma_s0 = gamma_s0[:, np.newaxis]

    #check empty arrays and fill with zeros if necessary
    if triplets.size == 0:
        triplets = np.zeros_like(singlets)
        chi_t = np.zeros_like(chi_s)
        gamma_t = np.zeros_like(gamma_s)
    if socs_complete.size == 0:
        socs_complete = np.zeros((singlets.shape[0], singlets.shape[1]**2))
    if h_ic.size == 0:
        h_ic = np.zeros((singlets.shape[0], singlets.shape[1]**2))   
    if socs_s0.size == 0:
        socs_s0 = np.zeros_like(triplets)

    #reshape socs_complete to have dimensions (number of geometries, number of initial states, number of final states)
    socs_complete = socs_complete.reshape((socs_complete.shape[0], singlets.shape[1], triplets.shape[1]))
    h_ic = h_ic.reshape((h_ic.shape[0], singlets.shape[1], singlets.shape[1]))

    # Emission Calculations

    energies = singlets if "s" in initial else triplets
    chis = chi_s if "s" in initial else chi_t
    gammas = gamma_s if "s" in initial else gamma_t

    initial_energy = energies - (gammas + chis) * alphast2
    initial_energy, energies, oscs, chis, gammas, socs_complete, h_ic, socs_s0 = sorting_parameters(initial_energy, energies, oscs, chis, gammas, socs_complete, h_ic, socs_s0)
    
    lambda_emission = lambda_solvent(chis, 0, alphaopt2, alphast2)
    final_energy_emission = - gamma_s0 * alphast2
    delta_emi = final_energy_emission - initial_energy + lambda_emission    
    delta_emi *= -1

    constante = (
        (refractive_index**2)
        * (E_CHARGE**2)
        / (2 * np.pi * HBAR_EV * MASS_E * (LIGHT_SPEED**3) * EPSILON_0)
    )
    if "t" in initial:
        constante *= 1 / 3

    delta_emi, initial_energy, oscs, lambda_emission, chis = select_columns(n_state, delta_emi, initial_energy, oscs, lambda_emission, chis)
    chis = chis[:,np.newaxis]
    sigma_int = np.std(delta_emi, axis=0)
    l_total = total_reorganization_energy(lambda_emission, kbt, sigma_int)
    espectro = constante * (delta_emi ** 2) * oscs
    
    tdm = nemo.tools.calc_tdm(oscs, delta_emi, espectro)
    x_axis = x_values(delta_emi, l_total)
    y_axis = espectro[:, np.newaxis] * nemo.tools.gauss(
        x_axis - delta_emi[:, np.newaxis], l_total[:, np.newaxis]
    )
    number_geoms = y_axis.shape[0]
    mean_y, error = rate_and_uncertainty(y_axis)
    emi_rate = np.mean(espectro, axis=0) / HBAR_EV
    with np.errstate(invalid='ignore', divide='ignore'):
        emi_error = np.sqrt(np.sum((espectro /HBAR_EV - emi_rate) ** 2, axis=0) / (number_geoms * (number_geoms - 1)))
    emi_error = np.nan_to_num(emi_error, nan=0.0)
    gap_emi = means(delta_emi, espectro, ensemble_average)
    mean_emi_coupling = 1000 * means(np.sqrt(espectro / 2 * np.pi), espectro, ensemble_average)
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

    #select columns corresponding to the initial state
    socs_complete = socs_complete[:, n_state, :]
    h_ic = h_ic[:, n_state, :]
    
    if "s" in initial:
        # Sm to Tn ISC
        final_energy_isc = triplets - (gamma_t + chi_t) * alphast2
        final_energy_isc, triplets, chi_t, gamma_t, socs_complete = sorting_parameters(final_energy_isc, triplets, chi_t, gamma_t, socs_complete)
        lambda_b_isc = lambda_solvent(chis, chi_t, alphaopt2, alphast2)
        delta_isc = final_energy_isc - initial_energy[:, np.newaxis] + lambda_b_isc
        sigma_int = np.std(delta_isc, axis=0)
        final = [
            i.split("_")[2].upper()
            for i in data.columns.values
            if "soc_" + initial.lower() + "_" in i
        ]
        # Sm to Sn IC
        final_energy_ic = singlets - (gamma_s + chi_s) * alphast2
        final_energy_ic, singlets, chi_s, gamma_s, h_ic = sorting_parameters(final_energy_ic, singlets, chi_s, gamma_s, h_ic)
        lambda_b_ic = lambda_solvent(chis, chi_s, alphaopt2, alphast2)
        delta_ic = final_energy_ic - initial_energy[:, np.newaxis] + lambda_b_ic
        #remove column corresponding to the initial state
        delta_ic = np.delete(delta_ic, n_state, axis=1)
        lambda_b_ic = np.delete(lambda_b_ic, n_state, axis=1)
        # add emission energy as the first column
        delta_ic = np.hstack((-1 * delta_emi[:, np.newaxis], delta_ic))
        lambda_b_ic = np.hstack((lambda_emission[:, np.newaxis], lambda_b_ic))
        sigma_int_ic = np.std(delta_ic, axis=0)
        #remove columns for states higher than n_state
        delta_ic = delta_ic[:, :n_state+1]
        lambda_b_ic = lambda_b_ic[:, :n_state+1]
        sigma_int_ic = sigma_int_ic[:n_state+1]
        h_ic = h_ic[:, :n_state+1]
        final = final + [f"S0"] + [f"S{j}" for j in range(1, 1 + delta_ic.shape[1]) if j != n_state+1]
    elif "t" in initial:
        # Tn to Sm ISC
        final_energy_isc = singlets - (gamma_s + chi_s) * alphast2
        final_energy_isc, singlets, chi_s, gamma_s, socs_complete = sorting_parameters(final_energy_isc, singlets, chi_s, gamma_s, socs_complete)
        lambda_b_isc = lambda_solvent(chis, chi_s, alphaopt2, alphast2)
        delta_isc = final_energy_isc - initial_energy[:, np.newaxis] + lambda_b_isc
        
        final = [
            i.split("_")[2].upper()
            for i in data.columns.values
            if "soc_" + initial.lower() + "_" in i and i.count("t") == 1
        ]
        # Tn to S0 ISC
        socs_s0 = socs_s0[:, n_state]
        socs_complete = np.hstack((socs_s0[:, np.newaxis], socs_complete))
        delta_isc = np.hstack((-1 * delta_emi[:, np.newaxis], delta_isc))
        lambda_b_isc = np.hstack((lambda_emission[:, np.newaxis], lambda_b_isc))
        sigma_int = np.std(delta_isc, axis=0)

    sigma = total_reorganization_energy(lambda_b_isc, kbt, sigma_int)
    y_axis = (
        (2 * np.pi / HBAR_EV) * (socs_complete**2) * nemo.tools.gauss(delta_isc, sigma)
    )
    #tau_D = (refractive_index**2 / eps)*1e-12
    #g_ad = 0#np.nan_to_num((2 * np.pi * (socs_complete**2) * tau_D)/ (HBAR_EV * lambda_b_isc))
    #y_axis =  y_axis / (1 + g_ad)
    if "s" in initial:
        sigma_ic = total_reorganization_energy(lambda_b_ic, kbt, sigma_int_ic)
        #check for 0 in sigma_ic
        y_axis_ic = (
            (2 * np.pi / HBAR_EV) * (h_ic) * nemo.tools.gauss(delta_ic, sigma_ic)
        )

        # hstack y and espectro
        y_axis = np.hstack((y_axis, y_axis_ic))
        sigma = np.hstack((sigma, sigma_ic))
        couplings = np.hstack((socs_complete, np.sqrt(h_ic)))
        gap = np.hstack((delta_isc,delta_ic))
        print('Shape',gap.shape, y_axis.shape)
    else:
        mean_soc = 1000 * means(socs_complete, y_axis, ensemble_average)[:, np.newaxis]
        gap = delta_isc
        couplings = socs_complete

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
                mean_emi_coupling,
                mean_sigma_emi,
                mean_part_emi,
            ]
        ]
    )
    mean_gap = means(gap, y_axis, ensemble_average)[:, np.newaxis]
    mean_sigma = means(sigma, y_axis, ensemble_average)[:, np.newaxis]
    mean_soc = 1000 * means(couplings, y_axis, ensemble_average)[:, np.newaxis] 
    #mean_g = means(g_ad, y_axis, ensemble_average)[:, np.newaxis]
    #print(np.round(mean_g, 3))
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
            "AvgCoupling(meV)",
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
    n_state = int(initial[1:])

    # excited state energies
    engs = fetch(data, [f"^e_{spin}"])

    #ground state susceptibility
    gamma_s0 = data['gamma_s0'].to_numpy()
    gamma_s0 = gamma_s0[:, np.newaxis]

    #excited state susceptibilities
    chis = fetch(data, [f"^chi_{spin}(?!0)"])
    gammas = fetch(data, [f"^gamma_{spin}(?!0)"])

    #oscillator strengths
    oscs = fetch(data, ["^osc_"])

    constante = (
        (np.pi * (E_CHARGE**2) * HBAR_EV)
        / (2 * refractive_index * MASS_E * LIGHT_SPEED * EPSILON_0)
        * 1e20
    )
    
    
    
    if initial == "s0":
        initial_energy = - gamma_s0 * alphast2
        final_energy = engs - gammas * alphast2 - chis * alphaopt2
        lambda_b = lambda_solvent(0, chis, alphaopt2, alphast2)
        deltae = final_energy - initial_energy
    else:
        initial_energy = engs - (gammas + chis) * alphast2
        initial_energy = initial_energy[:, n_state-1]
        chi_initial = chis[:, n_state-1]
        engs = engs[:, n_state:]
        chis = chis[:, n_state:]
        gammas = gammas[:, n_state:]

        final_energy = engs - gammas * alphast2 - chis * alphaopt2
        deltae = final_energy - initial_energy[:,np.newaxis]
        lambda_b = lambda_solvent(chi_initial[:, np.newaxis], chis, alphaopt2, alphast2)


    l_total = total_reorganization_energy(lambda_b, kbt, np.std(deltae, axis=0))
    x_axis = x_values(deltae, l_total)
    if nstates == -1:
        nstates = deltae.shape[1]
    # Add extra dimension to DE and Ltotal to match x shape
    deltae, l_total, oscs, lambda_b = another_dimension(
        nstates, deltae, l_total, oscs, lambda_b
    )
    y_axis = constante * oscs * nemo.tools.gauss(x_axis - deltae, l_total)
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
            initial, spin, n_state, constante*oscs, deltae, chis[:,:nstates], l_total
        )
        return abs_spec, breakdown
    return abs_spec


class Ensemble(object):
    def __init__(self, file, name=''):
        data = pd.read_csv(file)
        initial = data['ensemble'][0]
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