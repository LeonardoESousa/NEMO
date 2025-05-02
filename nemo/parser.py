import sys
import re
import numpy as np

##SOME CONSTANTS##############################################
EPSILON_0 = 8.854187817e-12  # F/m
HBAR_EV = 6.582119514e-16  # eV s
HBAR_J = 1.054571800e-34  # J s
MASS_E = 9.10938356e-31  # kg
LIGHT_SPEED = 299792458  # m/s
E_CHARGE = 1.60217662e-19  # C
BOLTZ_EV = 8.6173303e-5  # eV/K
AMU = 1.660539040e-27  # kg
###############################################################


##ERROR FUNCTION###############################################
def fatal_error(msg):
    print(msg)
    sys.exit()


##LIST OF KEYWORDS THAT SHOULD NOT BE READ#####################
def delist(elem):
    words = [
        "jobtype",
        "-----",
        "cis_n",
        "cis_s",
        "cis_t",
        "gui",
        "nto_",
        "soc",
        "sts_",
        "CIS_RELAXED_DENSITY",
        "solvent_method",
    ]
    for word in words:
        if word in elem.lower():
            return False
    return True


###############################################################


##CHECKS THE FREQUENCY LOG'S LEVEL OF THEORY###################
def busca_input(freqlog):
    charge_mult = None
    if ".out" in freqlog or ".log" in freqlog:
        search = False
    else:
        search = True
    rem = ""
    with open(freqlog, "r", encoding="utf-8") as freq_file:
        for line in freq_file:
            if "User input:" in line and not search:
                search = True
            elif search and delist(line):
                rem += line
            elif (
                "--------------------------------------------------------------" in line
                and search
                and rem != ""
            ):
                search = False
    rem = rem.split("$end")
    remove = ["$molecule", "$comment", "$pcm", "$solvent", "$end"]
    extra = ""
    for element in rem:
        if "$rem" in element:
            rem_section = element + "$end\n"
        elif "$molecule" in element:
            mol_section = element.split("\n")
            for mol in mol_section:
                mol = mol.split()
                if len(mol) == 2:
                    charge_mult = " ".join(mol)
        elif element.strip() and all(x not in element for x in remove):
            extra += element + "$end\n"
    return rem_section, charge_mult, extra


###############################################################


##GETS FREQUENCIES AND REDUCED MASSES##########################
def pega_freq(freqlog):
    freqs, masses = [], []
    with open(freqlog, "r", encoding="utf-8") as freq_file:
        for line in freq_file:
            if "Frequency:" in line:
                line = line.split()
                for j in range(1, len(line)):
                    if float(line[j]) in freqs:
                        pass
                    freqs.append(float(line[j]))
            elif "Red. Mass:" in line:
                line = line.split()
                for j in range(2, len(line)):
                    masses.append(float(line[j]))
    # conversion in angular frequency
    freqs = np.array(freqs) * (LIGHT_SPEED * 100 * 2 * np.pi)
    if len(freqs) == 0:
        fatal_error("No frequencies in the log file! Goodbye!")
    # conversion from amu to kg
    masses = np.asarray(masses) * AMU
    return freqs, masses


###############################################################


##GETS ATOMS AND LAST GEOMETRY IN FILE#########################
def pega_geom(freqlog):
    if ".out" in freqlog or ".log" in freqlog:
        start = False
        n_value = 0
        with open(freqlog, "r", encoding="utf-8") as freq_file:
            for line in freq_file:
                if "I     Atom           X                Y                Z" in line:
                    geometry = np.zeros((1, 3))
                    atomos = []
                    start = True
                elif start:
                    line = line.split()
                    try:
                        new_row = []
                        for j in range(2, len(line)):
                            new_row.append(float(line[j]))
                        atomos.append(line[1])
                        geometry = np.vstack((geometry, new_row))
                    except (ValueError, IndexError):
                        n_value += 1
                if n_value == 2:
                    start = False
                    n_value = 0

    else:
        geometry = np.zeros((1, 3))
        atomos = []
        with open(freqlog, "r", encoding="utf-8") as freq_file:
            for line in freq_file:
                line = line.split()
                try:
                    vetor = np.array([float(line[1]), float(line[2]), float(line[3])])
                    atomos.append(line[0])
                    geometry = np.vstack((geometry, vetor))
                except (ValueError, IndexError):
                    pass
    try:
        geometry = geometry[1:, :]
    except IndexError:
        fatal_error("No geometry in the log file! Goodbye!")
    return geometry, atomos


###############################################################


def pega_modos(G, freqlog):
    freqs, _ = pega_freq(freqlog)
    num_atoms = np.shape(G)[0]
    normal_modes = np.zeros((num_atoms, 3, len(freqs))) + np.nan
    mode = 0
    fetch = False
    atom = 0
    with open(freqlog, "r", encoding='utf-8') as f:
        for line in f:
            if "X      Y      Z" in line:
                fetch = True
            elif fetch:
                line = line.split()
                if 'TransDip' in line:
                    fetch = False
                    mode += 3
                    atom = 0
                else:
                    for j in range(1, len(line)):
                        normal_modes[atom, (j - 1) % 3, mode + (j - 1) // 3] = float(
                            line[j]
                        )
                    atom += 1
    return normal_modes

def parse_block(block, collect_corrections=False):
    """
    Parses a calculation block from a quantum chemistry log file, including composition.

    Parameters:
        block (str): The text block corresponding to one calculation.
        collect_corrections (bool): Whether to parse SS corrections.

    Returns:
        dict: Contains excited-state energies, spins, oscillator strengths, state indices,
              corrections, total energies, and composition.
    """
    data = {
        'energies': [],
        'spins': [],
        'oscillator': [],
        'indices': [],
        'correction': [],
        'correction2': [],
        'total_energy': [],
        'composition': [],
    }
    exc = False  # Flag for excited state section.
    corr = False # Flag for PCM correction section.
    current_comp = None  # To hold composition dict for each state.

    for line in block.splitlines():
        # Start a new excited state section.
        if "TDDFT/TDA Excitation Energies" in line or "TDDFT Excitation Energies" in line:
            data['energies'] = []
            data['spins'] = []
            data['oscillator'] = []
            data['indices'] = []
            data['composition'] = []
            exc = True
        elif exc:
            if "Excited state" in line:
                # Save previous state's composition before starting a new one
                if current_comp is not None:
                    data['composition'].append(current_comp)
                current_comp = {}  # Start a new composition dict

                parts = line.split()
                state_index = int(parts[2].replace(":", ""))
                data['indices'].append(state_index)

            elif "Total energy for state" in line:
                energy_val = float(line.split()[-2]) * 27.21139
                data['energies'].append(energy_val)

            elif "Multiplicity" in line:
                data['spins'].append(line.split()[-1])

            elif "Strength" in line:
                osc_val = float(line.split()[-1])
                data['oscillator'].append(osc_val)

            elif re.match(r'\s*X:\s*D\(', line):
                # Composition line: parse D, V, amplitude
                match = re.search(r'D\(\s*(\d+)\)\s*-->\s*V\(\s*(\d+)\)\s*amplitude\s*=\s*([-\d\.Ee]+)', line)
                if match:
                    d_idx = int(match.group(1))
                    v_idx = int(match.group(2))
                    amp = float(match.group(3))
                    current_comp[(d_idx, v_idx)] = amp

            elif "---------------------------------------------------" in line and len(data['energies']) > 0:
                # End of excited state section: append last composition
                if current_comp is not None:
                    data['composition'].append(current_comp)
                    current_comp = None
                exc = False

        # PCM correction info
        if collect_corrections:
            if "Excited-state properties with   relaxed density" in line:
                corr = True
            if corr:
                if "SS-PCM correction" in line:
                    val = float(line.split()[-2])
                    data['correction'].append(-1 * np.nan_to_num(val))
                elif "LR-PCM correction" in line:
                    val = float(line.split()[-2])
                    data['correction2'].append(-1 * np.nan_to_num(val))
                if "------------------------ END OF SUMMARY -----------------------" in line:
                    corr = False

        # Total energy in final basis set
        if "Total energy" in line and "=" in line:
            parts = line.split()
            data['total_energy'].append(float(parts[-1]) * 27.21139)

    # If the last state's composition wasn't added yet
    if current_comp is not None:
        data['composition'].append(current_comp)

    # Postprocessing
    spins = np.array(data['spins'])
    singlet_idx = np.where(spins == "Singlet")[0]
    triplet_idx = np.where(spins == "Triplet")[0]
    data['ind_s'] = np.array(data['indices'])[singlet_idx]
    data['ind_t'] = np.array(data['indices'])[triplet_idx]
    data['singlets'] = np.array(data['energies'])[singlet_idx]
    data['triplets'] = np.array(data['energies'])[triplet_idx]
    data['osc_singlets'] = np.array(data['oscillator'])[singlet_idx]
    data['comp_singlets'] = np.array(data['composition'])[singlet_idx]
    data['comp_triplets'] = np.array(data['composition'])[triplet_idx]
    del data['composition']
    if collect_corrections:
        data['ss_s'] = np.array(data['correction'])[singlet_idx] + np.array(data['correction2'])[singlet_idx]
        data['ss_t'] = np.array(data['correction'])[triplet_idx] + np.array(data['correction2'])[triplet_idx]
    data['len'] = len(data['energies']) // 2

    return data


def compute_similarity(comp1, comp2):
    """
    Compute cosine similarity between two compositions (dict of (D,V): amplitude).
    Only overlapping keys are compared.
    """
    keys = set(comp1.keys()) & set(comp2.keys())
    if not keys:
        return 0.0
    v1 = np.array([comp1[k] for k in keys])**2
    v2 = np.array([comp2[k] for k in keys])**2
    return np.dot(v1, v2)

def match_compositions(comp_singlets_1, comp_singlets_2):
    """
    Matches each state in comp_singlets_1 to a unique state in comp_singlets_2
    by highest composition similarity (greedy one-to-one match, always assigns).

    Parameters:
        comp_singlets_1: array/list of dicts for case 1
        comp_singlets_2: array/list of dicts for case 2

    Returns:
        list of indices: index i gives index in comp_singlets_2 matching state i in comp_singlets_1
    """
    available = list(range(len(comp_singlets_1)))  # indices of available states in comp_singlets_2
    match_indices = []

    for comp1 in comp_singlets_1:
        best_idx = available[0]  # fallback: pick first available
        best_sim = -1
        for j in available:
            sim = compute_similarity(comp1, comp_singlets_2[j])
            #print(f"Comparing {comp2} with {comp_singlets_1[j]}: similarity = {sim}")
            if sim > best_sim:
                best_sim = sim
                best_idx = j
        match_indices.append(best_idx)
        available.remove(best_idx)

    return match_indices


def pega_energias(file):
    """
    Extracts excited-state and total energy properties from a quantum chemistry log file.

    The log file is expected to contain two calculations separated by the marker "Have a nice day".
    The first calculation provides vacuum (reference) excited-state energies, and the second includes
    PCM corrections. The function returns sorted arrays of singlet and triplet energies (vacuum),
    their oscillator strengths and state indices (for singlets), the calculated solvent shifts, and
    the total energy (and its difference between the two calculations) in eV.

    Parameters:
        file (str): Path to the log file.

    Returns:
        tuple: (singlets, triplets, oscs, ind_s, ind_t, ss_s, ss_t, total_energy, energy_diff)
            where:
              - singlets (np.array): Vacuum energies for singlet excited states.
              - triplets (np.array): Vacuum energies for triplet excited states.
              - oscs (np.array): Oscillator strengths for singlet states.
              - ind_s (np.array): State indices for singlet states.
              - ind_t (np.array): State indices for triplet states.
              - ss_s (np.array): Solvent shifts for singlet states.
              - ss_t (np.array): Solvent shifts for triplet states.
              - total_energy (float): Total energy of S0.
              - s0_corr (float): Solvent correction to S0.
    """
    with open(file, "r", encoding="utf-8") as f:
        content = f.read()

    # Split the file into two calculation blocks.
    blocks = content.split("Have a nice day")
    if len(blocks) < 2:
        raise ValueError("Log file does not contain two calculations as expected.")

    # Parse the vacuum (first) and PCM-corrected (second) calculation blocks.
    vac_data = parse_block(blocks[0], collect_corrections=False)
    corr_data = parse_block(blocks[1], collect_corrections=True)
    min_len = min(vac_data['len'], corr_data['len'])
    match_singlets = match_compositions(vac_data['comp_singlets'][:min_len], corr_data['comp_singlets'][:min_len])
    match_triplets = match_compositions(vac_data['comp_triplets'][:min_len], corr_data['comp_triplets'][:min_len])
    
    
    singlets_vac = vac_data['singlets'][:min_len]
    triplets_vac = vac_data['triplets'][:min_len]
    singlets_pcm = corr_data['singlets'][match_singlets][:min_len]
    triplets_pcm = corr_data['triplets'][match_triplets][:min_len]
    
    oscs = vac_data['osc_singlets'][:min_len]
    ind_s = vac_data['ind_s'][:min_len]
    ind_t = vac_data['ind_t'][:min_len]
    s0_vac = vac_data['total_energy'][0]
    s0_pcm = corr_data['total_energy'][0]
    s0_corr = s0_vac - s0_pcm
    ss_s = corr_data['ss_s'][match_singlets][:min_len]
    ss_t = corr_data['ss_t'][match_triplets][:min_len]
    y_s = (singlets_vac - s0_vac) - (singlets_pcm - s0_pcm)
    y_t = (triplets_vac - s0_vac) - (triplets_pcm - s0_pcm)
    y_s[y_s < 0] = 0
    y_t[y_t < 0] = 0  
    return singlets_vac, triplets_vac, oscs, ind_s, ind_t, ss_s, ss_t, s0_vac, s0_corr, y_s, y_t

#########################################################################################


##GETS SOC BETWEEN Sn STATE AND TRIPLETS#################################################
def pega_soc_singlet(file, n_state, ind_s, ind_t):
    socs = []
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
def pega_soc_triplet(file, n_state, ind_s, ind_t):
    socs = []
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
def pega_soc_ground(file, n_state, ind_s, ind_t):
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
def pega_soc_triplet_triplet(file, n_state, ind_s, ind_t):
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
def avg_socs(files, tipo, n_state, ind_s, ind_t):
    col = None
    if tipo == "singlet":
        pega_soc = pega_soc_singlet
    elif tipo == "triplet":
        pega_soc = pega_soc_triplet
    elif tipo == "ground":
        pega_soc = pega_soc_ground
    elif tipo == "tts":
        pega_soc = pega_soc_triplet_triplet
    i = 0    
    for file in files:
        socs = pega_soc(file, n_state, ind_s[i,:], ind_t[i,:])
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
                    except (NameError,ValueError):
                        dipoles = np.zeros((1, len(line) - 1))
                        dipoles = np.vstack((dipoles, dipole_line))
                elif line[1] == str(ind[state]) and int(line[0]) < int(line[1]):
                    dipole_line = [float(line[0])]
                    extra_line = [float(line[i]) for i in range(2, len(line))]
                    dipole_line.extend(extra_line)
                    dipole_line = np.array([dipole_line])
                    try:
                        dipoles = np.vstack((dipoles, dipole_line))
                    except (NameError,ValueError):
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

##CALCULATES TRANSITION DIPOLE MOMENTS FOR Tn TO S0 TRANSITIONS##########################
def moment(file, ess, ets, e_s0, dipss, dipts, n_triplet, ind_s, ind_t):
    # Conversion factor between a.u. = e*bohr to SI
    conversion = 8.4783533e-30
    fake_t = np.where(np.sort(ind_t) == ind_t[n_triplet])[0][0]
    ess = np.array(ess)
    ets = np.array(ets)
    ess = np.insert(ess, 0, e_s0)
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
            part_1 = (socss0 / (e_s0 - ets)) * dipts[:, i]
            part_1 = np.sum(part_1)
            part_2 = (socst1 / (ets[n_triplet] - ess)) * dipss[:, i]
            part_2 = np.sum(part_2)
            complex_dipole = part_1 + part_2
            # append magnitude squared
            moments.append((complex_dipole * complex_dipole.conjugate()).real)

    moments = np.array(moments)
    moments = np.sum(moments) * (conversion**2)
    return moments

def phosph_osc(file, n_state, ind_s, ind_t, singlets, triplets, e_s0): 
    zero = ["0"]
    zero.extend(ind_s)
    total_moments = []
    ground_dipoles = pega_dipolos(
        file, zero, "Electron Dipole Moments of Ground State", 0
    )
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
        moments = moment(
            file,
            singlets,
            triplets,
            e_s0,
            ground_dipoles,
            triplet_dipoles,
            n_triplet,
            ind_s,
            ind_t,
        )
        total_moments.append(moments)
    total_moments = np.array(total_moments)
    term = E_CHARGE * (HBAR_J**2) / (triplets - e_s0)
    osc_strength = (2 * MASS_E) * total_moments / (3 * term)
    return osc_strength[np.newaxis, :]

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
            check = False
            for line in log_file:
                if frase in line:
                    oscs = []
                    dip = True
                    check = False
                elif dip and "--" not in line:
                    if 'States' not in line:
                        check = True
                    line = line.split()
                    if (line[0] == ind and int(line[1]) in ind_s) or (
                        line[1] == ind and int(line[0]) in ind_s
                    ):
                        oscs.append(float(line[-1]))
                elif check and "---" in line:
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
                    real_part = (
                        line[1].replace("(", "").replace(")", "").replace("i", "")
                    )
                    real_part = float(
                        real_part.replace("--", "+")
                        .replace("+-", "-")
                        .replace("-+", "-")
                        .replace("++", "+")
                    )
                    img_part = line[2] + line[3].replace("(", "").replace(
                        ")", ""
                    ).replace("i", "")
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
                    real_part = (
                        line[1].replace("(", "").replace(")", "").replace("i", "")
                    )
                    real_part = float(
                        real_part.replace("--", "+")
                        .replace("+-", "-")
                        .replace("-+", "-")
                        .replace("++", "+")
                    )
                    img_part = line[2] + line[3].replace("(", "").replace(
                        ")", ""
                    ).replace("i", "")
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
