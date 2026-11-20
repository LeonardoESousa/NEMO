import sys
import re
from math import fsum

import numpy as np
from scipy.optimize import linear_sum_assignment

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

def check_software(freqlog):
    """
    Checks if the frequency log file is from Gaussian or Q-Chem.
    Returns True if it is a Gaussian log file, False otherwise.
    """
    with open(freqlog, "r", encoding="utf-8") as freq_file:
        content = freq_file.read()
        if "Entering Gaussian System" in content:
            return 'Gaussian'
        elif "Q-Chem" in content:
            return 'QChem'
        else:
            fatal_error("Unknown log file format! Goodbye!")


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
def pega_freq_gaussian(freqlog):
    freqs, masses = [], []
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Frequencies -- " in line:
                line = line.split()
                for j in range(2, len(line)):
                    freqs.append(float(line[j]))
            elif "Red. masses --" in line:
                line = line.split()
                for j in range(3, len(line)):
                    masses.append(float(line[j]))
            elif "Thermochemistry" in line:
                break
    if len(freqs) == 0 or len(masses) == 0:
        fatal_error(
            "No frequencies or reduced masses found. Check your frequency file. Goodbye."
        )
    # conversion in angular frequency
    freqs = np.array(freqs) * (LIGHT_SPEED * 100 * 2 * np.pi)
    # conversion from amu to kg
    masses = np.asarray(masses) * AMU
    return freqs, masses

##GETS FREQUENCIES AND REDUCED MASSES##########################
def pega_freq_qchem(freqlog):
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

def pega_freq(freqlog):
    software = check_software(freqlog)
    if software == 'Gaussian':
        return pega_freq_gaussian(freqlog)
    elif software == 'QChem':
        return pega_freq_qchem(freqlog)


###############################################################

##GETS INPUT PARAMS FROM LOG FILES#############################
def get_input_params(freqlog):
    nproc, mem, header = "", "", ""
    with open(freqlog, "r",encoding='utf-8') as f:
        search = False
        for line in f:
            if "%nproc" in line.lower():
                line = line.split("=")
                nproc = line[-1].replace("\n", "")
            elif "%mem" in line.lower():
                line = line.split("=")
                mem = line[-1].replace("\n", "")
            elif "#" in line and not search and header == "":
                search = True
                header += line.lstrip().replace("\n", "")
            elif search and "----------" not in line:
                header += line.lstrip().replace("\n", "")
            elif search and "----------" in line:
                search = False
                break
    return nproc, mem, header


###############################################################


##CHECKS FREQ FILES############################################
def double_check_gaussian(freqlog):
    freqs, _ = pega_freq(freqlog)
    if np.any(freqs < 0):
        fatal_error(
            "Imaginary frequencies detected. Check your frequency file. Goodbye."
        )
    _, _, header = get_input_params(freqlog)
    header = header.lower()
    optfreqissue = False
    stationary = False
    if "opt" in header and "freq" in header and "iop(" in header:
        optfreqissue = True
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Stationary point found" in line:
                stationary = True
            elif all(s in line for s in ["Item", "Value", "Threshold", "Converged?"]):
                stationary = False
    if not stationary:
        print("*" * 50)
        print("WARNING: Non-optimized parameters detected in your frequency file.")
        print(
            "Even though the frequencies may be all real, your structure is still not fully optimized."
        )
        print("This may lead to inaccurate results.")
        if optfreqissue:
            print(
                "In your case, this may be due to running an opt freq calculation using a single input file and IOP options."
            )
            print(
                "Gaussian does not carry the IOP options to the frequency calculation when using a single input file."
            )
            print(
                "To avoid this issue, run the optimization and frequency calculations separately."
            )
        else:
            print("To learn more about this issue, check https://gaussian.com/faq3/ .")
        print("Proceed at your own risk.")
        print("*" * 50)
        print("\n")
###############################################################

def double_check(freqlog):
    software = check_software(freqlog)
    if software == 'Gaussian':
        double_check_gaussian(freqlog)
    elif software == 'QChem':
        # Q-Chem does not have the same issues as Gaussian, so no checks are needed.
        pass

##FETCHES CHARGE AND MULTIPLICITY FROM G16#####################
def get_cm_gaussian(freqlog):
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Charge" in line and "Multiplicity" in line:
                line = line.split()
                charge = line[2]
                mult = line[5]
                break
    return charge + " " + mult
###############################################################

def get_cm_qchem(freqlog):
    _, cm, _ = busca_input(freqlog)
    return cm

def get_cm(freqlog):
    software = check_software(freqlog)
    if software == 'Gaussian':
        return get_cm_gaussian(freqlog)
    elif software == 'QChem':
        return get_cm_qchem(freqlog)


def pega_geom(freqlog):
    software = check_software(freqlog)
    if software == 'Gaussian':
            return pega_geom_gaussian(freqlog)
    elif software == 'QChem':
            return pega_geom_qchem(freqlog)

def pega_geom_gaussian(freqlog):
    elements = (
        "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca "
        "Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr "
        "Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
        "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
        "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm "
        "Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
    ).split()
    if ".log" in freqlog:
        busca = " orientation:"
        fetch = False
        with open(freqlog, "r",encoding='utf-8') as f:
            for line in f:
                if busca in line and "Dipole" not in line:
                    geom = np.zeros((1, 3))
                    atomos = []
                    fetch = True
                elif fetch:
                    line = line.split()
                    try:
                        float(line[0])
                        NG = []
                        for j in range(3, len(line)):
                            NG.append(float(line[j]))
                        atomos.append(line[1])
                        geom = np.vstack((geom, NG))
                    except ValueError:
                        if len(atomos) > 0:
                            fetch = False
    else:
        geom = np.zeros((1, 3))
        atomos = []
        with open(freqlog, "r",encoding='utf-8') as f:
            for line in f:
                line = line.split()
                try:
                    vetor = np.array([float(line[1]), float(line[2]), float(line[3])])
                    atomos.append(line[0])
                    geom = np.vstack((geom, vetor))
                except (IndexError, ValueError):
                    pass
    geom = geom[1:, :]
    atomos = [elements[int(atom)] if atom.isdigit() else atom for atom in atomos]
    return geom, atomos

##GETS ATOMS AND LAST GEOMETRY IN FILE#########################
def pega_geom_qchem(freqlog):
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


##GETS NORMAL COORDINATES IN HIGH PRECISION####################
def pega_modosHP(G, freqlog):
    F, _ = pega_freq(freqlog)
    num_atoms = np.shape(G)[0]
    normal_modes = np.zeros((num_atoms, 3, len(F))) + np.nan
    mode = 0
    fetch = False
    last_mode = False
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Coord Atom Element:" in line:
                fetch = True
            elif fetch:
                line = line.split()
                if "Harmonic" in line[0]:
                    fetch=False
                    continue
                if len(line) < 6:
                    if last_mode:
                        coord = int(line[0]) - 1
                        atom = int(line[1]) - 1
                        for j in range(3, len(line)):
                            normal_modes[atom, coord, mode + j - 3] = float(line[j])
                    if not last_mode:
                        fetch = False
                        mode += 5
                        if len(line)<5: 
                            last_mode=True
                else:
                    coord = int(line[0]) - 1
                    atom = int(line[1]) - 1
                    for j in range(3, len(line)):
                        normal_modes[atom, coord, mode + j - 3] = float(line[j])
    return normal_modes


###############################################################


##GETS NORMAL COORDINATES IN REGULAR PRECISION#################
def pega_modosLP(G, freqlog):
    F, _ = pega_freq(freqlog)
    num_atoms = np.shape(G)[0]
    normal_modes = np.zeros((num_atoms, 3, len(F))) + np.nan
    mode = 0
    fetch = False
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Atom  AN      X      Y      Z" in line:
                fetch = True
            elif fetch:
                line = line.split()
                if len(line) < 5:
                    fetch = False
                    mode += 3
                else:
                    atom = int(line[0]) - 1
                    for j in range(2, len(line)):
                        normal_modes[atom, (j - 2) % 3, mode + (j - 2) // 3] = float(
                            line[j]
                        )
    return normal_modes


###############################################################


##DETECTS WHETHER HIGH PRECISION IS USED#######################
def pega_modos_gaussian(geom, freqlog):
    x = "LP"
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Coord Atom Element:" in line:
                x = "HP"
                break
    if x == "LP":
        return pega_modosLP(geom, freqlog)
    else:
        return pega_modosHP(geom, freqlog)


def pega_modos_qchem(G, freqlog):
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

def pega_modos(G, freqlog):
    software = check_software(freqlog)
    if software == 'Gaussian':
        return pega_modos_gaussian(G, freqlog)
    elif software == 'QChem':
        return pega_modos_qchem(G, freqlog)


def parse_block(block, collect_corrections=False):
    data = {
        'energies': [],
        'spins': [],
        'oscillator': [],
        'indices': [],
        'correction': [],
        'correction2': [],
        'total_energy': [],
        'composition': [],
        'transition_dipoles': [],
    }
    orientation = block.split('Standard Nuclear Orientation (Angstroms)', 1)
    geometry = orientation[1].split('Nuclear Repulsion Energy', 1)[0] if len(orientation) == 2 else ''
    rows = re.findall(r'^\s*\d+\s+([A-Z][a-z]?)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*$',
                      geometry, re.MULTILINE)
    data['atoms'] = [row[0] for row in rows]
    data['geometry'] = np.array([row[1:] for row in rows], dtype=float).reshape(-1, 3)
    exc = False
    corr = False
    current_comp = None

    for line in block.splitlines():
        # Excited-state section
        if "TDDFT/TDA Excitation Energies" in line or "TDDFT Excitation Energies" in line:
            data['energies'] = []
            data['spins'] = []
            data['oscillator'] = []
            data['indices'] = []
            data['composition'] = []
            data['transition_dipoles'] = []
            exc = True
            current_comp = None
            continue

        if exc:
            if "Excited state" in line:
                if current_comp is not None:
                    data['composition'].append(current_comp)

                current_comp = {}

                parts = line.split()
                state_index = int(parts[2].replace(":", ""))
                energy_val = float(line.split()[-1])

                data['indices'].append(state_index)
                data['energies'].append(energy_val)
                data['transition_dipoles'].append([np.nan, np.nan, np.nan])

            elif "Multiplicity" in line:
                data['spins'].append(line.split()[-1])

            elif "Strength" in line:
                data['oscillator'].append(float(line.split()[-1]))

            elif 'Trans. Mom.:' in line:
                data['transition_dipoles'][-1] = [float(line.split()[i]) for i in (2, 4, 6)]

            # TDA prints D(...); full TDDFT prefixes excitation amplitudes with X:.
            # Keep Y: de-excitations separate from the excitation composition.
            elif re.match(r'\s*(?:X:\s*)?D\(', line):
                match = re.search(
                    r'D\(\s*(\d+)\)\s*-->\s*V\(\s*(\d+)\)\s*amplitude\s*=\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)',
                    line
                )
                if match and current_comp is not None:
                    d_idx = int(match.group(1))
                    v_idx = int(match.group(2))
                    amp = float(match.group(3))
                    current_comp[(d_idx, v_idx)] = amp

            elif "---------------------------------------------------" in line and len(data['energies']) > 0:
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
                    data['correction'].append(-np.nan_to_num(val))

                elif "LR-PCM correction" in line:
                    val = float(line.split()[-2])
                    data['correction2'].append(-np.nan_to_num(val))

                elif "------------------------ END OF SUMMARY -----------------------" in line:
                    corr = False


        if "Total energy" in line and "=" in line:
            parts = line.split()
            data['total_energy'].append(float(parts[-1]) * 27.21139)

    if current_comp is not None:
        data['composition'].append(current_comp)

    spins = np.array(data['spins'])
    energies = np.array(data['energies'])
    oscillator = np.array(data['oscillator'])
    composition = np.array(data['composition'], dtype=object)

    singlet_idx = np.where(spins == "Singlet")[0]
    triplet_idx = np.where(spins == "Triplet")[0]

    data['singlets'] = energies[singlet_idx]
    data['triplets'] = energies[triplet_idx]
    data['osc_singlets'] = oscillator[singlet_idx]
    data['comp_singlets'] = composition[singlet_idx]
    data['comp_triplets'] = composition[triplet_idx]
    data['tdm_singlets'] = np.array(data.pop('transition_dipoles')).reshape(-1, 3)[singlet_idx]

    del data['composition']

    if collect_corrections:
        correction = np.array(data['correction'])
        correction2 = np.array(data['correction2'])

        data['ss_s'] = correction[singlet_idx] + correction2[singlet_idx]
        data['ss_t'] = correction[triplet_idx] + correction2[triplet_idx]

    data['len_s'] = len(data['singlets'])
    data['len_t'] = len(data['triplets'])

    return data


def compute_similarity(comp1, comp2):
    """
    Cosine similarity between two excitation compositions.

    Uses squared amplitudes as weights.
    Missing transitions are treated as zero.
    """
    keys = set(comp1.keys()) | set(comp2.keys())

    if not keys:
        return 0.0

    v1 = np.array([comp1.get(k, 0.0) ** 2 for k in keys])
    v2 = np.array([comp2.get(k, 0.0) ** 2 for k in keys])

    denom = np.linalg.norm(v1) * np.linalg.norm(v2)

    if denom == 0:
        return 0.0

    return np.dot(v1, v2) / denom


def composition_scores(comp_list_1, comp_list_2):
    """Similarity matrix in the printed occupied/virtual orbital labels."""
    return np.array([[compute_similarity(a, b) for b in comp_list_2]
                     for a in comp_list_1], dtype=float).reshape(
                         len(comp_list_1), len(comp_list_2))


def _dipole_rotation(vacuum, pcm):
    """Align vacuum vectors with the PCM nuclear frame."""
    x, y = vacuum['geometry'], pcm['geometry']
    if (x.shape != y.shape or len(x) < 3 or vacuum['atoms'] != pcm['atoms']
            or not np.isfinite(x).all() or not np.isfinite(y).all()):
        return None
    x, y = x - x.mean(axis=0), y - y.mean(axis=0)
    u, singular, vh = np.linalg.svd(x.T @ y)
    if singular[1] < 1e-10:
        return None  # A linear geometry does not fix rotation about its molecular axis.
    rotation = u @ np.diag([1., 1., np.linalg.det(u @ vh)]) @ vh
    if np.sqrt(np.mean(np.sum((x @ rotation - y)**2, axis=1))) > 1e-3:
        return None
    return rotation


def _dipole_scores(a, b, rotation, transition=False):
    """Compare aligned vectors, leaving unavailable comparisons as NaN."""
    if rotation is None:
        return np.full((len(a), len(b)), np.nan)
    a = a @ rotation
    aa, bb = np.sum(a*a, axis=1)[:, None], np.sum(b*b, axis=1)[None, :]
    dots = a @ b.T
    if transition:
        # A transition vector can change sign with the phase of the excited state.
        dots = np.where((aa > 1e-6) & (bb > 1e-6), np.abs(dots), 0.)
    else:
        dots = np.maximum(dots, 0.)
    floor = 1e-6 if transition else 0.
    return np.divide(2*dots, aa + bb, out=np.full(dots.shape, np.nan),
                     where=np.isfinite(aa + bb) & ((aa > floor) | (bb > floor)))


def _singlet_anchors(vacuum, pcm):
    """High-similarity transition dipoles for inferring the orbital map."""
    scores = np.nan_to_num(_dipole_scores(vacuum['tdm_singlets'], pcm['tdm_singlets'],
                                         _dipole_rotation(vacuum, pcm), transition=True))
    n_vac, n_pcm = scores.shape
    # Dummy columns leave weak optical matches out of orbital-map inference.
    rows, cols = linear_sum_assignment(-np.hstack((scores - 0.9, np.zeros((n_vac, n_vac)))))
    return [(i, j, scores[i, j]) for i, j in zip(rows, cols)
            if j < n_pcm and scores[i, j] > 0.9
            and any(vacuum['comp_singlets'][i].values()) and any(pcm['comp_singlets'][j].values())]


def _triplet_dipole_changes(block, data, pcm=False):
    """Read unrelaxed triplet minus ground-state permanent dipoles in Debye."""
    number = r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)'
    if pcm:
        ground = re.findall(r'Dipole Moment \(Debye\)\s*\n\s*X\s+' + number
                            + r'\s+Y\s+' + number + r'\s+Z\s+' + number, block)
        summary = block.partition('Excited-state properties with unrelaxed density')[2]
        summary = summary.partition('END OF SUMMARY')[0]
        states = re.split(r'Unrelaxed (?:singlet|triplet) excited state\s+(\d+):', summary)
        moments = {}
        for root, state in zip(states[1::2], states[2::2]):
            match = re.search(r'Dipole Mom\. \(Debye\):\s*' + number + r'\s+X\s+'
                              + number + r'\s+Y\s+' + number + r'\s+Z', state)
            if match:
                moments[int(root)] = match.groups()
        ground = np.array(ground[-1] if ground else [np.nan]*3, dtype=float)
        scale = 1.
    else:
        moments = {}
        for label in ('Ground State', 'Triplet Excited State'):
            _, header, table = block.rpartition('Electron Dipole Moments of ' + label)
            if not header:
                continue
            sections = re.split(r'^[ \t]*-{5,}[ \t]*$', table, maxsplit=3, flags=re.MULTILINE)
            if len(sections) == 4:
                rows = re.findall(r'^[ \t]*(\d+)[ \t]+' + r'[ \t]+'.join([number]*3)
                                  + r'[ \t]*$', sections[2], re.MULTILINE)
                moments.update({int(root): (x, y, z) for root, x, y, z in rows})
        ground = np.array(moments.get(0, [np.nan]*3), dtype=float)
        scale = -2.541746  # The vacuum table prints electron position moments in a.u.
    roots = [root for root, spin in zip(data['indices'], data['spins']) if spin == 'Triplet']
    excited = np.array([moments.get(root, [np.nan]*3) for root in roots], dtype=float)
    return scale * (excited.reshape(-1, 3) - ground)


def _orbital_relabeling(vacuum, pcm):
    """Fit one PCM->vacuum orbital permutation to the singlet anchors only.

    Alternate occupied/virtual assignments; a 0.03 penalty per changed label
    favors small permutations over fitting changes in orbital mixing.
    """
    anchors = _singlet_anchors(vacuum, pcm)
    if not anchors:
        return {}, {}
    pairs = [(vacuum['comp_singlets'][i], pcm['comp_singlets'][j]) for i, j, _ in anchors]
    orbitals = [sorted({key[axis] for pair in pairs for comp in pair for key in comp})
                for axis in (0, 1)]
    lookup = [{label: i for i, label in enumerate(labels)} for labels in orbitals]
    # Squared amplitudes indexed by anchor, occupied orbital and virtual orbital.
    a = np.zeros((len(pairs), len(orbitals[0]), len(orbitals[1])))
    b = np.zeros_like(a)
    for k, pair in enumerate(pairs):
        for target, comp in zip((a, b), pair):
            for (occupied, virtual), amplitude in comp.items():
                target[k, lookup[0][occupied], lookup[1][virtual]] = amplitude**2
    weights = np.array([score for _, _, score in anchors]) / (
        np.linalg.norm(a, axis=(1, 2)) * np.linalg.norm(b, axis=(1, 2)))
    identity = [np.arange(len(labels)) for labels in orbitals]
    # Marginal signatures can seed a joint occupied/virtual relabeling when
    # neither type of label initially agrees. Also try the original labels.
    marginal = []
    for axis, order in zip((2, 1), identity):
        cost = np.einsum('ki,kj,k->ij', b.sum(axis=axis), a.sum(axis=axis), weights)
        marginal.append(linear_sum_assignment(-(cost + 0.03*np.eye(len(order))))[1])
    best_score, best_orders = -np.inf, identity
    for occupied, virtual in (identity, marginal):
        for _ in range(20):
            old_occupied, old_virtual = occupied.copy(), virtual.copy()
            # Assignment rows are PCM labels; columns are vacuum labels.
            cost = np.einsum('kiv,kjv,k->ij', b[:, :, np.argsort(virtual)], a, weights)
            occupied = linear_sum_assignment(-(cost + 0.03*np.eye(len(occupied))))[1]
            cost = np.einsum('koi,koj,k->ij', b[:, np.argsort(occupied), :], a, weights)
            virtual = linear_sum_assignment(-(cost + 0.03*np.eye(len(virtual))))[1]
            if np.array_equal(occupied, old_occupied) and np.array_equal(virtual, old_virtual):
                break
        mapped = b[:, np.argsort(occupied), :][:, :, np.argsort(virtual)]
        changes = sum(np.count_nonzero(order != original)
                      for order, original in zip((occupied, virtual), identity))
        score = np.einsum('kov,kov,k->', a, mapped, weights) - 0.03*changes
        if score > best_score + 1e-12:
            best_score, best_orders = score, (occupied, virtual)
    return tuple({labels[i]: labels[j] for i, j in enumerate(order) if i != j}
                 for labels, order in zip(orbitals, best_orders))


def _relabel_compositions(compositions, occupied, virtual):
    return [{(occupied.get(i, i), virtual.get(a, a)): c for (i, a), c in comp.items()}
            for comp in compositions]


def match_compositions(comp_list_1, comp_list_2, return_scores=False, dipole_scores=None):
    """
    Globally matches states from comp_list_1 to comp_list_2
    using composition similarity, averaged with dipole similarity where available.

    With return_scores, also return the composition similarity matrix.
    """
    indices, matrix, _ = _match_compositions(comp_list_1, comp_list_2, dipole_scores)
    return (indices, matrix) if return_scores else indices


def _energy_scores(vacuum, pcm):
    """Compare excitation energies in units of a typical same-spin root spacing."""
    gaps = np.concatenate((np.diff(np.sort(vacuum)), np.diff(np.sort(pcm))))
    gaps = gaps[np.isfinite(gaps) & (gaps > 0)]
    if not gaps.size:
        return None
    scale = np.median(gaps)
    return scale / (scale + np.abs(vacuum[:, None] - pcm))


def _match_compositions(comp_list_1, comp_list_2, dipole_scores, energy_scores=None):
    n1 = len(comp_list_1)
    n2 = len(comp_list_2)

    sim_matrix = composition_scores(comp_list_1, comp_list_2)
    if n1 == 0 or n2 == 0:
        return np.array([], dtype=int), sim_matrix, 0.

    scores = sim_matrix.copy()
    counts = np.ones(sim_matrix.shape, dtype=int)
    for extra in (dipole_scores, energy_scores):
        if extra is not None:
            valid = np.isfinite(extra)
            scores[valid] += extra[valid]
            counts += valid
    scores /= counts
    rows, cols = linear_sum_assignment(-scores)

    match_indices = np.full(n1, -1, dtype=int)
    match_indices[rows] = cols

    return match_indices, sim_matrix, fsum(scores[rows, cols])


def pega_energias(file, return_matching=False):
    """Read energies/corrections; optionally append vacuum/PCM matching diagnostics."""
    with open(file, "r", encoding="utf-8") as f:
        content = f.read()

    blocks = content.split("Have a nice day")

    if len(blocks) < 3:
        vac_data = parse_block(blocks[0], collect_corrections=False)
        corr_data = parse_block(blocks[0], collect_corrections=True)
    else:
        vac_data = parse_block(blocks[0], collect_corrections=False)
        corr_data = parse_block(blocks[1], collect_corrections=True)

    min_len_s = min(vac_data['len_s'], corr_data['len_s'])
    min_len_t = min(vac_data['len_t'], corr_data['len_t'])

    dipole_s = dipole_t = energy_s = energy_t = None
    if len(blocks) >= 3:
        rotation = _dipole_rotation(vac_data, corr_data)
        dipole_s = _dipole_scores(vac_data['tdm_singlets'][:min_len_s], corr_data['tdm_singlets'],
                                 rotation, transition=True)
        dipole_t = _dipole_scores(_triplet_dipole_changes(blocks[0], vac_data)[:min_len_t],
                                 _triplet_dipole_changes(blocks[1], corr_data, pcm=True), rotation)
        energy_s = _energy_scores(vac_data['singlets'], corr_data['singlets'])
        energy_t = _energy_scores(vac_data['triplets'], corr_data['triplets'])

    def match_states(orbital_map):
        matches = []
        for key, count, dipoles, energies in (('comp_singlets', min_len_s, dipole_s, energy_s),
                                              ('comp_triplets', min_len_t, dipole_t, energy_t)):
            pcm_comp = corr_data[key]
            if any(orbital_map):
                pcm_comp = _relabel_compositions(pcm_comp, *orbital_map)
            indices, scores, total = _match_compositions(
                vac_data[key][:count], pcm_comp, dipoles,
                energies[:count] if energies is not None else None,
            )
            matches.append((indices, scores, total, pcm_comp))
        return matches

    matches = match_states(({}, {}))
    if len(blocks) >= 3:
        orbital_map = _orbital_relabeling(vac_data, corr_data)
        if any(orbital_map):
            relabeled = match_states(orbital_map)
            baseline_totals = tuple(result[2] for result in matches)
            relabeled_totals = tuple(result[2] for result in relabeled)
            # Apply the shared map only if it helps one spin manifold without hurting the other.
            if relabeled_totals != baseline_totals and all(
                new >= old for new, old in zip(relabeled_totals, baseline_totals)
            ):
                matches = relabeled
    (match_singlets, scores_s, _, comp_s), (match_triplets, scores_t, _, comp_t) = matches

    singlets_vac = vac_data['singlets'][:min_len_s]
    triplets_vac = vac_data['triplets'][:min_len_t]

    singlets_pcm = corr_data['singlets'][match_singlets]
    triplets_pcm = corr_data['triplets'][match_triplets]

    oscs = vac_data['osc_singlets'][:min_len_s]

    s0_vac = vac_data['total_energy'][0]
    s0_pcm = corr_data['total_energy'][0]
    s0_corr = s0_vac - s0_pcm

    ss_s = corr_data['ss_s'][match_singlets]
    ss_t = corr_data['ss_t'][match_triplets]


    y_s = (singlets_vac - singlets_pcm) + s0_corr
    y_t = (triplets_vac - triplets_pcm) + s0_corr

    if return_matching:
        matching = {'two_calculations': len(blocks) >= 3}
        for spin, count, indices, scores, raw_shift in (
            ('s', min_len_s, match_singlets, scores_s, y_s),
            ('t', min_len_t, match_triplets, scores_t, y_t),
        ):
            key = 'comp_singlets' if spin == 's' else 'comp_triplets'
            matching[spin] = {
                'pcm_roots': indices + 1,
                'scores': scores[np.arange(count), indices],
                'matrix': composition_scores(vac_data[key], comp_s if spin == 's' else comp_t),
                'retained': count,
                'raw_shift': raw_shift.copy(),
            }

    y_s[y_s < 0] = 0
    y_t[y_t < 0] = 0

    result = (
        singlets_vac,
        triplets_vac,
        oscs,
        ss_s,
        ss_t,
        s0_corr,
        y_s,
        y_t,
    )
    return result + (matching,) if return_matching else result

#########################################################################################


##GETS SOC BETWEEN Sn STATE AND TRIPLETS#################################################
def pega_soc_singlet(file, n_state):
    socs = []
    n_state += 1
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
    return socs[np.newaxis, :] * 0.12398 / 1000


#########################################################################################


##GETS SOC BETWEEN Tn STATE AND SINGLETS#################################################
def pega_soc_triplet(file, n_state):
    socs = []
    n_state += 1
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
    return socs[np.newaxis, :] * 0.12398 / 1000


#########################################################################################


##GETS SOC BETWEEN Tn STATE AND S0#######################################################
def pega_soc_ground(file, n_state):
    socs = []
    n_state += 1
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
    n_state += 1
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
def avg_socs(files, tipo, n_state, nstates=None):
    col = nstates
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
        if nstates is not None and socs.shape[1] < nstates:
            raise ValueError(
                f"{file}: missing SOC targets for {tipo} root {n_state + 1}; "
                f"expected {nstates}, found {socs.shape[1]}."
            )
        socs = socs[:, :col]
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
            if "Have a nice day" in line:
                break
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
def moment(file, ess, ets, e_s0, dipss, dipts, n_triplet):
    # Conversion factor between a.u. = e*bohr to SI
    conversion = 8.4783533e-30
    ess = np.array(ess)
    ets = np.array(ets)
    ess = np.insert(ess, 0, e_s0)
    # Energies may have been truncated to the states shared by the vacuum
    # and PCM calculations. Restrict both perturbative sums to those states.
    dipss = dipss[:len(ess), :]
    dipts = dipts[:len(ets), :]
    moments = []
    for mqn in ["1", "-1", "0"]:
        socst1 = soc_t1(file, mqn, n_triplet)
        socss0 = soc_s0(file, mqn)
        socss0 = socss0[:len(ets), :]
        socst1 = np.vstack((socss0[n_triplet, :], socst1))
        # Conjugate to get <S0|H|T1>
        socst1[0] = socst1[0].conjugate()
        # Now conjugate socst1
        socst1 = socst1.conjugate()
        socst1 = socst1[:len(ess), :]
        if (len(socss0) != len(ets) or len(dipts) != len(ets)
                or len(socst1) != len(ess) or len(dipss) != len(ess)):
            raise ValueError(
                f"{file}: incomplete SOC or dipole data for phosphorescence "
                f"of T{n_triplet + 1} ({len(ess) - 1} singlets, "
                f"{len(ets)} triplets required)."
            )
        if 0 in ets[n_triplet] - ess:
            return 0
        for i in [0, 1, 2]:
            part_1 = (socss0[:, 0] / (e_s0 - ets)) * dipts[:, i]
            part_1 = np.sum(part_1)
            part_2 = (socst1[:, 0] / (ets[n_triplet] - ess)) * dipss[:, i]
            part_2 = np.sum(part_2)
            complex_dipole = part_1 + part_2
            # append magnitude squared
            moments.append((complex_dipole * complex_dipole.conjugate()).real)

    moments = np.array(moments)
    moments = np.sum(moments) * (conversion**2)
    return moments

def _state_numbers(file, spin):
    """Resolve Q-Chem's combined dipole-table IDs in printed per-spin order."""
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        data = parse_block(log_file.read().split("Have a nice day")[0],
                           collect_corrections=False)
    return np.array(data["indices"])[np.array(data["spins"]) == spin]


def phosph_osc(file, n_state, singlets, triplets, e_s0=0):
    zero = ["0"]
    zero.extend(_state_numbers(file, "Singlet")[:len(singlets)])
    triplet_states = _state_numbers(file, "Triplet")[:n_state]
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
            file, triplet_states, "Electron Dipole Moments of Triplet Excited State", n_triplet
        )
        triplet_triplet_dipoles = pega_dipolos(
            file, triplet_states, "Transition Moments Between Triplet Excited States", n_triplet
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
        )
        total_moments.append(moments)
    total_moments = np.array(total_moments)
    term = E_CHARGE * (HBAR_J**2) / (triplets - e_s0)
    osc_strength = (2 * MASS_E) * total_moments / (3 * term)
    return osc_strength[np.newaxis, :]

##GETS TRANSITION DIPOLE MOMENTS#########################################################
def pega_oscs(files, initial, nstates=None):
    spin = initial[0].upper()
    num = int(initial[1:]) - 1
    mapa = {"S": "Singlet", "T": "Triplet"}
    phrase = "Transition Moments Between " + mapa[spin] + " Excited States"
    rows = []
    for file in files:
        state_numbers = _state_numbers(file, mapa[spin])
        retained = len(state_numbers) if nstates is None else nstates
        if num < 0 or num >= retained or len(state_numbers) < retained:
            raise ValueError(f"{file}: insufficient {mapa[spin]} roots for {initial} "
                             f"with {retained} retained states.")
        source = int(state_numbers[num])
        targets = state_numbers[num + 1:retained]
        strengths = {}
        reading = False
        with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
            for line in log_file:
                if "Have a nice day" in line:
                    break
                if phrase in line:
                    reading = True
                    strengths = {}
                    continue
                fields = line.split()
                if reading and len(fields) >= 6:
                    try:
                        first, second = int(fields[0]), int(fields[1])
                        strength = float(fields[-1])
                    except ValueError:
                        continue
                    if first == source:
                        strengths[second] = strength
                    elif second == source:
                        strengths[first] = strength
                elif reading and strengths and "---" in line:
                    reading = False
        missing = [int(t) for t in targets if int(t) not in strengths]
        if missing:
            raise ValueError(f"{file}: missing oscillator strengths from {initial} "
                             f"to dipole-table state IDs {missing}.")
        rows.append([strengths[int(t)] for t in targets])
    return np.array(rows, dtype=float)


#########################################################################################


##GETS SOCS BETWEEN S0 AND EACH TRIPLET SUBLEVEL#########################################
def soc_s0(file, mqn):
    socs = np.zeros((1))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        read = False
        for line in log_file:
            if "Have a nice day" in line:
                break
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
    return socs * 0.12398 / 1000


#########################################################################################


##GETS SOCS BETWEEN Sm AND EACH Tn SUBLEVEL##############################################
def soc_t1(file, mqn, n_triplet):
    socs = np.zeros((1))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        read = False
        for line in log_file:
            if "Have a nice day" in line:
                break
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
    return socs * 0.12398 / 1000

