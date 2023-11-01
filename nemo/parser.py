import sys
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


def pega_modos(file):
    with open(file, "r", encoding="utf-8") as freq_file:
        start = False
        coords = []
        count = 0
        for line in freq_file:
            if "X      Y      Z" in line:
                start = True
                continue
            if start:
                if "TransDip" in line:
                    start = False
                    line = line.split()
                    count += int((len(line) - 1) / 3)
                    try:
                        arranged_coords = np.hstack((arranged_coords, np.array(coords)))
                    except NameError:
                        arranged_coords = np.array(coords)
                    coords = []
                else:
                    line = line.split()
                    buffer = []
                    for i in range(1, len(line)):
                        buffer.append(float(line[i]))
                    coords.append(buffer)
    final_coords = np.zeros(
        (arranged_coords.shape[0] * 3, arranged_coords.shape[1] // 3)
    )
    for i in range(0, arranged_coords.shape[1] // 3):
        final_coords[:, i] = arranged_coords[:, 3 * i : 3 * i + 3].flatten()
    return final_coords


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
