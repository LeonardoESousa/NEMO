import sys
import numpy as np
import nemo.parser

##GETS ENERGIES AND OSCS FOR Sn AND Tn STATES##################################
def pega_energias(file):
    with open(file, "r", encoding="utf-8") as log_file:
        fetch_osc = False
        correction, correction2 = [], []
        energies, spins, oscs, double = [], [], [], []
        spin = 'Singlet'
        for line in log_file:
            if "Solving for EOMEE-CCSD A triplet states." in line:
                spin = 'Triplet'
            elif "State A: ccsd: 0/A" in line:
                fetch_osc = True
            elif "Solute Internal Energy" in line:
                sol_int = float(line.split()[5])
            elif "Total Free Energy" in line:
                total_free = float(line.split()[9])
            elif "Excitation energy" in line:
                energies.append(float(line.split()[8]))
                spins.append(spin)
            #elif "R2^2" in line:
            #    double.append(float(line.split()[8]))
            elif "Oscillator strength (a.u.):" in line and fetch_osc:
                oscs.append(float(line.split()[-1]))
                fetch_osc = False
            elif "Total energy in the final basis set" in line:
                line = line.split()
        correction, correction2 = np.zeros(len(energies)), np.zeros(len(energies))      #pega_correction(file, len(energies)/2)
        singlets = np.array(
            [energies[i] for i in range(len(energies)) if spins[i] == "Singlet"]
        )
        #double_s = np.array([double[i] for i in range(len(double)) if spins[i] == "Singlet"])
        ss_s = np.array(
            [
                correction[i] + correction2[i]
                for i in range(len(correction))
                if spins[i] == "Singlet"
            ]
        )
        oscs = np.array(
            [oscs[i] for i in range(len(energies)) if spins[i] == "Singlet"]
        )
        triplets = np.array(
            [energies[i] for i in range(len(energies)) if spins[i] == "Triplet"]
        )
        #double_t = np.array([double[i] for i in range(len(double)) if spins[i] == "Triplet"])
        ss_t = np.array(
            [
                correction[i] + correction2[i]
                for i in range(len(correction))
                if spins[i] == "Triplet"
            ]
        )
        return (
            singlets,
            triplets,
            oscs,
            ss_s,
            ss_t,
            (sol_int - total_free) * 27.2114,
        )

#########################################################################################

##GETS SOC BETWEEN Sn STATE AND TRIPLETS#################################################
def pega_soc_singlet(file, n_state):
    socs = []
    n_state += 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B, catch_C = False, False, False
        for line in log_file:
            if "State A:" in line:
                catch_A = "eomee_ccsd/rhfref/singlets:" in line and line.split()[-1] == f'{n_state}/A'
                catch_B, catch_C = False, False
            elif "State B: eomee_ccsd/rhfref/triplets:" in line:                    
                catch_B = True
            elif catch_A and catch_B and "Arithmetically averaged transition SO matrices" in line:
                catch_C = True
            elif catch_A and catch_B and catch_C and "SOCC = " in line:
                socs.append(float(line.split()[2]))
                catch_A, catch_B, catch_C = False, False, False
    socs = np.array(socs)
    return socs[np.newaxis, :] * 0.12398 / 1000

#########################################################################################


##GETS SOC BETWEEN Tn STATE AND SINGLETS#################################################
def pega_soc_triplet(file, n_state):
    socs = []
    n_state += 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B, catch_C = False, False, False
        for line in log_file:
            if "State A:" in line:
                catch_B = "eomee_ccsd/rhfref/singlets:" in line
                catch_A, catch_C = False, False
            elif "State B:" in line:
                catch_A = "eomee_ccsd/rhfref/triplets:" in line and line.split()[-1] == f'{n_state}/A'
            elif catch_A and catch_B and "Arithmetically averaged transition SO matrices" in line:
                catch_C = True
            elif catch_A and catch_B and catch_C and "SOCC = " in line:
                socs.append(float(line.split()[2]))
                catch_A, catch_B, catch_C = False, False, False
    socs = np.array(socs)
    return socs[np.newaxis, :] * 0.12398 / 1000


#########################################################################################

##GETS SOC BETWEEN Tn STATE AND S0#######################################################
def pega_soc_ground(file, n_state):
    socs = []
    n_state += 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B, catch_C = False, False, False
        for line in log_file:
            if "State A:" in line:
                catch_A = "State A: ccsd: 0/A" in line
                catch_B, catch_C = False, False
            elif "State B:" in line:
                catch_B = "eomee_ccsd/rhfref/triplets:" in line and line.split()[-1] == f'{n_state}/A'
            elif catch_A and catch_B and "Arithmetically averaged transition SO matrices" in line:
                catch_C = True
            elif catch_A and catch_B and catch_C and "SOCC = " in line:
                socs.append(float(line.split()[2]))
                catch_A, catch_B, catch_C = False, False, False
    socs = np.array(socs)
    return socs[np.newaxis, :] * 0.12398 / 1000

#########################################################################################

def pega_soc_triplet_triplet(file, n_state):
    #temporary fix since triplet-triplet socs are not used yet 
    triplets = pega_energias("Geometries/" + file)[1]
    return np.zeros((1, len(triplets)-1))

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
def pega_dipole_ground(file):
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch = False
        for line in log_file:
            if 'Dipole Moment (Debye)' in line:
                catch = True
            elif catch:
                line = line.split()
                dipole = [float(line[1]),float(line[3]),float(line[5])]
                catch = False
                break    
    return np.array(dipole)[np.newaxis, :] / 2.541746

def _eom_dipoles(file, n_singlets, n_triplets):
    """Build Hermitian effective dipoles from the two EOM transition directions."""
    dipoles = {
        "singlets": np.zeros((n_singlets + 1, n_singlets + 1, 3)),
        "triplets": np.zeros((n_triplets, n_triplets, 3)),
    }
    dipoles["singlets"][0, 0] = pega_dipole_ground(file)[0]
    spin, root = "singlets", 1
    first = second = None
    reading_transition = False
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        for line in log_file:
            if "Solving for EOMEE-CCSD A triplet states." in line:
                spin, root = "triplets", 0
            elif "Dipole moment (a.u.):" in line:
                if root < len(dipoles[spin]):
                    dipoles[spin][root, root] = _eom_dipole_components(line)
                root += 1
            elif "State A:" in line:
                first = _eom_dipole_state(line)
                second, reading_transition = None, False
            elif "State B:" in line:
                second = _eom_dipole_state(line)
            elif "Transition dipole moment (a.u.):" in line:
                reading_transition = first is not None and second is not None
            elif reading_transition and "A->B:" in line:
                forward = _eom_dipole_components(line)
            elif reading_transition and "B->A:" in line:
                first_spin, i = first
                second_spin, j = second
                if first_spin == second_spin and max(i, j) < len(dipoles[first_spin]):
                    moment = (forward + _eom_dipole_components(line)) / 2
                    dipoles[first_spin][i, j] = dipoles[first_spin][j, i] = moment
                reading_transition = False
    return dipoles["singlets"], dipoles["triplets"]


def _eom_dipole_components(line):
    coordinates = line.rsplit("(X", 1)[1].replace("Y", "").replace("Z", "").rstrip().rstrip(")")
    return np.array([float(value) for value in coordinates.split(",")])


def _eom_dipole_state(line):
    if "ccsd: 0/" in line:
        return "singlets", 0
    marker = "eomee_ccsd/rhfref/"
    if marker + "singlets:" in line or marker + "triplets:" in line:
        spin, label = line.split(marker, 1)[1].split(":")
        root = int(label.split("/")[0])
        return spin, root if spin == "singlets" else root - 1
    return None


def phosph_osc(file, n_state, singlets, triplets):
    if not len(triplets):
        return np.empty((1, 0))
    dipoles_s, dipoles_t = _eom_dipoles(file, len(singlets), len(triplets))
    couplings = nemo.parser.singlet_triplet_soc(file, len(singlets), len(triplets), soc_s0, soc_t1)
    return nemo.parser.spin_orbit_oscillators(singlets, triplets, couplings, dipoles_s, dipoles_t)


##GETS TRANSITION DIPOLE MOMENTS#########################################################
def pega_oscs(files, initial, nstates=None):
    spin = initial[0].upper()
    source = int(initial[1:])
    mapa = {"S": "singlets", "T": "triplets"}
    if nstates is None:
        energies = pega_energias("Geometries/" + files[0])
        nstates = len(energies[0 if spin == "S" else 1])
    targets = range(source + 1, nstates + 1)
    rows = []
    for file in files:
        strengths = {}
        first, second = None, None
        with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
            for line in log_file:
                if "State A:" in line:
                    first = int(line.split()[-1].split("/")[0]) if f"eomee_ccsd/rhfref/{mapa[spin]}:" in line else None
                    second = None
                elif "State B:" in line:
                    second = int(line.split()[-1].split("/")[0]) if f"eomee_ccsd/rhfref/{mapa[spin]}:" in line else None
                elif "Oscillator strength (a.u.):" in line:
                    if first == source and second in targets:
                        strengths[second] = float(line.split()[-1])
                    elif second == source and first in targets:
                        strengths[first] = float(line.split()[-1])
                elif '-----------------------' in line:
                    first, second = None, None
        missing = [root for root in targets if root not in strengths]
        if missing:
            raise ValueError(f"{file}: missing {initial.upper()} oscillator strengths to roots {missing}")
        rows.append([strengths[root] for root in targets])
    return np.array(rows, dtype=float)


#########################################################################################

##GETS SOCS BETWEEN S0 AND EACH TRIPLET SUBLEVEL#########################################
def soc_s0(file, mqn):
    if mqn == '-1':
        mqn = '(L-)'
    elif mqn == '0':
        mqn = '(L0)'
    else:
        mqn = '(L+)'        
    socs = np.zeros((1))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B = False, False
        for line in log_file:
            if "State A: ccsd: 0/A" in line:
                catch_A = True
                catch_B = False
            elif catch_A and "State B: eomee_ccsd/rhfref/triplets:" in line:# and f'{n_state}/A' in line:                    
                catch_B = True
            elif "--------------------------------------" in line:
                catch_A, catch_B = False, False
            #elif catch_A and "State B: eomee_ccsd/rhfref/triplets:" in line and f'{n_state}/A' not in line:                    
            #    catch_A = False
            elif catch_A and catch_B and 'Hso'+mqn in line:
                line = line.split()
                number = line[-1].replace('(','').replace(')','')
                real_part = float(number.split(',')[0])
                img_part = float(number.split(',')[1])
                complex_soc = real_part + img_part * 1j
                complex_soc = complex_soc.conjugate()
                socs = np.vstack((socs, np.array([complex_soc])))
                catch_A, catch_B = False, False
    socs = socs[1:, :]
    return socs * 0.12398 / 1000


#########################################################################################


##GETS SOCS BETWEEN Sm AND EACH Tn SUBLEVEL##############################################
def soc_t1(file, mqn, n_state):
    if mqn == '-1':
        mqn = '(L-)'
    elif mqn == '0':
        mqn = '(L0)'
    else:
        mqn = '(L+)' 
    socs = np.zeros((1))
    n_state += 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B = False, False
        for line in log_file:
            if "State B: eomee_ccsd/rhfref/triplets:" in line and line.split()[-1] == f'{n_state}/A':
                catch_A = True
            elif "State A: eomee_ccsd/rhfref/singlets:" in line:                    
                catch_B = True
                catch_A = False
            elif "--------------------------------------" in line:
                catch_A, catch_B = False, False    
            elif catch_A and catch_B and 'Hso'+mqn in line:
                line = line.split()
                number = line[-1].replace('(','').replace(')','')
                real_part = float(number.split(',')[0])
                img_part = float(number.split(',')[1])
                complex_soc = real_part + img_part * 1j
                complex_soc = complex_soc.conjugate()
                socs = np.vstack((socs, np.array([complex_soc])))
                catch_A, catch_B = False, False
    socs = socs[1:, :]
    return socs * 0.12398 / 1000

