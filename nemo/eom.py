import sys
import numpy as np
import nemo.parser

##SOME CONSTANTS##############################################
HBAR_J      = nemo.parser.HBAR_J # J s
MASS_E      = nemo.parser.MASS_E # kg
E_CHARGE    = nemo.parser.E_CHARGE # C
###############################################################

##GETS ENERGIES, OSCS, AND INDICES FOR Sn AND Tn STATES##################################
def pega_energias(file):
    with open(file, "r", encoding="utf-8") as log_file:
        fetch_osc = False
        correction, correction2 = [], []
        energies, spins, oscs, ind, double = [], [], [], [], []
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
                if spin == 'Singlet':
                    num = spins.count('Singlet')
                else:
                    num = spins.count('Triplet')   
                ind.append(num)
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
        ind_s = np.array([ind[i] for i in range(len(ind)) if spins[i] == "Singlet"])
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
        #double_s = double_s[order_s]
        #double_t = double_t[order_t]
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
def pega_soc_singlet(file, n_state, ind_s, ind_t):
    socs = []
    order_s = np.argsort(ind_s)
    order_t = np.argsort(ind_t)
    n_state = order_s[n_state] + 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B, catch_C = False, False, False
        for line in log_file:
            if "State A: eomee_ccsd/rhfref/singlets:" in line and f'{n_state}/A' in line:
                catch_A = True
            elif "State B: eomee_ccsd/rhfref/triplets:" in line:                    
                catch_B = True
            elif catch_A and catch_B and "Arithmetically averaged transition SO matrices" in line:
                catch_C = True
            elif catch_A and catch_B and catch_C and "SOCC = " in line:
                socs.append(float(line.split()[2]))
                catch_A, catch_B, catch_C = False, False, False
    socs = np.array(socs)
    socs = socs[order_t]
    return socs[np.newaxis, :] * 0.12398 / 1000

#########################################################################################


##GETS SOC BETWEEN Tn STATE AND SINGLETS#################################################
def pega_soc_triplet(file, n_state, ind_s, ind_t):
    socs = []
    order_s = np.argsort(ind_s)
    order_t = np.argsort(ind_t)
    n_state = order_s[n_state] + 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B, catch_C = False, False, False
        for line in log_file:
            if "State B: eomee_ccsd/rhfref/triplets:" in line and f'{n_state}/A' in line:
                catch_A = True
            elif "State A: eomee_ccsd/rhfref/singlets:" in line:                    
                catch_B = True
                catch_A = False
            elif catch_A and catch_B and "Arithmetically averaged transition SO matrices" in line:
                catch_C = True
            elif catch_A and catch_B and catch_C and "SOCC = " in line:
                socs.append(float(line.split()[2]))
                catch_A, catch_B, catch_C = False, False, False
    socs = np.array(socs)
    socs = socs[order_t]
    return socs[np.newaxis, :] * 0.12398 / 1000


#########################################################################################

##GETS SOC BETWEEN Tn STATE AND S0#######################################################
def pega_soc_ground(file, n_state, ind_s, ind_t):
    socs = []
    #_, _, _, ind_s, ind_t, _, _, _ = pega_energias("Geometries/" + file)
    #order_s = np.argsort(ind_s)
    #order_t = np.argsort(ind_t)
    #n_state = order_s[n_state] + 1
    n_state += 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B, catch_C = False, False, False
        for line in log_file:
            if "State A: ccsd: 0/A" in line:
                catch_A = True
                catch_B = False
            elif catch_A and "State B: eomee_ccsd/rhfref/triplets:" in line and f'{n_state}/A' in line:                    
                catch_B = True
            elif catch_A and "State B: eomee_ccsd/rhfref/triplets:" in line and f'{n_state}/A' not in line:                    
                catch_A = False
            elif catch_A and catch_B and "Arithmetically averaged transition SO matrices" in line:
                catch_C = True
            elif catch_A and catch_B and catch_C and "SOCC = " in line:
                socs.append(float(line.split()[2]))
                catch_A, catch_B, catch_C = False, False, False
    socs = np.array(socs)
    #socs = socs[order_t]
    return socs[np.newaxis, :] * 0.12398 / 1000

#########################################################################################

def pega_soc_triplet_triplet(file, n_state, ind_s, ind_t):
    #temporary fix since triplet-triplet socs are not used yet 
    return np.zeros((1, len(ind_t)-1))

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
def moment(file, ess, ets, dipss, dipts, n_triplet, ind_s, ind_t):
    # Conversion factor between a.u. = e*bohr to SI
    conversion = 8.4783533e-30
    fake_t = np.where(np.sort(ind_t) == ind_t[n_triplet])[0][0]
    ess = np.array(ess)
    ets = np.array(ets)
    ess = np.insert(ess, 0, 0)
    moments = []
    for mqn in ["1", "-1", "0"]:
        socst1 = soc_t1(file, mqn, fake_t, ind_s, ind_t)
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
    return np.array(dipole)[np.newaxis, :] 

def pega_dipole_ground_singlet(file):
    dipole = np.zeros((1, 3))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B = False, False
        for line in log_file:
            if 'State A: ccsd: 0/A' in line:
                catch_A = True
            elif catch_A and "State B: eomee_ccsd/rhfref/singlets:" in line:   
                catch_B = True
            elif '-----------------------------------------------------' in line:
                catch_A, catch_B = False, False    
            elif catch_A and catch_B and 'A->B:' in line:
                line = line.split()
                dipole = np.vstack((dipole,[float(line[3].replace(',','')),float(line[5].replace(',','')),float(line[7].replace(')',''))]))
                catch_A, catch_B = False, False
    dipole = dipole[1:, :]
    return dipole

def pega_dipole_triplets(file, ind, state):
    dipole = np.zeros((1, 3))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A = False
        for line in log_file:
            if 'Solving for EOMEE-CCSD A triplet states.' in line:
                catch_A = True
            elif catch_A and "Dipole moment (a.u.)" in line:   
                line = line.split()
                dipole = np.vstack((dipole,[float(line[5].replace(',','')),float(line[7].replace(',','')),float(line[9].replace(')',''))]))
                
    dipole = dipole[1:, :]
    dipole = dipole[ind[state], :]
    return np.array(dipole)[np.newaxis, :] 

def pega_dipole_triplet_triplet(file,ind,state):
    dipole = np.zeros((1, 4))
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch = False
        states = []
        for line in log_file:
            if '----------------------------------------' in line:
                states = []
                catch = False
            if 'State A: eomee_ccsd/rhfref/triplets:' in line and len(states) == 0:
                line = line.split()
                state_num = int(line[-1].replace('/A',''))
                states.append(state_num)
            elif 'State B: eomee_ccsd/rhfref/triplets:' in line and len(states) == 1:
                line = line.split()
                state_num = int(line[-1].replace('/A',''))
                states.append(state_num)
                if ind[state]+1 in states:
                    catch = True    
            elif catch and "A->B:" in line:
                if states[0] == ind[state]+1:
                    x = states[1]
                else:
                    x = states[0]       
                line = line.split()
                dipole = np.vstack((dipole,[x,float(line[3].replace(',','')),float(line[5].replace(',','')),float(line[7].replace(')',''))]))
                catch = False
                states = []
    dipole = dipole[1:, :]
    #sort by first column
    dipole = dipole[dipole[:,0].argsort()]
    return dipole[:,1:]

def phosph_osc(file, n_state, ind_s, ind_t, singlets, triplets):
    zero = ["0"]
    zero.extend(ind_s)
    total_moments = []
    ground_dipoles = pega_dipole_ground(file)
    ground_singlet_dipoles = pega_dipole_ground_singlet(file)  
    ground_dipoles = np.vstack((ground_dipoles, ground_singlet_dipoles))
    for n_triplet in range(n_state):
        triplet_dipoles = pega_dipole_triplets(file,ind_t,n_triplet)
        triplet_triplet_dipoles = pega_dipole_triplet_triplet(file,ind_t,n_triplet)
        triplet_dipoles = np.vstack((triplet_dipoles, triplet_triplet_dipoles))
        # Fixing the order
        order = np.arange(1, n_state)
        order = np.insert(order, n_triplet, 0)
        triplet_dipoles = triplet_dipoles[order, :]
        moments = moment(
            file,
            singlets,
            triplets,
            ground_dipoles,
            triplet_dipoles,
            n_triplet,
            ind_s,
            ind_t,
        )
        total_moments.append(moments)
    total_moments = np.array(total_moments)
    term = E_CHARGE * (HBAR_J**2) / triplets
    #print(total_moments.shape,term.shape)
    osc_strength = (2 * MASS_E) * total_moments / (3 * term)
    return osc_strength[np.newaxis, :]

##GETS TRANSITION DIPOLE MOMENTS#########################################################
def pega_oscs(files, indices, initial):
    spin = initial[0].upper()
    num = int(initial[1:]) - 1
    mapa = {"S": "singlets", "T": "triplets"}
    states = []
    for i, file in enumerate(files):
        oscs = []
        #ind = indices[i, num]
        #ind_s = indices[i, :]
        #location = np.where(ind_s == ind)[0][0]
        #ind_s = ind_s[location + 1 :]
        #ind = str(ind)
        with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
            #check_A, check_B = False, False
            for line in log_file:
                if f'State A: eomee_ccsd/rhfref/{mapa[spin]}:' in line:
                    line = line.split()
                    state_num = int(line[-1].replace('/A',''))
                    states.append(state_num)
                elif f'State B: eomee_ccsd/rhfref/{mapa[spin]}:' in line:
                    #check_B = True
                    line = line.split()
                    state_num = int(line[-1].replace('/A',''))
                    states.append(state_num)  
                elif '-----------------------' in line:
                    states = []
                elif len(states) == 2 and "Oscillator strength (a.u.):" in line:
                    if num+1 in states:
                        if states[0] == num+1:
                            x = states[1]
                        else:
                            x = states[0] 
                        oscs.append(float(line.split()[3]))       
                        #total_oscs = np.vstack((total_oscs, [x,float(line.split()[3])]))
                    #check_A, check_B = False, False      
                    states = []
        try:        
            total_oscs = np.vstack((total_oscs, oscs))
        except NameError:
            total_oscs = np.array(oscs)[np.newaxis,:]            
    # sort by first column
    #total_oscs = total_oscs[total_oscs[:, 0].argsort()]                
    
    return total_oscs


#########################################################################################

##GETS SOCS BETWEEN S0 AND EACH TRIPLET SUBLEVEL#########################################
def soc_s0(file, mqn, ind_t):
    if mqn == '-1':
        mqn = '(L-)'
    elif mqn == '0':
        mqn = '(L0)'
    else:
        mqn = '(L+)'        
    socs = np.zeros((1))
    #n_state += 1  # order_t[n_state] + 1
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
    indice = np.argsort(ind_t)
    socs = socs[indice, :]
    return socs * 0.12398 / 1000


#########################################################################################


##GETS SOCS BETWEEN Sm AND EACH Tn SUBLEVEL##############################################
def soc_t1(file, mqn, n_state, ind_s, ind_t):
    if mqn == '-1':
        mqn = '(L-)'
    elif mqn == '0':
        mqn = '(L0)'
    else:
        mqn = '(L+)' 
    socs = np.zeros((1))
    order_s = np.argsort(ind_s)
    order_t = np.argsort(ind_t)
    n_state = order_s[n_state] + 1
    with open("Geometries/" + file, "r", encoding="utf-8") as log_file:
        catch_A, catch_B = False, False
        for line in log_file:
            if "State B: eomee_ccsd/rhfref/triplets:" in line and f'{n_state}/A' in line:
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
    socs = socs[order_t]
    return socs * 0.12398 / 1000

