#!/usr/bin/env python3
import numpy as np
import os
import nemo.tools 

c     = nemo.tools.c
pi    = nemo.tools.pi
hbar  = nemo.tools.hbar
hbar2 = nemo.tools.hbar2
kb    = nemo.tools.kb
e     = nemo.tools.e
mass  = nemo.tools.mass

##RETURNS LIST OF LOG FILES WITH NORMAL TERMINATION######################################
def check_normal(files):
    normal = []
    add = False
    for file in files:
        with open('Geometries/'+file, 'r') as f:
            for line in f:
                if 'TDDFT/TDA Excitation Energies' in line or 'TDDFT Excitation Energies' in line:
                    exc = True
                elif 'Excited state' in line and exc:
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
def pega_energias(file,relaxed=True):
    if relaxed:
        ss = 'Excited-state properties with   relaxed density'
    else:
        ss = 'Excited-state properties with unrelaxed density'    
    with open(file, 'r') as f:
        exc = False
        corr = False
        corrected = []
        for line in f:
            if 'TDDFT/TDA Excitation Energies' in line or 'TDDFT Excitation Energies' in line:
                energies, spins, oscs, ind = [], [], [], []
                exc = True
            elif ss in line:
                corrected = []
                corr = True
            elif 'Excited state' in line and exc:
                energies.append(float(line.split()[7]))
                ind.append(int(line.split()[2].replace(':','')))
            elif 'Multiplicity' in line and exc:
                spins.append(line.split()[1])
            elif 'Strength' in line and exc:
                oscs.append(float(line.split()[2]))
            elif '---------------------------------------------------' in line and exc and len(energies) > 0:
                exc = False
            elif 'Total  1st-order corrected excitation energy' in line and corr:
                corrected.append(float(line.split()[6]))
            elif '------------------------ END OF SUMMARY -----------------------' in line and corr:
                corr = False      

        if len(corrected) > 0:
            energies = corrected

        singlets   = np.array([energies[i] for i in range(len(energies)) if spins[i] == 'Singlet'])
        ind_s      = np.array([ind[i] for i in range(len(ind)) if spins[i] == 'Singlet'])
        oscs       = np.array([oscs[i] for i in range(len(energies)) if spins[i] == 'Singlet'])
        triplets   = np.array([energies[i] for i in range(len(energies)) if spins[i] == 'Triplet'])
        ind_t      = np.array([ind[i] for i in range(len(ind)) if spins[i] == 'Triplet'])
        
        oscs       = np.array([x for _, x in zip(singlets, oscs)])
        ind_s      = np.array([x for _, x in zip(singlets, ind_s)])
        ind_t      = np.array([x for _, x in zip(triplets, ind_t)])

        order_s  = np.argsort(singlets)
        order_t  = np.argsort(triplets)
        singlets = np.sort(singlets)
        triplets = np.sort(triplets)
        oscs     = oscs[order_s]
        ind_s    = ind_s[order_s]
        ind_t    = ind_t[order_t]

        return singlets, triplets, oscs, ind_s, ind_t
#########################################################################################        


##GETS SOC BETWEEN Sn STATE AND TRIPLETS#################################################
def pega_soc_S(file,n_state):
    socs = []
    _, _, _, ind_s, ind_t = pega_energias('Geometries/'+file)
    order_s = np.argsort(ind_s)
    order_t = np.argsort(ind_t)
    n_state = order_s[n_state] + 1
    with open('Geometries/'+file, 'r') as f:
        catch = False
        for line in f:
            if "Total SOC between the S"+str(n_state)+" state and excited triplet states:" in line:
                catch = True
            elif catch and 'T' in line and '(' not in line:
                try:
                    socs.append(float(line.split()[1]))
                except:
                    catch = False
    socs = np.array(socs)
    socs = socs[order_t]
    return socs[np.newaxis,:]*0.12398/1000            
#########################################################################################

##GETS SOC BETWEEN Tn STATE AND SINGLETS#################################################      
def pega_soc_T(file,n_state):
    socs = []
    _, _, _, ind_s, ind_t = pega_energias('Geometries/'+file)
    order_s = np.argsort(ind_s)
    order_t = np.argsort(ind_t)
    n_state = order_t[n_state] + 1
    with open('Geometries/'+file, 'r') as f:
        catch = False
        for line in f:
            if "Total SOC between the S" in line and "state and excited triplet states:" in line:
                catch = True
            elif catch and  'T'+str(n_state)+' ' in line and '(' not in line:
                try:
                    socs.append(float(line.split()[1]))
                except:
                    catch = False
    socs = np.array(socs)
    socs = socs[order_s]
    return socs[np.newaxis,:]*0.12398/1000
#########################################################################################

##DECIDES WHICH FUNCTION TO USE IN ORDER TO GET SOCS#####################################
def avg_socs(tipo,n_state):
    files =  [i for i in os.listdir('Geometries') if '.log' in i]    
    files =  check_normal(files)
    files =  sorted(files, key=lambda pair: float(pair.split('-')[1]))
    if tipo == 'singlet':
        pega_soc = pega_soc_S
    elif tipo == 'triplet':
        pega_soc = pega_soc_T
    for file in files:
        socs = pega_soc(file,n_state)
        try:
            Socs  = np.vstack((Socs,socs))
        except:
            Socs  = socs        
    return Socs  
#########################################################################################

##GETS TRANSITION DIPOLE MOMENTS#########################################################
def pega_dipolos(file, ind,frase, state):
    mu = np.zeros((1,1))
    with open('Geometries/'+file, 'r') as f:
        dip = False
        for line in f:
            if frase in line:
                dip = True
            elif dip and '--' not in line:
                line = line.split()
                if line[0] == str(ind[state]):
                    a = [float(line[i]) for i in range(1,len(line))]
                    a = np.array([a])
                    try:
                        mu = np.vstack((mu,a))
                    except:
                        mu = np.zeros((1,len(line)-1))
                        mu = np.vstack((mu,a))
                elif line[1] == str(ind[state]) and int(line[0]) < int(line[1]):    
                    a = [float(line[0])]
                    b = [float(line[i]) for i in range(2,len(line))]
                    a.extend(b)
                    a = np.array([a])
                    try:
                        mu = np.vstack((mu,a))
                    except:
                        mu = np.zeros((1,len(line)-1))
                        mu = np.vstack((mu,a))
            elif np.shape(mu)[0] > 1 and '---' in line:
                dip = False   
    mu = mu[1:,:]
    muf = np.zeros((1,3))
    ind = np.array(ind).astype(float)
    col = np.shape(mu)[1]
    if col > 4:
        for i in ind:
            index = np.where(mu[:,0] == i)
            try:
                index = index[0][0]
                up = mu[index,1:-1]
                muf = np.vstack((muf,up))
            except:
                pass
        muf = muf[1:,:]
    else:
        muf = mu      
    return muf    
#########################################################################################

##GETS TRANSITION DIPOLE MOMENTS#########################################################
def pega_oscs(file, ind,spin, ind_s):
    frase    = "Transition Moments Between "+spin+" Excited States"
    location = np.where(ind_s == ind)[0][0]
    ind_s = ind_s[location+1:]
    order_s  = np.argsort(ind_s)
    ind = str(ind)
    oscs = []
    with open('Geometries/'+file, 'r') as f:
        dip = False
        for line in f:
            if frase in line:
                dip = True
            elif dip and '--' not in line:
                line = line.split()
                if (line[0] == ind and int(line[1]) in ind_s) or (line[1] == ind and int(line[0]) in ind_s):
                    oscs.append(float(line[5]))
            elif len(oscs) > 0 and '---' in line:
                dip = False 
    oscs = np.array(oscs)
    oscs = oscs[order_s]                         
    return oscs    
#########################################################################################

##GETS SOCS BETWEEN S0 AND EACH TRIPLET SUBLEVEL#########################################
def soc_s0(file,m, ind_t):
    socs = np.zeros((1,2))
    with open('Geometries/'+file, 'r') as f:
        read = False
        for line in f:
            if "SOC between the singlet ground state and excited triplet states (ms="+m in line:
               read = True 
            elif read:
                if "T" in line and "Total" not in line:
                    line = line.split()
                    if line[2] == '-':
                        sign = -1
                    else:
                        sign = 1    
                    c1   = float(line[1].replace('(','').replace(')','').replace('i',''))
                    c2   = sign*float(line[3].replace('(','').replace(')','').replace('i',''))
                    c    = np.array([[c1,c2]])
                    socs = np.vstack((socs,c))
                else:
                    read = False
    socs = socs[1:,:]
    indice = np.argsort(ind_t)
    socs = socs[indice,:]                
    return socs*0.12398/1000
#########################################################################################                    

##GETS SOCS BETWEEN Sm AND EACH Tn SUBLEVEL##############################################
def soc_t1(file,m,n_triplet,ind_s):
    socs = np.zeros((1,2))
    with open('Geometries/'+file, 'r') as f:
        read = False
        for line in f:
            if "SOC between the S" in line and "(ms="+m in line:
               read = True 
            elif read:
                if "T"+str(n_triplet+1)+'(' in line:
                    line = line.split()
                    if line[2] == '-':
                        sign = -1
                    else:
                        sign = 1    
                    c1 = line[1].replace('(','').replace(')','').replace('i','')
                    c2 = sign*line[3].replace('(','').replace(')','').replace('i','')
                    c1 = float(c1.replace('--',''))
                    c2 = float(c2.replace('--',''))
                    c = np.array([[c1,c2]])
                    socs = np.vstack((socs,c))
                    read = False
    socs = socs[1:,:]
    indice = np.argsort(ind_s)
    socs = socs[indice,:]                
    return socs*0.12398/1000
######################################################################################### 

##CALCULATES TRANSITION DIPOLE MOMENTS FOR Tn TO S0 TRANSITIONS##########################
def moment(file,ess,ets,dipss,dipts,n_triplet,ind_s,ind_t):
    #Conversion factor between a.u. = e*bohr to SI
    conversion = 8.4783533E-30 
    fake_t = np.where(np.sort(ind_t) == ind_t[n_triplet])[0][0]
    ess = np.array(ess)
    ets = np.array(ets)
    ess = np.insert(ess, 0, 0)
    Ms = []
    for m in ['1','-1','0']:
        socst1 = soc_t1(file,m,fake_t,ind_s)
        socss0 = soc_s0(file,m,ind_t) 
        socst1 = np.vstack((socss0[0,:],socst1))    
        socst1[:,1] *= -1
        ess = ess[:np.shape(socst1)[0]]
        for i in [0,1,2]:
            Ps = []
            for j in [0,1]:
                p1 = (socss0[:,j]/(0-ets))*dipts[:,i]
                p1 = np.sum(p1)
                p2 = (socst1[:,j]/(ets[n_triplet]-ess))*dipss[:,i]
                p2 = np.sum(p2)
                Ps.append((p1+p2))
            Ms.append(Ps[0]**2+Ps[1]**2)    
    
    Ms = np.array(Ms)
    Ms = np.sum(Ms)*(conversion**2)
    return Ms
#########################################################################################

##READS NUMBER OF EXCITED STATES FROM INPUT FILE#########################################
def read_cis(file):
    file = file[:-3]+'com'
    with open('Geometries/'+file, 'r') as f:
        for line in f:
            if 'cis_n_roots' in line.lower():
                line = line.split()
                for elem in line:
                    try:
                        n_state = int(elem)
                        break
                    except:
                        pass
    return n_state                
#########################################################################################      

##GETS ALL RELEVANT INFORMATION FROM LOG FILES###########################################
def analysis():         
    files =  [i for i in os.listdir('Geometries') if '.log' in i]    
    files = check_normal(files)
    files = sorted(files, key=lambda pair: float(pair.split('-')[1]))
    n_state = read_cis(files[0])
    Ms = np.zeros((1,n_state))

    for file in files:
        singlets, triplets, oscs, ind_s, ind_t = pega_energias('Geometries/'+file)       
        singlets, triplets, oscs, ind_s, ind_t = singlets[:n_state], triplets[:n_state], oscs[:n_state], ind_s[:n_state], ind_t[:n_state]     
        zero = ['0']
        zero.extend(ind_s)

        MMs = []
        MS0      = pega_dipolos(file, zero,"Electron Dipole Moments of Ground State",0)            
        MS0resto = pega_dipolos(file, zero,"Transition Moments Between Ground and Singlet Excited States",0) 
        MS0      = np.vstack((MS0,MS0resto))    
        for n_triplet in range(n_state):
            MT1      = pega_dipolos(file, ind_t,"Electron Dipole Moments of Triplet Excited State",n_triplet)            
            MT1resto = pega_dipolos(file, ind_t,"Transition Moments Between Triplet Excited States",n_triplet)           
            MT1      = np.vstack((MT1,MT1resto))
            #Fixing the order
            order = np.arange(1,n_state)
            order = np.insert(order,n_triplet,0)
            MT1   = MT1[order,:]
            ms    = moment(file,singlets,triplets,MS0,MT1,n_triplet,ind_s,ind_t)  
            MMs.append(ms)
        MMs = np.array(MMs)
        MMs = MMs[np.newaxis,:]
        Ms = np.vstack((Ms,MMs))

        singlets = np.array([singlets[:n_state]])
        triplets = np.array([triplets[:n_state]])
        oscs     = np.array([oscs[:n_state]]).astype(float)
        try:
            Singlets = np.vstack((Singlets,singlets))
            Triplets = np.vstack((Triplets,triplets))
            Oscs     = np.vstack((Oscs,oscs))
        except:
            Singlets = singlets
            Triplets = triplets
            Oscs     = oscs

    Ms = Ms[1:,:]
    term = e*(hbar2**2)/Triplets
    Os = (2*mass)*Ms/(3*term)
    return Os, Singlets, Triplets, Oscs
#########################################################################################


##CALCULATES ISC RATES FROM INITIAL STATE TO SEVERAL STATES OF OPPOSITE SPIN#############
def isc(initial):
    n_state = int(initial[1:]) -1
    kbT = nemo.tools.detect_sigma()
    if 's' in initial.lower():
        tipo = 'singlet'
        final = 'T'
    elif 't' in initial.lower():
        tipo = 'triplet'
        final = 'S'    
    _, Singlets, Triplets, _ = analysis()
    delta_s = np.mean(np.diff(Singlets,axis=1),axis=0)
    delta_t = np.mean(np.diff(Triplets,axis=1),axis=0)
    socs_complete = avg_socs(tipo,n_state)
    try:
        lambdas_list = np.loadtxt('lambdas.txt')
    except:
        nemo.tools.fatal_error('No lambdas.txt file found. Reorganization energies are required for this calculation! Goodbye!')
    with open('ISC_rates_{}_.txt'.format(initial.upper()), 'w') as f:
        f.write('#Intersystem Crossing Rates:\n')
        f.write('#Transition    Rate(s^-1)    Error(s^-1)   AvgGap(eV)  AvgSOC(meV)\n')
        for j in range(np.shape(socs_complete)[1]):
            try:
                lambdas = lambdas_list[j]
            except:
                break
            if tipo == 'singlet':
                delta =  Triplets[:,j] - Singlets[:,n_state]   #Tn (final) - Sm (initial)    
            elif tipo == 'triplet':
                delta = Singlets[:,j] - Triplets[:,n_state]    #Sm (final) - Tn (initial)
            socs = socs_complete[:,j]
            sigma = np.sqrt(2*lambdas*kbT + (kbT)**2)
            y = (2*np.pi/hbar)*(socs[:,np.newaxis]**2)*nemo.tools.gauss(delta[:,np.newaxis]+lambdas,0,sigma)
            N = len(Singlets)
            rate  = np.sum(y)/N 
            #Error estimate
            error = np.sqrt(np.sum((y-rate)**2)/(N*(N-1)))
            gap = np.mean(delta)
            mean_soc = 1000*np.mean(socs)
            f.write('{}->{}{}         {:5.2e}      {:5.2e}      {:+5.3f}      {:5.3f}\n'.format(initial.upper(),final,j+1,rate,error,gap,mean_soc))

        f.write('\n#Transition    AvgGap(eV)    Transition    AvgGap(eV)\n')
        for j in range(len(delta_s)):
            f.write('S{0:}->S{1:}         {2:+5.3f}        T{0:}->T{1:}         {3:+5.3f}\n'.format(j+1,j+2,delta_s[j],delta_t[j]))


    print('Results are written in the {} file'.format('ISC_rates_{}_.txt'.format(initial.upper())))        
#########################################################################################    