#!/usr/bin/env python3
import numpy as np
import os
from nemo.tools import *


##GETS ENERGIES, OSCS, AND INDICES FOR Sn AND Tn STATES##################################
def pega_energias(file):
    with open(file, 'r') as f:
        exc = False
        for line in f:
            if 'TDDFT/TDA Excitation Energies' in line or 'TDDFT Excitation Energies' in line:
                energies, spins, oscs, ind = [], [], [], []
                exc = True
            elif 'Excited state' in line and exc:
                energies.append(float(line.split()[7]))
                ind.append(line.split()[2].replace(':',''))
            elif 'Multiplicity' in line and exc:
                spins.append(line.split()[1])
            elif 'Strength' in line and exc:
                oscs.append(line.split()[2])
            elif '---------------------------------------------------' in line and exc and len(energies) > 0:
                exc = False
                
        singlets   = [energies[i] for i in range(len(energies)) if spins[i] == 'Singlet']
        ind_s      = [ind[i] for i in range(len(ind)) if spins[i] == 'Singlet']
        oscs       = [oscs[i] for i in range(len(energies)) if spins[i] == 'Singlet']
        triplets   = [energies[i] for i in range(len(energies)) if spins[i] == 'Triplet']
        ind_t      = [ind[i] for i in range(len(ind)) if spins[i] == 'Triplet']
        
        oscs       = [x for _, x in sorted(zip(singlets, oscs))]
        ind_s      = [x for _, x in sorted(zip(singlets, ind_s))]
        ind_t      = [x for _, x in sorted(zip(triplets, ind_t))]
        singlets = sorted(singlets)
        triplets = sorted(triplets)
        return singlets, triplets, oscs, ind_s, ind_t
#########################################################################################        

##GETS SOC BETWEEN Sn STATE AND TRIPLETS#################################################
def pega_soc_S(file,n_state):
    socs = []
    with open('Geometries/'+file, 'r') as f:
        catch = False
        for line in f:
            if "Total SOC between the S"+str(n_state)+" state and excited triplet states:" in line:
                catch = True
            elif catch and 'T' in line and '(' not in line: #len(socs) < n_state:
                try:
                    socs.append(float(line.split()[1]))
                except:
                    catch = False
    socs = np.array(socs)
    return socs[np.newaxis,:]*0.12398/1000            
#########################################################################################

##GETS SOC BETWEEN Tn STATE AND SINGLETS#################################################      
def pega_soc_T(file,n_state):
    socs = []
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
    return socs[np.newaxis,:]*0.12398/1000
#########################################################################################

##DECIDES WHICH FUNCTION TO USE UIN ORDER TO GET SOCS####################################
def avg_socs(tipo,n_state):
    files =  [i for i in os.listdir('Geometries') if '.log' in i]    
    files = sorted(files, key=lambda pair: float(pair.split('-')[1]))
    if tipo == 'singlet':
        pega_soc = pega_soc_S
    elif tipo == 'triplet':
        pega_soc = pega_soc_T
    for file in files:
        socs = pega_soc(file,n_state)
        try:
            Socs     = np.vstack((Socs,socs))
        except:
            Socs     = socs        
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
                if line[0] == ind[state]:
                    a = [float(line[i]) for i in range(1,len(line))]
                    a = np.array([a])
                    try:
                        mu = np.vstack((mu,a))
                    except:
                        mu = np.zeros((1,len(line)-1))
                        mu = np.vstack((mu,a))
                elif line[1] == ind[state] and int(line[0]) < int(line[1]):    
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

##GETS SOCS BETWEEN S0 AND EACH TRIPLET SUBLEVEL#########################################
def soc_s0(file,m):
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
    return socs[1:,:]*0.12398/1000
#########################################################################################                    

##GETS SOCS BETWEEN Sm AND EACH Tn SUBLEVEL##############################################
def soc_t1(file,m,n_triplet):
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
    return socs[1:,:]*0.12398/1000
######################################################################################### 

##CALCULATES TRANSITION DIPOLE MOMENTS FOR Tn TO S0 TRANSITIONS##########################
def moment(file,ess,ets,dipss,dipts,n_triplet):
    #Conversion factor between a.u. = e*bohr^-1 to SI
    conversion = 8.4783533-30 
    ess = np.array(ess)
    ets = np.array(ets)
    ess = np.insert(ess, 0, 0)
    Ms = []
    for m in ['1','-1','0']:
        socst1 = soc_t1(file,m,n_triplet)
        socss0 = soc_s0(file,m) 
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
                Ps.append((p1+p2)*conversion)
            Ms.append(Ps[0]**2+Ps[1]**2)    
    
    Ms = np.array(Ms)
    Ms = np.sum(Ms)
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
    files = sorted(files, key=lambda pair: float(pair.split('-')[1]))
    n_state = read_cis(files[0])
    Ms = np.zeros((1,n_state))

    for file in files:
        singlets, triplets, oscs, ind_s, ind_t = pega_energias('Geometries/'+file)            
        zero = ['0']
        zero.extend(ind_s)

        MMs = []
        for n_triplet in range(n_state):
            MS0      = pega_dipolos(file, zero,"Electron Dipole Moments of Ground State",0)            
            MS0resto = pega_dipolos(file, zero,"Transition Moments Between Ground and Singlet Excited States",0) 
            MS0      = np.vstack((MS0,MS0resto))
            MT1      = pega_dipolos(file, ind_t,"Electron Dipole Moments of Triplet Excited State",n_triplet)            
            MT1resto = pega_dipolos(file, ind_t,"Transition Moments Between Triplet Excited States",n_triplet)           
            MT1      = np.vstack((MT1,MT1resto))
            #Fixing the order
            MT1[[0,n_triplet]] = MT1[[n_triplet,0]]
            ms       = moment(file,singlets,triplets,MS0,MT1,n_triplet)  
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
    Ms /= (1/0.529177)*1e10
    term = e*(hbar**2)/Triplets
    Os = (2*mass)*(Ms**2)/(3*term)
    return Os, Singlets, Triplets, Oscs
#########################################################################################


##CALCULATES ISC RATES FROM INITIAL STATE TO SEVERAL STATES OF OPPOSITE SPIN#############
def isc(initial):
    n_state = int(initial[1:]) -1
    kbT = detect_sigma()
    if 's' in initial.lower():
        tipo = 'singlet'
        final = 'T'
    elif 't' in initial.lower():
        tipo = 'triplet'
        final = 'S'    
    _, Singlets, Triplets, _ = analysis()
    socs_complete = avg_socs(tipo,n_state+1)
    try:
        lambdas_list = np.loadtxt('lambdas.txt')
    except:
        fatal_error('No lambdas.txt file found. Reorganization energies are required for this calculation! Goodbye!')
    with open('ISC_rates_{}_.txt'.format(initial), 'w') as f:
        f.write('Intersystem Crossing Rates:\n')
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
            y  = []
            for i in range(np.shape(socs)[0]):
                contribution = (2*np.pi/hbar)*(socs[i]**2)*gauss(delta[i]+lambdas,0,sigma)
                y.append(contribution)
            y = np.array(y)
            N = len(Singlets)
            rate  = np.sum(y)/N 
            #Error estimate
            error = np.sqrt(np.sum((y-rate)**2)/(N*(N-1)))
            f.write('{} -> {}{} : {:5.2e} +/- {:5.2e} s^-1\n'.format(initial.upper(),final,j+1,rate,error))
#########################################################################################    