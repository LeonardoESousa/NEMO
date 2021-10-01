#!/usr/bin/env python3
import numpy as np
import os
import sys


n_state = int(sys.argv[1])

#n_triplet = int(sys.argv[2])

def pega_energias(file):
    with open(file, 'r') as f:
        energies, spins, corrected, oscs, ind = [], [], [], [], []
        exc, pcm, SO = False, True, False
        s0t1, s1t1 = False, False
        for line in f:
            if 'TDDFT/TDA Excitation Energies' in line:
                exc = True
            elif 'Excited state' in line and exc:
                energies.append(float(line.split()[7]))
                ind.append(line.split()[2].replace(':',''))
                corrected.append(float(line.split()[7]))
            elif 'Multiplicity' in line and exc:
                spins.append(line.split()[1])
            elif 'Strength' in line and exc:
                oscs.append(line.split()[2])
            elif '---------------------------------------------------' in line and exc and len(energies) > 0:
                exc = False
            #elif 'SUMMARY OF LR-PCM AND SS-PCM' in line and len(energies) > 0:
            #    pcm = True
            #elif '1st-order LR-PCM corrected excitation energy' in line and pcm: #'Total  1st-order corrected excitation energy'
            #    corrected.append(float(line.split()[6]))
            #elif '------------------------ END OF SUMMARY -----------------------' in line and pcm and len(corrected) > 0:
                #pcm = False
            #    break
                
        singlets   = [corrected[i] for i in range(len(corrected)) if spins[i] == 'Singlet']
        ind_s      = [ind[i] for i in range(len(ind)) if spins[i] == 'Singlet']
        oscs       = [oscs[i] for i in range(len(corrected)) if spins[i] == 'Singlet']
        n_singlets = [i for i in range(1,len(singlets)+1)]
        triplets   = [corrected[i] for i in range(len(corrected)) if spins[i] == 'Triplet']
        ind_t      = [ind[i] for i in range(len(ind)) if spins[i] == 'Triplet']
        n_triplets = [i for i in range(1,len(triplets)+1)]                        
        
        n_singlets = [x for y, x in sorted(zip(singlets, n_singlets))]
        oscs       = [x for y, x in sorted(zip(singlets, oscs))]
        ind_s      = [x for y, x in sorted(zip(singlets, ind_s))]
        n_triplets = [x for y, x in sorted(zip(triplets, n_triplets))]
        ind_t      = [x for y, x in sorted(zip(triplets, ind_t))]
        singlets = sorted(singlets)
        triplets = sorted(triplets)
        
        return singlets, triplets, oscs, ind_s, ind_t

def pega_dipolos(file, ind,frase, state):
    mu = np.zeros((1,1))
    with open(file, 'r') as f:
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
                    #print(np.shape(a))
                    #print(a)
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

def pega_soc(file):
    soc = np.nan
    socs = []
    with open(file, 'r') as f:
        catch = False
        for line in f:
            if "Total SOC between the S1 state and excited triplet states:" in line:
                catch = True
            elif catch and len(socs) < n_state:# and line.split()[0] == 'T1':
                socs.append(float(line.split()[1]))
                #catch = False
            elif len(socs) == n_state:
                catch = False
    socs = np.array(socs)
    return socs[np.newaxis,:]            
        
def soc_s0(file,m):
    socs = np.zeros((1,2))
    with open(file, 'r') as f:
        read = False
        for line in f:
            if "SOC between the singlet ground state and excited triplet states (ms="+m in line:
               read = True 
            elif read:
                if "T" in line and "Total" not in line:
                    line = line.split()
                    c1 = float(line[1].replace('(','').replace(')','').replace('i',''))
                    c2 = float(line[3].replace('(','').replace(')','').replace('i',''))
                    c = np.array([[c1,c2]])
                    socs = np.vstack((socs,c))
                else:
                    read = False
    return socs[1:,:]*0.12398/1000
                    
def soc_t1(file,m,n_triplet):
    socs = np.zeros((1,2))
    with open(file, 'r') as f:
        read = False
        for line in f:
            if "SOC between the S" in line and "(ms="+m in line:
               read = True 
            elif read:
                if "T"+str(n_triplet+1)+'(' in line:
                    line = line.split()
                    c1 = line[1].replace('(','').replace(')','').replace('i','')
                    c2 = line[3].replace('(','').replace(')','').replace('i','')
                    c1 = float(c1.replace('--',''))
                    c2 = float(c2.replace('--',''))
                    c = np.array([[c1,c2]])
                    socs = np.vstack((socs,c))
                    read = False
    return socs[1:,:]*0.12398/1000
        
def moment(file,ess,ets,dipss,dipts,n_triplet):
    conversion = 8.478353*10**(-30)
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
                #print(i,p1,p2)
                Ps.append((p1+p2)*conversion)
            Ms.append(Ps[0]**2+Ps[1]**2)    
    
    Ms = np.array(Ms)
    Ms = np.sum(Ms)
    return Ms

file = 'soc.out'
         
folders =  [x[0] for x in os.walk('.') if 'geom' in x[0] and x[0].count('/') < 2 and 'txt' not in x[0]]    
folders = sorted(folders, key=lambda pair: float(pair.split('_')[1]))

import sys
Ms = np.zeros((1,2))


for folder in folders:
    #print(folder)
    os.chdir(folder) 
    singlets, triplets, oscs, ind_s, ind_t = pega_energias(file)            
    zero = ['0']
    zero.extend(ind_s)

    MMs = []
    for n_triplet in range(2):
        MS0 = pega_dipolos(file, zero,"Electron Dipole Moments of Ground State",0)            
        MS0resto = pega_dipolos(file, zero,"Transition Moments Between Ground and Singlet Excited States",0) 
        MS0 = np.vstack((MS0,MS0resto))
        MT1 = pega_dipolos(file, ind_t,"Electron Dipole Moments of Triplet Excited State",n_triplet)            
        MT1resto = pega_dipolos(file, ind_t,"Transition Moments Between Triplet Excited States",n_triplet)           
        MT1 = np.vstack((MT1,MT1resto))
        MT1[[0,n_triplet]] = MT1[[n_triplet,0]]
        ms = moment(file,singlets,triplets,MS0,MT1,n_triplet)  
        MMs.append(ms)
    MMs = np.array(MMs)
    MMs = MMs[np.newaxis,:]
    Ms = np.vstack((Ms,MMs))
    #Ms.append(ms)
    
    socs = pega_soc(file)
    
    singlets = np.array([singlets[:n_state]])
    triplets = np.array([triplets[:n_state]])
    oscs     = np.array([oscs[:n_state]]).astype(float)
    try:
        Singlets = np.vstack((Singlets,singlets))
        Triplets = np.vstack((Triplets,triplets))
        Oscs     = np.vstack((Oscs,oscs))
        Socs     = np.vstack((Socs,socs))
    except:
        print('Entrei', folder)
        #print(np.shape(MT1))
        #print(MT1)
        #print(np.shape(MS0))
        #print(np.shape(MT1))
        print(np.shape(singlets))
        print(np.shape(triplets))
        print(np.shape(oscs))
        print(np.shape(socs))
        Singlets = singlets
        Triplets = triplets
        Oscs     = oscs
        Socs     = socs
    os.chdir('..')

Ms = Ms[1:,:]
np.savetxt('RTP.txt', Ms , delimiter='\t', fmt="%5.7e")     

engs = np.hstack((Singlets,Triplets))
np.savetxt("ENGs.txt", engs, delimiter='\t', fmt='%6.3f' )     
    
np.savetxt("OSCs.txt", Oscs, delimiter='\t', fmt='%6.3f' )    

np.savetxt("SOCs.txt", Socs, delimiter='\t', fmt='%6.3f' )
     
print(np.shape(engs))     

