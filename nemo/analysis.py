#!/usr/bin/env python3
import numpy as np
import os
import nemo.tools 
import warnings
import pandas as pd

c     = nemo.tools.c
pi    = nemo.tools.pi
hbar  = nemo.tools.hbar
hbar2 = nemo.tools.hbar2
kb    = nemo.tools.kb
e     = nemo.tools.e
mass  = nemo.tools.mass
epsilon0 = nemo.tools.epsilon0


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
def pega_energias(file):
    ss = 'Excited-state properties with   relaxed density'
    with open(file, 'r') as f:
        exc = False
        corr = False
        correction = []
        for line in f:
            if 'TDDFT/TDA Excitation Energies' in line or 'TDDFT Excitation Energies' in line:
                energies, spins, oscs, ind = [], [], [], []
                exc = True
            elif ss in line:
                corr = True
            elif 'Solute Internal Energy' in line:
                sol_int = float(line.split()[5])
            elif 'Total Free Energy' in line:  
                total_free = float(line.split()[9])
            elif 'Excited state' in line and exc:
                energies.append(float(line.split()[7]))
                ind.append(int(line.split()[2].replace(':','')))
            elif 'Multiplicity' in line and exc:
                spins.append(line.split()[1])
            elif 'Strength' in line and exc:
                oscs.append(float(line.split()[2]))
            elif '---------------------------------------------------' in line and exc and len(energies) > 0:
                exc = False
            elif 'SS-PCM correction' in line and corr:
                correction.append(-1*float(line.split()[3]))
            elif '------------------------ END OF SUMMARY -----------------------' in line and corr:
                corr = False      
            elif 'Total energy in the final basis set' in line:
                line = line.split()
                total_nopcm =  float(line[8])    
        if len(correction) == 0: #When run on logs that do not employ pcm
            correction = np.zeros(len(energies))
            sol_int = total_nopcm
            total_free = total_nopcm
        singlets   = np.array([energies[i] for i in range(len(energies)) if spins[i] == 'Singlet'])
        ss_s       = np.array([correction[i]   for i in range(len(correction))   if spins[i] == 'Singlet'])
        ind_s      = np.array([ind[i] for i in range(len(ind)) if spins[i] == 'Singlet'])
        oscs       = np.array([oscs[i] for i in range(len(energies)) if spins[i] == 'Singlet'])
        triplets   = np.array([energies[i] for i in range(len(energies)) if spins[i] == 'Triplet'])
        ss_t       = np.array([correction[i]   for i in range(len(correction))   if spins[i] == 'Triplet'])
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

        return singlets, triplets, oscs, ind_s, ind_t, ss_s, ss_t, (sol_int - total_free)*27.2114
#########################################################################################        


##GETS SOC BETWEEN Sn STATE AND TRIPLETS#################################################
def pega_soc_S(file,n_state):
    socs = []
    _, _, _, ind_s, ind_t, _, _, _ = pega_energias('Geometries/'+file)
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
    _, _, _, ind_s, ind_t, _, _, _ = pega_energias('Geometries/'+file)
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

##GETS SOC BETWEEN Tn STATE AND S0#######################################################      
def pega_soc_ground(file,n_state):
    socs = []
    #_, _, _, ind_s, ind_t, _, _, _ = pega_energias('Geometries/'+file)
    #order_s = np.argsort(ind_s)
    #order_t = np.argsort(ind_t)
    n_state += 1# order_t[n_state] + 1
    with open('Geometries/'+file, 'r') as f:
        catch = False
        for line in f:
            if "Total SOC between the singlet ground state and excited triplet states:" in line:
                catch = True
            elif catch and  'T'+str(n_state)+' ' in line and '(' not in line:
                try:
                    socs.append(float(line.split()[1]))
                except:
                    catch = False
            elif len(line.split()) < 2:
                    catch = False        
    socs = np.array(socs)
    #socs = socs[order_s]
    return socs[np.newaxis,:]*0.12398/1000
#########################################################################################

##GETS SOC BETWEEN Tn STATE AND SINGLETS#################################################      
def pega_soc_TT(file,n_state):
    socs = []
    #_, _, _, ind_s, ind_t, _, _, _ = pega_energias('Geometries/'+file)
    #order_s = np.argsort(ind_s)
    #order_t = np.argsort(ind_t)
    n_state += 1 #order_t[n_state] + 1
    with open('Geometries/'+file, 'r') as f:
        catch, catch2 = False, False
        for line in f:
            if "Total SOC between the T"+str(n_state)+" state and excited triplet states:" in line:
                catch2 = True
            elif "Total SOC between the T" in line and "state and excited triplet states:" in line:
                catch = True
            elif (catch and  'T'+str(n_state)+' ' in line and '(' not in line) or (catch2 and 'T' in line and '(' not in line):
                try:
                    socs.append(float(line.split()[1]))
                except:
                    catch, catch2 = False, False
            elif len(line.split()) < 2:
                catch, catch2 = False, False
    socs = np.array(socs)
    #socs = socs[order_s]
    return socs[np.newaxis,:]*0.12398/1000
#########################################################################################

##DECIDES WHICH FUNCTION TO USE IN ORDER TO GET SOCS#####################################
def avg_socs(files,tipo,n_state):
    if tipo == 'singlet':
        pega_soc = pega_soc_S
    elif tipo == 'triplet':
        pega_soc = pega_soc_T
    elif tipo == 'ground':
        pega_soc = pega_soc_ground 
    elif tipo == 'tts':
        pega_soc = pega_soc_TT        
    for file in files:
        socs = pega_soc(file,n_state)
        try:
            socs  = socs[:,:col]
            Socs  = np.vstack((Socs,socs))
        except:
            col   = socs.shape[1] 
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
def pega_oscs(files, IND,initial):
    spin = initial[0].upper()
    num  = int(initial[1:]) -1
    mapa = {'S':'Singlet','T':'Triplet'}
    frase    = "Transition Moments Between "+mapa[spin]+" Excited States"
    for i in range(len(files)):
        oscs  = []
        ind   = IND[i,num]
        ind_s = IND[i,:]
        location = np.where(ind_s == ind)[0][0]
        ind_s = ind_s[location+1:]
        ind   = str(ind)
        with open('Geometries/'+files[i], 'r') as f:
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
            try:
                Oscs = np.vstack((Oscs,np.array(oscs)[np.newaxis,:]))
            except:
                Oscs = np.array(oscs)[np.newaxis,:]         
    return Oscs    
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
        if 0 in ets[n_triplet]-ess:
            return 0
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

def phosph_osc(file,n_state,ind_s,ind_t,singlets,triplets):
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
    MMs  = np.array(MMs)
    term = e*(hbar2**2)/triplets
    Os = (2*mass)*MMs/(3*term)
    return Os[np.newaxis,:]

def get_osc_phosph(files,Singlets, Triplets, Ss_s, Ss_t, IND_S, IND_T):
    n_state = read_cis(files[0])
    Es  = Singlets #- (alphast2/alphaopt1)*Ss_s
    Et  = Triplets #- (alphast2/alphaopt1)*Ss_t #removed correction from phosph_osc calculation 
    for j in range(Singlets.shape[0]):
        tos = phosph_osc(files[j],n_state,IND_S[j,:],IND_T[j,:],Es[j,:],Et[j,:])        
        try:
            Os  = np.vstack((Os,tos))
        except:
            Os  = tos
    return Os

##GETS ALL RELEVANT INFORMATION FROM LOG FILES###########################################
def analysis(files):         
    n_state = read_cis(files[0])
    Numbers = []
    for file in files:
        singlets, triplets, oscs, ind_s, ind_t, ss_s, ss_t, gp = pega_energias('Geometries/'+file)       
        singlets = np.array([singlets[:n_state]])
        triplets = np.array([triplets[:n_state]])
        oscs     = np.array([oscs[:n_state]])
        ss_s     = np.array([ss_s[:n_state]])
        ss_t     = np.array([ss_t[:n_state]])
        ind_s    = np.array([ind_s[:n_state]])
        ind_t    = np.array([ind_t[:n_state]])
        gp       = np.array([gp])
        try:
            Singlets = np.vstack((Singlets,singlets))
            Triplets = np.vstack((Triplets,triplets))
            Oscs     = np.vstack((Oscs,oscs))
            Ss_s     = np.vstack((Ss_s,ss_s))
            Ss_t     = np.vstack((Ss_t,ss_t))
            IND_S    = np.vstack((IND_S,ind_s))
            IND_T    = np.vstack((IND_T,ind_t))
            GP       = np.append(GP,gp) 
        except:
            Singlets = singlets
            Triplets = triplets
            Oscs     = oscs
            Ss_s     = ss_s
            Ss_t     = ss_t    
            IND_S    = ind_s
            IND_T    = ind_t
            GP       = gp
        Numbers.append(int(file.split('-')[1]))
    Numbers = np.array(Numbers)[:,np.newaxis]    
    return Numbers, Singlets, Triplets, Oscs, Ss_s, Ss_t, GP, IND_S, IND_T
#########################################################################################

##PRINTS EMISSION SPECTRUM###############################################################
def printa_espectro_emi(initial,eps,nr,tdm,x,mean_y,error):
        mean_rate, error_rate = nemo.tools.calc_emi_rate(x, mean_y,error)
        primeira   = f"{'#Energy(ev)':4s} {'diff_rate':4s} {'error':4s} TDM={tdm:.3f} au\n"
        primeira  += f'# Total Rate {initial} -> S0: {mean_rate:5.2e} +/- {error_rate:5.2e} s^-1\n'
        primeira  += f'#Epsilon: {eps:.3f} nr: {nr:.3f}\n'
        arquivo    = nemo.tools.naming(f'differential_rate_{initial.upper()}.lx')
        with open(arquivo, 'w') as f:
            f.write(primeira)
            for i in range(0,len(x)):
                text = f"{x[i]:.6f} {mean_y[i]:.6e} {error[i]:.6e}\n"
                f.write(text)
        print(f'Spectrum printed in the {arquivo} file')        
#######################################################################################    

###CALCULATES WEIGHTED AVERAGES WHEN POSSIBLE##########################################
def means(y,weight,ensemble_mean=False):
    if ensemble_mean:
        try:
            mean = np.mean(y,axis=0)
        except:
            mean = np.mean(y)
    else:    
        try:
            mean = np.average(y,axis=0,weights=weight)
        except:
            mean = np.average(y,axis=0)
    return mean        
########################################################################################

###FORMATS RATES AND ERRORS IN THE SAME EXPONENT########################################
def format_rate(r,dr):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exp = np.nan_to_num(np.floor(np.log10(r)))
    try:
        exp[exp < -99] = -99
    except:
        exp = max(exp,-99)
    pre_r = r/10**exp
    pre_dr= dr/10**exp
    return pre_r, pre_dr, exp
#########################################################################################

###SAVES ENSEMBLE DATA#################################################################
def gather_data(initial,save=True):
    files =  [i for i in os.listdir('Geometries') if '.log' in i]    
    files = check_normal(files)
    files = sorted(files, key=lambda pair: float(pair.split('-')[1]))
    n_state     = int(initial[1:]) -1
    eps_i, nr_i = nemo.tools.get_nr()
    kbT         = nemo.tools.detect_sigma()
    if 's' in initial.lower():
        Numbers, Singlets, Triplets, Oscs, Ss_s, Ss_t, GP, IND_S, IND_T   = analysis(files) 
        if 's0' == initial.lower():
            label_oscs = [f'osc_s{i+1}' for i in range(Oscs.shape[1])]
        else:
            Oscs       = Oscs[:,n_state][:,np.newaxis]
            label_oscs = ['osc']
            noscs      = pega_oscs(files, IND_S,initial)
            label_oscs.extend([f'osc_s{n_state+2+i}' for i in range(noscs.shape[1])])
            Oscs       = np.hstack((Oscs,noscs))     
        try:
            socs_complete  = avg_socs(files,'singlet',n_state)
            header7 = ['soc_t'+str(i) for i in range(1,1+socs_complete.shape[1])]
        except:
            pass
    else:
        Numbers, Singlets, Triplets, _, Ss_s, Ss_t, GP, IND_S, IND_T = analysis(files)
        Oscs       = get_osc_phosph(files,Singlets, Triplets, Ss_s, Ss_t, IND_S, IND_T)
        Oscs       = Oscs[:,n_state][:,np.newaxis]
        label_oscs = ['osc']
        noscs      = pega_oscs(files, IND_T,initial)
        Oscs       = np.hstack((Oscs,noscs)) 
        label_oscs.extend([f'osc_t{n_state+2+i}' for i in range(noscs.shape[1])])
        try:
            socs_complete = np.hstack((avg_socs(files,'ground',n_state),avg_socs(files,'triplet',n_state),avg_socs(files,'tts',n_state)))
            indices  = [i+1 for i in range(Triplets.shape[1]) if i != n_state] #Removed Tn to Tn transfers
            header7  = ['soc_s0']
            header7.extend(['soc_s'+str(i) for i in range(1,1+Singlets.shape[1])])
            header7.extend(['soc_t'+str(i) for i in indices])   
        except:
            pass
    header = ['geometry']
    header.extend(['e_s'+str(i) for i in range(1,1+Singlets.shape[1])])
    header.extend(['e_t'+str(i) for i in range(1,1+Triplets.shape[1])])
    header.extend(['d_s'+str(i) for i in range(1,1+Ss_s.shape[1])])
    header.extend(['d_t'+str(i) for i in range(1,1+Ss_t.shape[1])])
    header.extend(['gp'])
    header.extend(label_oscs)
    try:
        header.extend(header7)
        data = np.hstack((Numbers,Singlets,Triplets,Ss_s,Ss_t,GP[:,np.newaxis],Oscs,socs_complete))
    except:
        data = np.hstack((Numbers,Singlets,Triplets,Ss_s,Ss_t,GP[:,np.newaxis],Oscs))
    arquivo = f'Ensemble_{initial.upper()}_.lx'
    data = pd.DataFrame(data,columns=header)
    data.insert(0,'ensemble',initial.upper())
    data.insert(0,'kbT',kbT)
    data.insert(0,'nr',nr_i)
    data.insert(0,'eps',eps_i)
    if save:
        data.to_csv(arquivo,index=False)
    return data
#######################################################################################

###PRINTS RATES AND EMISSION SPECTRUM##################################################
def export_results(data,emission,dielec):
    data    = data.copy()
    initial = data['Transition'][0].split('>')[0][:-1]
    printa_espectro_emi(initial,dielec[0],dielec[1],emission['TDM'][0],emission['Energy'].values,emission['Diffrate'].values,emission['Error'].values)
    pre_r, pre_dr, exp = format_rate(data['Rate(s^-1)'],data['Error(s^-1)'])
    rate    = [f'{pre_r[i]:5.2f}e{exp[i]:+03.0f}' for i in range(len(pre_r))]
    error   = [f'{pre_dr[i]:5.2f}e{exp[i]:+03.0f}' for i in range(len(pre_dr))]
    headers = [i for i in data.columns.values if i != 'Rate(s^-1)' and i != 'Error(s^-1)' and i != 'Transition' ]
    for header in headers:
        if header == 'Prob(%)' or header == 'AvgConc(%)':
            data[header] = data[header].map('{:.1f}'.format)
        else:
            data[header] = data[header].map('{:.3f}'.format)
    data['Rate(s^-1)']  = rate
    data['Error(s^-1)'] = error
    arquivo = nemo.tools.naming(f'rates_{initial}_.lx')  
    solvent = f'#Epsilon: {dielec[0]:.3f} nr: {dielec[1]:.3f}\n'
    with open(arquivo, 'w') as f:
        f.write(solvent+data.to_string(header=True, index=False))
    print(f'Rates are written in the {arquivo} file')  
#######################################################################################

###CALCULATES ISC AND EMISSION RATES & SPECTRA#########################################
def rates(initial,dielec,data=None,ensemble_average=False):
    if data is None:
        data        = gather_data(initial,save=True) 
        eps_i, nr_i = nemo.tools.get_nr()
        kbT         = nemo.tools.detect_sigma()
    else:
        eps_i  = data['eps'][0]
        nr_i   = data['nr'][0]
        kbT    = data['kbT'][0]
    eps, nr    = dielec[0], dielec[1]
    alphast1   = nemo.tools.get_alpha(eps_i)
    alphast2   = nemo.tools.get_alpha(eps)  
    alphaopt1  = nemo.tools.get_alpha(nr_i**2)
    alphaopt2  = nemo.tools.get_alpha(nr**2)
    n_state    = int(initial[1:]) -1
    initial    = initial.lower()       

    #Emission Calculations
    lambda_be  = (alphast2/alphast1 - alphaopt2/alphast1)*data['gp']
    Ltotal     = np.sqrt(2*lambda_be*kbT + kbT**2)
    delta_emi  = data['e_'+initial] - (alphast2/alphaopt1)*data['d_'+initial]
    if 's' in initial:
        constante = ((nr**2)*(e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
    elif 't' in initial:
        constante = (1/3)*((nr**2)*(e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))

    espectro  = (constante*((delta_emi-lambda_be)**2)*data['osc'])
    tdm       = nemo.tools.calc_tdm(data['osc'],data['e_'+initial],espectro)    
    left      = max(min(delta_emi-2*Ltotal),0.01)
    right     = max(delta_emi+2*Ltotal)    
    x         = np.linspace(left,right, int((right-left)/0.01))
    y         = espectro.to_numpy()[:,np.newaxis]*nemo.tools.gauss(x,delta_emi.to_numpy()[:,np.newaxis],Ltotal.to_numpy()[:,np.newaxis])
    N         = y.shape[0]
    mean_y    = np.sum(y,axis=0)/N 
    error     = np.sqrt(np.sum((y-mean_y)**2,axis=0)/(N*(N-1))) 
    emi_rate, emi_error = nemo.tools.calc_emi_rate(x, mean_y,error)     
    gap_emi        = means(delta_emi,espectro,ensemble_average)    
    mean_sigma_emi = means(Ltotal,espectro,ensemble_average) 
    mean_part_emi  = (100/N)/means(espectro/np.sum(espectro),espectro,ensemble_average) 
    emi            = np.hstack((x[:,np.newaxis],mean_y[:,np.newaxis],error[:,np.newaxis]))
    emi            = pd.DataFrame(emi,columns=['Energy','Diffrate','Error'])
    emi.insert(0,'TDM',tdm)

    #Checks number of logs
    if data is None:
        N    = data.shape[0]
        coms = nemo.tools.start_counter()
        if N != coms:
            print(f"There are {coms} inputs and just {N} log files. Something is not right! Computing the rates anyway...")
    
    #Intersystem Crossing Rates
    Singlets    =  data[[i for i in data.columns.values if 'e_s' in i]].to_numpy()
    Triplets    =  data[[i for i in data.columns.values if 'e_t' in i]].to_numpy()
    Ss_s        =  data[[i for i in data.columns.values if 'd_s' in i]].to_numpy()
    Ss_t        =  data[[i for i in data.columns.values if 'd_t' in i]].to_numpy()
    if 's' in initial:
        socs_complete =  data[[i for i in data.columns.values if 'soc_t' in i]].to_numpy()
        delta    = Triplets + np.repeat((alphast2/alphaopt1)*Ss_s[:,n_state][:,np.newaxis] - Singlets[:,n_state][:,np.newaxis],Triplets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_t   #Tn (final) - Sm (initial) + lambda_b
        lambda_b = (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_t
        final    = [i.upper()[4:] for i in data.columns.values if 'soc_t' in i]
        ##FOR WHEN IC IS AVAILABLE
        #socs_complete = np.hstack((socs_complete,0.0001*np.ones((Singlets.shape[0],Singlets.shape[1]-1))))
        #delta_ss = Singlets + np.repeat((alphast2/alphaopt1)*Ss_s[:,n_state][:,np.newaxis] - Singlets[:,n_state][:,np.newaxis],Singlets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_s    #Sm (final) - Sn (initial) + lambda_b
        #indices  = [i for i in range(Singlets.shape[1]) if i != n_state] #Removed Sn to Sn transfers
        #delta    = np.hstack((delta,delta_ss[:,indices]))
        #lambda_bt= (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_s
        #lambda_b = np.hstack((lambda_b,lambda_bt[:,indices]))
    elif 't' in initial: 
        #Tn to Sm ISC
        socs_complete =  data[[i for i in data.columns.values if 'soc_s' in i and 'soc_s0' not in i]].to_numpy()
        delta         = Singlets + np.repeat((alphast2/alphaopt1)*Ss_t[:,n_state][:,np.newaxis] - Triplets[:,n_state][:,np.newaxis],Singlets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_s    #Sm (final) - Tn (initial) + lambda_b
        lambda_b      = (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_s
        socs_s0       = data['soc_s0'].to_numpy()[:,np.newaxis]
        socs_complete = np.hstack((socs_s0,socs_complete))
        delta         = np.hstack((delta_emi.to_numpy()[:,np.newaxis],delta))
        lambda_b      = np.hstack((lambda_be.to_numpy()[:,np.newaxis],lambda_b))
        final         = [i.upper()[4:] for i in data.columns.values if 'soc_s' in i]
        #Tn to Tm ISC
        #delta_tt = Triplets + np.repeat((alphast2/alphaopt1)*Ss_t[:,n_state][:,np.newaxis] - Triplets[:,n_state][:,np.newaxis],Triplets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_t    #Tm (final) - Tn (initial) + lambda_b
        #indices  = [i for i in range(Triplets.shape[1]) if i != n_state] #Removed Tn to Tn transfers
        #delta    = np.hstack((delta,delta_tt[:,indices]))
        #lambda_bt= (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_t
        #lambda_b = np.hstack((lambda_b,lambda_bt[:,indices]))
        #final.extend([i.upper()[4:] for i in data.columns.values if 'soc_t' in i])

    sigma = np.sqrt(2*lambda_b*kbT + kbT**2)
    y     = (2*np.pi/hbar)*(socs_complete**2)*nemo.tools.gauss(delta,0,sigma)
    N     = y.shape[0]
    rate  = np.sum(y,axis=0)/N
    total = emi_rate + np.sum(rate)
    #Error estimate
    error = np.sqrt(np.sum((y-rate)**2,axis=0)/(N*(N-1)))  
    
    results    = np.array([[emi_rate,emi_error,100*emi_rate/total,gap_emi,np.nan,mean_sigma_emi,mean_part_emi]])
    mean_gap   = means(delta,y,ensemble_average)[:,np.newaxis]
    mean_soc   = 1000*means(socs_complete,y,ensemble_average)[:,np.newaxis]
    mean_sigma = means(sigma,y,ensemble_average)[:,np.newaxis]
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mean_part  = np.nan_to_num(100*rate/means(y,y,ensemble_average))
    rate       = rate[:,np.newaxis]
    error      = error[:,np.newaxis]
    labels = [f'{initial.upper()}->S0'] + [f'{initial.upper()}~>{j}' for j in final]
    labels = np.array(labels)
    results_isc = np.hstack((rate,error,100*rate/total,mean_gap,mean_soc,mean_sigma,mean_part[:,np.newaxis])) 
    results = np.vstack((results,results_isc))
    results = pd.DataFrame(results,columns=['Rate(s^-1)','Error(s^-1)','Prob(%)','AvgDE+L(eV)','AvgSOC(meV)','AvgSigma(eV)','AvgConc(%)'])
    results.insert(0,'Transition',labels)
    return results, emi       
#########################################################################################    

###COMPUTES ABSORPTION SPECTRA########################################################### 
def absorption(initial,dielec,data=None, save=False):
    if data is None:
        data        = gather_data(initial,save=True) 
        eps_i, nr_i = nemo.tools.get_nr()
        kbT         = nemo.tools.detect_sigma()
    else:
        eps_i  = data['eps'][0]
        nr_i   = data['nr'][0]
        kbT    = data['kbT'][0]
    eps, nr    = dielec[0], dielec[1]
    alphast1   = nemo.tools.get_alpha(eps_i)
    alphast2   = nemo.tools.get_alpha(eps)  
    alphaopt1  = nemo.tools.get_alpha(nr_i**2)
    alphaopt2  = nemo.tools.get_alpha(nr**2)
    initial    = initial.lower()
    constante  = (np.pi*(e**2)*hbar)/(2*nr*mass*c*epsilon0)*1e20
    if initial == 's0':
        engs = [i for i in data.columns if 'e_s' in i]
        ds   = [i for i in data.columns if 'd_s' in i]
        oscs = [i for i in data.columns if 'osc_s' in i]
        oscs = data[oscs].values
        engs = data[engs].values
        ds   = data[ds].values
        lambda_b = (alphast2/alphaopt1 - alphaopt2/alphaopt1)*ds
        DE   = engs - (alphaopt2/alphaopt1)*ds
    else:
        spin = initial[0]
        num  = int(initial[1:])
        engs = [i for i in data.columns if f'e_{spin}' in i]
        engs = [i for i in engs if int(i.split('_')[1][1:]) > num]
        ds   = [i for i in data.columns if f'd_{spin}' in i]
        ds   = [i for i in ds if int(i.split('_')[1][1:]) > num]
        oscs = [i for i in data.columns if 'osc_' in i]
        oscs = [i for i in oscs if int(i.split('_')[1][1:]) > num]
        base = data[f'e_{initial}'].values[:,np.newaxis]
        bs   = data[f'd_{initial}'].values[:,np.newaxis]
        engs = data[engs].values
        ds   = data[ds].values
        oscs = data[oscs].values
        lambda_b = (alphast2/alphaopt1 - alphaopt2/alphaopt1)*ds
        DE   = engs - (alphaopt2/alphaopt1)*ds - np.repeat(base -(alphast2/alphaopt1)*bs,engs.shape[1],axis=1)  
    lambda_b = lambda_b.flatten()
    DE       = DE.flatten()
    oscs     = oscs.flatten()
    Ltotal = np.sqrt(2*lambda_b*kbT + kbT**2)
    left   = max(min(DE-2*Ltotal),0.01)
    right  = max(DE+2*Ltotal)    
    x      = np.linspace(left,right,int((right-left)/0.01))
    y      = constante*oscs[:,np.newaxis]*nemo.tools.gauss(x,(DE+lambda_b)[:,np.newaxis],Ltotal[:,np.newaxis])
    N      = oscs.shape[0]
    mean_y = np.sum(y,axis=0)/N 
    #Error estimate
    sigma    = np.sqrt(np.sum((y-mean_y)**2,axis=0)/(N*(N-1))) 
    if save:
        arquivo  = nemo.tools.naming(f'cross_section_{initial.upper()}_.lx')
        primeira = f"{'#Energy(ev)':8s} {'cross_section(A^2)':8s} {'error(A^2)':8s}\n# Absorption from State: {initial.upper()}\n#Epsilon: {eps:.3f} nr: {nr:.3f}\n"
        with open(arquivo, 'w') as f:
            f.write(primeira)
            for i in range(0,len(x)):
                text = f"{x[i]:.6f} {mean_y[i]:.6e} {sigma[i]:.6e}\n"
                f.write(text)
        print(f'Spectrum printed in the {arquivo} file')
    return np.hstack((x[:,np.newaxis],mean_y[:,np.newaxis],sigma[:,np.newaxis]))                 
#########################################################################################