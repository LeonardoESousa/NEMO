#!/usr/bin/env python3
import numpy as np
import os
import nemo.tools 
import warnings

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
def pega_oscs(file, ind,spin, ind_s):
    mapa = {'1':'Singlet','3':'Triplet'}
    frase    = "Transition Moments Between "+mapa[spin]+" Excited States"
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

def get_osc_phosph(alphast2,alphaopt1,Singlets, Triplets, Ss_s, Ss_t, IND_S, IND_T):
    files =  [i for i in os.listdir('Geometries') if '.log' in i]    
    files = check_normal(files)
    files = sorted(files, key=lambda pair: float(pair.split('-')[1]))
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
def analysis(phosph=True):         
    files =  [i for i in os.listdir('Geometries') if '.log' in i]    
    files = check_normal(files)
    files = sorted(files, key=lambda pair: float(pair.split('-')[1]))
    n_state = read_cis(files[0])

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
    if phosph:
        return Singlets, Triplets, Oscs, Ss_s, Ss_t, GP, IND_S, IND_T
    else:    
        return Singlets, Triplets, Oscs, Ss_s, Ss_t, GP
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

###SAVES ENSEMBLE DATA#################################################################
def save_data(Singlets,Triplets,Ss_s,Ss_t, GP,socs_complete,oscs,espectro,y,initial):
    dados   = np.hstack((Singlets,Triplets,Ss_s,Ss_t, GP[:,np.newaxis],socs_complete,oscs,espectro,y))
    header1 = ['S'+str(i) for i in range(1,1+Singlets.shape[1])]
    header2 = ['T'+str(i) for i in range(1,1+Triplets.shape[1])]
    header3 = ['DS'+str(i) for i in range(1,1+Ss_s.shape[1])]
    header4 = ['DT'+str(i) for i in range(1,1+Ss_t.shape[1])]
    header5 = ['GP']
    header6 = ['SOC'+str(i) for i in range(1,1+socs_complete.shape[1])]
    header7 = ['OSC']
    header8 = ['ke']
    header9 = ['kisc'+str(i) for i in range(1,1+y.shape[1])]
    header  = ','.join(header1+header2+header3+header4+header5+header6+header7+header8+header9)
    arquivo = nemo.tools.naming(f'Ensemble_{initial}_.lx')
    np.savetxt(arquivo,dados,fmt='+%.4e',header=header, delimiter=',')
#######################################################################################    

###CALCULATES WEIGHTED AVERAGES WHEN POSSIBLE##########################################
def means(y,weigh):
    try:
        mean = np.average(y,axis=0,weights=weigh)
    except:
        mean = np.average(y,axis=0)
    return mean        
########################################################################################

###FORMATS RATES AND ERRORS IN THE SAME EXPONENT########################################
def format_rate(r,dr):
    exp = np.nan_to_num(np.floor(np.log10(r)))
    pre_r = r/10**exp
    pre_dr= dr/10**exp
    return pre_r, pre_dr, exp
#########################################################################################

##CALCULATES ISC RATES FROM INITIAL STATE TO SEVERAL STATES OF OPPOSITE SPIN#############
def rates(initial,dielec):
    eps, nr     = dielec[0], dielec[1]
    eps_i, nr_i = nemo.tools.get_nr()
    alphast1    = nemo.tools.get_alpha(eps_i)
    alphast2    = nemo.tools.get_alpha(eps)  
    alphaopt1   = nemo.tools.get_alpha(nr_i**2)
    alphaopt2   = nemo.tools.get_alpha(nr**2)
    n_state     = int(initial[1:]) -1
    kbT         = nemo.tools.detect_sigma()
    coms        = nemo.tools.start_counter()
    
    #Emission Calculations
    if 's' in initial.lower():
        Singlets, Triplets, Oscs, Ss_s, Ss_t, GP   = analysis(phosph=False)
        tipo      = 'singlet'
        final     = 'T'
        lambda_b  = (alphast2/alphast1 - alphaopt2/alphast1)*GP 
        delta_emi = Singlets[:,n_state] - (alphast2/alphaopt1)*Ss_s[:,n_state]
        constante = ((nr**2)*(e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
        espectro  = (constante*((delta_emi-lambda_b)**2)*Oscs[:,n_state])
        tdm       = nemo.tools.calc_tdm(Oscs[:,n_state],Singlets[:,n_state],espectro)
    elif 't' in initial.lower():
        Singlets, Triplets, _, Ss_s, Ss_t, GP, IND_S, IND_T = analysis(phosph=True)
        Oscs      = get_osc_phosph(alphast2,alphaopt1,Singlets, Triplets, Ss_s, Ss_t, IND_S, IND_T)
        tipo      = 'triplet'
        final     = 'S'
        lambda_b  = (alphast2/alphast1 - alphaopt2/alphast1)*GP 
        delta_emi = Triplets[:,n_state] - (alphast2/alphaopt1)*Ss_t[:,n_state]
        constante = (1/3)*((nr**2)*(e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
        espectro  = (constante*((delta_emi-lambda_b)**2)*Oscs[:,n_state])
        tdm       = nemo.tools.calc_tdm(Oscs[:,n_state],Triplets[:,n_state],espectro)
    Ltotal    = np.sqrt(2*lambda_b*kbT + kbT**2)
    left      = max(min(delta_emi-2*Ltotal),0.01)
    right     = max(delta_emi+2*Ltotal)    
    x         = np.linspace(left,right, int((right-left)/0.01))
    y         = espectro[:,np.newaxis]*nemo.tools.gauss(x,delta_emi[:,np.newaxis],Ltotal[:,np.newaxis])
    N         = y.shape[0]
    mean_y    = np.sum(y,axis=0)/N 
    error     = np.sqrt(np.sum((y-mean_y)**2,axis=0)/(N*(N-1))) 
    emi_rate, emi_error = nemo.tools.calc_emi_rate(x, mean_y,error)     
    gap_emi        = np.average(delta_emi,weights=espectro)
    mean_sigma_emi = np.average(Ltotal,weights=espectro)
    mean_part_emi  = (100/N)/np.average(espectro/np.sum(espectro),weights=espectro)
    printa_espectro_emi(initial,eps,nr,tdm,x,mean_y,error)

    
    try:
        socs_complete = avg_socs(tipo,n_state)
    except:
        socs_complete = np.zeros(Singlets.shape)
        print('ISC rates are not available!')
    #Checks number of logs
    N = min(Singlets.shape[0],Triplets.shape[0],socs_complete.shape[0])
    if N != coms:
        print(f"There are {coms} inputs and just {N} log files. Something is not right! Computing the rates anyway...")
    
    #Intersystem Crossing Rates
    if tipo == 'singlet':
        delta    = Triplets + np.repeat((alphast2/alphaopt1)*Ss_s[:,n_state][:,np.newaxis] - Singlets[:,n_state][:,np.newaxis],Triplets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_t   #Tn (final) - Sm (initial) + lambda_b
        lambda_b = (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_t
    elif tipo == 'triplet':
        delta    = Singlets + np.repeat((alphast2/alphaopt1)*Ss_t[:,n_state][:,np.newaxis] - Triplets[:,n_state][:,np.newaxis],Singlets.shape[1],axis=1) - (alphaopt2/alphaopt1)*Ss_s    #Sm (final) - Tn (initial) + lambda_b
        lambda_b = (alphast2/alphaopt1 - alphaopt2/alphaopt1)*Ss_s
    sigma = np.sqrt(2*lambda_b*kbT + kbT**2)
    y     = (2*np.pi/hbar)*(socs_complete**2)*nemo.tools.gauss(delta,0,sigma)
    N     = y.shape[0]
    rate  = np.sum(y,axis=0)/N
    rate[rate< 1e-99] = 0 
    total = emi_rate + np.sum(rate)
    #Error estimate
    error = np.sqrt(np.sum((y-rate)**2,axis=0)/(N*(N-1)))
    error[error< 1e-99] = 0 
    save_data(Singlets,Triplets,Ss_s,Ss_t,GP,socs_complete,Oscs[:,n_state][:,np.newaxis],espectro[:,np.newaxis],y,initial.upper())
    arquivo = nemo.tools.naming(f'rates_{initial.upper()}_.lx')    
    with open(arquivo, 'w') as f:
        pre_r, pre_dr, exp = format_rate(emi_rate,emi_error)
        f.write(f'#Epsilon: {eps:.3f} nr: {nr:.3f}\n')
        f.write('#Transition    Rate(s^-1)    Error(s^-1)   Prob(%)   AvgDE+L(eV)  AvgSOC(meV)  AvgSigma(eV)   AvgConc(%)\n')     
        f.write(f'{initial.upper()}->S0         {pre_r:5.2f}e{exp:+03.0f}      {pre_dr:5.2f}e{exp:+03.0f}      {100*emi_rate/total:5.1f}         {gap_emi:+5.3f}       {"-":5}         {mean_sigma_emi:5.3f}        {mean_part_emi:5.1f}%\n')

        gap        = means(delta,y)
        mean_soc   = 1000*means(socs_complete,y)
        mean_sigma = means(sigma,y)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mean_part  = np.nan_to_num(100*rate/means(y,y))
        pre_r, pre_dr, exp = format_rate(rate,error)
        for j in range(delta.shape[1]):
            f.write(f'{initial.upper()}->{final}{j+1}         {pre_r[j]:5.2f}e{exp[j]:+03.0f}      {pre_dr[j]:5.2f}e{exp[j]:+03.0f}      {100*rate[j]/total:5.1f}         {gap[j]:+5.3f}       {mean_soc[j]:5.3f}         {mean_sigma[j]:5.3f}        {mean_part[j]:5.1f}%\n')

        #Internal conversion avg deltaE+L
        delta_s = np.mean(np.diff(Singlets - (alphast2/alphaopt1)*Ss_s,axis=1) + (alphast2/alphaopt1 -alphaopt2/alphaopt1)*Ss_s[:,1:],axis=0)
        delta_t = np.mean(np.diff(Triplets - (alphast2/alphaopt1)*Ss_t,axis=1) + (alphast2/alphaopt1 -alphaopt2/alphaopt1)*Ss_t[:,1:],axis=0)
        f.write('\n#Estimated Internal Conversion Driving Energies:\n')
        f.write('#Transition    AvgDE+L(eV)    Transition    AvgDE+L(eV)\n')
        for j in range(len(delta_s)):
            f.write(f'S{j+1}->S{j+2}         {delta_s[j]:+5.3f}         T{j+1}->T{j+2}         {delta_t[j]:+5.3f}\n')


    print(f'Rates are written in the {arquivo} file')        
#########################################################################################    
