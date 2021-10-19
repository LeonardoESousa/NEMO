#!/usr/bin/env python3
import numpy as np
import os
import sys
from scipy.stats import norm

##SOME CONSTANTS##############################################
epsilon0 = 8.854187817e-12   #F/m
hbar = 6.582119514e-16       #eV s
hbar2 = 1.054571800e-34      #J s
mass = 9.10938356e-31        #kg
c = 299792458                #m/s
e = 1.60217662e-19           #C
kb = 8.6173303e-5            #eV/K
amu = 1.660539040e-27        #kg
pi = np.pi
###############################################################

##ERROR FUNCTION###############################################
def fatal_error(msg):
    print(msg)
    sys.exit()
###############################################################

##GETS FREQUENCIES AND REDUCED MASSES##########################
def pega_freq(freqlog):
    F, M = [], []
    with open(freqlog, 'r') as f:
        for line in f:
            if "Frequency:" in line:
                line = line.split()
                for j in range(1,len(line)):
                    if float(line[j]) in F:
                        pass
                    F.append(float(line[j]))
            elif "Red. Mass:" in line:
                line = line.split()
                for j in range(2,len(line)):
                    M.append(float(line[j]))
    #conversion in angular frequency
    F = np.array(F)*(c*100*2*pi) 
    try:
        f = F[0]
    except:
        fatal_error("No frequencies in the log file! Goodbye!")
    #conversion from amu to kg
    M = np.asarray(M)*amu
    return F, M
###############################################################

##GETS ATOMS AND LAST GEOMETRY IN FILE#########################
def pega_geom(freqlog):
    if ".out" in freqlog:
        busca = "Nuclear Orientation"
        n = -1
        with open(freqlog, 'r') as f:
            for line in f:
                if busca in line and 'Dipole' not in line:
                    n = 0
                    G = np.zeros((1,3))
                    atomos = []
                elif n >= 0 and n < 1:
                    n += 1
                elif n >= 1 and "----------------------------------------------------------------" not in line:    
                    line = line.split()
                    NG = []
                    for j in range(2,len(line)):
                        NG.append(float(line[j]))
                    atomos.append(line[1])
                    G = np.vstack((G,NG))       
                    n += 1  
                elif "----------------------------------------------------------------" in line and n>1:
                    n = -1       
    else:
        G = np.zeros((1,3))
        atomos = []
        with open(freqlog, 'r') as f:
            for line in f:
                line = line.split()
                try:
                    vetor = np.array([float(line[1]),float(line[2]), float(line[3])])
                    atomos.append(line[0])
                    G = np.vstack((G,vetor))
                except:
                    pass
    try:
        G = G[1:,:]                 
    except:
        fatal_error("No geometry in the log file! Goodbye!")
    return G, atomos
###############################################################

##GETS NORMAL COORDINATES IN REGULAR PRECISION#################
def pega_modos(G,freqlog):
    F, M = pega_freq(freqlog)
    C = []
    n = -1
    num_atom = np.shape(G)[0]
    with open(freqlog, 'r') as f:
        for line in f:
            if n < 0 or n >= num_atom:
                if "X      Y      Z" in line:
                    n = 0
                else:
                    pass
            elif n >= 0 and n < num_atom:
                line = line.split()
                for j in range(1,len(line)):
                    C.append(float(line[j]))
                n += 1  
                
    num_modos = len(F)
    
    l = 0
    p = 0
    NNC = np.zeros((num_atom,1))
    while l < num_modos:
        NC = np.zeros((1,3))
        k =0
        while k < num_atom:     
            B = np.asarray(C[3*(l+3*k)+p:3*(l+3*k)+3+p])
            NC = np.vstack((NC,B))
            k += 1      
        NNC = np.hstack((NNC,NC[1:,:]))
        l += 1
        if l%3 == 0 and l != 0:
            p = p + (num_atom-1)*9  
    NNC = NNC[:,1:] #matriz com as coordenadas normais de cada modo
    D = np.zeros((3*num_atom,1))
    for i in range(0,len(F)):
        normal = NNC[:,3*i:3*i+3].flatten()
        normal = np.expand_dims(normal,axis=1)
        D = np.hstack((D,normal))
    D = D[:,1:]
    MM = np.zeros((1,len(F)))
    M = np.expand_dims(M,axis=0)
    for i in range(0,3*num_atom):
        MM = np.vstack((MM,M))
    M = MM[1:,:]
    return D
###############################################################

##WRITES ATOMS AND XYZ COORDS TO FILE##########################
def write_input(atomos,G,header,bottom,file):
    with open(file, 'w') as f:
        f.write(header)
        for i in range(0,len(atomos)):
            texto = "{:2s}  {:.14f}  {:.14f}  {:.14f}\n".format(atomos[i],G[i,0],G[i,1],G[i,2])
            f.write(texto)
        f.write(bottom+'\n')
###############################################################

##CHECKS FOR EXISTING GEOMETRIES###############################
def start_counter():
    files = [file for file in os.listdir('Geometries') if ".com" in file and "Geometr" in file]
    return len(files)
###############################################################

##SAMPLES GEOMETRIES###########################################
def sample_geom(freqlog, num_geoms, T, header):
    F, M = pega_freq(freqlog)
    if F[0] < 0:
        fatal_error("Imaginary frequency! Goodbye!")
    try:
        os.mkdir('Geometries')
    except:
        pass        
    counter = start_counter()
    G, atomos = pega_geom(freqlog)
    write_input(atomos,G,'#','','geom.lx')
    print("The optimized geometry that is used is saved in the opt_geom.lx file!")
    NNC = pega_modos(G,freqlog)
    num_atom = np.shape(G)[0]   
    print("\nGenerating geometries...\n")
    with open('Magnitudes_{:.0f}K_.lx'.format(T), 'a') as file:
        for n in range(1,num_geoms+1):
            A = np.zeros((3*num_atom,1))
            numbers = []
            for i in range(0,len(F)):
                scale = np.sqrt(hbar2/(2*M[i]*F[i]*np.tanh(hbar*F[i]/(2*kb*T))))
                normal = norm(scale=scale,loc=0)
                #Displacements in  Å
                q = normal.rvs()*1e10
                numbers.append(q)
                A += q*(np.expand_dims(NNC[:,i],axis=1))
            numbers = np.round(np.array(numbers)[np.newaxis,:],4)
            np.savetxt(file, numbers, delimiter='\t', fmt='%s')
            A = np.reshape(A,(num_atom,3))
            Gfinal = A + G  
            write_input(atomos,Gfinal,header,'$end',"Geometries/Geometry-"+str(n+counter)+"-.com")
            progress = 100*n/num_geoms
            text = "{:2.1f}%".format(progress)
            print(' ', text, "of the geometries done.",end="\r", flush=True)
    
    print("\n\nDone! Ready to run.")   
###############################################################    

            
##COLLECTS RESULTS############################################## 
def gather_data(opc):
    files = [file for file in os.listdir('Geometries') if ".log" in file and "Geometr" in file ]    
    files = sorted(files, key=lambda file: float(file.split("-")[1])) 
    from nemo.analysis import analysis
    Os, Singlets, Triplets, Oscs = analysis()
    num = np.shape(Singlets)[1]
    with open("Samples.lx", 'w') as f:
        for i in range(np.shape(Singlets)[0]):
            f.write("Geometry "+str(i)+":  Vertical transition (eV) Oscillator strength Broadening Factor (eV) Spin \n")
            for j in range(num):
                f.write("Excited State {}:\t{}\t{:.5e}\t{}\t{}\n".format(j+1,Singlets[i,j],Oscs[i,j],opc,'Singlet'))        
            for j in range(num):
                f.write("Excited State {}:\t{}\t{:.5e}\t{}\t{}\n".format(j+1,Triplets[i,j],Os[i,j],opc,'Triplet'))
############################################################### 


##NORMALIZED GAUSSIAN##########################################
def gauss(x,v,s):
    y =  (1/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-v)/s)**2)
    return y
###############################################################


##COMPUTES AVG TRANSITION DIPOLE MOMENT########################
def calc_tdm(O,V):
    #Energy terms converted to J
    term = e*(hbar**2)/V
    dipoles = np.sqrt(3*term*O/(2*mass))
    #Conversion in au
    dipoles *= (1/0.529177)*1e10
    return np.mean(dipoles)
###############################################################

##PREVENTS OVERWRITING#########################################
def naming(arquivo):
    new_arquivo = arquivo
    if arquivo in os.listdir('.'):
        duplo = True
        vers = 2
        while duplo:
            new_arquivo = str(vers)+arquivo
            if new_arquivo in os.listdir('.'):
                vers += 1
            else:
                duplo = False
    return new_arquivo        
###############################################################

##CASK FOR THE RELEVANT STATE##################################
def ask_states(frase):
    estados = input(frase)
    try:
        estados = int(estados)
    except:
        fatal_error("It must be an integer! Goodbye!")
    return estados
###############################################################

##COMPUTES SPECTRA############################################# 
def spectra(tipo, num_ex, nr):
    if tipo == "abs":
        spin = 'Singlet'
        num_ex = range(0,num_ex+1)
        num_ex = list(map(int,num_ex))
        constante = (np.pi*(e**2)*hbar)/(2*nr*mass*c*epsilon0)*10**(20)
    elif tipo == 'fluor':
        spin = 'Singlet'
        num_ex = [num_ex]
        constante = ((nr**2)*(e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
    elif tipo == 'phosph':
        spin = 'Triplet'
        num_ex = [num_ex]
        constante = (1/3)*((nr**2)*(e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
    V, O, S = [], [], []
    N = 0
    with open("Samples.lx", 'r') as f:
        for line in f:
            if "Geometry" in line:
                N += 1
            elif "Excited State" in line and int(line.split()[2][:-1]) in num_ex and spin in line:
                line = line.split()
                V.append(float(line[3]))
                O.append(float(line[4]))
                S.append(float(line[5]))
    coms = start_counter()
    if len(V) == 0 or len(O) == 0:
        fatal_error("You need to run steps 1 and 2 first! Goodbye!")
    elif len(V) != coms*len(num_ex):
        print("Number of log files is less than the number of inputs. Something is not right! Computing the spectrum anyway...")
    V = np.array(V)
    O = np.array(O)
    S = np.array(S)
    if tipo == 'abs':
        espectro = (constante*O)
    else:
        espectro = (constante*(V**2)*O)
        tdm = calc_tdm(O,V)
    x  = np.linspace(min(V)-3*max(S), max(V)+ 3*max(S), 200)
    y  = np.zeros((1,len(x)))
    if tipo == 'abs':
        arquivo = 'cross_section.lx'
        primeira = "{:8s} {:8s} {:8s}\n".format("#Energy(ev)", "cross_section(A^2)", "error")
    else:
        arquivo = tipo+'_differential_rate.lx'
        primeira = "{:4s} {:4s} {:4s} TDM={:.3f} au\n".format("#Energy(ev)", "diff_rate", "error",tdm)
    arquivo = naming(arquivo)
    for i in range(0,len(espectro)):
        contribution = espectro[i]*gauss(x,V[i],S[i])
        y  = np.vstack((y,contribution[np.newaxis,:]))

    y = y[1:,:]
    mean_y =   np.sum(y,axis=0)/N 
    #Error estimate
    sigma  =   np.sqrt(np.sum((y-mean_y)**2,axis=0)/(N*(N-1))) 
    
    if tipo == 'fluor' or tipo == 'phosph':
        #Emission rate calculations
        mean_rate, error_rate = calc_emi_rate(x, mean_y,sigma) 
        segunda = '# Total Rate {}{} -> S0: {:5.2e} +/- {:5.2e} s^-1\n'.format(spin[0],num_ex[0],mean_rate,error_rate)
    else:
        segunda = ''

    print(N, "geometries considered.")     
    with open(arquivo, 'w') as f:
        f.write(primeira)
        f.write(segunda)
        for i in range(0,len(x)):
            text = "{:.6f} {:.6e} {:.6e}\n".format(x[i],mean_y[i], sigma[i])
            f.write(text)
    print('Spectrum printed in the {} file'.format(arquivo))                
############################################################### 

##CHECKS THE FREQUENCY LOG'S LEVEL OF THEORY###################
def busca_input(freqlog):
    spec = 'ABSSPCT'
    root = '1'
    with open(freqlog, 'r') as f:
        search = False
        molec  = False
        for line in f:
            if 'User input:' in line:
                rem = ''    
                search = True
            elif search and 'jobtype' not in line.lower() and '$molecule' not in line.lower() and '-----' not in line and 'cis_' not in line.lower():
                rem += line
            elif 'CIS_STATE_DERIV' in line.upper():
                spec = 'EMISPCT'
                root = line.split()[-1]
            elif search and '$molecule' in line.lower():
                molec = True
                search = False
            elif molec:
                line = line.split()
                if len(line) == 2:
                    cm = ' '.join(line)
                elif '$end' in line:
                    molec = False
                    search = True
            elif '--------------------------------------------------------------' in line and search and rem != '':
                search = False
    if spec == 'EMISPCT':
        from nemo.analysis import pega_energias
        _, _, _, ind_s, ind_t = pega_energias(freqlog)
        if root in ind_s:
            spec = 'FLUORSPCT'
        elif root in ind_t:
            spec = 'PHOSPHSPCT'
    return rem, cm, spec                
###############################################################

##CHECKS PROGRESS##############################################
def andamento():
    try:
        coms = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.com' in file]
        logs = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.log' in file]
        factor = 1
        with open('Geometries/'+coms[0], 'r') as f:
            for line in f:
                if 'Link1' in line:
                    factor = 2
        count = 0
        error = 0 
        for file in logs:
            with open('Geometries/'+file, 'r') as f:
                for line in f:
                    if "Have a nice day" in line:
                        count += 1
                    elif "fatal error" in line:
                        error += 1    
        print("\n\nThere are", int(count/factor), "successfully completed calculations out of", len(coms), "inputs")
        if error > 0:
            print("There are {} failed jobs. If you used option 2, check the nohup.out file for details.".format(error))                
        print(np.round(100*(count+error)/(factor*len(coms)),1), "% of the calculations have been run.")
    except:
        print('No files found! Check the folder!')                
###############################################################


##FETCHES  FILES###############################################
def fetch_file(frase,ends):
    files = []
    for file in [i for i in os.listdir('.')]:
        for end in ends:
            if end in file:
                 files.append(file)
    if len(files) == 0:
        fatal_error("No {} file found. Goodbye!".format(frase))
    freqlog = 'nada0022'    
    for file in files:
        print("\n"+file)
        resp = input('Is this the {} file? y ou n?\n'.format(frase))
        if resp.lower() == 'y':
            freqlog = file
            break
    if freqlog == 'nada0022':
        fatal_error("No {} file found. Goodbye!".format(frase))
    return freqlog  
###############################################################  
   
##RUNS TASK MANAGER############################################
def batch():
    script = fetch_file('batch script?',['.sh'])    
    limite = input("Maximum number of batches to be submitted simultaneously?\n")
    nproc  = input('Number of processors for each individual job\n')
    num    = input('Number of jobs in each batch\n')
    try:
        limite = int(limite)
        int(nproc)
        int(num)
    except:
        fatal_error("It must be an integer. Goodbye!")
    
    import subprocess
    folder = os.path.dirname(os.path.realpath(__file__)) 
    with open('limit.lx','w') as f:
        f.write(str(limite))
    subprocess.Popen(['nohup', 'python3', folder+'/batch_lx.py', script, nproc, num, '&'])
###############################################################

##FINDS SUITABLE VALUE FOR STD#################################    
def detect_sigma():
    try:
        files = [i for i in os.listdir('.') if 'Magnitudes' in i and '.lx' in i]
        file  = files[0]
        temp = float(file.split('_')[1].strip('K'))
        sigma =  np.round(kb*temp,3)
    except:
        sigma = 0.000
    return sigma
###############################################################    

##CHECKS SPECTRUM TYPE#########################################
def get_spec():
    coms = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.com' in file]
    with open('Geometries/'+coms[0],'r') as f:
        for line in f:
            if 'ABSSPCT' in line:
                tipo = 'absorption'
                break
            elif 'EMISPCT' in line:
                tipo = 'emission'
                break
            elif 'FLUORSPCT' in line:
                tipo = 'fluorescence'
            elif 'PHOSPHSPCT' in line:
                tipo = 'phosphorescence'    
    return tipo       
###############################################################

##FETCHES REFRACTIVE INDEX##################################### 
def get_nr():
    nr = 1
    coms = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.com' in file]
    with open('Geometries/'+coms[0],'r') as f:
        for line in f:
            if 'opticaldielectric' in line.lower():
                line = line.split()
                for elem in line:
                    try:
                        nr = np.sqrt(float(elem))
                        break
                    except:
                        pass
    return nr                
###############################################################

##QUERY FUNCTION###############################################
def default(a,frase):
    b = input(frase)
    if b == '':
        return a
    else:
        return b    
###############################################################

##STOP SUBMISSION OF JOBS######################################
def abort_batch():
    choice = input('Are you sure you want to prevent new jobs from being submitted? y or n?\n')
    if choice == 'y':
        try:
            os.remove('limit.lx')
            print('Done!')
        except:
            print('Could not find the files. Maybe you are in the wrong folder.')
    else:
        print('OK, nevermind')
###############################################################

##DELETES CHK FILES############################################
def delchk(input,term):
    num = input.split('-')[1]
    if term == 1:
        a = ''
    elif term == 2:
        a = '2'
    try:        
        os.remove('step{}_{}.chk'.format(a,num))
    except:
        pass      
###############################################################

##CHECKS WHETHER JOBS ARE DONE#################################
def watcher(files,counter):
    rodando = files.copy()
    done = []
    for input in rodando: 
        term = 0
        error = False
        try:
            with open(input[:-3]+'log', 'r') as f:
                for line in f:
                    if 'Have a nice day' in line:
                        term += 1
                        if counter == 2:
                            delchk(input,term)
                    elif 'fatal error' in line:
                        error = True
                        print('The following job returned an error: {}'.format(input))
                        print('Please check the file for any syntax errors.')        
            if term == counter or error:
                done.append(input)
        except:
            pass 
    for elem in done:
        del rodando[rodando.index(elem)]                                
    return rodando
###############################################################


##GETS SPECTRA#################################################
def search_spectra():
    Abs, Emi = 'None', 'None'
    candidates = [i for i in os.listdir('.') if '.lx' in i]
    for candidate in candidates:
        with open(candidate, 'r') as f:
            for line in f:
                if 'cross_section' in line:
                    Abs = candidate
                elif 'diff_rate' in line:     
                    Emi = candidate
                break
    return Abs, Emi
###############################################################

##CALCULATES FLUORESCENCE LIFETIME IN S########################
def calc_emi_rate(xd,yd,dyd):
    #Integrates the emission spectrum
    IntEmi = np.trapz(yd,xd)
    taxa   = (1/hbar)*IntEmi
    error  = (1/hbar)*np.sqrt(np.trapz((dyd**2),xd))
    return taxa, error 
###############################################################