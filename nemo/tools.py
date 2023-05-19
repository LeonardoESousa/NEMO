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
def pega_freq_gauss(freqlog):
    F, M = [], []
    with open(freqlog, 'r') as f:
        for line in f:
            if "Frequencies --" in line:
                line = line.split()
                for j in range(2,len(line)):
                    if float(line[j]) in F:
                        pass
                    F.append(float(line[j]))
            elif "Red. masses --" in line:
                line = line.split()
                for j in range(3,len(line)):
                    M.append(float(line[j]))
            elif 'Thermochemistry' in line:
                break        
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
            texto = f"{atomos[i]:2s}  {G[i,0]:.14f}  {G[i,1]:.14f}  {G[i,2]:.14f}\n"
            f.write(texto)
        f.write(bottom+'\n')
###############################################################

##CHECKS FOR EXISTING GEOMETRIES###############################
def start_counter():
    files = [file for file in os.listdir('Geometries') if ".com" in file and "Geometr" in file]
    return len(files)
###############################################################

##SAMPLES GEOMETRIES###########################################
def sample_geometries(freqlog,num_geoms,T, limit=np.inf):
    G, atomos = pega_geom(freqlog)
    F, M      = pega_freq(freqlog)
    F[F < 0] *= -1
    NNC       = pega_modos(G,freqlog)
    mask = F < limit*(c*100*2*pi)
    F = F[mask]
    NNC = NNC[:,mask]
    num_atom  = np.shape(G)[0]
    A = np.zeros((3*num_atom,num_geoms))
    for i in range(0,len(F)):
        scale = np.sqrt(hbar2/(2*M[i]*F[i]*np.tanh(hbar*F[i]/(2*kb*T))))
        normal = norm(scale=scale,loc=0)
        #Displacements in  Å
        q = normal.rvs(size=num_geoms)*1e10
        try:
            numbers = np.hstack((numbers,q[:,np.newaxis]))
        except:
            numbers = q[:,np.newaxis]
        A += np.outer(NNC[:,i],q)
    for n in range(np.shape(A)[1]):
        A1 = np.reshape(A[:,n],(num_atom,3))
        try:
            Gfinal = np.hstack((Gfinal,A1 + G))
        except:
            Gfinal = A1 + G     
    numbers = np.round(numbers,4)
    return numbers, atomos, Gfinal
###############################################################

##MAKES ENSEMBLE###############################################
def make_ensemble(freqlog, num_geoms, T, header, bottom):
    try:
        os.mkdir('Geometries')
    except:
        pass        
    counter = start_counter()   
    print("\nGenerating geometries...\n")
    numbers, atomos, A = sample_geometries(freqlog,num_geoms,T)
    with open(f'Magnitudes_{T:.0f}K_.lx', 'a') as file:
        np.savetxt(file, numbers, delimiter='\t', fmt='%s')
    for n in range(0,np.shape(A)[1],3):
        Gfinal = A[:,n:n+3]  
        write_input(atomos,Gfinal,header,bottom,"Geometries/Geometry-"+str((n+3)//3+counter)+"-.com")
        progress = 100*((n+3)//3)/num_geoms
        text = f"{progress:2.1f}%"
        print(' ', text, "of the geometries done.",end="\r", flush=True)
    print("\n\nDone! Ready to run.")   
################################################################

##NORMALIZED GAUSSIAN##########################################
def gauss(x,v,s):
    y =  (1/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-v)/s)**2)
    return y
###############################################################


##COMPUTES AVG TRANSITION DIPOLE MOMENT########################
def calc_tdm(O,V,pesos):
    #Energy terms converted to J
    term = e*(hbar2**2)/V
    dipoles = np.sqrt(3*term*O/(2*mass))
    #Conversion in au
    dipoles *= 1.179474389E29
    return np.average(dipoles,weights=pesos)
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
        int(estados[1:])
    except:
        fatal_error("It must be S or T and an integer! Goodbye!")
    if estados[0].upper() != 'S' and estados[0].upper() != 'T':
        fatal_error("It must be S or T and an integer! Goodbye!")
    return estados.upper()
###############################################################

def get_alpha(eps):
    return (eps-1)/(eps+1)

##LIST OF KEYWORDS THAT SHOULD NOT BE READ#####################
def delist(elem):
    words = ['jobtype', '-----', 'cis_n', 'cis_s', 'cis_t', 'gui', 'nto_', 'soc', 'sts_', 'CIS_RELAXED_DENSITY', 'solvent_method']
    for w in words:
        if w in elem.lower():
            return False
    return True        
###############################################################

##CHECKS THE FREQUENCY LOG'S LEVEL OF THEORY###################
def busca_input(freqlog):
    search = True
    with open(freqlog, 'r') as f:
        for line in f:
            if 'A Quantum Leap Into The Future Of Chemistry' in line:
                search = False
                break           
    spec = 'ABSSPCT'
    rem  = ''
    with open(freqlog, 'r') as f:
        for line in f:
            if 'User input:' in line and not search:    
                search = True    
            elif search and delist(line):
                rem += line
            elif 'CIS_STATE_DERIV' in line.upper():
                spec = 'EMISPCT'
            elif '--------------------------------------------------------------' in line and search and rem != '':
                search = False
    rem = rem.split('$end')
    for r in rem:
        if '$rem' in r:
            rem = r+'$end\n'
        elif '$molecule' in r:
            mol = r.split('\n')
            for m in mol:
                m = m.split()
                if len(m) == 2:
                    cm = ' '.join(m)
                    break                 
    return rem, cm, spec                
###############################################################

##CHECKS PROGRESS##############################################
def andamento():
    try:
        coms = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.com' in file and '.com_' not in file]
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
            print(f"There are {error} failed jobs. If you used option 2, check the nohup.out file for details.")                
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
        fatal_error(f"No {frase} file found. Goodbye!")
    freqlog = 'nada0022'    
    for file in files:
        print("\n"+file)
        resp = input(f'Is this the {frase} file? y ou n?\n')
        if resp.lower() == 'y':
            freqlog = file
            break
    if freqlog == 'nada0022':
        fatal_error(f"No {frase} file found. Goodbye!")
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
        print('WARNING: Magnitudes.lx file is absent! Temperature is being set to 300 K!')
        sigma = 0.026
    return sigma
###############################################################    

##FETCHES REFRACTIVE INDEX##################################### 
def get_nr():
    nr = 1
    epsilon = 1
    coms = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.com' in file]
    with open('Geometries/'+coms[0],'r') as f:
        for line in f:
            if 'opticaldielectric' in line.lower():
                nr = np.sqrt(float(line.split()[1]))    
            elif 'dielectric' in line.lower() and 'optical' not in line.lower():
                epsilon = float(line.split()[1])
    return epsilon, nr                
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

##CHECKS WHETHER JOBS ARE DONE#################################
def watcher(files,counter,first):
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
                    elif ('fatal error' in line or 'failed standard') in line and not first:
                        error = True
                        print(f'The following job returned an error: {input}')
                        print('Please check the file for any syntax errors.') 
                    elif ('fatal error' in line or 'failed standard') in line and first:
                        os.remove(input[:-3]+'log')          
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

