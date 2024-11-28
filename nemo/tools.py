#!/usr/bin/env python3
import os
import subprocess
import sys
import time
import requests
import pkg_resources
from subprocess import Popen
import numpy as np
import pandas as pd
from scipy.stats import norm
from joblib import Parallel, delayed
import lx.tools
import lx.parser
import nemo.parser

LIGHT_SPEED = nemo.parser.LIGHT_SPEED
HBAR_EV = nemo.parser.HBAR_EV
HBAR_J = nemo.parser.HBAR_J
BOLTZ_EV = nemo.parser.BOLTZ_EV
E_CHARGE = nemo.parser.E_CHARGE
MASS_E = nemo.parser.MASS_E
EPSILON_0 = nemo.parser.EPSILON_0

###############################################################

##WRITES ATOMS AND XYZ COORDS TO FILE##########################
def write_input(atomos, geometry, header, bottom, file):
    with open(file, "w", encoding="utf-8") as input_file:
        input_file.write(header)
        for i, atomo in enumerate(atomos):
            texto = f"{atomo:2s}  {geometry[i,0]:.7f}  {geometry[i,1]:.7f}  {geometry[i,2]:.7f}\n"
            input_file.write(texto)
        input_file.write(bottom + "\n")


###############################################################


##CHECKS FOR EXISTING GEOMETRIES###############################
def start_counter():
    files = [
        file
        for file in os.listdir("Geometries")
        if ".com" in file and "Geometr" in file
    ]
    return len(files)


###############################################################


def sample_single_geometry(args):
    geom, atomos, old, scales, normal_coord, warning = args
    rejected_geoms = 0
    ok = False
    
    while not ok:
        start_geom = geom.copy()
        qs = [norm(scale=scale, loc=0).rvs(size=1) for scale in scales]
        qs = np.array(qs)
        start_geom += np.sum(qs.reshape(1, 1, -1) * normal_coord, axis=2)
        new = lx.tools.adjacency(start_geom, atomos)
        if 0.5 * np.sum(np.abs(old - new)) < 1 or not warning:
            ok = True
            return (start_geom, qs.T, rejected_geoms)
        else:
            rejected_geoms += 1

def sample_geometries(freqlog, num_geoms, temp, limit=np.inf, warning=True, show_progress=False):
    geom, atomos = nemo.parser.pega_geom(freqlog)
    old = lx.tools.adjacency(geom, atomos)
    freqs, masses = nemo.parser.pega_freq(freqlog)
    normal_coord = nemo.parser.pega_modos(geom, freqlog)

    if not warning:
        freqs[freqs < 0] *= -1
        mask = freqs < limit * (LIGHT_SPEED * 100 * 2 * np.pi)
        freqs = freqs[mask]
        masses = masses[mask]
        normal_coord = normal_coord[:, :, mask]

    scales = 1e10 * np.sqrt(
        HBAR_J / (2 * masses * freqs * np.tanh(HBAR_EV * freqs / (2 * BOLTZ_EV * temp)))
    )

    args = [(geom, atomos, old, scales, normal_coord, warning) for _ in range(num_geoms)]

    # Use joblib to parallelize the geometry generation
    results = Parallel(n_jobs=-1, verbose=show_progress)(
        delayed(sample_single_geometry)(arg) for arg in args
    )

    structures = np.zeros((geom.shape[0], geom.shape[1], num_geoms))
    numbers = np.zeros((num_geoms, len(scales)))

    progress, rejected = 0, 0

    for j, output in enumerate(results):
        if output is None:
            continue
        structures[:, :, j] = output[0]
        numbers[j] = output[1].flatten()
        rejected += output[-1]
        progress += 1

    if show_progress:
        print(f"\nAccepted Geometries: {progress} Rejected Geometries: {rejected}")

    numbers = np.round(numbers, 4)
    return numbers, atomos, structures


###############################################################


##MAKES ENSEMBLE###############################################
def make_ensemble(freqlog, num_geoms, temperature, header, bottom):
    try:
        os.mkdir("Geometries")
    except FileExistsError:
        pass
    counter = nemo.tools.start_counter()
    print("\nGenerating geometries...\n")
    numbers, atomos, A = sample_geometries(freqlog, num_geoms, temperature,warning=False, show_progress=True)
    F, M = nemo.parser.pega_freq(freqlog)
    # convert numbers to dataframe
    numbers = pd.DataFrame(
        numbers, columns=[f"mode_{i+1}" for i in range(np.shape(numbers)[1])]
    )
    # check if file exists
    if os.path.isfile(f"Magnitudes_{temperature:.0f}K_.lx"):
        data = pd.read_csv(f"Magnitudes_{temperature:.0f}K_.lx")
        # get only columns with mode_ in the name
        data = data.filter(regex="mode_")
        # remove nan values
        data = data.dropna()
        # join data and numbers on axis 0
        numbers = pd.concat([data, numbers], axis=0, ignore_index=True)
    # concatenate frequencies and masses to numbers
    numbers = pd.concat(
        [pd.DataFrame(F, columns=["freq"]), pd.DataFrame(M, columns=["mass"]), numbers],
        axis=1,
    )
    numbers.to_csv(f"Magnitudes_{temperature:.0f}K_.lx", index=False)
    for n in range(np.shape(A)[2]):
        gfinal = A[:, :, n]
        nemo.tools.write_input(
            atomos,
            gfinal,
            header,
            bottom,
            f"Geometries/Geometry-{n+1+counter}-.com",
        )
    print("\n\nDone! Ready to run.")


################################################################

def check_dielectric(eps,nr):
    if eps < 1 or nr**2 > eps:
        nemo.parser.fatal_error("Dielectric constant must be higher than 1 and the refractive index squared must be lower than the static dielectric constant! Goodbye!")

def add_header(rem, num_ex, soc, static, refrac):
    rem = rem.lower().strip()
    method = rem.split()
    try:
        method = method[method.index('method')+1]
    except ValueError:
        method = 'td-dft'
    if method == 'eom-ccsd':
        header =(f"$rem\n"
            f"ee_singlets             {num_ex}\n"
            f"ee_triplets             {num_ex}\n"
            f"cc_trans_prop           2\n"
            f"calc_soc                {soc}\n"
            f"solvent_method          PCM\n"
            f"EOM_DAVIDSON_MAXVECTORS 300\n"
            f"EOM_DAVIDSON_MAX_ITER 300\n"
            f"$end\n"
            f"\n"
            f"$trans_prop\n"
            f"state_list\n"
            f"ref\n"
            f"ee_singlets 0 0\n"
            f"end_list\n"
            f"calc dipole linmom soc opdm_norm\n\n"
            f"state_list\n"
            f"ref\n"
            f"ee_triplets 0 0\n"
            f"end_list\n"
            f"calc dipole linmom soc opdm_norm\n\n"
            f"state_list\n"
            f"ee_singlets 0 0\n"
            f"ee_singlets 0 0\n"
            f"end_list\n"
            f"calc dipole linmom opdm_norm\n\n"
            f"state_list\n"
            f"ee_singlets 0 0\n"
            f"ee_triplets 0 0\n"
            f"end_list\n"
            f"calc dipole linmom soc opdm_norm\n\n"
            f"state_list\n"
            f"ee_triplets 0 0\n"
            f"ee_triplets 0 0\n"
            f"end_list\n"
            f"calc dipole linmom  opdm_norm\n"
            f"$end\n")
 
    else:
        header =(f"$rem\n"
            f"cis_n_roots             {num_ex}\n"
            f"cis_singlets            true\n"
            f"cis_triplets            true\n"
            f"calc_soc                {soc}\n"
            f"STS_MOM                 true\n"
            f"CIS_RELAXED_DENSITY     TRUE\n"
            f"solvent_method          PCM\n"
            f"MAX_CIS_CYCLES          200\n"
            f"MAX_SCF_CYCLES          200\n"
            f"$end\n")
    header += ("\n$pcm\n"
            "theory                  IEFPCM\n"
            "ChargeSeparation        Marcus\n"
            "StateSpecific           Perturb\n"
            "$end\n")
    header += (f"\n$solvent\n"
            f"Dielectric              {static}\n"
            f"OpticalDielectric       {refrac**2}\n"
            f"$end\n\n")   
    #remove $end from rem
    rem = rem.replace("$end", "")
    header = header.replace("$rem", rem)
    return header    

def setup_ensemble():
    freqlog = fetch_file("frequency", [".out", ".log"])
    print(f"\n\nFrequency log file: {freqlog}")
    with open(freqlog, "r", encoding="utf-8") as frequency_file:
        for line in frequency_file:
            if "Entering Gaussian System" in line:
                gaussian = True
            else:
                gaussian = False
            break
    if gaussian:
        print("You are using a Gaussian log file.")
        template = fetch_file("QChem template", [".in"])
        charge_multiplicity = lx.parser.get_cm(freqlog)
        rem, _, extra = nemo.parser.busca_input(template)
    else:
        template = fetch_file("QChem template", [".in"])
        rem, _, extra = nemo.parser.busca_input(template)
        _, charge_multiplicity, _ = nemo.parser.busca_input(freqlog)
    print(f"QChem template file: {template}")
    print("\nThe configurations to be used are:\n")
    rem += extra + "\n"
    print(rem)
    static = input("Solvent's static dielectric constant?\n")
    refrac = input("Solvent's refractive index?\n")
    try:
        static = float(static)
        refrac = float(refrac)
    except ValueError:
        nemo.parser.fatal_error(
            "Dielectric constant and refractive index must be numbers!"
        )
    check_dielectric(static,refrac)
    num_ex = input("How many excited states?\n")
    try:
        num_ex = int(num_ex)
    except ValueError:
        nemo.parser.fatal_error("This must be a number! Better luck next time!")
    abs_only = input("Are you interested in absorption spectra ONLY? (y or n)\n")
    if abs_only.lower() == "y":
        print(
            ("Ok, calculations will only be suitable for absorption "
             "or fluorescence spectrum simulations!\n")
        )
        header = add_header(rem, num_ex, 'false', static, refrac) 
    else:
        print(
            "Ok, calculations will be suitable for all spectra and ISC rate estimates!\n"
        )
        header = add_header(rem, num_ex, 'true', static, refrac)
    header += f"$molecule\n{charge_multiplicity}\n"
    num_geoms = int(input("How many geometries to be sampled?\n"))
    temperature = float(input("Temperature in Kelvin?\n"))
    if temperature <= 0:
        nemo.parser.fatal_error("Have you heard about absolute zero? Goodbye!")
    if gaussian:
        lx.tools.make_ensemble(freqlog, num_geoms, temperature, header, "$end\n")
    else:
        make_ensemble(freqlog, num_geoms, temperature, header, "$end\n")



##NORMALIZED GAUSSIAN##########################################
def gauss(x_value, mean, std):
    y_value = (1 / (np.sqrt(2 * np.pi) * std)) * np.exp(-0.5 * ((x_value - mean) / std) ** 2)
    return y_value


###############################################################


##COMPUTES AVG TRANSITION DIPOLE MOMENT########################
def calc_tdm(osc, delta_e, pesos):
    # Energy terms converted to J
    term = E_CHARGE * (HBAR_J**2) / delta_e
    dipoles = np.sqrt(3 * term * osc / (2 * MASS_E))
    # Conversion in au
    dipoles *= 1.179474389e29
    return np.average(dipoles, weights=pesos)


###############################################################


##PREVENTS OVERWRITING#########################################
def naming(arquivo):
    new_arquivo = arquivo
    if arquivo in os.listdir("."):
        duplo = True
        vers = 2
        while duplo:
            new_arquivo = str(vers) + arquivo
            if new_arquivo in os.listdir("."):
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
    except ValueError:
        nemo.parser.fatal_error("It must be S or T and an integer! Goodbye!")
    if estados[0].upper() != "S" and estados[0].upper() != "T":
        nemo.parser.fatal_error("It must be S or T and an integer! Goodbye!")
    return estados.upper()


###############################################################


def get_alpha(eps):
    return max((eps - 1) / (eps + 1),1e-10)


##FETCHES  FILES###############################################
def fetch_file(frase, ends):
    for file in [i for i in os.listdir(".")]:
        for end in ends:
            if end in file:
                return file
    nemo.parser.fatal_error(f"No {frase} file found. Goodbye!")
###############################################################


##RUNS TASK MANAGER############################################
def batch():
    script = fetch_file("batch.sh", ["batch.sh"])
    limite = input("Maximum number of jobs to be submitted simultaneously?\n")
    nproc = input("Number of threads for each individual calculation?\n")
    num = input("Number of calculations in each job?\n")
    try:
        limite = int(limite)
        int(nproc)
        int(num)
    except ValueError:
        nemo.parser.fatal_error("It must be an integer. Goodbye!")

    with open("limit.lx", "w",encoding="utf-8") as limit_file:
        limit_file.write(str(limite))
    Popen(
        ["nohup", "nemo_batch_run", script, nproc, num, "&"]
    )

###############################################################


##FINDS SUITABLE VALUE FOR STD#################################
def detect_sigma():
    try:
        files = [i for i in os.listdir(".") if "Magnitudes" in i and ".lx" in i]
        file = files[0]
        temp = float(file.split("_")[1].strip("K"))
        sigma = np.round(BOLTZ_EV * temp, 3)
    except IndexError:
        print(
            "WARNING: Magnitudes.lx file is absent! Temperature is being set to 300 K!"
        )
        sigma = 0.026
    return sigma


###############################################################

def fetch_nr(file):
    refractive_index = None
    epsilon = None
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if "opticaldielectric" in line.lower():
                refractive_index = np.sqrt(float(line.split()[1]))
            elif "dielectric" in line.lower() and "optical" not in line.lower():
                epsilon = float(line.split()[1])
            if refractive_index is not None and epsilon is not None:
                return epsilon, refractive_index
    return epsilon, refractive_index        

def susceptibility_check(file, tuning=False):
    # Fetch energy levels and other data
    es, et, _, _, _, ss_s, ss_t, _ = nemo.parser.pega_energias(file)
    _, nr = fetch_nr(file)
    
    # Calculate alpha and susceptibility chi values
    alpha = (nr**2 - 1) / (nr**2 + 1)
    chi_s = ss_s / alpha
    chi_t = ss_t / alpha
    
    if tuning:
        return es[0], chi_s[0]
    else:
        chi_symbol = '\u03C7(eV)'
        # Print header with aligned columns
        print(fr"{'State':<6} {'E_vac(eV)':<12} {chi_symbol:<10}")
        
        # Print singlet states
        for i, (e, chi) in enumerate(zip(es, chi_s), start=1):
            print(f"S{i:<5} {e:<12.3f} {chi:<10.3f}")

        # Print triplet states
        for i, (e, chi) in enumerate(zip(et, chi_t), start=1):
            print(f"T{i:<5} {e:<12.3f} {chi:<10.3f}")



##FETCHES REFRACTIVE INDEX#####################################
def get_nr():
    coms = [
        file
        for file in os.listdir("Geometries")
        if "Geometr" in file and ".com" in file
    ]
    epsilon, refractive_index = fetch_nr('Geometries/'+coms[0])
    if epsilon is None:
        return 1, 1
    else:
        return epsilon, refractive_index


###############################################################


##QUERY FUNCTION###############################################
def default(default_value, frase):
    new_value = input(frase)
    if new_value == "":
        return default_value
    else:
        return new_value


###############################################################


##STOP SUBMISSION OF JOBS######################################
def abort_batch():
    choice = input(
        "Are you sure you want to prevent new jobs from being submitted? y or n?\n"
    )
    if choice == "y":
        try:
            os.remove("limit.lx")
            print("Done!")
        except FileNotFoundError:
            print("Could not find the files. Maybe you are in the wrong folder.")
    else:
        print("OK, nevermind")


###############################################################


class Watcher:
    def __init__(self, folder, key="Geometr"):
        self.folder = folder
        self.key = key
        self.files = [i[:-4] for i in os.listdir(folder) if i.endswith('.com') and key in i]
        self.files = sorted(self.files, key=lambda pair: float(pair.split("-")[1]))
        self.number_inputs = len(self.files)
        self.done = []
        self.license_error = []
        self.error = []
        self.running = []
        self.running_batches = 0

    def check(self):
        list_to_check = self.files.copy()
        for input_file in list_to_check:
            try:
                with open(self.folder + "/" + input_file + ".log", "r",encoding="utf-8") as log_file:
                    for line in log_file:
                        if "Have a nice day" in line:
                            self.done.append(input_file)
                            del self.files[self.files.index(input_file)]
                            break
                        elif "fatal error" in line:
                            self.error.append(input_file)
                            del self.files[self.files.index(input_file)]
                            break
                        elif "failed standard" in line:
                            self.license_error.append(input_file)
                            del self.files[self.files.index(input_file)]
                            break
            except FileNotFoundError:
                pass

    def report(self):
        self.check()
        print('\n\n')
        print(f'There are {len(self.done)} successfully completed calculations out of {self.number_inputs} inputs.')
        print(f'{100 * len(self.done) / self.number_inputs:.1f}% of the calculations have been run.')
        if len(self.error) > 0:
            print(f"There are {len(self.error)} failed jobs.")
            print('These are: ', self.error)
        if len(self.license_error) > 0:
            print(f"There are {len(self.license_error)} failed jobs due to license error.")
            print('These are: ', self.license_error)

    def limit(self):
        if self.key == "Geometr":
            try:
                return np.loadtxt("../limit.lx",encoding='utf-8')
            except (OSError,FileNotFoundError):
                sys.exit()
        else:
            try:
                return np.loadtxt("limit.lx",encoding='utf-8')
            except (OSError,FileNotFoundError):
                sys.exit()

    def keep_going(self,num):
        if len(self.running) / num < self.limit():
            return False
        return True

    def clean_failed(self):
        for failed in self.error + self.license_error:
            os.remove(self.folder + "/" + failed + ".log")
        self.files += self.error + self.license_error
        self.error = []
        self.license_error = []

    def run(self, batch_file, nproc, num):
        total_threads = int(nproc) * int(num)
        self.check()
        self.clean_failed()
        inputs = self.files.copy()
        while len(inputs) > 0:
            next_inputs = inputs[:int(num)]
            num_proc = int(total_threads / len(next_inputs))
            command = ''
            for input_file in next_inputs:
                command += f"qchem -nt {num_proc} {input_file}.com {input_file}.log &\n"
                command += 'sleep 5\n'
                self.running.append(input_file)
                inputs.remove(input_file)
            command += "wait"
            with open(f"cmd_{self.running_batches}_.sh", "w",encoding='utf-8') as cmd:
                cmd.write(command)
            sub = subprocess.call(["bash", batch_file, f"cmd_{self.running_batches}_.sh"])
            self.running_batches += 1
            keep = self.keep_going(num)
            while keep:
                time.sleep(20)
                self.check()
                concluded = self.done + self.error + self.license_error
                self.running = [elem for elem in self.running if elem not in concluded]
                keep = self.keep_going(num)    

    def hold_watch(self):
        while len(self.files) > 0:
            time.sleep(20)
            self.check()

##CHECKS PROGRESS##############################################
def andamento():
    the_watcher = Watcher('Geometries')
    the_watcher.report()

###############################################################

def check_for_updates(package_name):
    try:
        # Get the currently installed version
        installed_version = pkg_resources.get_distribution(package_name).version
        
        # Fetch the latest version from PyPI
        response = requests.get(f'https://pypi.org/pypi/{package_name}/json')
        response.raise_for_status()
        latest_version = response.json()['info']['version']

        # Compare versions
        if installed_version != latest_version:
            print(f"ATTENTION: Update available! {package_name} {installed_version} -> {latest_version}")
            print("Run `pip install --upgrade {}` to update.".format(package_name))

    except Exception as e:
        print(f"An error occurred while checking for updates: {e}")

##RUNS W TUNING################################################
def empirical_tuning():
    geomlog = fetch_file("input or log", [".com", ".log"])
    rem, _, extra = nemo.parser.busca_input(geomlog)
    print(f"QChem template file: {geomlog}")
    rem += extra + "\n"
    #iterate over lines of rem
    for line in rem.split("\n"):
        if "method" in line.lower() or 'exchange' in line.lower():
            functional = line.split()[-1]
        if "basis" in line.lower():
            basis = line.split()[-1]
        if  'mem_total' in line.lower():
            mem = line.split()[-1]   
    omega1 = "0.1"
    passo = "0.025"
    relax = 'yes'
    print(f"Total memory: {mem}")
    print(f"Functional: {functional}")
    print(f"Basis: {basis}")
    print(f"Initial Omega: {omega1} bohr^-1")
    print(f"Step: {passo} bohr^-1")
    print(f'Optimize at each step? {relax}')
    change = input("Are you satisfied with these parameters? y or n?\n")
    if change == "n":
        omega1 = default(
            omega1,
            f"Initial omega is {omega1} bohr^-1. If ok, Enter. Otherwise, type it.\n",
        )
        passo = default(
            passo,
            f"Initial step is {passo} bohr^-1. If ok, Enter. Otherwise, type it.\n",
        )
        relax = default(
            relax,
            f"Optimize at each step? {relax}. If ok, Enter. Otherwise, type n.\n",
        )
    script = fetch_file("batch script", ["batch.sh"])
    nproc = input("Number of threads for each calculation\n")
    e_exp = input("Experimental vacuum energy and uncertainty in eV? (space separated)\n")
    chi_exp = input("Experimental susceptibility and uncertainty in eV? (space separated)?\n")
    
    with open("limit.lx", "w",encoding="utf-8") as f:
        f.write("10")
    subprocess.Popen(
        [
            "nohup",
            "nemo_tuning",
            geomlog,
            functional,
            basis,
            nproc,
            omega1,
            passo,
            relax,
            script,
            e_exp,
            chi_exp,
            mem,
            "&",
        ]
    )

###############################################################
