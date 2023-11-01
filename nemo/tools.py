#!/usr/bin/env python3
import os
import subprocess
import sys
import time
from subprocess import Popen
import numpy as np
from scipy.stats import norm
import lx.tools
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


##SAMPLES GEOMETRIES###########################################
def sample_geometries(freqlog, num_geoms, temperature, limit=np.inf):
    geometry, atomos = nemo.parser.pega_geom(freqlog)
    freqs, masses = nemo.parser.pega_freq(freqlog)
    freqs[freqs < 0] *= -1
    normal_coordinates = nemo.parser.pega_modos(freqlog)
    mask = freqs < limit * (LIGHT_SPEED * 100 * 2 * np.pi)
    freqs = freqs[mask]
    normal_coordinates = normal_coordinates[:, mask]
    num_atom = np.shape(geometry)[0]
    structure = np.zeros((3 * num_atom, num_geoms))
    for i, freq in enumerate(freqs):
        scale = np.sqrt(HBAR_J / (2 * masses[i] * freq *
                             np.tanh(HBAR_EV * freq / (2 * BOLTZ_EV * temperature))))
        normal = norm(scale=scale, loc=0)
        # Displacements in  Å
        displacement = normal.rvs(size=num_geoms) * 1e10
        try:
            numbers = np.hstack((numbers, displacement[:, np.newaxis]))
        # If it is the first iteration, numbers is not defined yet
        except NameError:
            numbers = displacement[:, np.newaxis]
        structure += np.outer(normal_coordinates[:, i], displacement)
    for i in range(np.shape(structure)[1]):
        reshaped_structure = np.reshape(structure[:, i], (num_atom, 3))
        try:
            final_geometry = np.hstack((final_geometry, reshaped_structure + geometry))
        except NameError:
            final_geometry = reshaped_structure + geometry
    numbers = np.round(numbers, 4)
    return numbers, atomos, final_geometry


###############################################################


##MAKES ENSEMBLE###############################################
def make_ensemble(freqlog, num_geoms, temperature, header, bottom):
    try:
        os.mkdir("Geometries")
    except FileExistsError:
        pass
    counter = start_counter()
    print("\nGenerating geometries...\n")
    numbers, atomos, final_geometry = sample_geometries(freqlog, num_geoms, temperature)
    with open(f"Magnitudes_{temperature:.0f}K_.lx", "a", encoding="utf-8") as file:
        np.savetxt(file, numbers, delimiter="\t", fmt="%s")
    for i in range(0, np.shape(final_geometry)[1], 3):
        gfinal = final_geometry[:, i : i + 3]
        write_input(
            atomos,
            gfinal,
            header,
            bottom,
            "Geometries/Geometry-" + str((i + 3) // 3 + counter) + "-.com",
        )
        progress = 100 * ((i + 3) // 3) / num_geoms
        text = f"{progress:2.1f}%"
        print(" ", text, "of the geometries done.", end="\r", flush=True)
    print("\n\nDone! Ready to run.")


################################################################

def check_dielectric(eps,nr):
    if eps <= 1 or nr**2 > eps:
        nemo.parser.fatal_error("Dielectric constant must be higher than 1 and the refractive index squared must be lower than the static dielectric constant! Goodbye!")

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
        charge_multiplicity = lx.tools.get_cm(freqlog)
        rem, _, extra = nemo.parser.busca_input(template)
    else:
        template = fetch_file("QChem template", [".in"])
        rem, _, extra = nemo.parser.busca_input(template)
        _, charge_multiplicity, _ = nemo.parser.busca_input(freqlog)
    print(f"QChem template file: {template}")
    print("\nThe configurations to be used are:\n")
    rem += extra + "\n"
    print(rem)
    rem += ("\n$pcm\n"
            "theory                  IEFPCM\n"
            "ChargeSeparation        Marcus\n"
            "StateSpecific           Perturb\n"
            "$end\n")
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
    rem += (f"\n$solvent\n"
            f"Dielectric              {static}\n"
            f"OpticalDielectric       {refrac**2}\n"
            f"$end\n\n")
    num_ex = input("How many excited states?\n")
    try:
        num_ex = int(num_ex)
    except ValueError:
        nemo.parser.fatal_error("This must be a number! Better luck next time!")
    header =(f"$rem\n"
            f"cis_n_roots             {num_ex}\n"
            f"cis_singlets            true\n"
            f"cis_triplets            true\n"
            f"STS_MOM                 true\n"
            f"CIS_RELAXED_DENSITY     TRUE\n"
            f"solvent_method          PCM\n"
            f"MAX_CIS_CYCLES          200\n"
            f"MAX_SCF_CYCLES          200\n")
    abs_only = input("Are you interested in absorption spectra ONLY? (y or n)\n")
    if abs_only.lower() == "y":
        print(
            ("Ok, calculations will only be suitable for absorption "
             "or fluorescence spectrum simulations!\n")
        )
        header += "calc_soc                false\n"
    else:
        print(
            "Ok, calculations will be suitable for all spectra and ISC rate estimates!\n"
        )
        header += "calc_soc             true\n"
    header = rem.replace("$rem", header)
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
    return (eps - 1) / (eps + 1)


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


##FETCHES REFRACTIVE INDEX#####################################
def get_nr():
    refractive_index = 1
    epsilon = 1
    coms = [
        file
        for file in os.listdir("Geometries")
        if "Geometr" in file and ".com" in file
    ]
    with open("Geometries/" + coms[0], "r", encoding="utf-8") as com_file:
        for line in com_file:
            if "opticaldielectric" in line.lower():
                refractive_index = np.sqrt(float(line.split()[1]))
            elif "dielectric" in line.lower() and "optical" not in line.lower():
                epsilon = float(line.split()[1])
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
    def __init__(self, folder):
        self.folder = folder
        self.files = [i[:-4] for i in os.listdir(folder) if i.endswith('.com') and "Geometr" in i]
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
        try:
            return np.loadtxt("../limit.lx",encoding='utf-8')
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

##CHECKS PROGRESS##############################################
def andamento():
    the_watcher = Watcher('Geometries')
    the_watcher.report()

###############################################################

##CALCULATES FLUORESCENCE LIFETIME IN S########################
def calc_emi_rate(energy, diff_rate, diff_rate_error):
    # Integrates the emission spectrum
    int_emi = np.trapz(diff_rate, energy)
    taxa = (1 / HBAR_EV) * int_emi
    error = (1 / HBAR_EV) * np.sqrt(np.trapz((diff_rate_error**2), energy))
    return taxa, error


###############################################################
