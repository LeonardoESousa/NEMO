#!/usr/bin/env python3
"""
This module contains the main function for the NEMO program, 
which is a tool for simulating photophysical processes in organic materials.
"""
import sys
import lx.tools
import nemo.tools
from nemo.analysis import rates, export_results, gather_data, absorption

def interface():
    """
    This function displays a menu of options for the user to choose from
    and prompts the user for input.
    The user's input is then used to call the appropriate function
    to perform the selected operation.
    """
    print(
        """
    ███▄▄▄▄      ▄████████    ▄▄▄▄███▄▄▄▄    ▄██████▄ 
    ███▀▀▀██▄   ███    ███  ▄██▀▀▀███▀▀▀██▄ ███    ███
    ███   ███   ███    █▀   ███   ███   ███ ███    ███
    ███   ███  ▄███▄▄▄      ███   ███   ███ ███    ███
    ███   ███ ▀▀███▀▀▀      ███   ███   ███ ███    ███
    ███   ███   ███    █▄   ███   ███   ███ ███    ███
    ███   ███   ███    ███  ███   ███   ███ ███    ███
     ▀█   █▀    ██████████   ▀█   ███   █▀   ▀██████▀ 
    ------------------Photophysics--------------------
    \n"""
    )
    print("Choose your option:\n")
    print("ENSEMBLE SETUP:")
    print("\t1 - Generate the inputs for the nuclear ensemble calculation")
    print("\t2 - Run the ensemble calculations")
    print("\t3 - Check the progress of the calculations")
    print("\t4 - Abort my calculations")
    print("ABSORPTION:")
    print("\t5 - Generate the absorption spectrum")
    print("EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):")
    print("\t6 - Estimate rates and compute emission spectrum")
    print("ENSEMBLE DATA:")
    print("\t7 - Gather ensemble data only")
    print("EXCITON ANALYSIS:")
    print(
        "\t8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths"
    )
    print("EXTRA FEATURES:")
    print(
        "\t9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)"
    )
    operation = input()
    if operation == "1":
        freqlog = nemo.tools.fetch_file("frequency", [".out", ".log"])
        print(f"\n\nFrequency log file: {freqlog}")
        with open(freqlog, "r", encoding="utf-8") as frequency_file:
            for line in frequency_file:
                if "Entering Gaussian System" in line:
                    gauss = True
                else:
                    gauss = False
                break
        if gauss:
            print("You are using a Gaussian log file.")
            template = nemo.tools.fetch_file("QChem template", [".in"])
            charge_multiplicity = lx.tools.get_cm(freqlog)
            rem, _, extra = nemo.tools.busca_input(template)
        else:
            template = nemo.tools.fetch_file("QChem template", [".in"])
            rem, _, extra = nemo.tools.busca_input(template)
            _, charge_multiplicity, _ = nemo.tools.busca_input(freqlog)
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
            nemo.tools.fatal_error(
                "Dielectric constant and refractive index must be numbers!"
            )
        rem += (f"\n$solvent\n"
                f"Dielectric              {static}\n"
                f"OpticalDielectric       {refrac**2}\n"
                f"$end\n\n")
        num_ex = input("How many excited states?\n")
        go_ahead = False
        while not go_ahead:
            try:
                num_ex = int(num_ex)
                go_ahead = True
            except ValueError:
                print("This must be a number! Try again!\n")
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
            nemo.tools.fatal_error("Have you heard about absolute zero? Goodbye!")
        if gauss:
            lx.tools.make_ensemble(freqlog, num_geoms, temperature, header, "$end\n")
        else:
            nemo.tools.make_ensemble(freqlog, num_geoms, temperature, header, "$end\n")
    elif operation == "2":
        nemo.tools.batch()
    elif operation == "3":
        nemo.tools.andamento()
    elif operation == "4":
        nemo.tools.abort_batch()
    elif operation == "5":
        epsilon, refractive_index = nemo.tools.get_nr()
        print("The spectrum will be run with the following parameters:\n")
        print(f"Solvent dielectric constant: {epsilon:.3f}")
        print(f"Solvent refractive index: {refractive_index:.3f}\n")
        change = input("Are you satisfied with these parameters? y or n?\n")
        if change.lower() == "n":
            epsilon = nemo.tools.default(
                epsilon,
                (f"Solvent dielectric constant is {epsilon:.3f}. "
                f"If ok, Enter. Otherwise, type value.\n"),
            )
            refractive_index = nemo.tools.default(
                refractive_index,
                (f"Refractive index is {refractive_index:.3f}. "
                 "If ok, Enter. Otherwise, type value.\n")
            )
            try:
                epsilon = float(epsilon)
                refractive_index = float(refractive_index)
            except ValueError:
                nemo.tools.fatal_error(
                    "Dielectric constant and refractive index must be numbers. Bye!"
                )
        state = input(
            ("What is the initial state (S0, S1, T1, S2 ...)? "
             "Accepts comma separated values Ex: T1,T2\n")
        )
        states = state.split(",")
        for state in states:
            absorption(state, [epsilon, refractive_index], save=True)
    elif operation == "6":
        epsilon, refractive_index = nemo.tools.get_nr()
        print("The rates will be calculated with the following parameters:\n")
        print(f"Solvent dielectric constant: {epsilon:.3f}")
        print(f"Solvent refractive index: {refractive_index:.3f}\n")
        change = input("Are you satisfied with these parameters? y or n?\n")
        if change.lower() == "n":
            epsilon = nemo.tools.default(
                epsilon,
                (f"Solvent dielectric constant is {epsilon:.3f}. "
                 "If ok, Enter. Otherwise, type value.\n")
            )
            refractive_index = nemo.tools.default(
                refractive_index,
                (f"Refractive index is {refractive_index:.3f}. "
                 f"If ok, Enter. Otherwise, type value.\n")
            )
            try:
                epsilon = float(epsilon)
                refractive_index = float(refractive_index)
            except ValueError:
                nemo.tools.fatal_error(
                    "Dielectric constant and refractive index must be numbers. Bye!"
                )
        state = input(
            "What is the initial state (S1, T1, S2 ...)? Accepts comma separated values Ex: T1,T2\n"
        )

        states = state.split(",")
        for state in states:
            res, emi = rates(state, [epsilon, refractive_index])
            export_results(res, emi, [epsilon, refractive_index])
    elif operation == "7":
        state = input(
            ("What is the initial state (S0, S1, T1, S2 ...)? "
             "Accepts comma separated values Ex: T1,T2\n")
        )
        states = state.split(",")
        for state in states:
            gather_data(state, save=True)
    elif operation == "8":
        lx.tools.ld()
    elif operation == "9":
        lx.tools.omega_tuning()
    else:
        nemo.tools.fatal_error("It must be one of the options... Goodbye!")


def main():
    """
    This function is the entry point for the NEMO program. It takes command line arguments
    and either prints out the geometry of a molecule or launches the program's interface.
    """
    if len(sys.argv) > 1:
        try:
            freqlog = sys.argv[1]
            geometry, atomos = nemo.tools.pega_geom(freqlog)
            for i, atomo in enumerate(atomos):
                print(
                    f"{atomo:2s}  {geometry[i, 0]:.8f}  {geometry[i, 1]:.8f}  {geometry[i, 2]:.8f}"
                )
        except IndexError:
            interface()
    else:
        interface()


if __name__ == "__main__":
    sys.exit(main())
