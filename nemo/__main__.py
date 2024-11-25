#!/usr/bin/env python3
"""
This module contains the main function for the NEMO program, 
which is a tool for simulating photophysical processes in organic materials.
"""
import sys
import argparse
import lx.tools
import nemo.tools
import nemo.parser
from nemo import __version__ as nemo_version
from nemo.analysis import rates, export_results, gather_data, absorption
# pylint: disable=unbalanced-tuple-unpacking

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
    print(f"Version: {nemo_version.__version__}\n")
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
    print('\n')
    nemo.tools.check_for_updates('nemophoto')
    operation = input()
    if operation == "1":
        nemo.tools.setup_ensemble()
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
                nemo.parser.fatal_error(
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
                nemo.parser.fatal_error(
                    "Dielectric constant and refractive index must be numbers. Bye!"
                )
        nemo.tools.check_dielectric(epsilon,refractive_index)        
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
        nemo.tools.empirical_tuning()
        #nemo.parser.fatal_error("It must be one of the options... Goodbye!")


def main():
    """
    This function is the entry point for the NEMO program. It takes command line arguments
    and either prints out the geometry of a molecule or launches the program's interface.
    """
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Nemo script with susceptibility check option.")
    
    # Add the `-c` flag for susceptibility check with a file argument
    parser.add_argument('-c', '--check', type=str, help="Run susceptibility check on the specified file.")

    parser.add_argument('-g', '--geom', type=str, help="Gets geometry from a log file.")    
    # Parse arguments
    args = parser.parse_args()
    
    # If `-c` is provided, call the susceptibility_check function
    if args.check:
        nemo.tools.susceptibility_check(args.check)
        sys.exit(0)

    elif args.geom:
        try:
            geometry, atomos = nemo.parser.pega_geom(args.geom)
            print(len(atomos))
            print("")
            for i, atomo in enumerate(atomos):
                print(
                    f"{atomo:2s}  {geometry[i, 0]:.7f}  {geometry[i, 1]:.7f}  {geometry[i, 2]:.7f}"
                )
            sys.exit(0)

        except FileNotFoundError:
            interface()
    else:
        interface()


if __name__ == "__main__":
    sys.exit(main())
