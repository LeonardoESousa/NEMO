# NEMO - Photophysics with the Nuclear Ensemble Method

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![license](https://img.shields.io/github/license/LeonardoESousa/NEMO?style=plastic)]()
[![down](https://img.shields.io/pypi/dm/nemophoto)]()
[![maint](https://img.shields.io/maintenance/yes/2021)]()
[![commit](https://img.shields.io/github/last-commit/LeonardoESousa/NEMO?style=plastic)]()


Absorption, fluorescence and phosphorescence spectrum simulations. Fluorescence, phosphorescence and intersystem crossing (ISC) rate calculations. Foerster radius and singlet exciton diffusion length estimates. Interfaces with the QChem package. 

Table of Contents
=================
<!--ts-->
* [Cite as:](#cite-as)
* [What does this program do?](#what-does-this-program-do)
* [What is necessary to use it?](#what-is-necessary-to-use-it)
* [How to install it?](#how-to-install-it)
* [How to use it?](#how-to-use-it)
   
<!--te-->

## Cite as:

> Leonardo Evaristo de Sousa and Piotr de Silva
Journal of Chemical Theory and Computation 2021 17 (9), 5816-5824
DOI: 10.1021/acs.jctc.1c00476


## What does this program do?

1.  Spectrum simulation:
    - Absorption, Fluorescence and Phosphorescence spectrum simulations using TD(A)-DFT.
    - Calculations include vibrational contributions to the spectra. Optionally, it may also include solvent effects with PCM.
    - Calculates fluorescence and phosphorescence rates.
2.  Intersystem Crossing Rates:
    - Calculates ISC rates from a given singlet state to several triplet states and vice-versa 
3.  Exciton properties:   
    - Calculates the Förster radius for transfers between two molecules of equal or different type.
    - Calculates fluorescence lifetimes.
    - Calculates singlet exciton diffusion lengths.
4.  Extra features:
    - Extracts last geometry from QChem log file.
    
## What is necessary to use it?

 -  The program requires that the QChem quantum chemistry software be installed, as it interfaces with it.

 -  The first step for running spectrum calculations is providing a QChem log file for a frequency calculation in the S0, S1 or T1 state, if the goal is computing an absorption, fluorescence or phosphorescence spectrum, respectively. All frequencies must be real.  

 -  For calculating ISC rates from the Sn state to different triplet states, a QChem frequency calculation at the Sn state must be provided. In addition, reorganization energies for the Sn -> Tm (m = 1,2,3...) transfers of interest must be included.

 -  Similarly, for reverse ISC rates from the Tn state to different triplet states, a QChem frequency calculation at the Tn state must be provided. In addition, reorganization energies for the Tn -> Sm (m = 1,2,3...) transfers of interest must be included. 
 
 -  To obtain the estimates of Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths, it is necessary to first perform both absorption and fluorescence spectra calculations for the molecule of interest.

## How to install it?

Run:

`pip install nemophoto`

Alternatively, clone the repository to your computer. Inside the NEMO folder, run:

`pip install .`

Once installed, you should be able to run the program from any folder by just using the `nemo` command.

## How to use it?

1. Initial steps:
    - Create a folder for your project. Add the log file for the frequency calculation to your folder. A frequency calculation in the S0 state is suitable for computing an absorption spectrum. For fluorescence spectra and/or ISC rates calculations from Sn states to triplet states, a Sn frequency calculation is expected. Finally, for phosphorescence spectra and/or rISC rates calculations from Tn states to singlet states, a Tn frequency calculation is expected.   
    - Run the `nemo` command. Choose option 1 and follow the instructions to select the parameters of the calculation.
    - Add a bash script file to the folder. This file depends on which batch system you use. Examples of this file for users of slurm or task spooler (ts) are presented in [here](https://github.com/LeonardoESousa/NEMO/tree/main/batch_examples)).
    - Run the `nemo` command again, choose option 2 and follow the instructions. Once the calculations are running, you may use option 3 to check the progress or option 4 to abort.

2. For spectrum simulations and fluorescence or phosphorescence rate calculations:
    - Once all calculations from step 1 are done, run the `nemo` command and choose option 5. Follow the instructions to set the parameters and the spectrum will be generated. The corresponding fluorescence or phosphorescence rate will be shown in the spectrum file.

3. For ISC rates:
    - Create a file named lambdas.txt inside your project folder. In this file you should write the reorganization energies in eV for the processes you are interested. Suppose you are calculating ISC rates from the T1 state to several singlet states. Each line in the file should correspond to the reorganziation energy for the T1 -> Sn transfer (n = 1,2,3...). The first line will correspond to T1 -> S1, the second to T1 -> S2, and so on. Include as many as you want. So long as there are enough excited states in the ensemble calculations, they will be used to calculate rates. Check the example [here](https://github.com/LeonardoESousa/NEMO/tree/main/batch_examples)).
    - Once all calculations from step 1 are done and the reorganization energies have been set, run the `nemo` command and choose option 6. Follow the instructions to set the parameters and the intersystem crossing rates will be generated.

4. For exciton properties:
    - For exciton properties, you must first calculate the fluorescence and absorption spectra of the donor and acceptor molecules of interest to you. Copy both spectra to a folder and inside this folder run the `nemo` command. Choose option 6. Follow the instructions to set the calculation parameters. A file will be generated with all the information. Importantly, diffusion length estimates are only sensible if donor and acceptor molecules are of the same kind. These estimations follow from the procedures described in: "de Sousa, L. E., Bueno, F. T., e Silva, G. M., da Silva Filho, D. A., & de Oliveira Neto, P. H. (2019). Fast predictions of exciton diffusion length in organic materials. Journal of Materials Chemistry C, 7(14), 4066-4071." 


