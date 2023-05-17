# NEMO - Photophysics with the Nuclear Ensemble Method

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![license](https://img.shields.io/github/license/LeonardoESousa/NEMO?style=plastic)]()
[![down](https://img.shields.io/pypi/dm/nemophoto)]()
[![maint](https://img.shields.io/maintenance/yes/2023)]()
[![commit](https://img.shields.io/github/last-commit/LeonardoESousa/NEMO?style=plastic)]()


 Fluorescence, phosphorescence and intersystem crossing (ISC) rate calculations. Absorption, fluorescence and phosphorescence spectrum simulations. Förster radius and singlet exciton diffusion length estimates. Interfaces with the QChem package. 

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

1.  Photophysics with TD(A)-DFT:
    - Calculates fluorescence and phosphorescence rates from an excited state.
    - Calculates ISC rates from a given singlet state to several triplet states and vice-versa.
    - Absorption, Fluorescence and Phosphorescence spectrum simulations.
    - Rate calculations and spectra include vibrational contributions and state specific solvation effects.
2.  Exciton properties:   
    - Calculates the Förster radius for transfers between two molecules of equal or different type.
    - Estimates singlet exciton diffusion lengths.

    
## What is necessary to use it?

 -  The program requires that the QChem quantum chemistry software be installed, since it interfaces with it.

 -  The first step for running spectrum calculations is providing either a QChem or Gaussian log file for a frequency calculation in the S0, S1 or T1 state, if the goal is computing an absorption, fluorescence or phosphorescence spectrum, respectively. All frequencies must be real.  

 -  For calculating ISC rates from the Sn state to different triplet states, a QChem/Gaussian frequency calculation at the Sn state must be provided.

 -  Similarly, for reverse ISC rates from the Tn state to different triplet states, a QChem/Gaussian frequency calculation at the Tn state must be provided.
 
 -  To obtain the estimates of Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths, it is necessary to first perform both absorption and fluorescence spectra calculations for the molecule of interest.

## How to install it?

Run:

`pip install nemophoto`

To get the latest commit, run:

`pip install git+https://github.com/LeonardoESousa/NEMO`

Alternatively, clone the repository to your computer. Inside the NEMO folder, run:

`pip install .`

Once installed, you should be able to run the program from any folder by just using the `nemo` command.

## How to use it?

Here is a quick guide on how to use the software. For a detailed tutorial, click [here](https://github.com/LeonardoESousa/NEMO/blob/main/Tutorial/Tutorial.md).

1. Initial steps:
    - Create a folder for your project. Add the log file for the frequency calculation to your folder. If the frequency calculation was run with Gaussian, you must also provide a QChem input file containing some settings you wish to apply in the ensemble calculations (e.g. functional, basis set, omega value, charge and multiplicity etc). An example of such file (td.in) is provided [here](https://github.com/LeonardoESousa/NEMO/tree/main/batch_examples).
    - A frequency calculation in the S0 state is suitable for computing an absorption spectrum. For fluorescence spectra and/or ISC rates calculations from Sn states to triplet states, a Sn frequency calculation is expected. Finally, for phosphorescence spectra and/or rISC rates calculations from Tn states to singlet states, a Tn frequency calculation is expected.  
    - Run the `nemo` command. Choose option 1 and follow the instructions to select the parameters of the calculation. This includes the dielectric constant and refractive index of the medium. This information will be used to obtain state-specific solvent corrections to the TD(A)-DFT energies.  
    - Add a bash script file to the folder. This file depends on which batch system you use. Examples of this file for users of slurm or task spooler (ts) are presented [here](https://github.com/LeonardoESousa/NEMO/tree/main/batch_examples)).
    - Run the `nemo` command again, choose option 2 and follow the instructions. Alternatively, just run all calculations created in the Geometries folder. Once the calculations are running, you may use option 3 to check the progress or option 4 to abort.

2. For absorption spectrum simulations:
    - Once all calculations from step 1 are done, run the `nemo` command and choose option 5. Follow the instructions to set the parameters and the spectrum will be generated. 

3. For photophysical rates:
    - Once all calculations from step 1 are done, run the `nemo` command and choose option 6. Follow the instructions to set the parameters. Three files will be generated: an Ensemble file, with data from the ensemble of geometries; a differential_rate file, with the emission spectrum; a rates file, with all available rates (-> denote radiative transitions and ~> denote ISC transitions). 

4. IMPORTANT:
    - You may choose to calculate spectra and rates with different solvent dielectric constant and refractive index than the ones selected as input in Option 1. To do so, NEMO will resort to the extrapolation procedure described in ADD PAPER to adjust the results to the new solvent. 

5. For exciton properties:
    - For exciton properties, you must first calculate the fluorescence and absorption spectra of the donor and acceptor molecules of interest to you. Copy both spectra to a folder and inside this folder run the `nemo` command. Choose option 7. Follow the instructions to set the calculation parameters. A file will be generated with all the information. Importantly, diffusion length estimates are only sensible if donor and acceptor molecules are of the same kind. These estimations follow from the procedures described in: "de Sousa, L. E., Bueno, F. T., e Silva, G. M., da Silva Filho, D. A., & de Oliveira Neto, P. H. (2019). Fast predictions of exciton diffusion length in organic materials. Journal of Materials Chemistry C, 7(14), 4066-4071." 


For better visualization of results, consider using  [Nemoview](https://github.com/LeonardoESousa/nemoview). 
