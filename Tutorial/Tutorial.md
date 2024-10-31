# **NEMO** Tutorial Version 1.2.0

<img src="logo.png" alt="Alt Text" width="100%">


# Table of Contents

- [Introduction](#introduction)
- [Installation of the **NEMO** package](#installation-of-the-nemo-package)
- [Preliminary Steps](#preliminary-steps)
- [Generating Ensembles](#generating-ensembles)
- [Running Ensemble Calculations](#running-ensemble-calculations)
- [Obtaining Results](#obtaining-results)
- [Extra Features](#extra-features)

# Introduction

**NEMO** is a package devoted to predicting the photophysics of molecules. This means being able to simulate absorption and emission spectra; estimating transition rates for radiative and non-radiative processes; and providing insight into these phenomena. **NEMO** comes from the realization that for properly predicting the behavior of excited molecules, it is necessary to go beyond the paradigm of running single-point excited state calculations on optimized molecular structures. 

To overcome this limitation, **NEMO** resorts to the nuclear ensemble method. Taking an optimized structure at a particular electronic state as starting point and computing its normal modes' frequencies we can generate an ensemble of molecular conformations that correspond to configurations a molecule may attain as a result of vibrations during the lifetime of its electronic state. Running excited state calculations on these various conformations allows us to obtain a fuller picture of the photphysics of a given molecule with vibrational and solvent effects being accounted for. In this way we move towards more accurate simulations that may help shed light on experimental results as well as assisting in materials discovery by exploring the optoelectronic potential of candidate molecules before synthesis.

The theoretical aspects behind **NEMO** can be found in the following papers:


> 1- [de Sousa, Leonardo Evaristo, and Piotr de Silva. "Unified framework for photophysical rate calculations in tadf molecules." Journal of Chemical Theory and Computation 17.9 (2021): 5816-5824](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00476)

>  2 - [de Sousa, Leonardo Evaristo, and Piotr de Silva. "Photophysics of Solvated Molecules: Computational Protocol Combining Nuclear Ensemble and Nonequilibrium State-Specific Solvation Methods." The Journal of Physical Chemistry A (2023)](https://pubs.acs.org/doi/full/10.1021/acs.jpca.3c03533).




# Installation of the **NEMO** package

The easiest way to install is to use pip:

`pip install nemophoto`

This will install the lates released version.

To install the version with the latest commit, run:

`pip install git+https://github.com/LeonardoESousa/NEMO`

Alternatively, clone the repository to your computer. 

`git clone https://github.com/LeonardoESousa/NEMO`

A new folder named **NEMO** will appear. Move to this folder and install using pip:

```
cd NEMO
pip install .
```

Once installed, you should be able to run the program from any folder in your computed by just using the `nemo` command.


# Preliminary Steps

**NEMO** requires output files for DFT frequency calculations as input in order to generate ensembles. Frequencies must all be real, which means that these calculations have to performed on optimized geometries. The kind of frequency calculation required depends on the kind of photophysical process you are interested in. 

1. Absorption spectrum from the ground state: in this case, ground state $(S_0)$ frequencies are needed to generate an $S_0$ ensemble.

2. Fluorescence and/or ISC rates from state $S_n$ $(n = 1,2, \cdots)$: in this case,  $S_n$ frequencies are needed to generate an $S_n$ ensemble. Unless you are interested in molecules with anti-Kasha behavior, you would go for the $S_1$ state.

3. Phosphorescence and/or rISC rates from state $T_n$ $(n = 1,2, \cdots)$: here,  $T_n$ frequencies are needed to generate a $T_n$ ensemble. Most often, you would be interested in the $T_1$ state, considering that fast relaxation from higher lying triplet states is usually expected, but calculations from any triplet state is possible. 

To get a full picture of the photophysics of your molecule you will need to run ensembles for the three cases mentioned above.

## Best Practices

### Gaussian 09/16 vs QChem

Frequency calculations used as input for **NEMO** can be run with either QChem or Gaussian 09/16. The actual ensemble calculations can only be run with QChem, though Qchem versions >= 5.0 are supported. 

### Method Selection

**NEMO** is typically run with TD(A)-DFT. However, ensemble calculation run with EOM-CCSD are also available as this method can be useful for studying molecules for which DFT is known to fail, such as for molecules with inverted singlet-triplet gaps. Note that EOM-CCSD calculations are much more computationally expensive that TD(A)-DFT ones, which restrict them to samll molecules only.

### Functional Selection

For TADF molecules, it is often the case that charge-transfer (CT) states are relevant. Proper description of CT states in DFT requires the use of range-separated functionals such as $\omega\text{B97XD}$ or $\text{LRC-}\omega\text{PBE}$. For better results, consider performing the non-empirically tuning of the $\omega$ parameter. This can be done within **NEMO**.  See details in the section on [Functional Tuning](#functional-tuning). The next best alternative would be to use the M062X functional.

### The Tamn-Dancoff Approximation (TDA)

It has been shown that employment of the Tamn-Dancoff approximation (TDA) may improve the accuracy of triplet state energies by mitigating triplet instability issues, as well
as providing accurate results in spectrum simulations. Therefore it is a good idea to use this approximation when running ensemble calculations that aim at obtaining transition rates involving triplet states. Note that TDA use is optional.

- In Gaussian, TDA can be invoked by substituting the ``TD`` keyword used for excited state calculations with ``TDA``.

- In QChem, TDA can be invoked by including the line ``RPA    False`` in the ``$rem`` section of the input file. 

## Examples of Input Files

Here are some examples of geometry optimization and frequency input files for QChem and Gaussian.  

### Gaussian 09/16

Some observations:

1. In the following, the $\text{LC-}\omega\text{hPBE}$ functional is used with the 6-31G(d,p) basis set.
2. The long-range parameter is set to 0.1843 $\text{bohr}^{-1}$. This is done by including the keywords `iop(3/108=0184300000) iop(3/107=0184300000)`. If these keywords are not added, the default $\omega$ value is used.
3. In the excited state calculations, TDA is being invoked. If you do not wish to use it, substitute the keyword `TDA` with `TD`.
4. It is recommended to use the `tight` keyword in the optimization step, especially for molecules with very low frequency modes. However, this ins not mandatory. 
5. The `noraman` keyword is used in the frequency jobs to speed up calculations. 
6. In the examples below, the optimization and frequency steps are separated in two consecutive jobs. If you make changes in the inputs for one job, make sure to make the same modification to the second job. 

#### **IMPORTANT!**

 Even though optimization and frequency calculations can be performed in a single Gaussian job by using the `opt freq` keywords simultaneously, in the examples below we separate the optimization and frequency parts in two jobs. The reason is that when `iop` keywords such as the one needed to change the value of the $\omega$ parameter are present, Gaussian applies them to the first job, i.e. optimization, but does not include it in the second job, i.e. frequency. As a result, frequencies end up being calculated with slightly different functionals than the one used for optimization. To make matters worse, this disparity does not necessarily translate into the presence of imaginary frequencies, so it is easy for this issue to go unnoticed. However, once the frequencies have been computed, Gaussian performs a final check on the geometry, which may still fail regardless of the absence of imaginary frequencies. Importantly, failing this final check does not throw an error, meaning that you may find yourself with a frequency output file with no imaginary frequencies and a "NORMAL TERMINATION OF GAUSSIAN" message that nevertheless does not correspond to an equilibrium geometry. The reasons for this are detailed [here](https://gaussian.com/faq3/). In view of this issue, **NEMO** performs a check on the frequency file when an ensemble is generated to see whether any of these problems is present.     

### Examples 

#### Ground State $(S_0)$

```
%nproc=40
%mem=100GB
%chk=s0.chk
# lc-whpbe/6-31G(d,p) iop(3/108=0184300000) iop(3/107=0184300000) opt=tight 

TITLE

0 1
#ADD XYZ COORDINATES HERE

--Link1--
%nproc=40
%mem=100GB
%chk=s0.chk
# lc-whpbe/6-31G(d,p) iop(3/108=0184300000) iop(3/107=0184300000) freq=noraman geom=allcheck
```


#### Singlet Excited State $(S_1)$ 

```
%nproc=40
%mem=100GB
%chk=s1.chk
# lc-whpbe/6-31G(d,p) iop(3/108=0184300000) iop(3/107=0184300000) opt=tight TDA

TITLE

0 1
#ADD XYZ COORDINATES HERE

--Link1--
%nproc=40
%mem=100GB
%chk=s1.chk
# lc-whpbe/6-31G(d,p) iop(3/108=0184300000) iop(3/107=0184300000) TDA freq=noraman geom=allcheck
```


#### Triplet Excited State $(T_1)$

```
%nproc=40
%mem=100GB
%chk=t1.chk
# lc-whpbe/6-31G(d,p) iop(3/108=0184300000) iop(3/107=0184300000) opt=tight TDA=(Triplets)

TITLE

0 1
#ADD XYZ COORDINATES HERE

--Link1--
%nproc=40
%mem=100GB
%chk=t1.chk
# lc-whpbe/6-31G(d,p) iop(3/108=0184300000) iop(3/107=0184300000) TDA=(Triplets) freq=noraman geom=allcheck
```

### QChem

Some observations:

1. In the following, the $\text{LRC-}\omega\text{PBE}$ functional is used with the 6-31G(d,p) basis set.
2. The long-range parameter is set to 0.1843 $\text{bohr}^{-1}$. This is done by including the keywords `omega      184` in the `$rem` section of the calculations. If this keyword is not added, the default $\omega$ value is used. You should also remove it if your functional is not long-range corrected.
3. In the excited state calculations, TDA is being invoked. If you do not wish to use it, substitute the keyword `RPA     False` with `RPA     True`. 
4. The input files here contain two jobs, an optimization followed by a frequency calculation. If you make changes in the inputs for one job, make sure to make the same modification to the second job.

#### Ground State $(S_0)$

```
$rem
JOBTYPE         opt
METHOD          LRC-wPBE
OMEGA           184
BASIS           6-31G(d,p)
MEM_TOTAL       100000
MEM_STATIC      2000
$end

$molecule
0 1
#ADD XYZ COORDINATES HERE
$end

@@@

$molecule
read
$end

$rem
JOBTYPE         freq
METHOD          LRC-wPBE
OMEGA           184
BASIS           6-31G(d,p)
MEM_TOTAL       100000
MEM_STATIC      2000
$end
```


#### Singlet Excited State $(S_1)$ 

```
$rem
JOBTYPE         opt
METHOD          LRC-wPBE
OMEGA           184
BASIS           6-31G(d,p)
MEM_TOTAL       100000
MEM_STATIC      2000
CIS_N_ROOTS     3
CIS_SINGLETS    true 
CIS_TRIPLETS    false
CIS_STATE_DERIV 1
RPA             false
$end

$molecule
0 1
#ADD XYZ COORDINATES HERE
$end

@@@

$molecule
read
$end

$rem
JOBTYPE         freq
METHOD          LRC-wPBE
OMEGA           184
BASIS           6-31G(d,p)
MEM_TOTAL       100000
MEM_STATIC      2000
CIS_N_ROOTS     3
CIS_SINGLETS    true 
CIS_TRIPLETS    false
CIS_STATE_DERIV 1
RPA             false
$end
```


#### Triplet Excited State $(T_1)$ 

```
$rem
JOBTYPE         opt
METHOD          LRC-wPBE
OMEGA           184
BASIS           6-31G(d,p)
MEM_TOTAL       100000
MEM_STATIC      2000
CIS_N_ROOTS     3
CIS_SINGLETS    false 
CIS_TRIPLETS    true
CIS_STATE_DERIV 1
RPA             false
$end

$molecule
0 1
#ADD XYZ COORDINATES HERE
$end

@@@

$molecule
read
$end

$rem
JOBTYPE         freq
METHOD          LRC-wPBE
OMEGA           184
BASIS           6-31G(d,p)
MEM_TOTAL       100000
MEM_STATIC      2000
CIS_N_ROOTS     3
CIS_SINGLETS    false 
CIS_TRIPLETS    true
CIS_STATE_DERIV 1
RPA             false
$end
```

## Dealing With Imaginary Frequencies

If an imaginary frequency is found (represented by a negative frequency in the output file), the structure does not correspond to a minimum in the energy potential surface. You may tighten the optimization criteria or start from a different initial structure to try to get rid of the imaginary frequencies. 

# Generating Ensembles

Once the optimization and frequencies are computed, it is time to generate ensembles.

We start by creating a folder and pasting the frequency output file in it. Make sure your output file ends with either `.log` or `.out`, as **NEMO** will look for it. Let's say we are generating an $S_1$ ensemble and we have a frequency file run with either Gaussian 09/16 or QChem called **freqS1.log**. We create a folder named **EnsembleS1** and move our frequency file there. Note that you may name the folder whichever way you prefer. It is a good idea, though, to identify to which electronic state the ensemble belongs.

Now we need to add to the folder a QChem template file. This file contains the `$rem` section of a QChem input file from where functional, basis set and other information are going to be retrieved. Other sections may also be included if needed. There is no need to include solvent related keywords and subsections as these are automatically added later on. The template file may have any name, just make sure its name ends with `.in`. In this example we name it `template.in`.  

## TD(A)-DFT

For TDA-DFT calculations, the `template.in` file will look something like this:

```
$rem
mem_total             100000
mem_static            2000
method                LRC-wPBE
omega                 184
basis                 6-31G(d,p)
RPA                   false
$end
```

## EOM-CCSD

For EOM-CCSD calculations, the `template.in` file will look something like this:

```
$rem
METHOD eom-ccsd
BASIS cc-pVDZ
MEM_TOTAL 200000
MEM_STATIC 2000
CC_MEMORY 50000
$end
```

At this point you folder should contain the following files:

```
    EnsembleS1
    ├── freqS1.log
    └── template.in
```

Now run the `nemo` command inside the **EnsembleS1** folder. The following menu will show up on the screen. Select option **1 - Generate the inputs for the nuclear ensemble calculation** by typing 1 and pressing Enter. 

```
[EnsembleS1]$ nemo
███▄▄▄▄      ▄████████    ▄▄▄▄███▄▄▄▄    ▄██████▄  
███▀▀▀██▄   ███    ███  ▄██▀▀▀███▀▀▀██▄ ███    ███ 
███   ███   ███    █▀   ███   ███   ███ ███    ███ 
███   ███  ▄███▄▄▄      ███   ███   ███ ███    ███ 
███   ███ ▀▀███▀▀▀      ███   ███   ███ ███    ███ 
███   ███   ███    █▄   ███   ███   ███ ███    ███ 
███   ███   ███    ███  ███   ███   ███ ███    ███ 
 ▀█   █▀    ██████████   ▀█   ███   █▀   ▀██████▀  
------------------Photophysics--------------------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
EXTRA FEATURES:
        9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)
1


Frequencty log file: freqS1.log
You are using a Gaussian log file.
QChem template file: template.in

The configurations to be used are:

$rem
mem_total             100000
mem_static            2000
method                LRC-wPBE
omega                 184
basis                 6-31G(d,p)
RPA                   false
$end
```

Here we see that **NEMO** has recognized the frequency file and the template file. Then it prints the configurations that you have provided so we may double check before proceeding.  

Next we must provide two solvent properties: the stactic dielectric constant $(\epsilon)$ and the refractive index $(n_r)$. In this example we are going with toluene, for which $\epsilon = 2.38$ and $n_r = 1.497$.

```
Solvent's static dielectric constant?
2.38 
Solvent's refractive index?
1.497
```

Now we must select the number of excited states included in the calculation. For ISC rate calculations, 5 excited states is usually more than enough. For absorption spectra, we may want to go with a higher number. 

```
How many excited states?
5
Are you interested in absorption spectra ONLY? (y or n)
n
```

If you are interested in absorption spectra only, type "y". This is the case, for instance, when generating $S_0$ ensembles, which are useful only for obtaining absorption spectra. If you reply with "y", spin-orbit couplings are not going to be computed. As a result, all calculations will run faster, but it will not be possible to compute ISC or phosphorescence rates.  

Next we select the number of geometries to be sampled in the ensemble. There is no clear cut answer as to how many geometries are needed. Typically, a number of 500 geometries is able to generate results with reasonably sized uncertainties. If in the end we are not happy with accuracy, we may add more geometries to the ensemble until the sampling error is low enough.  

```
How many geometries to be sampled?
500 
```

Finally, we select the temperature to be used in the sampling procedure, which is usually set at 300 K.

```
Temperature in Kelvin?
300

Generating geometries...

  100.0% of the geometries done. 2 geometries rejected

Done! Ready to run.
```

At this point, the folder should contain the following:

```
    EnsembleS1
    ├── freqS1.log
    ├── template.in
    ├── Magnitudes_300K_.lx
    └── Geometries
```

The folder named **Geometries** will contain 500 QChem input files corresponding to the sampled geometries. The file named **Magnitudes_300K_.lx** contains the displacements along each vibrational mode that produce each sampled geometry. There is no need to do anything with it, just keep it there. 

That's it. We are ready to run the calculations.

## A Note on Rejected Geometries

**NEMO** computes the adjacency matrix of the optimized molecule and compares it with the adjacency matrix of each geometry sampled in the ensemble. If the matrices do not match, i.e. if one or more bonds are broken, the sampled geometry is rejected. This helps reduce issues with unphysical geometries and negative transition energies.

# Running Ensemble Calculations

The next step consists of running all the individual calculations present in the **Geometries** folder. This can be done by whatever method. However, **NEMO** comes with a built in batch scheme. Usage will vary depending on the user's queue system. Here we present an example for use with slurm:  

First, we need to create a submission script, which we name **nemo.sh**. The following is an example that first sets up some options and loads the appropriate modules. Adjust it according to the cluster you use. The essential line is the `bash $1` line.   

```
#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --time=1-0
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --partition=xeon40
module load QChem/5.2-multicore
export $QCLOCALSCR=/scratch/

#This is the key line:
bash $1

rm -rf /scratch/*
rm slurm*out 
```

For conveninence, you may put the **nemo.sh** file in your home directory so you may use it for all your **NEMO** projects.

Now we create a file named **batch.sh** in the **EnsembleS1** directory. Assuming the **nemo.sh** file is in the user's home directory, this can be done by running:

```
echo "sbatch ~/nemo.sh \$1" >> batch.sh
```

The folder structure should look like the following :

```
    EnsembleS1
    ├── freqS1.log
    ├── template.in
    ├── Magnitudes_300K_.lx
    ├── batch.sh
    └── Geometries
```

Now, to run the ensemble calculations, run the `nemo` command and select  option **2 - Run the ensemble calculations**:

```
[EnsembleS1]$ nemo
███▄▄▄▄      ▄████████    ▄▄▄▄███▄▄▄▄    ▄██████▄  
███▀▀▀██▄   ███    ███  ▄██▀▀▀███▀▀▀██▄ ███    ███ 
███   ███   ███    █▀   ███   ███   ███ ███    ███ 
███   ███  ▄███▄▄▄      ███   ███   ███ ███    ███ 
███   ███ ▀▀███▀▀▀      ███   ███   ███ ███    ███ 
███   ███   ███    █▄   ███   ███   ███ ███    ███ 
███   ███   ███    ███  ███   ███   ███ ███    ███ 
 ▀█   █▀    ██████████   ▀█   ███   █▀   ▀██████▀  
------------------Photophysics--------------------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
EXTRA FEATURES:
        9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)
2
```

Now we are prompted for the maximum number of jobs to be submitted at the same time. Clusters may have limitations on the number of jobs on the queue an user may have, so here we can select an appropriate number, let's say 10. Each job may contain more than one calculation to be run simultaneously. This may be useful because spin-orbit coupling calculations are not parallelized in QChem, so when calculations reach this point, a single thread ends up being even though you requested many more. This has a tendency of making cluster administrators unhappy, so bundling several calculations into a single job helps mitigate this issue. Finally, each individual calculation is run with $T$ threads. If the total number of CPU cores available is $C$ and we are including $n$ calculations in a single job, then the number of threads for each individual job should be $T = C/n$.  

For the sake of example, let's say we have a total of 24 cores available and we are bundling 2 calculations together. This leaves us with 12 threads for each calculation. This is implemented as follows  

```
Maximum number of jobs to be submitted simultaneously?
10
Number of threads for each individual calculation?
12
Number of calculations in each job?
2 
```

If everything is set correctly, all calculations in the **Geometries** folder will eventually be run. When a job finishes, a new one will be submitted until all are completed. 

After running option 2, the folder structure will look as follows:

```
    EnsembleS1
    ├── freqS1.log
    ├── template.in
    ├── Magnitudes_300K_.lx
    ├── batch.sh
    ├── limit.lx
    ├── nohup.out
    └── Geometries
```

The `nohup.out` file will register the submitted jobs. More importantly, the `limit.lx` file contains a number. This is the maximum number of jobs to be submitted simultaneously that was selected earlier. We may change this number at any point as needed.

## Following the Progress

We may check the progress of the calculations by running the `nemo` command from the **EnsembleS1** folder and choosing option **3 - Check the progress of the calculations**. It should look as follows:

```
[EnsembleS1]$ nemo
███▄▄▄▄      ▄████████    ▄▄▄▄███▄▄▄▄    ▄██████▄  
███▀▀▀██▄   ███    ███  ▄██▀▀▀███▀▀▀██▄ ███    ███ 
███   ███   ███    █▀   ███   ███   ███ ███    ███ 
███   ███  ▄███▄▄▄      ███   ███   ███ ███    ███ 
███   ███ ▀▀███▀▀▀      ███   ███   ███ ███    ███ 
███   ███   ███    █▄   ███   ███   ███ ███    ███ 
███   ███   ███    ███  ███   ███   ███ ███    ███ 
 ▀█   █▀    ██████████   ▀█   ███   █▀   ▀██████▀  
------------------Photophysics--------------------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
EXTRA FEATURES:
        9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)
3

There are 500 successfully completed calculations out of 500 inputs
100.0 % of the calculations have been run.
```

A message is also displayed if failed jobs are detected. Failed jobs must be investigated individually for their cause. 

## Aborting a Run

To abort an ensemble calculation, run the `nemo` command inside the folder **EnsembleS1** and choose option **4 - Abort my calculations**. 

Alternatively, you may just delete the `limit.lx` file from the folder. This will kill the script controlling job submission. Note that this will not kill jobs already submitted, it will only prevent new jobs from being added to the queue. 


# Obtaining Results

Once all calculations are run, there are three ways to get results: 

1. The recommended option is using a separate visualization tool called [**NEMOview**](https://github.com/LeonardoESousa/nemoview). With **NEMOview** you will be able to interactively plot spectra, visualize rates, generate levels diagrams and more. 

2. Use options 5 and 6 in the **NEMO** menu.

3. Use the python API.    


## A Note on Solvent Selection

Regardless of the method used to compile results, you will be prompted to select the solvent parameters (dielectric constant and refractive index) for your rates and/or spectra. These solvent parameters may be different than the ones used in the ensemble generation step. In this case, **NEMO** will employ an extrapolation procedure to obtain solvent corrections for the new dielectric constant and refractive index. Details on how this is performed can be found in [this paper](https://pubs.acs.org/doi/full/10.1021/acs.jpca.3c03533).

## Using **NEMOview**

Go to the [**NEMOview**](https://github.com/LeonardoESousa/nemoview) page and follow the installation instructions there. Once the program is installed, the next step is to run option **7 - Gather ensemble data only** in the **EnsembleS1** folder. 

```
[EnsembleS1]$ nemo
███▄▄▄▄      ▄████████    ▄▄▄▄███▄▄▄▄    ▄██████▄  
███▀▀▀██▄   ███    ███  ▄██▀▀▀███▀▀▀██▄ ███    ███ 
███   ███   ███    █▀   ███   ███   ███ ███    ███ 
███   ███  ▄███▄▄▄      ███   ███   ███ ███    ███ 
███   ███ ▀▀███▀▀▀      ███   ███   ███ ███    ███ 
███   ███   ███    █▄   ███   ███   ███ ███    ███ 
███   ███   ███    ███  ███   ███   ███ ███    ███ 
 ▀█   █▀    ██████████   ▀█   ███   █▀   ▀██████▀  
------------------Photophysics--------------------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
EXTRA FEATURES:
        9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)
7
What is the initial state (S0, S1, T1, S2 ...)? Accepts comma separated values Ex: T1,T2
s1
```

Here we are prompted for the initial state. In our example, we are running an $S_1$ ensemble so $S_1$ is the initial state. This will generate a file called `Ensemble_S1_.lx`, which we may load into [**NEMOview**](https://github.com/LeonardoESousa/nemoview) for visualization. Several ensembles belonging to different molecules and/or different states can be uploaded to [**NEMOview**](https://github.com/LeonardoESousa/nemoview) simultaneously for comparisons. Please, check the [**NEMOview**](https://github.com/LeonardoESousa/nemoview) tutorial.


## Using the **NEMO** Menu


### Absorption Spectra

For absorption spectrum simulations: Run the `nemo` command and choose option **5 - Generate the absorption spectrum**.

```
[EnsembleS1]$ nemo
███▄▄▄▄      ▄████████    ▄▄▄▄███▄▄▄▄    ▄██████▄  
███▀▀▀██▄   ███    ███  ▄██▀▀▀███▀▀▀██▄ ███    ███ 
███   ███   ███    █▀   ███   ███   ███ ███    ███ 
███   ███  ▄███▄▄▄      ███   ███   ███ ███    ███ 
███   ███ ▀▀███▀▀▀      ███   ███   ███ ███    ███ 
███   ███   ███    █▄   ███   ███   ███ ███    ███ 
███   ███   ███    ███  ███   ███   ███ ███    ███ 
 ▀█   █▀    ██████████   ▀█   ███   █▀   ▀██████▀  
------------------Photophysics--------------------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
EXTRA FEATURES:
        9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)
5
The spectrum will be run with the following parameters:

Solvent dielectric constant: 2.380
Solvent refractive index: 1.497

Are you satisfied with these parameters? y or n?
y
What is the initial state (S0, S1, T1, S2 ...)? Accepts comma separated values Ex: T1,T2
s1
Spectrum printed in the cross_section_S1_.lx file
```

If you are not satisfied with the solvent parameters, you may change them. The parameters suggested by **NEMO** correspond to those present in the ensemble calculations. 

The spectrum will be printed in a file named `cross_section_s1.lx`. Note that since we set $S_1$ as the initial state, this corresponds to the absorption spectrum from the $S_1$ state. The `cross_section_s1.lx`. In this file, the first column corresponds to energies (in eV). The following columns correspond to the absorption cross-section (in $Å^2$) of each individual state (in this example, $S_2$ to $S_5$). The `Total` column gives the total absorption cross-section including all states and the last columns provides the uncertainty to the `Total` cross-section. 

## Emission Spectra and Rates

To obtain emission spectra and rates for both emission and ISC,  run the `nemo` command and choose option **6 - Estimate rates and compute emission spectrum**. Follow the instructions to set the parameters. Three files will be generated: an Ensemble file, with data from the ensemble of geometries; a differential_rate file, with the emission spectrum; a rates file, with all available rates (-> denote radiative transitions and ~> denote ISC transitions).

```
[EnsembleS1]$ nemo
███▄▄▄▄      ▄████████    ▄▄▄▄███▄▄▄▄    ▄██████▄  
███▀▀▀██▄   ███    ███  ▄██▀▀▀███▀▀▀██▄ ███    ███ 
███   ███   ███    █▀   ███   ███   ███ ███    ███ 
███   ███  ▄███▄▄▄      ███   ███   ███ ███    ███ 
███   ███ ▀▀███▀▀▀      ███   ███   ███ ███    ███ 
███   ███   ███    █▄   ███   ███   ███ ███    ███ 
███   ███   ███    ███  ███   ███   ███ ███    ███ 
 ▀█   █▀    ██████████   ▀█   ███   █▀   ▀██████▀  
------------------Photophysics--------------------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
EXTRA FEATURES:
        9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)
6
The rates will be calculated with the following parameters:

Solvent dielectric constant: 2.380
Solvent refractive index: 1.497

Are you satisfied with these parameters? y or n?
y
What is the initial state (S1, T1, S2 ...)? Accepts comma separated values Ex: T1,T2
s1
Spectrum printed in the differential_rate_S1.lx file
Rates are written in the rates_S1_.lx file
```

As in the case of absorption spectra, you may select new values for the solvent parameters. Then you must select the initial state to be considered, in our case the $S_1$ state. Finally, three files will be produced: `Ensemble_S1_.lx`,`differential_rate_S1.lx` and `rates_S1_.lx`. The first one is the same we obtain when running option 7 as mentioned above. The second contains three columns: energy (eV), the differential emission rate (dimensionless) and the uncertainty associated with the differential emission rate (also dimensionless). The total emission rate will also be printed. Note that for singlet ensembles ($S_1$, $S_2$ etc), emission means fluorescence. For triplet ensembles ($T_1$, $T_2$ etc) we would have phosphorescence. 

The third file, `rates_S1_.lx`, contains the following:

```
#Epsilon: 2.380 nr: 1.497
Transition Rate(s^-1) Error(s^-1) Prob(%) AvgDE+L(eV) AvgSOC(meV) AvgSigma(eV) AvgConc(%)
    S1->S0   1.83e+07    0.36e+07     2.5       2.410         nan        0.038       35.1
    S1~>T1   6.83e+08    1.19e+08    94.8       0.019       0.418        0.035        6.2
    S1~>T2   1.93e+07    1.55e+07     2.7      -0.005       0.235        0.029        0.3
    S1~>T3   5.16e+03    5.10e+03     0.0       0.108       0.218        0.027        0.2
    S1~>T4   8.47e+00    8.47e+00     0.0       0.142       0.271        0.026        0.2
    S1~>T5   3.30e-16    3.30e-16     0.0       0.271       0.203        0.027        0.2
```

Here we have rate estimates for several processes. In the `Transition` column, processes containing `->` refer to radiative transitions, i.e. emission. Those with `~>` refer to non-radiative transitions, i.e. ISC. The `Rate` columns gives rate estimates in $s^{-1}$. Associated uncertainties are presented in the `Error` column. The `Prob` columns gives the yield of each transition (Note that internal conversion rates are not yet available, so they are not accounted for). The remaining columns correspond to average values of the energy gap  (`AvgDE+L`), spin-orbit coupling (`AvgSOC`) and broadening factor (`AvgSigma`). Finally, the `AvgConc` column provides the fraction of conformations in the ensemble that actively contribute to generate the estimated rate. 


## Using the API

The first step in this case is to use option 7 on the **NEMO** menu to generate the ensemble files, just as shown in the [Using **NEMOview**](#using-nemoview) section. 

A jupyter notebook with instructions on how to proceed can be found [here](https://github.com/LeonardoESousa/nemoview/tree/main/Tutorial/API.ipynb).


# Extra Features

## Functional Tuning 

**Note**: Functional tuning is currently only available for Gaussian 09/16 users.

The tuning of long-range corrected functionals consists in finding the value for the long-range parameter $\omega$ that minimizes an objective function $J(\omega)$ given by:

$J^2(\omega) = (E_{HOMO}^{0} + (E_{SCF}^{+1} -E_{SCF}^{0}))^2\\ \hspace{10.8mm} + (E_{HOMO}^{-1} + (E_{SCF}^{0} -E_{SCF}^{-1}))^2 $

In this expression, $E_{HOMO}$ is the energy of the highest occupied molecular orbital and $E_{SCF}$ is the self-consistent ground-state energy. The superscripts refer to the molecular charge from which these energies must be calculated. Details can be found in:

> Stein, Tamar, Leeor Kronik, and Roi Baer. "Reliable prediction of charge transfer excitations in molecular complexes using time-dependent density functional theory." Journal of the American Chemical Society 131.8 (2009): 2818-2820.

To use this feature, first create a new folder and move there. 

```
mkdir tuning
cd tuning
```

Now we must add a Gaussian input file containing our molecular structure, charge and multiplicity, functional and basis set as well as the `%nproc` and `%mem` sections. Let's call it `geom.com`. 

```
%nproc=40
%mem=100GB
# lc-whpbe/6-31G(d,p)

TITLE

0 1
#ADD XYZ COORDINATES HERE

```

Now we need our `batch.sh` file. For a slurm system, it will read

```
sbatch g16.sh $1
```

The `g16.sh` file includes the slurm commands, loads the Gaussian 16 module and should look something like the following

```
#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --time=1-0
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --partition=xeon40

module load Gaussian/16
export GAUSS_SCRDIR=/scratch/


bash $1

rm -rf /scratch/*
rm slurm*.out
```

The key line here is the `bash $1` line. Everything else should be adjusted to your particular circumstances.

At this point the folder should contain 

```
    tuning
    ├── geom.com
    ├── g16.sh
    └── batch.sh

```

Finally, run the nemo command, choose option **9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)** and follow the instructions.

```
[omega] nemo
███▄▄▄▄      ▄████████    ▄▄▄▄███▄▄▄▄    ▄██████▄ 
███▀▀▀██▄   ███    ███  ▄██▀▀▀███▀▀▀██▄ ███    ███
███   ███   ███    █▀   ███   ███   ███ ███    ███
███   ███  ▄███▄▄▄      ███   ███   ███ ███    ███
███   ███ ▀▀███▀▀▀      ███   ███   ███ ███    ███
███   ███   ███    █▄   ███   ███   ███ ███    ███
███   ███   ███    ███  ███   ███   ███ ███    ███
 ▀█   █▀    ██████████   ▀█   ███   █▀   ▀██████▀
------------------Photophysics--------------------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
EXTRA FEATURES:
        9 - Perform tuning of long range corrected functional (Gaussian 09/16 only)
9

geom.com
Is this the input or log file? y ou n?
y
This is the configuration taken from the file:

Functional/basis: lc-whpbe/6-31G(d,p)
%nproc=40
%mem=100GB
Initial Omega: 0.1 bohr^-1
Step: 0.025 bohr^-1
Optimize at each step: yes
Are you satisfied with these parameters? y or n?
y
```

Functional, basis set, nproc and mem are taken from the input file provided. The remaining parameters are the default options, which should be good enough for any application. Note the option `Optimize at each step:`. By default it is set to `yes`, which means that for each $\omega$ value, the structure will be reoptimized. This means that the initial structure provided in the `geom.com` file does not need to be an equilibrium structure. On the other hand, if you decide not to optimize at each step, then the initial structure should already be optimized. 

Now we get prompted for the version of Gaussian we want to use. We are also asked whether we want to minimize the submission of jobs. If yes, the calculations that can be done in parallel will be instead performed in series in a single job. This is useful when the time jobs spend waiting in the queue is too long.

```
g16 or g09?
g16
Minimize submission of jobs: y/n
y
nohup: ignoring input and appending output to ‘nohup.out’
```

At this point, input files will be generated and jobs submitted. A file named `omega.lx` will be created, from where you may follow the progress. Eventually, it will look something like this:

```
#w(10^4 bohr^-1)    J(eV)
01000               0.8557
01250               0.5571
01500               0.3003
01750               0.0789
01812               0.0294
01843               0.0079
01874               0.0216
01875               0.0221
01905               0.0444
01936               0.0675
02000               0.1141

#Done! Optimized value: 01843
```

After the optimization procedure is concluded, a file named `tuned_w.com` will be created. This file is a Gaussian input file for which the `iop` keywords have been included to set up the functional with optimized $\omega$ value. If you ran the optimization with `Optimize at each step: yes`, then the stucture in the input file will be the optimized geometry using the optimized $\omega$ value. Otherwise, it will be the structure provided in the `geom.com` file. 

### Aborting a Functional Tuning

To abort a functional tuning procedure, just delete the `limit.lx` file from the folder. This way, no further jobs will be submitted.


