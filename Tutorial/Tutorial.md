# Introduction

In the following, we will learn how to use NEMO package (Leonardo Evaristo de Sousa and Piotr de Silva, Journal of Chemical Theory and Computation 2021 17 (9), 5816-5824 DOI: 10.1021/acs.jctc.1c00476)

This step by step tutorial will use the example of Heptazine (Hz)
Firsr, one will need to compute the frequencies of the optimized structures at the S0, S1, and T1 states using **QChem** ot **Gaussian**
Then, an ensemble of geometries will be created using the **NEMO** package.

Finally, one will visualize the results using nemoview.

# Installation of the NEMO package

Run:

`pip install nemophoto`

To get the latest commit, run:

`pip install git+https://github.com/LeonardoESousa/NEMO`

Alternatively, clone the repository to your computer. Inside the NEMO folder, run:

`pip install .`

Once installed, you should be able to run the program from any folder by just using the `nemo` command.

# Ground-state optimization

**NEMO** is interfaced with QChem, but it also accepts Gaussian frequency calculations as input.

Here is an example of optimization and frequencies calculation of the ground state of the molecule using qchem. The inputfile is **Hz_optfreqS0.com** and the outputfile is **Hz_optfreqS0.out**

```
$rem
GEOM_OPT_PRINT  6 !Print all the information of the optimization process (optional)
JOBTYPE         opt
METHOD  wb97xD
BASIS   cc-pVDZ
MEM_TOTAL       4000
MEM_STATIC      100
$end

$molecule
0 1
N
C  1 1.390948
N  2 1.317471 1 119.100767
C  3 1.315069 2 116.705881 1 -0.004215 0
H  4 1.075573 3 115.805871 2 179.983742 0
N  4 1.315046 3 128.387403 2 0.010584 0
C  6 1.317493 4 116.706392 3 -0.011785 0
N  7 1.317470 6 121.799814 4 179.997928 0
C  8 1.315070 7 116.706019 6 -179.995454 0
H  9 1.075573 8 115.805882 7 -179.990374 0
N  9 1.315045 8 128.387319 7 -0.002529 0
C  11 1.317494 9 116.706427 8 0.006074 0
N  12 1.317470 11 121.799668 9 179.993285 0
C  13 1.315070 12 116.705933 11 -179.996228 0
H  14 1.075572 13 115.805919 12 -179.993637 0
N  14 1.315045 13 128.387316 12 -0.002210 0
$end

@@@

$molecule
read
$end

$rem
JOBTYPE         freq
METHOD  wb97xD
BASIS   cc-pVDZ
MEM_TOTAL       4000
MEM_STATIC      100
$end
```

If an imaginary frequency is found, the molecule is unstable and the following computations might be erronous. Hence, a strategy is to tighten the convergence criteria and re-run the calculation. Or rerun the optimization from the saddle point.

Once the optimization and the frequencies are computed, it is time to generate the ensemble.

## Generating the ensemble

We start by creating a folder and pasting the file *Hz_optfreqS0.out* in it. Let's name the folder **EnsembleS0**
Now go to the **EnsembleS0** folder and run the `nemo` command
    
```
        cd EnsembleS0
        nemo

```

The following menu will show up in the screen. Select the first option *Generate the inputs for the nuclear ensemble calculation* and follow the instructions.

```
[EnsembleS0]$ nemo
#     # ####### #     # #######
##    # #       ##   ## #     #
# #   # #       # # # # #     #
#  #  # #####   #  #  # #     #
#   # # #       #     # #     #
#    ## #       #     # #     #
#     # ####### #     # #######
----------Photophysics---------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations #(deletes the limit.lx file in the folder you are. This stops the submission of further jobs. It does not kill jobs already on the queue)
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
1

Hz_optfreqS0.out #File contains structure and vibrational modes
Is this the frequency file? y ou n? #Press "n" and the program will propose another file
y

The suggested configurations for you are: #Check the calculation parameters : method, basis set and so on

$rem
GEOM_OPT_PRINT  6
METHOD  wb97xD
BASIS   cc-pVDZ
MEM_TOTAL       4000
MEM_STATIC      100
$end

Are you satisfied with these parameters? y or n?
y
Solvent's static dielectric constant?
2.38 #For toluene
Solvent's refractive index?
1.4 #For toluene
How many excited states?
10
Prepare input for absorption or fluorescence spectrum only? (y or n)  press "n" only for excited-state geometries. Allow to compute SOC rates /!\ SOC not possible with ADC(2) method 
Y 
Ok, calculations will only be suitable for absorption or fluorescence spectrum simulations!

How many geometries to be sampled?
500 #To have a representative ensemble, at least 200
Temperature in Kelvin?
300

Generating geometries...

  100.0% of the geometries done.

Done! Ready to run.
```

The folder **Geometries** is generated with all the geometries required. These are generated according to the molecule's vibrational modes. 
In addition to the folder **Geometries**, the file Magnitudes_300K_.lx is generated. It contains the magnitudes of the displacements along the vibrational modes and allows for an ensemble to me rebuilt.

Check if the input file **Geometry-1-.com** looks like :

```
[EnsembleS0]$ cat Geometries/Geometry-1-.com

$comment
ABSSPCT
$end

$rem
cis_n_roots             10
cis_singlets            true
cis_triplets            true
calc_soc                false
STS_MOM                 true
CIS_RELAXED_DENSITY     TRUE
solvent_method          PCM
GEOM_OPT_PRINT  6
METHOD  wb97xD
BASIS   cc-pVDZ
MEM_TOTAL       4000
MEM_STATIC      100
$end

$pcm
theory                  IEFPCM
ChargeSeparation        Marcus
StateSpecific           Perturb
$end

$solvent
Dielectric              2.38
OpticalDielectric       1.9599999999999997 !Square of 1.4
$end

$molecule
0 1
N   0.03971718136936  -0.04954075376651  -0.09305162034207
C   -0.60242722734286  -1.22998167855624  -0.03840524045065
N   0.08748486595726  -2.26651486282793  0.15508198357482
C   1.39480492258405  -2.21674883250099  0.01480361474069
H   1.99322406536358  -3.10087112748250  0.07320854372390
N   2.14168805921999  -1.02737120411578  -0.09049653528144
C   1.37350101239272  0.01864790878427  0.01099262644877
N   1.97799997771598  1.21738510587816  0.04235194464137
C   1.21101778083085  2.25865379842196  0.02131885328317
H   1.70919782331199  3.18060859112108  0.15211877100936
N   -0.10550847045592  2.36094988141728  0.00161851202900
C   -0.80469867294343  1.19326561222245  -0.04495861619107
N   -2.11861436012403  1.14065212096364  0.06200387308776
C   -2.62996485441531  -0.12475124453871  0.01684246292799
H   -3.79497911209860  -0.08861776682951  0.20606551546139
N   -1.96641949693936  -1.28848651869247  -0.09233838212287
$end

```

## Run the ensemble calculation

The next step is running all the inidividual calculations in the ensemble. This can be done by whatever method. However, **NEMO** comes with a built in batch scheme. Usage will vary depending on the user's queue system. Here we present an example for use with slurm:  

First, one needs to create the submission script, which we call here **~/nemo.sh** . The following is an example that loads the appropriate modules. The essential line is the ```bash $1``` line.   

```#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --time=1-0
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --partition=xeon24
module purge
module use /home/energy/modules/modules/all
module --ignore-cache load "binutils/2.31.1-GCCcore-8.2.0"
module load iomkl
module load QChem/5.2-multicore
export $QCLOCALSCR=/scratch/lajour
bash $1
rm -rf /scratch/lajour/*
rm slurm*out 
```

Then, one needs to create the **batch.sh** file in the **EnsembleS0** directory. Assuming the **nemo.sh** file is in the user's home directory, this can be done as:

```
        echo "sbatch ~/nemo.sh \$1" >> EnsembleS0/batch.sh
```

The file structure should look like the following :

```    
    EnsembleS0
    ├── batch.sh
    └── Hz_optfreqS0.out
    └── Geometries
```

Now, to run the ensemble calculation, run the `nemo` command and select he second option:

```
[EnsembleS0]$ nemo
#     # ####### #     # #######
##    # #       ##   ## #     #
# #   # #       # # # # #     #
#  #  # #####   #  #  # #     #
#   # # #       #     # #     #
#    ## #       #     # #     #
#     # ####### #     # #######
----------Photophysics---------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations (deletes the limit.lx file in the folder you are. This stops the submission of further jobs. It does not kill jobs already on the queue)
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
2 # Run the ensemble calculations


batch.sh
Is this the batch script? file? y ou n?
y
Maximum number of batches to be submitted simultaneously?
10 # Number of jobs in the queue at one time
Number of processors for each individual job
12 # Depends on the partition uses: The multiplication with the “ jobs in each batch has to be equal to the number of processors in the partition
Number of jobs in each batch. Here, 24 processors were available, so 2 job were run simultaneously with 12 procs each
2 # Number of jobs in each batch file

```
If everything is set correctly, all calculations in the Geometries folder will be run. You may check the progress of the calculations by running the `nemo` command from the EnsembleS0 folder and choosing option 3. It should look as follows:

```
#     # ####### #     # #######
##    # #       ##   ## #     #
# #   # #       # # # # #     #
#  #  # #####   #  #  # #     #
#   # # #       #     # #     #
#    ## #       #     # #     #
#     # ####### #     # #######
----------Photophysics---------

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
3


There are 283 successfully completed calculations out of 500 inputs
There are 3 failed jobs. If you used option 2, check the nohup.out file for details.
57.2 % of the calculations have been run.
```

Failed jobs may require extra keywords increasing, for instance, the number of SCF cycles. 


# Excited state optimization

Here is an example of optimization and frequencies calculation of the first singlet and triplet excited state of the molecule using **Qchem**. The input files are **Hz_optfreqS1.com** and **Hz_optfreqT1.com** and the output files are **Hz_optfreqS1.out** **Hz_optfreqT1.out**

```
$rem
GEOM_OPT_PRINT  6
JOBTYPE         opt
METHOD  wb97xD
BASIS   cc-pVDZ
CIS_N_ROOTS     5
CIS_SINGLETS    true !False for triplet optimization
CIS_TRIPLETS    false !True for triplet optimization
CIS_STATE_DERIV 1
MEM_TOTAL       4000
MEM_STATIC      100
GEOM_OPT_TOL_GRADIENT 150 !Tighened the convergence criteria because negatives frequencies were first computed.
GEOM_OPT_TOL_DISPLACEMENT       600
GEOM_OPT_TOL_ENERGY     50
$end

$molecule
0 1
N          0.00006        0.00006        0.15645
C         -0.62938        1.26435        0.01273
N          0.11879        2.35949       -0.03235
C          1.45189        2.18957        0.00053
H          2.05506        3.09906        0.00514
N          2.12743        1.02776       -0.03146
C          1.40950       -0.08726        0.01379
N          1.98427       -1.28242       -0.03253
C          1.17057       -2.35197       -0.00089
H          1.65699       -3.32891        0.00264
N         -0.17328       -2.35610       -0.03200
C         -0.78055       -1.17699        0.01450
N         -2.10289       -1.07699       -0.03129
C         -2.62259        0.16238        0.00018
H         -3.71162        0.22968        0.00460
N         -1.95396        1.32817       -0.03288
$end

@@@

$molecule
read
$end

$rem
JOBTYPE         freq
METHOD  wb97xD
BASIS   cc-pVDZ
CIS_STATE_DERIV     1
CIS_N_ROOTS     5
CIS_SINGLETS    true !False for triplet optimization
CIS_TRIPLETS    false !True for triplet optimization
MEM_TOTAL       4000
MEM_STATIC      100
MAX_SCF_CYCLES 200
MAX_CIS_CYCLES 200
$end

```

## Generating the ensemble

One needs to create the **EnsembleS1** (**EnsembleT1**) directory and to paste the file *Hz_optfreqS1.out* (*Hz_optfreqT1.out*)in it.

Then, one needs to create the **nemo.sh** file in the **EnsembleS1** (**EnsembleT1**)  directory. This can be done by running the following cell.

```
        mkdir EnsembleS1/
        cp Hz_optfreqS1.out EnsembleS1/
        echo "sbatch ~/nemo.sh \$1" >> EnsembleS1/batch.sh
```
 or
```
        mkdir EnsembleT1/
        cp Hz_optfreqS1.out EnsembleT1/
        echo "sbatch ~/nemo.sh \$1" >> EnsembleT1/batch.sh       
```

The folder structure should look like the following :

```    
    EnsembleS1
    ├── batch.sh
    └── Hz_optfreqS1.out
or    
    EnsembleT1
    ├── batch.sh
    └── Hz_optfreqT1.out
```

Finally, go to the **EnsembleS1** (**EnsembleT1**) folder and generate the ensemble with the following steps using the `nemo` command with the first option *Generate the inputs for the nuclear ensemble calculation*
    
```    
        cd EnsembleS1
        nemo
```

Differently from the ground-state ensemble, for the question "Prepare input for absorption or fluorescence spectrum only? (y or n)" press "n". This will enable the calculation of spin-orbit couplings necessary for ISC and phosphorescence calculations.       

```
[EnsembleS1]$ nemo
#     # ####### #     # #######
##    # #       ##   ## #     #
# #   # #       # # # # #     #
#  #  # #####   #  #  # #     #
#   # # #       #     # #     #
#    ## #       #     # #     #
#     # ####### #     # #######
----------Photophysics---------

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
1

Hz_optfreqT1.out
Is this the frequency file? y ou n?
y

The suggested configurations for you are:

$rem
METHOD  wb97xD
BASIS   cc-pVDZ
MEM_TOTAL       4000
MEM_STATIC      100
$end

Are you satisfied with these parameters? y or n?
y
Solvent's static dielectric constant?
2.38
Solvent's refractive index?
1.4
How many excited states?
5
Prepare input for absorption or fluorescence spectrum only? (y or n)
n #Will compute emission spectra AND SOC constants
Ok, calculations will be suitable for all spectra and ISC rate estimates!

How many geometries to be sampled?
500
Temperature in Kelvin?
300

Generating geometries...

  100.0% of the geometries done.

Done! Ready to run.
```

The folder **Geometries** is generated with all the geometries required. These are generated along the vibrational modes. 

Check if the input file **Geometry-01-.com** looks like the following cell, if so run the ensemble calculation as described in the previous section.

```
[EnsembleS1]$ cat Geometries/Geometry-01-.com

$comment
EMISPCT
$end

$rem
cis_n_roots             5
cis_singlets            true
cis_triplets            true
calc_soc                true
STS_MOM                 true
CIS_RELAXED_DENSITY     TRUE
solvent_method          PCM
METHOD  wb97xD
BASIS   cc-pVDZ
MEM_TOTAL       4000
MEM_STATIC      100
$end

$pcm
theory                  IEFPCM
ChargeSeparation        Marcus
StateSpecific           Perturb
$end

$solvent
Dielectric              2.38
OpticalDielectric       1.9599999999999997
$end

$molecule
0 1
N   -0.01011961841039  0.01824575984271  0.21051878677385
C   1.24735574010264  0.75023265696078  -0.03288757696271
N   1.09139146015891  2.09709104326643  -0.06542076101527
N   2.37242521442645  0.09379384367060  -0.13265607243157
C   2.37566672791496  -1.14830750383109  0.01744774668435
H   3.31585768898542  -1.87876384180200  -0.01434545362253
C   -0.14462254987935  2.66447437623288  -0.03203222413041
H   -0.16146613685169  3.78870640870462  -0.02704544490080
N   1.28561147502865  -1.99275787331816  0.06349582722178
C   0.06176384120582  -1.38848408496539  0.08939217806750
N   -1.13329317963977  -2.09146318228738  -0.04583829059091
C   -2.22929197379487  -1.45777439689018  -0.17528818348807
H   -3.06464026558267  -2.07798031725186  -0.27846665437763
N   -2.33318945246678  -0.07533322512782  0.00175307139228
C   -1.23970913567289  0.55699589133919  0.05118818767545
N   -1.33974510675090  1.98137815058730  0.06565404616007
$end
```

# Obtaining results

Once all calculations are run, there are basically two options to get results. The first is generating files with the spectra and rates. The second is using a separate  visualization tool called [Nemoview](https://github.com/LeonardoESousa/nemoview).

For the first method, there are two cases:

1. For absorption spectrum simulations:
    - Run the `nemo` command and choose option 5. Follow the instructions to set the parameters and the spectrum will be generated. 

2 . For photophysical rates:
    - Run the `nemo` command and choose option 6. Follow the instructions to set the parameters. Three files will be generated: an Ensemble file, with data from the ensemble of geometries; a differential_rate file, with the emission spectrum; a rates file, with all available rates (-> denote radiative transitions and ~> denote ISC transitions).

For the second method, which is the one we recommend, the first step is running the `nemo` command with option 7 (Gather ensemble data only)

```
    cd EnsembleS0
    nemo
```

choose option 7 and S0 state. As a result, the file Ensemble_S0_.lx is generated
Do the same for the S1 and T1 states and the files Ensemble_S1_.lx and Ensemble_T1_.lx are generated.

```
[EnsembleS0]$ nemo
#     # ####### #     # #######
##    # #       ##   ## #     #
# #   # #       # # # # #     #
#  #  # #####   #  #  # #     #
#   # # #       #     # #     #
#    ## #       #     # #     #
#     # ####### #     # #######
----------Photophysics---------

Choose your option:

ENSEMBLE SETUP:
        1 - Generate the inputs for the nuclear ensemble calculation
        2 - Run the ensemble calculations
        3 - Check the progress of the calculations
        4 - Abort my calculations (deletes the limit.lx file in the folder you are. This stops the submission of further jobs. It does not kill jobs already on the queue)
ABSORPTION:
        5 - Generate the absorption spectrum
EXCITED STATE PROPERTIES (FLUORESCENCE, PHOSPHORESCENCE, ISC):
        6 - Estimate rates and compute emission spectrum
ENSEMBLE DATA:
        7 - Gather ensemble data only
EXCITON ANALYSIS:
        8 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
7
What is the initial state (S0, S1, T1, S2 ...)? Accepts comma separated values Ex: T1,T2
S0 # I visualize the results obtianed from the ground-state geometry
```

Copy the files **Ensemble_S0_.lx**, **Ensemble_S1_.lx**, **Ensemble_T1_.lx** to your work directory on your computer and follow the tutorial on the [Nemoview](https://github.com/LeonardoESousa/nemoview) page.
