## Version 1.3.1

1 - Bug fixes.

## Version 1.3.0

1 - Added EOM-CCSD capabilities do NEMO

2 - Introduced a tool for quickly checking energies and susceptibilities in QChem log files (nemo -c file.log)

3 - Added a tool for retrieving geometries from qchem files (nemo -g file.log)

4 - ptSS + ptLR as standard state specific correction.

5 - Paralellized geometry generation.


## Version 1.2.0

1 - Introduced Ensemble class to make it easier to get ensemble results without NEMOview.

2 - Accepts calculations in gas phase.

3 - Bug fixes.

## Version 1.1.0

1 - During ensemble generation, **NEMO** now computes the adjacency matrix of the optimized molecule and compares it with the adjacency matrix of each geometry sampled in the ensemble. If the matrices do not match, i.e. if one or more bonds are broken, the sampled geometry is rejected. This helps reduce issues with unphysical geometries and negative transition energies.

2 - Bug fixes.