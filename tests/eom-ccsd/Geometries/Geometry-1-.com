$rem
method eom-ccsd
basis cc-pvdz
mem_total 200000
mem_static 2000
cc_memory 50000
ee_singlets             5
ee_triplets             5
cc_trans_prop           2
calc_soc                true
solvent_method          PCM
EOM_DAVIDSON_MAXVECTORS 300
EOM_DAVIDSON_MAX_ITER   300
$end

$trans_prop
state_list
ref
ee_singlets 0 0
end_list
calc dipole linmom soc opdm_norm

state_list
ref
ee_triplets 0 0
end_list
calc dipole linmom soc opdm_norm

state_list
ee_singlets 0 0
ee_singlets 0 0
end_list
calc dipole linmom opdm_norm

state_list
ee_singlets 0 0
ee_triplets 0 0
end_list
calc dipole linmom soc opdm_norm

state_list
ee_triplets 0 0
ee_triplets 0 0
end_list
calc dipole linmom  opdm_norm
$end

$pcm
theory                  IEFPCM
ChargeSeparation        Marcus
StateSpecific           Perturb
$end

$solvent
Dielectric              3.0
OpticalDielectric       1.9599999999999997
$end

$molecule
0 1
6   1.9611652  1.1450516  -0.1450743
6   2.7837312  0.0186978  -0.0796412
6   1.9255275  -1.0868623  0.1000639
6   0.5478679  -0.7423380  0.0705901
6   0.5165601  0.6535761  0.0224950
6   -0.5752637  1.5701480  0.2150331
6   -1.9040875  1.2887039  -0.0771032
6   -2.5568943  0.0318012  -0.1028485
6   -1.9126616  -1.2859716  0.0116610
6   -0.5955469  -1.6502646  -0.0020400
1   2.2570402  2.2881252  -0.1089511
1   3.8666613  0.0410531  -0.1636328
1   2.3286167  -1.9977637  0.1103980
1   -0.3873593  2.7528231  0.3976268
1   -2.4124550  2.2855743  -0.0381378
1   -3.6786730  -0.0192040  -0.2771689
1   -2.5965133  -2.0828390  0.0936294
1   -0.5049468  -2.5774982  -0.1580103

$end

