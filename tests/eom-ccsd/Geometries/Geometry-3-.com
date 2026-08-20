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
6   1.9736077  1.1884370  0.1573985
6   2.7707657  -0.0217177  0.0352793
6   1.9877736  -1.0589501  -0.0419632
6   0.5374085  -0.7039273  -0.1731269
6   0.5197759  0.6087987  0.0226234
6   -0.5849115  1.6235355  -0.0990227
6   -1.9608044  1.2612440  -0.1072756
6   -2.5248356  0.0310425  0.0419776
6   -1.9384745  -1.3008460  0.1641581
6   -0.6315023  -1.6226426  -0.1304039
1   2.2426950  2.0684524  0.2469028
1   4.1100745  -0.2536418  0.2725875
1   2.4262553  -2.1143421  -0.2111393
1   -0.1759199  2.6819300  -0.0325469
1   -2.6971616  2.1005399  0.0687878
1   -3.5778971  0.1618191  0.3535643
1   -2.6358149  -2.1521953  0.5284505
1   -0.3057512  -2.6371293  0.2639295

$end

