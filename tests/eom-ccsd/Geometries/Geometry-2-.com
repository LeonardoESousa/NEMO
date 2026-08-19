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
6   1.9394166  1.0857192  0.0906671
6   2.7773842  0.0229088  0.0078427
6   1.8996456  -1.1813541  -0.1080158
6   0.5659947  -0.6750579  -0.0197417
6   0.5831181  0.7406569  0.0446518
6   -0.5650408  1.6236073  -0.0018002
6   -1.9256581  1.2576998  -0.0115367
6   -2.5643651  -0.0370245  0.0019250
6   -1.9261635  -1.3216079  -0.0241625
6   -0.5610379  -1.5181277  0.0990197
1   2.2168800  2.1784940  0.2426179
1   3.8708371  0.0285842  -0.3835765
1   2.2209261  -2.2882193  -0.2083451
1   -0.4080133  2.6593313  -0.4032034
1   -2.7682503  2.1810081  -0.1757321
1   -3.5766289  -0.0275177  -0.2808954
1   -2.8143195  -2.0517670  0.1041367
1   -0.3364163  -2.5926637  0.2180251

$end

