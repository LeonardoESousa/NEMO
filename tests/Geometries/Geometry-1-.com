$rem
mem_total             100000
mem_static            2000
method                lrc-wpbe
omega                 236
basis                 6-31+g(d,p)
rpa                     false
cis_n_roots             5
cis_singlets            true
cis_triplets            true
calc_soc                true
STS_MOM                 true
MAX_CIS_CYCLES          200
MAX_SCF_CYCLES          200
$end


$molecule
0 1
6   1.9565000  1.1416199  0.0344570
6   2.7802161  0.0175052  -0.1641292
6   1.9201116  -1.1533898  -0.0542209
6   0.5359966  -0.7258667  0.1155509
6   0.5436831  0.7033633  0.1133370
6   -0.5578864  1.6849502  0.0726208
6   -1.9107093  1.2484829  0.0129635
6   -2.5690027  -0.0140426  -0.2617069
6   -1.9444340  -1.2824696  -0.0060940
6   -0.5715485  -1.6439192  0.1361374
1   2.2843582  2.2445118  -0.1287414
1   4.0008227  -0.0825135  -0.2906663
1   2.4982547  -2.1710093  0.0749154
1   -0.4494857  2.8690147  0.0861159
1   -2.6765596  2.1146768  0.1893838
1   -3.6650495  0.1408791  -0.0911736
1   -2.6242902  -2.2364740  0.0039728
1   -0.4596092  -2.6793569  0.1634059

$end


@@@

$molecule
read
$end

$rem
mem_total             100000
mem_static            2000
method                lrc-wpbe
omega                 236
basis                 6-31+g(d,p)
rpa                     false
cis_n_roots             5
cis_singlets            true
cis_triplets            true
CIS_RELAXED_DENSITY     TRUE
solvent_method          PCM
MAX_CIS_CYCLES          200
MAX_SCF_CYCLES          200
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


@@@

$molecule
read
$end

$rem
mem_total             100000
mem_static            2000
method                lrc-wpbe
omega                 236
basis                 6-31+g(d,p)
rpa                     false
cis_n_roots             5
cis_singlets            true
cis_triplets            false
MAX_CIS_CYCLES          200
MAX_SCF_CYCLES          200
CIS_DER_NUMSTATE        6
CALC_NAC                true
$end

$derivative_coupling
   comment
    0 1 2 3 4 5
$end


