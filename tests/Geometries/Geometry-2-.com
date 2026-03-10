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
6   1.9507019  1.1518843  0.0882933
6   2.7425470  0.0065766  -0.1270286
6   1.9034884  -1.0944585  0.0415173
6   0.4702916  -0.7329234  -0.0401467
6   0.5530476  0.6617104  0.1928025
6   -0.5419833  1.5319212  -0.0408423
6   -1.9438949  1.3044254  -0.0814773
6   -2.5982561  0.0183557  0.0095644
6   -1.8284338  -1.1989435  0.0130459
6   -0.5142947  -1.6341486  -0.0102233
1   2.2597577  2.1328661  -0.1289746
1   3.7320170  -0.0732893  -0.1588954
1   2.3022110  -2.1069076  -0.1104044
1   -0.2756423  2.6414068  -0.1165420
1   -2.6235249  1.9973782  -0.0794899
1   -3.7013795  -0.0070448  0.0527718
1   -2.6489380  -2.0991928  0.0951124
1   -0.2499671  -2.6407620  -0.0559816

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


