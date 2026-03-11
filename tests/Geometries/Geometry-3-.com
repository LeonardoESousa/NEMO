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
6   1.9488172  1.1281006  0.0244125
6   2.7777608  0.0258653  -0.0964186
6   1.9187278  -1.0638598  -0.1071018
6   0.5040414  -0.7027053  0.0720193
6   0.5935080  0.6356204  0.0313895
6   -0.5264892  1.6294414  0.1893091
6   -1.9149407  1.2816270  -0.0488260
6   -2.5726368  -0.0365157  -0.2126380
6   -1.9401207  -1.2544829  -0.0075567
6   -0.5908022  -1.6839286  0.1969405
1   2.1190285  2.1706557  -0.2310780
1   3.9467737  0.1543386  -0.1750892
1   2.0716159  -2.0435541  -0.0540367
1   -0.4213297  2.6887182  0.1933019
1   -2.3762304  2.2917135  -0.1251856
1   -3.6337117  0.0328138  -0.2694317
1   -2.6611076  -2.0730243  0.0360945
1   -0.2579012  -2.6637740  0.1847277

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


