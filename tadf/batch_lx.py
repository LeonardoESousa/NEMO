#!/usr/bin/env python3
import os
import time
import subprocess
import numpy as np
import sys
import shutil
from lx.tools import *
        

##CHECKS WHETHER THE JOB IS TWO STEP###########################
def set_factor(file):
    factor = 1
    with open(file, 'r') as f:
        for line in f:
            if 'Link1' in line:
                factor = 2
                break
    return factor
###############################################################

#SETS NUMBER OF SIMULTANEOUS JOBS##############################
def limite():
    numero = np.loadtxt('../limit.lx')
    return numero
###############################################################


##MAIN LOOP####################################################
try:
    batch_file = sys.argv[1]
    scripts = [i for i in os.listdir('.') if '.sh' in i]
    for file in scripts:
        shutil.copy(file,'Geometries')

    os.chdir('Geometries')
    inputs = [i for i in os.listdir('.') if 'Geometr' in i and '.com' in i]
    inputs = sorted(inputs, key=lambda pair: float(pair.split('-')[1]))

    factor = set_factor(inputs[0])

    inputs = watcher(inputs,factor)
    if len(inputs) == 0:
        print('No jobs left to run! Goodbye!')
        sys.exit()
    rodando = []
    for input in inputs:
        rodando = watcher(rodando,factor)
        nlim = limite()
        a = subprocess.call(['bash',batch_file, input])
        rodando.append(input)
        while len(rodando) >= nlim:
            time.sleep(20)
            rodando = watcher(rodando,factor)
            nlim = limite()
except:
    print('Something went wrong! Abort.')           
###############################################################    