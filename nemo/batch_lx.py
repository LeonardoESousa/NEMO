#!/usr/bin/env python3
import os
import time
import subprocess
import numpy as np
import sys
import shutil
import nemo.tools        

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
    nproc      = int(sys.argv[2])
    num        = int(sys.argv[3])
    command    = 'qchem -nt MMMM'
    scripts    = [i for i in os.listdir('.') if '.sh' in i]
    for file in scripts:
        shutil.copy(file,'Geometries')
    os.chdir('Geometries')
    inputs = [i for i in os.listdir('.') if 'Geometr' in i and '.com' in i]
    inputs = sorted(inputs, key=lambda pair: float(pair.split('-')[1]))
    factor = set_factor(inputs[0])
    inputs = nemo.tools.watcher(inputs,factor,True)
    if len(inputs) == 0:
        print('No jobs left to run! Goodbye!')
        sys.exit()
    rodando = []
    queue, batch_num = 0, 0
    newcommand = ''
    leftover = len(inputs)%num
    for input in inputs:
        rodando = nemo.tools.watcher(rodando,factor,False)
        nlim    = limite()
        newcommand += f'{command} {input} {input[:-3]}log &\n'
        queue += 1
        if queue == num or (queue == leftover and batch_num >= len(inputs) - leftover):
            newcommand += 'wait'
            NN = num*nproc/newcommand.count('qchem')
            newcommand = newcommand.replace('MMMM',f'{NN:.0f}')
            with open(f'cmd_{batch_num}_.sh', 'w') as q:
                q.write(newcommand)
            a = subprocess.call(['bash',batch_file, f'cmd_{batch_num}_.sh'])
            queue = 0
            newcommand = ''
        rodando.append(input)
        batch_num += 1
        while len(rodando)/num >= nlim:
            time.sleep(20)
            rodando = nemo.tools.watcher(rodando,factor,False)
            nlim = limite()
except:
    print('Something went wrong! Abort.')           
###############################################################    