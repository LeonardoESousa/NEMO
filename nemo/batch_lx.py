#!/usr/bin/env python3
import os
import time
import subprocess
import sys
import shutil
import numpy as np
import nemo.tools

# SETS NUMBER OF SIMULTANEOUS JOBS##############################
def limite():
    numero = np.loadtxt("../limit.lx",encoding='utf-8')
    return numero


###############################################################


##MAIN LOOP####################################################
try:
    BATCH_FILE = sys.argv[1]
    NPROC = int(sys.argv[2])
    NUM = int(sys.argv[3])
    COMMAND = "qchem -nt MMMM"
    scripts = [i for i in os.listdir(".") if ".sh" in i]
    for file in scripts:
        shutil.copy(file, "Geometries")
    os.chdir("Geometries")
    inputs = [i for i in os.listdir(".") if "Geometr" in i and ".com" in i]
    inputs = sorted(inputs, key=lambda pair: float(pair.split("-")[1]))
    FACTOR = 1
    inputs = nemo.tools.watcher(inputs, FACTOR, True)
    if len(inputs) == 0:
        print("No jobs left to run! Goodbye!")
        sys.exit()
    rodando = []
    QUEUE, batch_num = 0, 0
    NEWCOMMAND = ""
    leftover = len(inputs) % NUM
    for input_file in inputs:
        rodando = nemo.tools.watcher(rodando, FACTOR, False)
        nlim = limite()
        NEWCOMMAND += f"{COMMAND} {input_file} {input_file[:-3]}log &\n"
        QUEUE += 1
        if QUEUE == NUM or (QUEUE == leftover and batch_num >= len(inputs) - leftover):
            NEWCOMMAND += "wait"
            NN = NUM * NPROC / NEWCOMMAND.count("qchem")
            NEWCOMMAND = NEWCOMMAND.replace("MMMM", f"{NN:.0f}")
            with open(f"cmd_{batch_num}_.sh", "w",encoding='utf-8') as q:
                q.write(NEWCOMMAND)
            A = subprocess.call(["bash", BATCH_FILE, f"cmd_{batch_num}_.sh"])
            QUEUE = 0
            NEWCOMMAND = ""
        rodando.append(input_file)
        batch_num += 1
        while len(rodando) / NUM >= nlim:
            time.sleep(20)
            rodando = nemo.tools.watcher(rodando, FACTOR, False)
            nlim = limite()
except IndexError:
    print("Something went wrong! Abort.")
###############################################################
