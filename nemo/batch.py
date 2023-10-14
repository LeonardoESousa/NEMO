#!/usr/bin/env python3
import os
import sys
import shutil
import nemo.tools

BATCH_FILE = sys.argv[1]
NPROC = int(sys.argv[2])
NUM = int(sys.argv[3])

def run_batch():
    scripts = [i for i in os.listdir(".") if ".sh" in i]
    for file in scripts:
        shutil.copy(file, "Geometries")
    os.chdir("Geometries")
    the_watcher = nemo.tools.Watcher('.')
    the_watcher.run(BATCH_FILE, NPROC, NUM)
