import os
import time
import subprocess 

def ler():
    with open('raiz.sh', 'r') as f:
        batch = f.read()
    return batch

batch = ler()

def watcher():
    p = subprocess.check_output(['squeue','-u','ledso'])
    p = p.decode('UTF-8').rstrip()
    p = p.split()
    num = p.count('ledso')
    return num

with open('falta.txt', 'r') as f:
    count, ident = 1, 1
    for line in f:
        if count <= 6:
            batch = batch.replace('wait', line[:-1]+' \nwait')
            count += 1
        else:
            with open('batch'+str(ident)+'.sh', 'w') as g:
                g.write(batch)
            num = watcher()
            while num >= 20:
                time.sleep(20)
                num = watcher()
            subprocess.call(['sbatch','batch'+str(ident)+'.sh'])
            print('mandou', 'batch'+str(ident)+'.sh')
            batch = ler()
            batch = batch.replace('wait', line[:-1]+' \nwait')
            count = 2
            ident +=  1

if count > 2:
    with open('batch'+str(ident)+'.sh', 'w') as g:
        g.write(batch)
    subprocess.call(['sbatch','batch'+str(ident)+'.sh'])
    print('mandou', 'batch'+str(ident)+'.sh')



