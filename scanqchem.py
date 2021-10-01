import os
import numpy as np
import subprocess
import os
import time

geoms = [i for i in os.listdir('Geometrias') if 'Geometria' in i]
geoms = sorted(geoms, key=lambda pair: float(pair.split('-')[1]))

def pega_geom(file):
    geom = ''
    with open(file) as f:
        for line in f:
            try:
                float(line.split()[1])
                float(line.split()[2])
                float(line.split()[3])
                geom += line
            except:
                pass
    return geom

def watcher():
    p = subprocess.check_output(['squeue','-u','ledso'])
    p = p.decode('UTF-8').rstrip()
    p = p.split()
    num = p.count('ledso')
    return num

with open('../qchem_template.txt', 'r') as f:
    text = f.read()

with open('../qc.sh', 'r') as g:
    batch = g.read()

for geom in geoms:
    #print(geom)
    #num = watcher()
    #while num >= 20:
    #    time.sleep(20)
    #    num = watcher()
    a = geom.split('-')[1]
    new_geom = pega_geom('Geometrias/'+geom)
    newtext = text.replace('DDDD',new_geom)
    try:
        os.mkdir("geom_"+a)
    except:
        pass
    roda = True
    try:
        with open("geom_"+a+"/soc.out",'r') as g:
            for line in g:
                if 'nice day' in line:
                    roda = False
                    break
    except:
        pass
    if roda == True:
        with open('falta.txt','a') as f:
            f.write('qchem -nt 4 geom_'+a+'/soc.in geom_'+a+'/soc.out &\n')    
        print(geom)
        with open("geom_"+a+"/soc.in",'w') as g:
            g.write(newtext)
        #with open("geom_"+a+'/qc.sh','w') as g:
        #    g.write(batch)
    #    os.chdir("geom_"+a)
    #    subprocess.call(['sbatch', 'qc.sh', 'soc.in'])
    #    os.chdir("..")

