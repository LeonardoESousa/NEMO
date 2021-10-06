from tadf.tools import *
from tadf.analysis import *
import numpy as np

np.set_printoptions(suppress=True,formatter={'float': '{: 6.6f}'.format})

freqlog = 'freq.out'

#G, atomos = pega_geom('freq.out')
#f, m = pega_freq(freqlog)

#NC = pega_modosLP(G, freqlog)

#busca_input(freqlog)

file = 'Geometry-1-.log'
#socs = soc_t1(file,'1',0)

#print(socs/(0.12398/1000))

singlets, triplets, oscs, ind_s, ind_t = pega_energias(file)
zero = ['0']
zero.extend(ind_s)
n_triplet = 1
MS0       = pega_dipolos(file, zero,"Electron Dipole Moments of Ground State",0)            
MS0resto  = pega_dipolos(file, zero,"Transition Moments Between Ground and Singlet Excited States",0) 
MS0       = np.vstack((MS0,MS0resto))
MT1       = pega_dipolos(file, ind_t,"Electron Dipole Moments of Triplet Excited State",n_triplet)            
MT1resto  = pega_dipolos(file, ind_t,"Transition Moments Between Triplet Excited States",n_triplet)           
MT1       = np.vstack((MT1,MT1resto))            
MT1[[0,n_triplet]] = MT1[[n_triplet,0]]


print(np.shape(MT1))
#ms       = moment(file,singlets,triplets,MS0,MT1,n_triplet)



