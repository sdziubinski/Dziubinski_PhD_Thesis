import numpy as np
import fnmatch
import os
import sys

mainpath = '/Users/sean/Desktop/Work/Geant4/elossf/build/'
folder = sys.argv[1]

n_threads = len(fnmatch.filter(os.listdir(mainpath+folder), 'pos*.dat'))

photons = []
for thread in range(n_threads):
    filename = folder+"/photon{0}.dat".format(thread)
    photons.append(np.loadtxt(filename, dtype='int'))

positions = []
for thread in range(n_threads):
    filename = folder+"/pos{0}.dat".format(thread)
    positions.append(np.loadtxt(filename, dtype='int'))

energies = []
for thread in range(n_threads):
    filename = folder+"/energy{0}.dat".format(thread)
    energies.append(np.loadtxt(filename))
    
t_photon = np.vstack(photons)

with open(folder+"/t_photon.dat", 'w') as filehandle:
    for line in t_photon:
        #print(line)
        for col in line:
            #print(str(col))
            filehandle.write(str(col)+" ")
        filehandle.write("\n")
filehandle.close()

t_pos = np.vstack(positions)

with open(folder+"/t_position.dat", 'w') as filehandle2:
    for line in t_pos:
        for col in line:
            filehandle2.write(str(col)+" ")
        filehandle2.write("\n")
filehandle2.close()

#t_energy = np.vstack(energies)

#with open(folder+"/t_energy.dat", 'w') as filehandle3:
#    for line in t_energy:
#        for col in line:
#            filehandle3.write(str(col)+" ")
#        filehandle3.write("\n")
#filehandle3.close()