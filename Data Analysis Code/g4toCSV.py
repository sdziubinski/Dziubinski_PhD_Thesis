#Author: Sean Dziubinski
#Date Created: 2/23/24
#Last Edited: 2/23/24

#Objectives:
#1. Convert r8520.photon, and r8520.position into csv files

#--------****--------Imports-------****--------#
import os
import numpy as np
import pandas as pd

#--------****------Directories-----****--------#
#Base directory
directory = "/mnt/c/Users/dziubins/Home/Research/Energy Loss Optical Scintillation System/"
#data directory
data_directory = os.path.join(directory, "Data/results_Se_Z34_A81_Q34_140MeVu/")
#python files directory
python_directory = os.path.join(directory, "Scripts/Python/")
#file names and directories
photon_fname = "r8520.photon"
position_fname = "r8520.position"
photon_directory = os.path.join(data_directory, photon_fname)
position_directory = os.path.join(data_directory, position_fname)

#grab photon data as a list of strings for each line in file
f_photons = open(photon_directory, "r")
photons = f_photons.read()
f_photons.close()

X = []
Y = []
#grab position data
with open(position_directory, 'r') as file:
    # Iterate through each line in the file
    for line in file:
        # Split the line into individual coordinates
        x, y = map(float, line.split())
        X.append(x)
        Y.append(y)



#Create list of lists. The elements in the main list are each particle/run's data. The data is a list of # of photons counted for each PMT.
photon_list = photons.split(" ")
s_0 = 0
particle = []
for s in range(120,len(photon_list),120):
    particle.append([float(ele) for ele in photon_list[s_0:s]])
    s_0 = s

#particle is a list of input data, each simulation is one entry to the model, each simulation has 120 data points
    
#create a pandas dataframe
df = pd.DataFrame(particle)
df2 = pd.DataFrame(data=zip(X,Y),columns=['X','Y'])

#save dataframe to csv file
file_name = input("Enter Photon file name: ")
file_name2 = input("Enter Position file name: ")
df.to_csv(file_name)
df2.to_csv(file_name2)
