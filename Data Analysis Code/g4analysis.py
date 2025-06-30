#Author: Sean Dziubinski
#Date Created: 2/20/24
#Last Edited: 2/20/24

#Make GeantAnalysis2.py into a function to be imported as a module
import os
import numpy as np
import matplotlib.pyplot as plt
import geoarray
from scipy.optimize import curve_fit

# Define the Gaussian function 
def gaussian(x, amplitude, mean, std_dev):
    return amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2)) 

# Define conversion from PMT # to distance (mm)
def pmt_to_dist(x, board):
    # change PMT # to distance
    if board == "bottom":
        offset = 15 #mm
        width = -300 #mm
        sign = 1
    elif board == "right":
        offset = 7.5 #mm
        width = -150 #mm
        sign = 1
    elif board == "top":
        offset = -15 #mm
        width = 300 #mm
        sign = -1
    else:
        offset = -7.5 #mm
        width = 150 #mm
        sign = -1
    gap = 60 #mm
    return width + offset + sign*x*gap

#main
def analyze(photon_file, position_file, num):

    #------****-----Photon Counting-----****------#
    f_photons = open(photon_file, "r")
    photons = f_photons.read()
    f_photons.close()

    list_photons = photons.split(" ")
    S1, S2, S3, S4, particle = geoarray.make_geometrical_array(list_photons, num)
    bottom = particle[num][0:10]
    max_bottom = bottom.index(np.max(bottom))
    print(max_bottom)
    max_bottom_mm = pmt_to_dist(max_bottom, "bottom")
    right = particle[num][10:15]
    top = particle[num][15:25]
    left = particle[num][25:30]

    bottom_x = list(range(0,10,1))
    right_x = list(range(0,5,1))
    top_x = list(range(0,10,1))
    left_x = list(range(0,5,1))

    #perform curve fitting for gaussian fucntion
    bottom_mm = pmt_to_dist(np.array(bottom_x), "bottom")
    top_mm = pmt_to_dist(np.array(top_x), "top")
    right_mm = pmt_to_dist(np.array(right_x), "right")
    left_mm = pmt_to_dist(np.array(left_x), "left")
    initial_guess_x = [np.max(bottom),bottom_mm[bottom.index(np.max(bottom))] , 50]
    parameters, _ = curve_fit(gaussian, bottom_mm, bottom, p0=initial_guess_x, maxfev=5000)
    paramtop, _ = curve_fit(gaussian, top_mm, top, p0=initial_guess_x, maxfev=5000)

    per_diff = (paramtop[0] - parameters[0]) / ((parameters[0] + paramtop[0])/2)
    initial_guess_y = [np.max(top), per_diff*75 , 100]
    paramright, _ = curve_fit(gaussian, right_mm, right, p0=initial_guess_y, maxfev=5000)
    paramleft, _ = curve_fit(gaussian, left_mm, left, p0=initial_guess_y, maxfev=5000)

    #create the x and y data to represent the fit
    fit_x = np.arange(-300,301,1)
    fit_y = fit_x/2
    fit_bottom = gaussian(fit_x, *parameters)
    fit_right = gaussian(fit_y, *paramright)
    fit_top = gaussian(fit_x, *paramtop)
    fit_left = gaussian(fit_y, *paramleft)

    #-------****------Positioning-----****--------#
    X=[]
    Y=[]
    with open(position_file, 'r') as file:
        # Iterate through each line in the file
        for line in file:
            # Split the line into individual coordinates
            x, y = map(float, line.split())
            X.append(x)
            Y.append(y)

    bottom_mean = parameters[1]
    top_mean = paramtop[1]
    wb = 1/abs(parameters[0])
    wt = 1/abs(paramtop[0])
    right_mean = paramright[1]
    left_mean = paramleft[1]
    wr = 1/abs(paramright[0])
    wl = 1/abs(paramleft[0])

    ave_x = round((wb*bottom_mean + wt*top_mean)/(wb+wt))
    est_x = 0
    if (wb>wt):
        est_x = round(bottom_mean, 2)   
    else:
        est_x = round(top_mean, 2)
    inc_x = X[num]
    ave_y = round((wl*left_mean + wr*right_mean)/(wl+wr))
    est_y = 0
    if (wl>wr):
        est_y = round(left_mean, 2)   
    elif (wl<wr):
        est_y = round(right_mean, 2)
    inc_y = Y[num]

    return  inc_x,  est_x, max_bottom_mm