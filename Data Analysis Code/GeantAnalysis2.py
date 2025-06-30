#Author: Sean Dziubinski
#Date Created: 2/20/24
#Last Edited: 2/20/24

#Objectives:
#1. Plot color map of PMT counts
#2. Plot PMT counts on each side as a histogram
#3. Fit photon counts histogram with a gaussian and extract mean
#4. Plot incident position vs measured position

#--------****--------Imports-------****--------#
import os
import numpy as np
import matplotlib.pyplot as plt
import geoarray
from scipy.optimize import curve_fit

#--------****------Directories-----****--------#
#Base directory
directory = "/mnt/c/Users/dziubins/Home/Research/Energy Loss Optical Scintillation System/"
#data directory
data_directory = os.path.join(directory, "Data/results_Se_Z34_A81_Q34_140MeVu/")
#python files directory
python_directory = os.path.join(directory, "Scripts/Python/")
#file names and directories
photon_fname = "photon.txt"
position_fname = "r8520.position"
photon_directory = os.path.join(data_directory, photon_fname)
position_directory = os.path.join(data_directory, position_fname)

#-------****--------Functions--------****-------#

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


#------****-----Photon Counting-----****------#

#grab photon data as a list of strings for each line in file
f_photons = open(photon_directory, "r")
photons = f_photons.read()
f_photons.close()

#ask what particle number you want to look at
num = int(input("Particle number to look at: "))

#take string data and put it into a 2D numpy array of floats that has the same shape as the optical readout geometry
#plotting S1, S2, S3, S4 gives a color plot of the photon counting in a geometrical view
list_photons = photons.split(" ")
S1, S2, S3, S4, particle = geoarray.make_geometrical_array(list_photons, num)
bottom = particle[num][0:10]
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
print(bottom_mm[bottom.index(np.max(bottom))])
initial_guess_x = [np.max(bottom),bottom_mm[bottom.index(np.max(bottom))] , 50]
parameters, _ = curve_fit(gaussian, bottom_mm, bottom, p0=initial_guess_x)
paramtop, _ = curve_fit(gaussian, top_mm, top, p0=initial_guess_x)

per_diff = (paramtop[0] - parameters[0]) / ((parameters[0] + paramtop[0])/2)
initial_guess_y = [np.max(top), per_diff*75 , 100]
paramright, _ = curve_fit(gaussian, right_mm, right, p0=initial_guess_y)
paramleft, _ = curve_fit(gaussian, left_mm, left, p0=initial_guess_y)

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
with open(position_directory, 'r') as file:
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


#-------****--------Outputs-------****--------#
#print("Fitted Bottom Mean: ", parameters[1])
#print("Fitted Top Mean: ", paramtop[1])
# X Position
print("Weighted Averaged X Value: ", round((wb*bottom_mean + wt*top_mean)/(wb+wt)), " mm")
if (wb>wt):
    print("Estimated X Position: ", round(bottom_mean, 2))
else:
    print("Estimated X Position: ", round(top_mean, 2))
print("Incident X Position: ", X[num], " mm")
# Y Position
print("Weighted Averaged Y Value: ", round((wl*left_mean + wr*right_mean)/(wl+wr)), " mm")
if (wl>wr):
    print("Estimated Y Position: ", round(left_mean, 2))
elif (wl<wr):
    print("Estimated Y Position: ", round(right_mean, 2))
print("Estimated Y Position by Using Top and Bottom: ", per_diff*75, " mm")
print("Incident Y Position: ", Y[num], " mm")

#print(list_photons)
floats = [float(ele) for ele in list_photons[0:2400]]
print("Maximum photoelectrons collected: ", np.max(floats))
print("Minimum photoelectrons collected: ", np.min(floats))
print("average photoelectrons collected: ", np.mean(floats))
print("stdev photoelectrons collected: ", np.std(floats))


#-------****-------Plotting-------****--------#
fig, [[ax1, ax2],[ax3, ax4]] = plt.subplots(nrows=2, ncols=2)


ax1.set_title("Sector 1")
ax1.imshow(S1)
ax2.set_title("Sector 2")
ax2.imshow(S2)
ax3.set_title("Sector 3")
ax3.imshow(S3)
ax4.set_title("Sector 4")
ax4.imshow(S4)

#plot each side of the sector as a bar graph instead of color plot
fig2, [[ax5, ax6],[ax7, ax8]] = plt.subplots(nrows=2, ncols=2)
fig2.tight_layout()
ax5.set_title("Bottom S1")
ax5.set_xlim(-300,300)
ax5.set_xlabel("X Position (mm)")
ax5.bar(bottom_mm, bottom, width=40)
ax6.set_title("Right S1")
ax6.set_xlim(-150,150)
ax6.set_xlabel("X Position (mm)")
ax6.bar(right_mm, right, width=20)
ax7.set_title("Top S1")
ax7.set_xlim(-300,300)
ax7.set_xlabel("X Position (mm)")
ax7.bar(top_mm, top, width=40)
ax8.set_title("Left S1")
ax8.set_xlim(-150,150)
ax8.set_xlabel("X Position (mm)")
ax8.bar(left_mm, left, width=20)

#plot the gaussian fitting
ax5.scatter(bottom_mm, bottom, color="yellow")
ax5.plot(fit_x, fit_bottom, color="red")
ax6.scatter(right_mm, right, color="yellow")
ax6.plot(fit_y, fit_right, color="red")
ax7.scatter(top_mm, top, color="yellow")
ax7.plot(fit_x, fit_top, color="red")
ax8.scatter(left_mm, left, color="yellow")
ax8.plot(fit_y, fit_left, color="red")

plt.figure(3)
plt.hist(floats, bins=100)

plt.show(block=True)







