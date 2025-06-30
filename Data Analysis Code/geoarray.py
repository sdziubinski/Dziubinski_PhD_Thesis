import numpy as np

def make_geometrical_array(photon_list, n):
    #Create list of lists. The elements in the main list are each particle/run's data. The data is a list of # of photons counted for each PMT.
    s_0 = 0
    particle = []
    for s in range(120,len(photon_list),120):
        particle.append([float(ele) for ele in photon_list[s_0:s]])
        s_0 = s

    #template for ELOSS optical readout geometry
    rows, cols = (7, 12)
    blank_map = [[0]*cols]*rows

    #For each particle, seperate data into each sector.      
    S1_ion = particle[n][0:30]
    S2_ion = particle[n][30:60]
    S3_ion = particle[n][60:90]
    S4_ion = particle[n][90:120]
    #Create arrays to represent each sector's readout/PMT geometry
    S1 = np.array(blank_map)
    S2 = np.array(blank_map)
    S3 = np.array(blank_map)
    S4 = np.array(blank_map)
    #Plot the bottom 10 PMT array
    S1[-1] = [0] + S1_ion[0:10] + [0]
    S2[-1] = [0] + S2_ion[0:10] + [0]
    S3[-1] = [0] + S3_ion[0:10] + [0]
    S4[-1] = [0] + S4_ion[0:10] + [0]
    #Plot the right 5 PMT array
    for r in range(1,rows-1):
        right1 = S1_ion[10:15]
        right2 = S2_ion[10:15]
        right3 = S3_ion[10:15]
        right4 = S4_ion[10:15]
        right1.reverse()
        right2.reverse()
        right3.reverse()
        right4.reverse()
        S1[r][11] = right1[r-1]
        S2[r][11] = right2[r-1]
        S3[r][11] = right3[r-1]
        S4[r][11] = right4[r-1]
    #Plot the top 10 PMT array
    top1 = S1_ion[15:25]
    top2 = S2_ion[15:25]
    top3 = S3_ion[15:25]
    top4 = S4_ion[15:25]
    top1.reverse()
    top2.reverse()
    top3.reverse()
    top4.reverse()
    S1[0] = [0] + top1 + [0]
    S2[0] = [0] + top2 + [0]
    S3[0] = [0] + top3 + [0]
    S4[0] = [0] + top4 + [0]
    #Plot the left 5 PMT array
    for r in range(1,rows-1):
        S1[r][0] = S1_ion[25:30][r-1]
        S2[r][0] = S2_ion[25:30][r-1]
        S3[r][0] = S3_ion[25:30][r-1]
        S4[r][0] = S4_ion[25:30][r-1]
    return S1, S2, S3, S4, particle
   