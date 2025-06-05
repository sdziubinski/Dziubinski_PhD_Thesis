# Dziubinski_PhD_Thesis
This repository affectively serves as the appendix to Sean Dziubinski's PhD thesis at Michigan State University. Contains all code and data relevant to the thesis work.

---

## Geant4
Contains the most recent Geant4 simulation of ELOSS. Can't run code from inside this directory. To use on a computer besides Sean's, update the CMAKE and G4INSTALL paths (among other things).

## Arduino Code
Contains the Arduino code used to control the LED mover motor.

## ELOSS Channel Mapping
Contains any useful information for mapping the individual channels of ELOSS. The excel file provides the HV channel, Cable and pin ID, PMT location, PMT# and DAQ channel for each PMT. Can be used as reference for testing connectivity. Picture references are provided for HV and DAQ channel mapping (currently incomplete). Please update the version number in the excel file name everytime it is changed. Format: ELOSS_ChMap_MMDDYYYYv#.xlsx

## Data Analysis Code
Contains a few Jupyter notebooks used to analyze PMT Calibration data, evaluate different conditions applied to the G4 output photon counts before using the DNN to correct for the position dependent photon collection efficiency.