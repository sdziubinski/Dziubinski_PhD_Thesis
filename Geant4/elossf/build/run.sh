#!/bin/bash
echo -n "Start simulation"
rm -f r8520*.*
# isotope As
Z_iso=50
A_iso=108
n=1000
folder="sean_108Sn_80MeV_13ph"
./exe_eloss $Z_iso $A_iso $n
printf "\n"
mkdir $folder
mv pos*.dat $folder/.
mv energy*.dat $folder/.
mv photon*.dat $folder/.
python3 combine_G4outputs.py $folder

