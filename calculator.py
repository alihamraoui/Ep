#!/usr/bin/env python
# coding: utf-8

import math
import sys
import numpy as np
from scipy.constants import e, epsilon_0

###################################################################################
######### CHARMM Parameters
###################################################################################
r0=0.9572 #  minimal distance in Å
kb=450.000 #  kcal.mol-1.Å-2
theteq= 104.5200 # minimal angle in degree (°)
kt= 55.000 #kcal.mol-1.rad-2
q_O= -0.834 
q_H= 0.417
F=332.0716
eps_HH= 0.0460 #kcal.mol-1 
eps_OH= 0.0836 #kcal.mol-1 
eps_OO= 0.1521 #kcal.mol-1
rmin_HH= 0.4490 #Å
rmin_OH= 1.9927 #Å
rmin_OO= 3.5364 #Å
cutoff=14 #Å
L_BOX=30
HALF_BOX=L_BOX/2
###################################################################################
###################################################################################

try:
    filename= sys.argv[1]
    print ("opening {} file :-)".format(filename))
except:
    exit("Please specify an output PDB file!")

#this function extracts the Cartesian coordinates of atomes of each molecule in PDB file.
def get_coors(mol):
    coors = []
    with open(filename, "r") as f:
                for line in f:
                    if line.startswith('ATOM')                     and int(line[23:27])== mol:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        a = line[13:14]
                        coors.append({"atome":a,"x":x, "y": y, "z":z})
    return coors

#this function calculate the distances (dij) between Atomes of each water molecule
def calc_dist( i, j,mol):
        coors= get_coors(mol)
        # i & j are aom numbers.
        a,xi, yi, zi = coors[i].values()
        a,xj, yj, zj = coors[j].values()
        # Calc distance 
        dx = xj - xi
        dy = yj - yi
        dz = zj - zi        
        r2 = dx*dx + dy*dy + dz*dz
        r= math.sqrt(r2)
        return r, dx, dy, dz

# this function Calculate the angle between 3 atomes (ijk)
def calc_angle(i,j,k,mol):
        # i, j & k are atom numbers.
        # d1(d1x, d1y, d1z) and d2(d2x, d2y, d2z) are two Vector make an angle theta
        r1, d1x, d1y, d1z =calc_dist(i,j,mol)
        r2, d2x, d2y, d2z =calc_dist(i,k,mol)
        # norme of d1 & d2 vectors
        nor_d1= math.sqrt(d1x*d1x + d1y*d1y + d1z*d1z)
        nor_d2= math.sqrt(d2x*d2x + d2y*d2y + d2z*d2z)
        # scalaire Prodact
        ij_jk= d1x*d2x + d1y*d2y + d1z*d2z
        # Calc angle theta
        thet=math.acos((ij_jk / (nor_d1 * nor_d2)))
        return thet #get angle in radian

def bonded_energy(mol):
    # Initialize E.
    bonds = 0.0
    ang = 0.0
    # Get Distances OH2-H1 (r1) & OH2-H2 (r2)
    r1, dx1, dy1, dz1 = calc_dist(0,1,mol)
    r2, dx2, dy2, dz2 = calc_dist(0,2,mol)
    # Get theta angle between OH2-H1 & OH2-H2
    thet= calc_angle(0,1,2,mol)
    # Calculate Potential Energy in Kcal/mol
    bonds= kb *(r1 - r0)**2 + kb *(r2 - r0)**2
    angle= kt*(thet-math.radians(theteq))**2
    return bonds + angle 

### this function calculate the unbonded interaction according to cutuff and PBC
def non_bonded_energy( coors1, coors2):
        qi=0
        qj=0
        elec=0.0
        vdw=0.0
        # i & j are aom numbers.
        a1,xi, yi, zi = coors1.values()
        a2,xj, yj, zj = coors2.values()
        # # Calc dx & dy according to mimimum image convention (PBC).
        dx = xj - xi
        if dx > HALF_BOX:
            dx -= L_BOX
        if dx <= -HALF_BOX:
            dx += L_BOX
        dy = yj - yi
        if dy > HALF_BOX:
            dy -= L_BOX
        if dy <= -HALF_BOX:
            dy += L_BOX
        dz = zj - zi 
        if dz > HALF_BOX:
            dz -= L_BOX
        if dz <= -HALF_BOX:
            dz += L_BOX
        r2 = dx*dx + dy*dy + dz*dz
        r= math.sqrt(r2)
        #Cutoff
        if r == 0 or r > cutoff:
            return 0
        ###### 
        if a1 == 'O' : 
            if a2 == 'O' : 
                rmin= rmin_OO 
                epsilon= eps_OO
                qi = q_O
                qj = q_O 
            else : 
                rmin= rmin_OH  
                epsilon= eps_OH
                qi = q_O
                qj = q_H
        elif a1 == 'H' : 
            if a2 == 'O' : 
                rmin= rmin_OH 
                epsilon= eps_OH
                qi = q_H
                qj = q_O 
            else : 
                rmin= rmin_HH  
                epsilon= eps_HH
                qi = q_H
                qj = q_H    
        
        ##calculate Coulomb energy
        elec = ((qi * qj )/r)*F
        ##calculate van der Waals interaction energy
        vdw= epsilon * (np.power(rmin / r, 12) - 2 * np.power(rmin / r, 6))
        ## 
        return   elec+ vdw 


def get_PDB (pdb):
    coordAtomes=[]
    with open (pdb, 'r') as pdb:
        for line in pdb:
            if line.startswith('ATOM'):
                nb_h2o=int(line[23:27])
    for mol in range (1,nb_h2o+1):
        coordAtomes.append(get_coors(mol))
    return coordAtomes


def main():
    Pot_energy=0
    non_bonded=0
    bonded=0
    current_Atomes=[]
    coordAtomes = get_PDB(filename)
    nb_h2o=len(coordAtomes)
    ###### msg
    print("Number of molecule in PDB file is {}".format(nb_h2o))
    print ("processing ...")
    ###### calculate the bonded interactions
    for mol in range (1,nb_h2o+1) :
        bonded+= bonded_energy(mol)
    ##### calculate the non_bonded interactions
    for i in range (0,nb_h2o):
        for j in range (i+1,nb_h2o):
            current_Atomes=[coordAtomes[i],coordAtomes[j]]
            for i in range (0,3):
                for j in range(0,3):
                    non_bonded+= non_bonded_energy(current_Atomes[0][i],current_Atomes[1][j])              
    #### potential energy
    Pot_energy =  bonded + non_bonded 

################### output
    sp='_________'
    d = {'bonded': [ bonded , 'Kcal/mol'],
    "non bonded": [non_bonded , 'Kcal/mol'],
    "potentiel": [ Pot_energy, 'Kcal/mol']}
    print ("{:<14} {:<22} {:<12}".format('Energy','Value','Unit'))
    print ("{:<14} {:<22} {:<12}".format(sp,sp,sp))
    for k, v in d.items():
        lang, perc= v
        print ("{:<14} {:<22} {:<12}".format(k, lang, perc))
##################

main()

