###################################################################################
#########  CHARMM parameters
###################################################################################
r0=0.9572 # Å
kb=450.000 #  kcal.mol-1.Å-2
theteq= 104.5200 #°
kt= 55.000 #kcal.mol-1.rad-2
#################################################################

import math
import sys

coors = []
# Open PDB file & get cartesian coordinates 
with open(sys.argv[1], "r") as f:
                for line in f:
                    if line.startswith('ATOM'):
                        iCode = line[25]
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coors.append({"x":x, "y": y, "z":z})
#Calculate distances between atomes
def calc_dist( i, j):
        # i & j are atom numbers.
        xi, yi, zi = coors[i].values()
        xj, yj, zj = coors[j].values()
        # Calc distance 
        dx = xj - xi
        dy = yj - yi
        dz = zj - zi        
        r2 = dx*dx + dy*dy + dz*dz
        r= math.sqrt(r2)
        return r, dx, dy, dz
    
def calc_angle( i,j,k):
        # i & j and k are atom numbers.
        # d1(d1x, d1y, d1z) and d2(d2x, d2y, d2z) are two Vector make an theta angle
        r1, d1x, d1y, d1z =calc_dist(i,j)
        r2, d2x, d2y, d2z =calc_dist(i,k)
        #norm of vectors d1 & d2
        nor_d1= math.sqrt(d1x*d1x + d1y*d1y + d1z*d1z)
        nor_d2= math.sqrt(d2x*d2x + d2y*d2y + d2z*d2z)
        #scalaire Prodact of d1 X d2
        ij_jk= d1x*d2x + d1y*d2y + d1z*d2z
        # Calc angle theta
        thet=math.acos((ij_jk / (nor_d1 * nor_d2)))
        return math.degrees(thet) #get angle in degree

def calc_E():
    # Initialize E.
    E = 0.0
    # Get Distances OH2-H1 (r1) & OH2-H2 (r2)
    r1, dx1, dy1, dz1 = calc_dist(0,1)
    r2, dx2, dy2, dz2 = calc_dist(0,2)
    # Get theta angle between OH2-H1 & OH2-H2
    thet= calc_angle(0,1,2)
    # Calculate Potential Energy 
    E = kb *(r1 - r0)**2 + kb *(r2 - r0)**2
    E+= kt *( thet-theteq)**2
    return E
    
print("l'energie potentielle de ce systeme est de : {} Kcal/mol".format(calc_E()))
