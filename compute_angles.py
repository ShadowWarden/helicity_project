# compute_angles.py
#
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
# 
# Compute angles given three vectors
# 

import numpy as np
from astropy.io import fits
import montecarlo_isometric_headers as head

def compute_angles(v1,v2,v3):
    """
    INPUT: v1, v2, v3 (1x3 Row vectors)
    OUTPUT: alpha, beta, gamma (angle between v3 and v2, v1 and v2, and 
    v3 and v1 respectively) in radians
    """
    cosalpha = np.dot(v3-v1,v2-v1)/(np.linalg.norm(v3-v1)*np.linalg.norm(v2-v1))
    cosbeta = np.dot(v1-v3,v2-v3)/(np.linalg.norm(v1-v3)*np.linalg.norm(v2-v3))
    cosgamma = np.dot(v1-v2,v3-v2)/(np.linalg.norm(v1-v2)*np.linalg.norm(v3-v2))
    return np.sort(np.array([np.arccos(cosalpha), np.arccos(cosbeta), np.arccos(cosgamma)]))[1]

DEGTORAD = 0.01745

hdul_jpd = fits.open("data/gtselected_exp_map_50_60.fits")
Jpd = (hdul_jpd[0].data)[0]

Jpd = Jpd/np.sum(Jpd)

Jpd = np.ones(np.shape(Jpd))/np.sum(np.ones(np.shape(Jpd)))

#hdul1 = fits.open("data/FermiData_PASS8_30-40GeV.fits")
#hdul2 = fits.open("data/FermiData_PASS8_40-50GeV.fits")
#hdul3 = fits.open("data/FermiData_PASS8_50-60GeV.fits")

# Known high energy source file

#data = hdul1[1].data
#
#theta1_all = data.field(4)[-(data.field(4)) > 70]*DEGTORAD
#phi1_all = data.field(3)[-(data.field(4)) > 70]*DEGTORAD
#
#data = hdul2[1].data
#
#theta2_all = data.field(4)[-(data.field(4)) > 70]*DEGTORAD
#phi2_all = data.field(3)[-(data.field(4)) > 70]*DEGTORAD
#
#data = hdul3[1].data
#
#theta3 = data.field(4)[-(data.field(4)) > 80]*DEGTORAD
#phi3 = data.field(3)[-(data.field(4)) > 80]*DEGTORAD
#

theta1,phi1,theta2,phi2,theta3_all,phi3_all = head.genRays(5000,3000,1000,Jpd)

theta1_all = theta1[abs(theta1) > 70*DEGTORAD]
phi1_all = phi1[abs(theta1) > 70*DEGTORAD]
theta2_all = theta2[abs(theta2) > 70*DEGTORAD]
phi2_all = phi2[abs(theta2) > 70*DEGTORAD]
theta3 = theta3_all[abs(theta3_all) > 80*DEGTORAD]
phi3 = phi3_all[abs(theta3_all) > 80*DEGTORAD]

flag = 0

R = np.linspace(3,10,10)
Q = np.zeros(10)
std = np.zeros(10)

for l in range(10):
    flag = 0
    angles = np.zeros(100000)
    Qsub = np.zeros(len(theta3))
    print("Running iteration for R =",R[l])
    for i in range(len(theta3)):
        v1f = np.zeros(3)
        v2f = np.zeros(3)
        v3f = np.zeros(3)

        target = np.cos(R[l]*DEGTORAD)
        v3 = np.array([np.cos(theta3[i])*np.sin(phi3[i]),np.cos(theta3[i])*np.cos(phi3[i]),np.sin(theta3[i])])
        theta1 = theta1_all[np.dot(np.transpose(v3),np.array([np.cos(theta1_all)*np.sin(phi1_all),np.cos(theta1_all)*np.cos(phi1_all),np.sin(theta1_all)])) - target >= 0]
        phi1 = phi1_all[np.dot(np.transpose(v3),np.array([np.cos(theta1_all)*np.sin(phi1_all),np.cos(theta1_all)*np.cos(phi1_all),np.sin(theta1_all)])) - target >= 0]
        theta2 = theta2_all[np.dot(np.transpose(v3),np.array([np.cos(theta2_all)*np.sin(phi2_all),np.cos(theta2_all)*np.cos(phi2_all),np.sin(theta2_all)])) - target >= 0]
        phi2 = phi2_all[np.dot(np.transpose(v3),np.array([np.cos(theta2_all)*np.sin(phi2_all),np.cos(theta2_all)*np.cos(phi2_all),np.sin(theta2_all)])) - target >= 0]
        for j in range(len(theta2)):
            v2 = np.array([np.cos(theta2[j])*np.sin(phi2[j]),np.cos(theta2[j])*np.cos(phi2[j]),np.sin(theta2[j])])
            for k in range(len(theta1)):
                if(flag == 99999):
                    break
                v1 = np.array([np.cos(theta1[k])*np.sin(phi1[k]),np.cos(theta1[k])*np.cos(phi1[k]),np.sin(theta1[k])])
                angles[flag] = compute_angles(v1,v2,v3)/np.pi*180
                if((angles[flag] > 1) & (angles[flag] < 11)):
                    v1f += v1
                    v2f += v2
                    v3f = v3 
                flag += 1
                if(flag % 1000 == 0):
                    print("Computed Angle",flag)
        Qsub[i] = np.dot(np.cross(v1f,v2f),v3f)/len(theta1)/len(theta2)
    angles = angles[:flag]
    Q[l] = np.sum(Qsub)
    std[l] = np.std(Qsub)

Q /= len(theta3)
std /= np.sqrt(len(theta3))

print(Q)
