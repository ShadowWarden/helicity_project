# computeQ.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Compute Q for a small patch in the sky
#

import numpy as np
from astropy.io import fits

DEGTORAD = 0.01745


hdul1 = fits.open("data/gtselected_P8_10-20GeV.fits")
hdul2 = fits.open("data/gtselected_P8_20-30GeV.fits")
hdul3 = fits.open("data/gtselected_P8_50-60GeV.fits")

i = np.linspace(0,16,11)

data = hdul1[1].data

E1_all = data.field(0)
theta1_all = data.field(4)
phi1_all = data.field(3)

data = hdul2[1].data

E2_all = data.field(0)
theta2_all = data.field(4)
phi2_all = data.field(3)

data = hdul3[1].data

E3_all = data.field(0)[(data.field(4)) < -80]
theta3_all = data.field(4)[(data.field(4)) < -80]
phi3_all = data.field(3)[(data.field(4)) < -80]

R = np.linspace(0.0,0.25,10)
Q = np.zeros(len(R))
std = np.zeros(len(R))

for j in range(len(R)):
    print("Computing loop for R=",R[j])
    Qsub = np.zeros(len(theta3_all))
    for i in range(len(theta3_all)): 
#        print("Running iteration",i)
        target = np.cos(np.arctan(R[j]))

        theta3 = theta3_all[i]
        phi3 = phi3_all[i]
        v3 = np.array([np.cos(theta3)*np.sin(phi3),np.cos(theta3)*np.cos(phi3),np.sin(theta3)])
        v3 = np.transpose(v3)

        theta1 = theta1_all[np.dot(np.transpose(v3),np.array([np.cos(theta1_all)*np.sin(phi1_all),np.cos(theta1_all)*np.cos(phi1_all),np.sin(theta1_all)])) - target >= 0]
        phi1 = phi1_all[np.dot(np.transpose(v3),np.array([np.cos(theta1_all)*np.sin(phi1_all),np.cos(theta1_all)*np.cos(phi1_all),np.sin(theta1_all)])) - target >= 0]
        v1 = np.array([np.cos(theta1)*np.sin(phi1),np.cos(theta1)*np.cos(phi1),np.sin(theta1)])
        if(len(theta1) != 0):
            len1 = len(theta1)
        else:
            len1 = 1.0
        v1 = np.transpose(v1)
        v1 = np.sum(v1,axis=0)
        
        theta2 = theta2_all[np.dot(np.transpose(v3),np.array([np.cos(theta2_all)*np.sin(phi2_all),np.cos(theta2_all)*np.cos(phi2_all),np.sin(theta2_all)])) - target >= 0]
        phi2 = phi2_all[np.dot(np.transpose(v3),np.array([np.cos(theta2_all)*np.sin(phi2_all),np.cos(theta2_all)*np.cos(phi2_all),np.sin(theta2_all)])) - target >= 0]
        v2 = np.array([np.cos(theta2)*np.sin(phi2),np.cos(theta2)*np.cos(phi2),np.sin(theta2)])
        if(len(theta2) != 0):
            len2 = len(theta2)
        else:
            len2 = 1.0
        v2 = np.transpose(v2)
        v2 = np.sum(v2,axis=0)
        
        Qsub[i] = np.dot(np.cross(v1,v2),v3)/len2/len1 
        Q[j] += np.dot(np.cross(v1,v2),v3)/len2/len1
    std[j] = np.std(Qsub) 

Q /= len(theta3_all)
std /= np.sqrt(len(theta3_all))
print(Q)
