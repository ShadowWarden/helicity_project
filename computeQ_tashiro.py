# computeQ.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Compute Q for a small patch in the sky
#

import numpy as np
from astropy.io import fits

DEGTORAD = 0.01745


hdul1 = fits.open("data/FermiData_PASS8_10-40GeV.fits")
hdul2 = fits.open("data/FermiData_PASS8_40-50GeV.fits")
hdul3 = fits.open("data/FermiData_PASS8_50-100GeV.fits")

i = np.linspace(0,16,11)

# Known high energy source file
hdul_he = fits.open("data/gll_psch_v07.fit")
data_he = hdul_he[1].data
data_he_lon = data_he.field(4)[abs(data_he.field(4))>70]*DEGTORAD
data_he_lat = data_he.field(3)[abs(data_he.field(4))>70]*DEGTORAD

data = hdul1[1].data

E1_all = data.field(0)
theta1_all = data.field(4)*DEGTORAD
phi1_all = data.field(3)*DEGTORAD

data = hdul2[1].data

E2_all = data.field(0)
theta2_all = data.field(4)*DEGTORAD
phi2_all = data.field(3)*DEGTORAD

data = hdul3[1].data

E3_all = data.field(0)[abs(data.field(4)) > 70]
theta3_all = data.field(4)[abs(data.field(4)) > 70]*DEGTORAD
phi3_all = data.field(3)[abs(data.field(4)) > 70]*DEGTORAD

# Mask around known sources to 2 degrees
print("Masking Data around 1FHL sources")
for i in range(len(data_he_lon)):
    v = np.array([np.cos(data_he_lat[i])*np.cos(data_he_lon[i]),np.cos(data_he_lat[i])*np.sin(data_he_lon[i]),np.sin(data_he_lat[i])])
    # Mask to 2 degrees
    target = np.cos(2*DEGTORAD)

    theta3_new = theta3_all[np.dot(v,np.array([np.cos(theta3_all)*np.sin(phi3_all),np.cos(theta3_all)*np.cos(phi3_all),np.sin(theta3_all)])) - target <= 0]
    phi3_new = phi3_all[np.dot(v,np.array([np.cos(theta3_all)*np.sin(phi3_all),np.cos(theta3_all)*np.cos(phi3_all),np.sin(theta3_all)])) - target <= 0]
    theta2_new = theta2_all[np.dot(v,np.array([np.cos(theta2_all)*np.sin(phi2_all),np.cos(theta2_all)*np.cos(phi2_all),np.sin(theta2_all)])) - target <= 0]
    phi2_new = phi2_all[np.dot(v,np.array([np.cos(theta2_all)*np.sin(phi2_all),np.cos(theta2_all)*np.cos(phi2_all),np.sin(theta2_all)])) - target <= 0]
    theta1_new = theta1_all[np.dot(v,np.array([np.cos(theta1_all)*np.sin(phi1_all),np.cos(theta1_all)*np.cos(phi1_all),np.sin(theta1_all)])) - target <= 0]
    phi1_new = phi1_all[np.dot(v,np.array([np.cos(theta1_all)*np.sin(phi1_all),np.cos(theta1_all)*np.cos(phi1_all),np.sin(theta1_all)])) - target <= 0]

    theta3_all = theta3_new
    phi3_all = phi3_new
    theta2_all = theta2_new
    phi2_all = phi2_new
    theta1_all = theta1_new
    phi1_all = phi1_new



print("Completed Masking")
R = np.linspace(0.0,20.0,10)
Q = np.zeros(len(R))
std = np.zeros(len(R))

for j in range(len(R)):
    print("Computing loop for R=",R[j])
    Qsub = np.zeros(len(theta3_all))
    for i in range(len(theta3_all)): 
#        print("Running iteration",i)
        target = np.cos(R[j]*DEGTORAD)

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
