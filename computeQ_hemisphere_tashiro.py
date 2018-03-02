# computeQ.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Compute Q for a small patch in the sky
#

import numpy as np
from astropy.io import fits

import computeQ_headers as head

DEGTORAD = 0.01745

# Known high energy source file
hdul_he = fits.open("data/gll_psch_v07.fit")
data_he = hdul_he[1].data
data_he_lon = data_he.field(4)[abs(data_he.field(4))>70]*DEGTORAD
data_he_lat = data_he.field(3)[abs(data_he.field(4))>70]*DEGTORAD

interm = ["data/FermiData_PASS8_20-30GeV.fits",\
        "data/FermiData_PASS8_30-40GeV.fits",\
        "data/FermiData_PASS8_40-50GeV.fits"
        ]

hdul1 = fits.open("data/FermiData_PASS8_10-20GeV.fits")
hdul3 = fits.open("data/FermiData_PASS8_50-60GeV.fits")

Qnorth = np.zeros(len(interm))
Qsouth = np.zeros(len(interm))
Qtotal = np.zeros(len(interm))
stdnorth = np.zeros(len(interm))
stdsouth = np.zeros(len(interm))
stdtotal = np.zeros(len(interm))

for i in range(len(interm)):
    hdul2 = fits.open(interm[i])

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

    for j in range(3):
        if(j == 1):
            print("Running Qnorth computation loop for",interm[i])
            E3_all = data.field(0)[(data.field(4)) > 70]
            theta3_all = data.field(4)[(data.field(4)) > 70]*DEGTORAD
            phi3_all = data.field(3)[(data.field(4)) > 70]*DEGTORAD
        elif(j == 2):
            print("Running Qsouth computation loop for",interm[i])
            E3_all = data.field(0)[-(data.field(4)) > 70]
            theta3_all = data.field(4)[-(data.field(4)) > 70]*DEGTORAD
            phi3_all = data.field(3)[-(data.field(4)) > 70]*DEGTORAD
        else:
            print("Running Qtotal computation loop for",interm[i])
        # Mask around known sources to 2 degrees
        print("Masking Data around 1FHL sources")
        theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all = head.mask(theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all,data_he_lat,data_he_lon)
        print("Completed Masking")

        R = np.linspace(12.5,20.0,4)
        Q = np.zeros(len(R))
        std = np.zeros(len(R))

        Q,std = head.computeQ(theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all,R)
        if(j == 0):
            Qtotal[i] = np.mean(Q)
            stdtotal[i] = np.sqrt(np.mean(std**2))
        if(j == 1):
            Qnorth[i] = np.mean(Q)
            stdnorth[i] = np.sqrt(np.mean(std**2))
        if(j == 2):
            Qsouth[i] = np.mean(Q)
            stdsouth[i] = np.sqrt(np.mean(std**2))

