# computeQ.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Compute Q for a small patch in the sky
#

import numpy as np
from astropy.io import fits

DEGTORAD = 0.01745

def computeQ(fname,Theta_min,theta_max,dtheta,Phi_min,phi_max,dphi):
    print("Computing Q for data from",fname)
    hdul = fits.open(fname)

<<<<<<< HEAD
hdul1 = fits.open("data/gtselected_P8_10-20GeV.fits")
hdul2 = fits.open("data/gtselected_P8_20-30GeV.fits")
hdul3 = fits.open("data/gtselected_P8_50-60GeV.fits")

dtheta = 10
dphi = 20

Theta_min = -90
Phi_min = 0

theta_max = 90
phi_max = 360

theta_min = np.linspace(Theta_min,theta_max,int((theta_max-Theta_min)/dtheta)) 
phi_min = np.linspace(Phi_min,phi_max,int((phi_max-Phi_min)/dphi))

Q = np.zeros([int(len(theta_min)),int(len(phi_min))])

data = hdul1[1].data

E1_all = data.field(0)
theta1_all = data.field(4)
phi1_all = data.field(3)

data = hdul2[1].data

E2_all = data.field(0)
theta2_all = data.field(4)
phi2_all = data.field(3)

data = hdul3[1].data

E3_all = data.field(0)[abs(data.field(4)) > 70]
theta3_all = data.field(4)[abs(data.field(4)) > 70]
phi3_all = data.field(3)[abs(data.field(4)) > 70]
=======
    theta_min = np.linspace(Theta_min,theta_max,int((theta_max-Theta_min)/dtheta)) 
    phi_min = np.linspace(Phi_min,phi_max,int((phi_max-Phi_min)/dphi))

    Q = np.zeros([int(len(theta_min)),int(len(phi_min))])

    data = hdul[1].data

    E_all = data.field(0)
        
    theta_all = data.field(4)

    phi_all = data.field(3)

    time_all = data.field(9)
>>>>>>> 17068feca86ddb080d146cc22dc28f534e520410

# Find indices where theta_min < theta < theta_min + dtheta
# and phi_min < phi < phi_min + dphi
#y = np.where(abs(theta_all) - theta_min < dtheta)
#x = np.where(abs(phi_all[y]) - phi_min < dphi)

# Compute new vectors that only hold relevent information

<<<<<<< HEAD
for i in range(len(theta_min)):
    for j in range(len(phi_min)):
        print("Running iteration",[i,j])
        
        theta1 = theta1_all[((theta1_all > theta_min[i]) & (theta1_all < theta_min[i]+dtheta)) & ((phi1_all > phi_min[j]) & (phi1_all < phi_min[j]+dphi))]
        phi1 = phi1_all[((theta1_all > theta_min[i]) & (theta1_all < theta_min[i]+dtheta)) & ((phi1_all > phi_min[j]) & (phi1_all < phi_min[j]+dphi))]
        v1 = np.array([np.cos(theta1)*np.sin(phi1),np.cos(theta1)*np.cos(phi1),np.sin(theta1)])
        v1 = np.transpose(v1)
        v1 = np.sum(v1,axis=0)

        theta2 = theta2_all[((theta2_all > theta_min[i]) & (theta2_all < theta_min[i]+dtheta)) & ((phi2_all > phi_min[j]) & (phi2_all < phi_min[j]+dphi))]
        phi2 = phi2_all[((theta2_all > theta_min[i]) & (theta2_all < theta_min[i]+dtheta)) & ((phi2_all > phi_min[j]) & (phi2_all < phi_min[j]+dphi))]
        v2 = np.array([np.cos(theta2)*np.sin(phi2),np.cos(theta2)*np.cos(phi2),np.sin(theta2)])
        v2 = np.transpose(v2)
        v2 = np.sum(v2,axis=0)

        theta3 = theta3_all[((theta3_all > theta_min[i]) & (theta3_all < theta_min[i]+dtheta)) & ((phi3_all > phi_min[j]) & (phi3_all < phi_min[j]+dphi))]
        phi3 = phi3_all[((theta3_all > theta_min[i]) & (theta3_all < theta_min[i]+dtheta)) & ((phi3_all > phi_min[j]) & (phi3_all < phi_min[j]+dphi))]
        v3 = np.array([np.cos(theta3)*np.sin(phi3),np.cos(theta3)*np.cos(phi3),np.sin(theta3)])
        v3 = np.transpose(v3)
        v3 = np.sum(v3,axis=0)

        if(len(theta3) != 0):
            Q[i,j] += np.dot(np.cross(v1,v2),v3)/len(theta3)
=======
    for i in range(len(theta_min)):
        for j in range(len(phi_min)):
            print("Running iteration",[i,j])
            E = E_all[(theta_all <= theta_min[i] + dtheta) & (theta_all >= theta_min[i]) & (phi_all <= phi_min[j] + dphi) & (phi_all >= phi_min[j])]
            Efinal = E[(E >= 10e3) & (E <= 40e3)]
            N = np.histogram(Efinal,bins=3,range=(10e3,40e3))


            E1 = Efinal[(Efinal >= N[1][0]) & (Efinal <= N[1][1])]
            E2 = Efinal[(Efinal >= N[1][1]) & (Efinal <= N[1][2])]
            E3 = Efinal[(Efinal >= N[1][2]) & (Efinal <= N[1][3])]
        
            theta = theta_all[(theta_all <= theta_min[i] + dtheta) & (theta_all >= theta_min[i]) & (phi_all <= phi_min[j] + dphi) & (phi_all >= phi_min[j])]
            phi = phi_all[(theta_all <= theta_min[i] + dtheta) & (theta_all >= theta_min[i]) & (phi_all <= phi_min[j] + dphi) & (phi_all >= phi_min[j])]
            time = time_all[(theta_all <= theta_min[i] + dtheta) & (theta_all >= theta_min[i]) & (phi_all <= phi_min[j] + dphi) & (phi_all >= phi_min[j])]

            theta_final = theta[(E >= 10e3) & (E <= 40e3)]
            phi_final = phi[(E >= 10e3) & (E <= 40e3)]
            time_final = time[(E >= 10e3) & (E <= 40e3)]



        # Bin energies from 10 GeV to 60 GeV. Do 3 bins for now

        # Demarcate energies by bin
            theta1 = (theta_final[(Efinal >= N[1][0]) & (Efinal <= N[1][1])])*DEGTORAD
            phi1 = phi_final[(Efinal >= N[1][0]) & (Efinal <= N[1][1])]*DEGTORAD
            time1 = time_final[(Efinal >= N[1][0]) & (Efinal <= N[1][1])]
            v1 = np.array([np.cos(theta1)*np.sin(phi1),np.cos(theta1)*np.cos(phi1),np.sin(theta1)])
            v1 = np.transpose(v1)
            v1 = np.sum(v1,axis=0)

            theta2 = (theta_final[(Efinal >= N[1][1]) & (Efinal <= N[1][2])])*DEGTORAD
            phi2 = phi_final[(Efinal >= N[1][1]) & (Efinal <= N[1][2])]*DEGTORAD
            time2 = time_final[(Efinal >= N[1][1]) & (Efinal <= N[1][2])]
            v2 = np.array([np.cos(theta2)*np.sin(phi2),np.cos(theta2)*np.cos(phi2),np.sin(theta2)])
            v2 = np.transpose(v2)
            v2 = np.sum(v2,axis=0)

            theta3 = (theta_final[(Efinal >= N[1][2]) & (Efinal <= N[1][3])])*DEGTORAD
            phi3 = phi_final[(Efinal >= N[1][2]) & (Efinal <= N[1][3])]*DEGTORAD
            time3 = time_final[(Efinal >= N[1][2]) & (Efinal <= N[1][3])]
            v3 = np.array([np.cos(theta3)*np.sin(phi3),np.cos(theta3)*np.cos(phi3),np.sin(theta3)])
            v3 = np.transpose(v3)
            v3 = np.sum(v3,axis=0)
            
            Qlocal = np.dot(np.cross(v1,v2),v3)

            Q[i,j] += Qlocal
    return Q

files = ['data/lat_photon_weekly_w493_p302_v001.fits'\
        ,'data/lat_photon_weekly_w494_p302_v001.fits'\
        ,'data/lat_photon_weekly_w495_p302_v001.fits'\
        ,'data/lat_photon_weekly_w496_p302_v001.fits'\
        ,'data/lat_photon_weekly_w497_p302_v001.fits'\
        ,'data/lat_photon_weekly_w498_p302_v001.fits'\
        ,'data/lat_photon_weekly_w499_p302_v001.fits'\
        ,'data/lat_photon_weekly_w500_p302_v001.fits']

Q = computeQ(files[0],-50,50,10,0,360,36)

for i in range(1,5):
    Qt = computeQ(files[i],-50,50,10,0,360,36)
    Q += Qt
>>>>>>> 17068feca86ddb080d146cc22dc28f534e520410
