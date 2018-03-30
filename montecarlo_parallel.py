# montecarlo_parallel.py
# 
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
# 
# Parellelized montecarlo designed to run on CU's Summit supercomputer
#

import numpy as np
import montecarlo_isometric_headers as head
from astropy.io import fits

Nruns = 50
Nrays1 = 6055
Nrays2 = 3000
Nrays3 = 1000
DEGTORAD = np.pi/180

print("**********Montecarlo_Parallel***********")
print("Simulation Paramters")
print("Nruns=",Nruns)
print("N(E_1)=",Nrays1,",N(E_2)=",Nrays2,",N(E_3)=",Nrays3)

hdul_he = fits.open("data/gll_psch_v07.fit")
data_he = hdul_he[1].data
data_he_lon = data_he.field(4)[abs(data_he.field(4))> 0]*DEGTORAD
data_he_lat = data_he.field(3)[abs(data_he.field(4))> 0]*DEGTORAD

hdul_jpd = fits.open("data/gtselected_exp_map_50_60.fits")

Jpd = (hdul_jpd[0].data)[0]

print(Jpd.sum())

Jpd1 = Jpd / Jpd.sum()

R = np.linspace(2,20,10)

data = []

Q = np.zeros([Nruns,len(R)])
Qnorth = np.zeros([Nruns,len(R)])
Qsouth = np.zeros([Nruns,len(R)])
std = np.zeros([Nruns,len(R)])
stdnorth = np.zeros([Nruns,len(R)])
stdsouth = np.zeros([Nruns,len(R)])

for j in range(Nruns):
   
    theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all = head.genRays(Nrays1,Nrays2,Nrays3,Jpd1)
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
 
    phi3 = phi3_all[abs(theta3_all) > 80*DEGTORAD]
    theta3 = theta3_all[abs(theta3_all) > 80*DEGTORAD]

    Q[j],std[j] = head.computeQ(R,theta1_all,phi1_all,theta2_all,phi2_all,theta3,phi3)
   
    phi3 = phi3_all[(theta3_all) > 80*DEGTORAD]
    theta3 = theta3_all[(theta3_all) > 80*DEGTORAD]

    Qnorth[j],stdnorth[j] = head.computeQ(R,theta1_all,phi1_all,theta2_all,phi2_all,theta3,phi3)
 
    phi3 = phi3_all[-(theta3_all) > 80*DEGTORAD]
    theta3 = theta3_all[-(theta3_all) > 80*DEGTORAD]

    Qsouth[j],stdsouth[j] = head.computeQ(R,theta1_all,phi1_all,theta2_all,phi2_all,theta3,phi3) 


Qmean = np.zeros(len(R))
stdmean = np.zeros(len(R))
QmeanN = np.zeros(len(R))
stdmeanN = np.zeros(len(R))
QmeanS = np.zeros(len(R))
stdmeanS = np.zeros(len(R))
for j in range(len(R)):
    Qmean[j] = np.mean(Q[:,j])
    stdmean[j] = np.sqrt(np.mean(std[:,j]**2))
    QmeanN[j] = np.mean(Qnorth[:,j])
    stdmeanN[j] = np.sqrt(np.mean(stdnorth[:,j]**2))
    QmeanS[j] = np.mean(Qsouth[:,j])
    stdmeanS[j] = np.sqrt(np.mean(stdsouth[:,j]**2))

