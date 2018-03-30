# montecarlo_isometric.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Generates isometric rays on (cos(\theta),\phi)
#

import numpy as np
import matplotlib.pyplot as plt

DEGTORAD = np.pi/180

th = np.arccos(np.linspace(1,-1,1800))
phi = np.linspace(0,2*np.pi,3600)
dphi = phi[1] - phi[0]

def genRays(Nrays1,Nrays2,Nrays3,Jpd):
    """
    INPUT: Number of rays in bin 1, bin 2, bin 3 and the joint probability density.
    OUTPUT: theta1_all, phi1_all, theta2_all, phi2_all, theta3_all, phi3_all
    """
    phi1_all = np.zeros(Nrays1)
    theta1_all = np.zeros(Nrays1)

    phi2_all = np.zeros(Nrays2)
    theta2_all = np.zeros(Nrays2)

    phi3_all = np.zeros(Nrays3)
    theta3_all = np.zeros(Nrays3)

    # Compute the marginal probability for phi.
    phi_marginal = np.zeros(len(phi))
    for i in range(len(phi_marginal)):
        phi_marginal[i] = sum(Jpd[:,i])
    phi_marginal /= phi_marginal.sum()
    print(sum(phi_marginal))

    print("Generating rays for E1.")
    for i in range(Nrays1):
        phi1_all[i] = np.random.choice(phi,p=phi_marginal)
        ptheta_g_phi = np.zeros(len(th))
        ii = int(phi1_all[i]/dphi)
        for j in range(len(th)):
            ptheta_g_phi[j] = Jpd[j,ii]/phi_marginal[ii]
        ptheta_g_phi = ptheta_g_phi/ptheta_g_phi.sum()
        theta1_all[i] = np.random.choice(th,p=ptheta_g_phi)

        # Convert from theta to latitude
    theta1_all = (np.pi/2.-theta1_all)
    print("Generating rays for E2.")
    for i in range(Nrays2):
        phi2_all[i] = np.random.choice(phi,p=phi_marginal)
        ptheta_g_phi = np.zeros(len(th))
        ii = int(phi2_all[i]/dphi)
        for j in range(len(th)):
            ptheta_g_phi[j] = Jpd[j,ii]/phi_marginal[ii]
        ptheta_g_phi = ptheta_g_phi/ptheta_g_phi.sum()
        theta2_all[i] = np.random.choice(th,p=ptheta_g_phi)

    theta2_all = (np.pi/2.-theta2_all)

    print("Generating rays for E3.")
    for i in range(Nrays3):
        phi3_all[i] = np.random.choice(phi,p=phi_marginal)
        ptheta_g_phi = np.zeros(len(th))
        ii = int(phi3_all[i]/dphi)
        for j in range(len(th)):
            ptheta_g_phi[j] = Jpd[j,ii]/phi_marginal[ii]
        ptheta_g_phi = ptheta_g_phi/ptheta_g_phi.sum()
        theta3_all[i] = np.random.choice(th,p=ptheta_g_phi)

    theta3_all = (np.pi/2.-theta3_all)

    return theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all

R = np.linspace(0.0,20.,10)

def computeQ(R,theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all): 
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
            v1 = np.transpose(v1)
            if(len(theta1) != 0):
                len1 = len(theta1)
            else:
                len1 = 1.0
            v1 = np.sum(v1,axis=0)
            
            theta2 = theta2_all[np.dot(np.transpose(v3),np.array([np.cos(theta2_all)*np.sin(phi2_all),np.cos(theta2_all)*np.cos(phi2_all),np.sin(theta2_all)])) - target >= 0]
            phi2 = phi2_all[np.dot(np.transpose(v3),np.array([np.cos(theta2_all)*np.sin(phi2_all),np.cos(theta2_all)*np.cos(phi2_all),np.sin(theta2_all)])) - target >= 0]
            v2 = np.array([np.cos(theta2)*np.sin(phi2),np.cos(theta2)*np.cos(phi2),np.sin(theta2)])
            v2 = np.transpose(v2)
            if(len(theta2) != 0):
                len2 = len(theta2)
            else:
                len2 = 1.0
            v2 = np.sum(v2,axis=0)
            
            Qsub[i] = np.dot(np.cross(v1,v2),v3)/len2/len1 
            Q[j] += np.dot(np.cross(v1,v2),v3)/len2/len1
        std[j] = np.std(Qsub) 
    if(len(theta3_all) != 0):
        Q /= len(theta3_all)
        std /= np.sqrt(len(theta3_all))
    else:
        Q = np.zeros(len(R))
        std = np.zeros(len(R))
    return Q, std

#Nruns = 1000
#Q = np.zeros([Nruns,len(R)])
#std = np.zeros([Nruns,len(R)])

#for j in range(Nruns):
#print("Running Loop",j)
#    theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all = genRays(Nrays1,Nrays2,Nrays3)
#
#    phi3_all = phi3_all[abs(theta3_all) > 70*DEGTORAD]
#    theta3_all = theta3_all[abs(theta3_all) > 70*DEGTORAD]
#
#    Q[j],std[j] = computeQ(R,theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all)
#
#Qmean = np.zeros(len(R))
#stdmean = np.zeros(len(R))
#for j in range(len(R)):
#    Qmean[j] = np.mean(Q[:,j])
#    stdmean[j] = np.sqrt(np.mean(std[:,j]**2)/(Nruns-1.))
