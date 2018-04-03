# montecarlo_isometric.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Generates isometric rays on (cos(\theta),\phi)
#

import numpy as np
import matplotlib.pyplot as plt

DEGTORAD = np.pi/180
Nth = 1800
Nphi = 3600

th = np.arccos(np.linspace(1,-1,Nth))
phi = np.linspace(0,2*np.pi,Nphi)
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

def mask(theta1_all,phi1_all,path1FHL):
    theta1_all = theta1_all*DEGTORAD
    phi1_all = phi1_all*DEGTORAD

    hdul_he = fits.open(path1FHL)
    data_he = hdul_he[1].data
    data_he_lon = data_he.field(4)[abs(data_he.field(4)) > 0]*DEGTORAD
    data_he_lat = data_he.field(3)[abs(data_he.field(4)) > 0]*DEGTORAD
    for i in range(len(data_he_lon)):
        v = np.array([np.cos(data_he_lat[i])*np.cos(data_he_lon[i]),np.cos(data_he_lat[i])*np.sin(data_he_lon[i]),np.sin(data_he_lat[i])])
        # Mask to 2 degrees
        target = np.cos(2*DEGTORAD)

        theta1_new = theta1_all[np.dot(v,np.array([np.cos(theta1_all)*np.sin(phi1_all),np.cos(theta1_all)*np.cos(phi1_all),np.sin(theta1_all)])) - target <= 0]
        phi1_new = phi1_all[np.dot(v,np.array([np.cos(theta1_all)*np.sin(phi1_all),np.cos(theta1_all)*np.cos(phi1_all),np.sin(theta1_all)])) - target <= 0]

        theta1_all = theta1_new
        phi1_all = phi1_new

    return theta1_all,phi1_all

def flip(N):
    """ Convert Phi axis from 0->360 to 180->180 """
    n = np.shape(N)[1]
    N2 = np.zeros(np.shape(N))

    N2[:,:180] = N[:,180:]
    N2[:,180:] = N[:,:180]

    return N2

def derivative(N):
    n = np.shape(N)[0]
    m = np.shape(N)[1]

    diff = np.zeros(n-2)
    diffarray = np.zeros(m)

    dth = 180./n

    theta = np.linspace(0,180,n)
    theta = theta[1:-1]

    print(n,m)

    for i in range(n-2):
        costhinv = (np.cos(theta[i]*DEGTORAD))**-1
        for j in range(m):
            diffarray[j] = abs((N[i+2,j] - N[i,j])/(2*dth))

        diff[i] = np.mean(diffarray)
    return diff

def Flux(N,Th):
    nth = int(Th/180.*Nth)
    PhiN = np.sum(N[:nth])/(4*np.pi*(1-np.sin(Th*DEGTORAD)))
    PhiS = np.sum(N[-nth:])/(4*np.pi*(1-np.sin(Th*DEGTORAD)))
    PhiT = (PhiN+PhiS)
    PhiRest = np.sum(N[nth:-nth])/(4*np.pi*np.sin(Th*DEGTORAD))
    return np.array([PhiT,PhiRest])

def Correct_Galactic(N,Th,theta_all,phi_all,Jpd):
    PhiT, PhiRest = Flux(N,Th)
    
    nth = int(Th/180.*Nth)
   
    thi = np.linspace(0,180.,Nth)

    ii = np.where(theta_all > -(thi[nth]-90.0))
    theta_all = np.delete(theta_all,ii)
    phi_all = np.delete(phi_all,ii)

    ii = np.where(-theta_all > (thi[-nth]-90.0))
    theta_all = np.delete(theta_all,ii)
    phi_all = np.delete(phi_all,ii)

    print(PhiRest)

    for i in range(int(PhiRest)):
        Ptheta_all = np.ones(len(theta_all))/len(theta_all)
        theta_choice = np.random.choice(theta_all,p=Ptheta_all)
        ii = np.where(theta_all == theta_choice)
        theta_all = np.delete(theta_all,ii)
        phi_all = np.delete(phi_all,ii)
        
    return theta_all*DEGTORAD, phi_all*DEGTORAD, PhiT+PhiRest
