# montecarlo_isometric.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Generates isometric rays on (cos(\theta),\phi)
#

import numpy as np
import matplotlib.pyplot as plt
import subprocess
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

DEGTORAD = np.pi/180

costh = np.linspace(-1,1,100001)
phi = np.linspace(0,2*np.pi,100001)

Pcosth = np.ones(len(costh))/len(costh)
Pphi = np.ones(len(phi))/len(phi)

Nrays1 = 705
Nrays2 = 400
Nrays3 = 100

R = np.linspace(2.0,20.,10)
# Nruns per processor
Nruns = 5

result = subprocess.run(['git','rev-parse', 'HEAD'], stdout=subprocess.PIPE)

if(rank == 0):
    print("**************Monte-Carlo Simulation**************")
    print("Omkar H. Ramachandran, Laboratory for Atmospheric and Space Physics")
    print("Revision Hash:",result.stdout)
    print("Parameters:")
    print("Nruns:",Nruns*size,"\nNrays1:",Nrays1,"\nNrays2:",Nrays2,"\nNrays3:",Nrays3)


def genRays(Nrays1,Nrays2,Nrays3):
    phi1_all = np.zeros(Nrays1)
    theta1_all = np.zeros(Nrays1)

    phi2_all = np.zeros(Nrays2)
    theta2_all = np.zeros(Nrays2)

    phi3_all = np.zeros(Nrays3)
    theta3_all = np.zeros(Nrays3)
    for i in range(Nrays1):
        phi1_all[i] = np.random.choice(phi,p=Pphi)
        theta1_all[i] = np.arccos(np.random.choice(costh,p=Pcosth))

        # Convert from theta to latitude
    theta1_all = (np.pi/2.-theta1_all)
    for i in range(Nrays2):
        phi2_all[i] = np.random.choice(phi,p=Pphi)
        theta2_all[i] = np.arccos(np.random.choice(costh,p=Pcosth))

    theta2_all = (np.pi/2.-theta2_all)

    for i in range(Nrays3):
        phi3_all[i] = np.random.choice(phi,p=Pphi)
        theta3_all[i] = np.arccos(np.random.choice(costh,p=Pcosth))

    theta3_all = (np.pi/2.-theta3_all)

    return theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all


def computeQ(R,theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all): 
    Q = np.zeros(len(R))
    std = np.zeros(len(R))

    for j in range(len(R)):
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

    Q /= len(theta3_all)
    std /= np.sqrt(len(theta3_all))
    return Q, std

Q = np.zeros([Nruns,len(R)])
std = np.zeros([Nruns,len(R)])
print("Processor",rank,"is starting data build run")
for j in range(Nruns):
    theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all = genRays(Nrays1,Nrays2,Nrays3)

    phi3_all = phi3_all[abs(theta3_all) > 70*DEGTORAD]
    theta3_all = theta3_all[abs(theta3_all) > 70*DEGTORAD]

    Q[j],std[j] = computeQ(R,theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all)
print("Processor",rank,": Computation complete!")
Qmean = np.zeros(len(R))
stdmean = np.zeros(len(R))
for j in range(len(R)):
    Qmean[j] = np.mean(Q[:,j])
    stdmean[j] = np.sqrt(np.mean(std[:,j]**2)/(Nruns-1.))

print("Gather at rank 0")
if(rank != 0):
    comm.send(np.array([Qmean,stdmean]),dest=0)
else:
    flag = 0
    for i in range(1,size):
        flag += 1
        A = comm.recv(source=i)
        Qmean = Qmean*(flag)/(flag+1)+A[0]/(flag+1)
        stdmean = np.sqrt(stdmean**2*flag/(flag+1)+A[1]**2/(flag+1))
    A = np.array([Qmean,stdmean])
    np.savetxt("Q.txt",A)

