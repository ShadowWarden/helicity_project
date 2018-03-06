# montecarlo_parallel.py
# 
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
# 
# Parellelized montecarlo designed to run on CU's Summit supercomputer
#

import numpy as np
from mpi4py import MPI
import montecarlo_isometric_headers as head

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

Nruns = 10000
Nrays1 = 2055
Nrays2 = 2000
Nrays3 = 2000

print("**********Montecarlo_Parallel***********")
print("Simulation Paramters")
print("Nruns=",Nruns)
print("N(E_1)=",Nrays1,",N(E_2)=",Nrays2,",N(E_3)=",Nrays3)

R = np.linspace(2,20,10)

data = []

Q = np.zeros([Nruns,len(R)])
std = np.zeros([Nruns,len(R)])

print("Processor",rank,"ist starting data build run")

for j in range(Nruns):
    theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all = head.genRays(Nrays1,Nrays2,Nrays3)

    phi3_all = phi3_all[abs(theta3_all) > 70*DEGTORAD]
    theta3_all = theta3_all[abs(theta3_all) > 70*DEGTORAD]

    Q[j],std[j] = head.computeQ(R,theta1_all,phi1_all,theta2_all,phi2_all,theta3_all,phi3_all)
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
