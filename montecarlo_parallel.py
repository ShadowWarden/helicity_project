# montecarlo_parallel.py
# 
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
# 
# Parellelized montecarlo designed to run on CU's Summit supercomputer
#

import numpy as np
from mpi4py import MPI

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

#if(rank == size-1):
#    len_data = Nruns % size
#    data = np.zeros([len_data,len(R)])
#    for i in range(len_data):
#        genRays(Nrays1,Nrays2,Nrays3)
