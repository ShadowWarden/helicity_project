# motion.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Motion of a charged particle in a helical magnetic field
#

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

q = 1.6e-19
qminus = -1.6e-19
k=9e9

def force(q,B,v,traj1,traj2):
    dist = np.linalg.norm(traj1-traj2)
    vec_norm = (traj1-traj2)/dist
    print(q*np.cross(v,B),k*q*qminus/dist**2*vec_norm)
    return q*np.cross(v,B)-k*q*qminus/dist**2*vec_norm

Zmax = 100.0
Nz = 1001
Nt = 1001
dt = 1e-6

fig = plt.figure()
ax = fig.gca(projection='3d')

z = np.linspace(0,Zmax,Nz)

traj = np.zeros([Nt,3])
traj_elec = np.zeros([Nt,3])
traj_elec[0] = np.array([1e-12,1e-12,0.0])
vel = np.zeros([Nt,3])
vel_elec = np.zeros([Nt,3])
vel_elec[0] = np.array([0,0,1])
vel[0] = np.array([0,0,1])

for i in range(1,Nt):
    B = np.array([np.cos(traj[i-1][2]),np.sin(traj[i-1][2]),0])
    vel[i] = vel[i-1] + dt*(force(q,B,vel[i-1],traj[i-1],traj_elec[i-1]))
    traj[i] = traj[i-1] + dt*(vel[i-1]) 
    B = np.array([np.cos(traj_elec[i-1][2]),np.sin(traj_elec[i-1][2]),0])
    vel_elec[i] = vel_elec[i-1] + dt*(force(qminus,B,vel_elec[i-1],traj_elec[i-1],traj[i-1]))
    traj_elec[i] = traj_elec[i-1] + dt*(vel_elec[i-1]) 


ax.plot(traj[:,0],traj[:,1],traj[:,2],color='r',label='positron')
ax.plot(traj_elec[:,0],traj_elec[:,1],traj_elec[:,2],color='b',label='electron')

plt.show()
