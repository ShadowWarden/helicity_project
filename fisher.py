# fisher.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Symbolic fisher matrix computation
#

import numpy as np
import sympy as syp
import montecarlo_isometric_headers as mchead

Npts = 5 

p = []
f = []

def triple_vector(p1,p2,p3):
    vect = p3[0]*(p1[1]*p2[2]-p2[1]*p1[2])- p3[1]*(p1[0]*p2[2]-p2[0]*p1[2])\
            +p3[2]*(p1[0]*p2[1]-p2[0]*p1[1])
    return vect

for i in range(Npts):
    p.append(syp.symbols('p%dx p%dy p%dz' % (i+1,i+1,i+1)))

for i in range(Npts-2):
    for j in range(i+1,Npts-1):
        f.append(triple_vector(p[i],p[j],p[-1]))

M = syp.zeros(Npts)
sigma = syp.symbols('sigma')
for i in range(Npts):
    for j in range(Npts):
        for k in range(len(f)):
            p1 = syp.symbols('p%dx p%dy p%dz' % (i+1,i+1,i+1))
            p2 = syp.symbols('p%dx p%dy p%dz' % (j+1,j+1,j+1))
            print(p1,p2)
            V1 = [syp.diff(f[k],p1[0]),syp.diff(f[k],p1[1]),syp.diff(f[k],p1[2])]
            V2 = [syp.diff(f[k],p2[0]),syp.diff(f[k],p2[1]),syp.diff(f[k],p2[2])]
            M[i,j] += (V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2])

theta = np.zeros(Npts)
phi = np.zeros(Npts)

theta[0],phi[0],theta[1],phi[1],theta[2],phi[2] = mchead.genRays(1,1,1)

v = np.zeros([Npts,3])

for i in range(Npts):
    v[i] = np.array([np.random.uniform(low=-1.0,high=1.0),np.random.uniform(low=-1.0,high=1.0),np.random.uniform(low=-1.0,high=1.0)])
    v[i] /= np.linalg.norm(v[i])

for i in range(Npts):
    for j in range(3):
        M = M.subs({p[i][j] : v[i][j]})

print(M)
print(M.inv())
