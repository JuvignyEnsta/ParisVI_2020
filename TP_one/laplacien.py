from math import pi
import numpy as np
import time
import matplotlib.pyplot as plt

nx,ny = (1000,1000)
xcoords = np.linspace(-5,5,nx)
ycoords = np.linspace(-5,5,ny)

xv,yv = np.meshgrid(xcoords,ycoords)

u = np.cos(1.4*(xv**2+yv**2))/(xv**2+yv**2+0.0001)

h = plt.contourf(xv,yv,u)
plt.show()

v = np.zeros(u.shape)

t1 = time.time()
for i in range(1,u.shape[0]-1):
    for j in range(1,u.shape[1]-1):
        v[i,j] = u[i-1,j]+u[i+1,j]+u[i,j-1]+u[i,j+1]-4*u[i,j]
t2 = time.time()
print(f"Temps calcul laplacien : {t2-t1} secondes")

h2 = plt.contourf(xv,yv,v)
plt.show()

v2 = np.zeros(u.shape)

t1 = time.time()
v2[1:-1,1:-1] = u[0:-2,1:-1]+u[2:,1:-1]+u[1:-1,0:-2]+u[1:-1,2:]-4*u[1:-1,1:-1]
t2 = time.time()
print(f"Temps calcul laplacien vectorise : {t2-t1} secondes")

h3 = plt.contourf(xv,yv,v2)
plt.show()

diff = v2 - v
h3 = plt.contourf(xv,yv,diff)
plt.show()

#t1 = time.time()
###for i in range(1,u.shape[0]-1):
##    v[i] = u[i+1]+u[i-1]-2*u[i]
#t2 = time.time()
