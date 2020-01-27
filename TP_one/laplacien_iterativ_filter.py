from math import pi
import numpy as np
import time
import matplotlib.pyplot as plt

nx,ny = (1000,1000)
xcoords = np.linspace(-5,5,nx)
ycoords = np.linspace(-5,5,ny)

xv,yv = np.meshgrid(xcoords,ycoords)

u = np.cos(1.4*(xv**2+yv**2))
h = plt.contourf(xv,yv,u)
plt.show()

t1 = time.time()
for i in range(1,u.shape[0]-1):
    for j in range(1,u.shape[1]-1):
        u[i,j] = 0.25*(u[i-1,j]+u[i+1,j]+u[i,j-1]+u[i,j+1])
t2 = time.time()
print(f"Temps calcul moyenne : {t2-t1} secondes")

h2 = plt.contourf(xv,yv,u)
plt.show()

u2 = np.cos(1.4*(xv**2+yv**2))
t1 = time.time()
u2[1:-1,1:-1] = 0.25*(u2[0:-2,1:-1]+u2[2:,1:-1]+u2[1:-1,0:-2]+u2[1:-1,2:])
t2 = time.time()
print(f"Temps calcul moyenne : {t2-t1} secondes")

h3 = plt.contourf(xv,yv,u)
plt.show()

diff = u2 - u
h4 = plt.contourf(xv,yv,diff)
plt.show()

u = np.cos(1.4*(xv**2+yv**2))
v2 = np.zeros(u.shape)

t1 = time.time()
v2[1:-1,1:-1] = 0.25*(u[0:-2,1:-1]+u[2:,1:-1]+u[1:-1,0:-2]+u[1:-1,2:])
t2 = time.time()
print(f"Temps calcul moyenne : {t2-t1} secondes")

diff = u2 - v2
h5 = plt.contourf(xv,yv,diff)
plt.show()
