import numpy as np
import matplotlib.pyplot as plt

A=1
def potential(x,y,z):
    return A*(-x**2-y**2+2*z**2)

def eletric_field(x,y):

    return np.array([2*A*x,2*A*y,-4*A*z])

Lx = 5
Ly = 5
Lz = 5
N = 30
x = np.linspace(-Lx,Lx,N)
y = np.linspace(-Ly,Ly,N)
z = np.linspace(-Lz,Lz,N)
rx,rz = np.meshgrid(x,z)

V = np.zeros((N,N),float)
for i in range(len(rx.flat)):
    V.flat[i] = potential(rx.flat[i],0,rz.flat[i])

Ez,Ex = np.gradient(-V)

plt.figure(figsize=(10,10))
plt.contourf(rx,rz,V)
plt.quiver(rx,rz,Ex,Ez)
plt.show()