from re import M
import numpy as np
import matplotlib.pyplot as plt

V0_d2=9.65
B0=9.65e1
q=1
m=1

x_0=20
v_0=25
z_0=0

w_0=q*B0/m
w_z2=2*q*V0_d2/m

w_p=(w_0+np.sqrt(w_0**2-2*w_z2))/2
w_m=(w_0-np.sqrt(w_0**2-2*w_z2))/2

A_p=(v_0 + w_m*x_0)/(w_m - w_p)
A_m=-(v_0 + w_m*x_0)/(w_m - w_p)

phi_m=0
phi_p=0

def f(t):
    x=A_p*np.cos(-(w_p*t + phi_p)) + A_m*np.cos(-(w_m*t + phi_m))
    y=A_p*np.sin(-(w_p*t + phi_p)) + A_m*np.sin(-(w_m*t + phi_m))
    z=z_0*np.cos(np.sqrt(w_z2)*t)

    return x,y,z

t=np.linspace(0,0.1,1000)

x,y,z=f(t)

plt.plot(x,y)
plt.show()
