import pandas as pd
import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Open binary file
data= pa.cx_mat()
data.load("./data/U_double.bin", pa.arma_binary)
U= np.array(data,dtype=np.clongdouble)
print(U[0,0])
print(np.shape(U))

U_0= U[:,0]
#reshape U_0
h=0.005
dt=2.5e-5
N=int(1/h-2)
U_0= np.reshape(U_0,(N,N))


P= np.conjugate(U)*U
P_cube= np.real(np.reshape(P,(N,N,len(U[0,:]))))

V_data= pa.cx_mat()
V_data.load("./data/V.bin", pa.arma_binary)
V= np.array(V_data,dtype=np.clongdouble)






x_points = np.linspace(0, 1, N)
y_points = np.linspace(0, 1, N)
x,y = np.meshgrid(x_points,y_points)
t_points= np.arange(0,1+dt,dt)
#adjust colorcale to max value of P_cube
plt.contourf(P_cube[:,:,0], 100, cmap='jet')
#adjust colorcale to max value of P_cube

plt.show()

for i in range(len(P_cube[0,0,:])):
    plt.plot(t_points[i],np.sum(P_cube[:,:,i]),'o')
plt.show()

# Some settings
fontsize = 12
t_min = t_points[0]
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

# Create figure
fig = plt.figure()
ax = plt.gca()

# Create a colour scale normalization according to the max z value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(P_cube[:,:,0]))

# Plot the first frame
img = ax.imshow(P_cube[:,:,0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("jet"), norm=norm)

# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("P", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(P_cube[:,:,i]))
    img.set_norm(norm)

    # Update z data
    img.set_data((P_cube[:,:,i]))

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(P_cube[0,0,:]), 2), repeat=True, blit=0)

# Run the animation!
plt.show() 
