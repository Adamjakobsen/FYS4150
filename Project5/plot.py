
import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

#Defining some parameters
h=0.005
dt=2.5e-5
N=int(1/h-1)

x_points = np.linspace(0, 1, N)
y_points = np.linspace(0, 1, N)
x,y = np.meshgrid(x_points,y_points)
t_points= np.arange(0,1+dt,dt)


def get_data(problem):
    # Open binary file
    data= pa.cx_mat()
    data.load("./data/U_input_"+problem+ ".bin", pa.arma_binary)
    U= np.array(data,dtype=np.clongdouble)
    print(np.shape(U))
    return U




def plot_V():
    V_data= pa.cx_mat()
    V_data.load("./data/V.bin", pa.arma_binary)
    V= np.array(V_data,dtype=np.clongdouble)
    plt.contourf(np.real(V))
    plt.savefig("./fig/"+"V.pdf")
    plt.close()


def plot_colormaps(U,problem):

    #Wave function t=0, t=0.001 and t=0.002
    U_0= U[:,0]
    U_1= U[:,int(0.001/dt)]
    U_2= U[:,-1]
    #reshape
    U_0= np.reshape(U_0,(N,N))
    U_1= np.reshape(U_1,(N,N))
    U_2= np.reshape(U_2,(N,N))
    #Probability function at t=0, t=0.001 and t=0.002
    P_0= np.conjugate(U_0)*U_0
    P_1= np.conjugate(U_1)*U_1
    P_2= np.conjugate(U_2)*U_2

    #Plotting contourf subplots real(U0) and im(U0)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('Real and imaginary part of the wave function at t=0')
    ax1.contourf(x,y,np.real(U_0))
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    

    ax1.set_title("Real part")
    ax2.contourf(x,y,np.imag(U_0))
    ax2.set_title("Imaginary part")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    fig.colorbar(ax1.contourf(x,y,np.real(U_0)),ax=ax2)
    plt.savefig("./fig/"+"U0_"+problem+".pdf")
    plt.close()

    #Plotting contourf subplots real(U1) and im(U1)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('Real and imaginary part of the wave function at t=0.001')
    ax1.contourf(x,y,np.real(U_1))
    ax1.set_title("Real part")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax2.contourf(x,y,np.imag(U_1))
    ax2.set_title("Imaginary part")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    fig.colorbar(ax1.contourf(x,y,np.real(U_1)),ax=ax2)
    plt.savefig("./fig/"+"U1_"+problem+".pdf")
    plt.close()

    #Plotting contourf subplots real(U2) and im(U2)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('Real and imaginary part of the wave function at t=0.002')
    ax1.contourf(x,y,np.real(U_2))
    ax1.set_title("Real part")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    
    ax2.contourf(x,y,np.imag(U_2))
    ax2.set_title("Imaginary part")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    fig.colorbar(ax1.contourf(x,y,np.real(U_2)),ax=ax2)
    plt.savefig("./fig/"+"U2_"+problem+".pdf")
    plt.close()

    #Plotting contourf subplots P0, P1 and P2
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.suptitle('Probability function at t=0, t=0.001 and t=0.002')
    ax1.contourf(x,y,P_0)
    ax1.set_title("t=0")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax2.contourf(x,y,P_1)
    ax2.set_title("t=0.001")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax3.contourf(x,y,P_2)
    ax3.set_title("t=0.002")
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    fig.colorbar(ax1.contourf(x,y,P_0),ax=ax3)
    plt.savefig("./fig/"+"P_"+problem+".pdf")
    plt.close()
    









def plot_detector(P_cube,problem):

    idx_detector= int(0.8/h)
    plt.plot(y_points,P_cube[idx_detector,:,-1])
    plt.savefig("./fig/"+"detector_"+problem+ ".pdf")
    plt.close()

def plot_deviation(P_cube,problem):
    for i in range(len(P_cube[0,0,:])):
        plt.plot(t_points[i],np.abs(np.sum(P_cube[:,:,i])-1),'or')

    plt.xlabel("t[s]")
    plt.ylabel("Deviation of total probability")
    plt.savefig("./fig/"+"deviation_"+problem+".pdf")
    plt.close()

def make_animation(P_cube):
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
        img.set_data((P_cube[:,:,i].T))

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(P_cube[0,0,:]), 2), repeat=True, blit=0)

    # Run the animation!
    anim.save("./animation.mp4",writer="ffmpeg", fps=30, dpi=300)
    plt.show()


def main():

    problem = sys.argv[1]
    U= get_data(problem)

    
    P= np.conjugate(U)*U
    P_cube= np.real(np.reshape(P,(N,N,len(U[0,:]))))

    if problem=="P7_1" or problem=="P7_2":

        plot_deviation(P_cube,problem)
    
    if problem=="P8":
        
        plot_colormaps(U,problem)
    
    if problem=="P9_single" or problem=="P9_double" or problem=="P9_triple":
        plot_detector(P_cube,problem)
        
    



if __name__ == "__main__":
    main()

