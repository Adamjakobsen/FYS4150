import numpy as np

data = np.loadtxt('../positions_rk4.txt')
x = data[:,0]
y = data[:,1]
z = data[:,2]


def animation_3d():
    """
    3d animation of data using matplotlib.animation



    x: x data
    y: y data
    z: z data
    """
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import animation

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d([-1.0, 1.0])
    ax.set_xlabel('X')
    ax.set_ylim3d([-1.0, 1.0])
    ax.set_ylabel('Y')
    ax.set_zlim3d([-1.0, 1.0])
    ax.set_zlabel('Z')
    ax.set_title('3D Test')

    ax.plot(x, y, z, label='parametric curve')
    ax.legend()

    plt.show()
    



