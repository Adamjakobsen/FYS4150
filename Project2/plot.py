import numpy as np
import matplotlib.pyplot as plt
import sys

problem = sys.argv[1]

if problem == "P5":

    data1 = np.loadtxt("NOperationsP5Dense.txt", float)
    data2 = np.loadtxt("NOperationsP5Tridiag.txt", float)

    x = np.zeros(len(data1))
    y1 = np.zeros(len(data1))
    y2 = np.zeros(len(data2))

    for i in range(len(data1)):
        x[i] = data1[i][0]
        y1[i] = data1[i][1]
        y2[i] = data2[i][1]

    plt.plot(np.log(x), np.log(y1), label="Dense")
    plt.plot(np.log(x), np.log(y2), label="Tridiagonal")

    plt.xlabel("log(N)")
    plt.ylabel("log(Number of iterations)")
    plt.title("Number of iterations as function of N")
    plt.legend()
    plt.savefig("P5_fig_N200.pdf")


elif problem == "P6":

    # Get data from P6eigenboys.txt, first column is v1,second column is v2, third column is v3
    data = np.loadtxt("P6eigenboys.txt", float)

    v1 = data[:, 0]
    v2 = data[:, 1]
    v3 = data[:, 2]

    # add 0 to start of v1, v2, v3
    v1 = np.insert(v1, 0, 0)
    v2 = np.insert(v2, 0, 0)
    v3 = np.insert(v3, 0, 0)

    # add 0 to end of v1, v2, v3
    v1 = np.append(v1, 0)
    v2 = np.append(v2, 0)
    v3 = np.append(v3, 0)

    x = np.linspace(0, 1, len(v1))

    # Plot the data
    plt.plot(x, v1, label="v1")
    plt.plot(x, v2, label="v2")
    plt.plot(x, v3, label="v3")
    plt.xlabel(r"$\^x$")
    plt.ylabel("v")
    plt.title("Solutions to the eigenvalue problem")
    plt.legend()
    plt.savefig("P6_fig" + f"_{len(x)}" + ".pdf")
