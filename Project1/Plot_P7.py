import matplotlib.pyplot as plt
import numpy as np
import sys


def Analytical(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)


X = np.linspace(0, 1, 1000)
U = Analytical(X)

plt.plot(X, U, label="Analytical")

for i in range(len(sys.argv) - 1):
    data = np.loadtxt("General_solution_" + sys.argv[i + 1] + ".txt")
    x = data[:, 0]
    v = data[:, 1]

    plt.plot(x, v, linestyle="dashed", label=f"n_steps={sys.argv[i+1]}")


plt.xlabel("x")
plt.ylabel("Solutions")
plt.legend()
plt.savefig("P7_fig.pdf")
plt.show()
