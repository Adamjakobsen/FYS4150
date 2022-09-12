import numpy as np
import matplotlib.pyplot as plt
import sys

n = ["10", "100", "1000"]


def Exact(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)


fig, ax = plt.subplots(2, 1)
for i in n:
    # Loading Data file
    data2 = np.loadtxt("General_solution_" + i + ".txt")
    x = data2[:, 0]
    v = data2[:, 1]

    u = Exact(x)
    Error_log = np.log10(np.abs((u[1:-1] - v[1:-1]) / u[1:-1]))
    Error_log_abs = np.log10(np.abs((u[1:-1] - v[1:-1])))

    ax[0].plot(x[1:-1], Error_log, "o", label=f"n_Steps={i}")

    ax[1].plot(x[1:-1], Error_log_abs, "o", label=f"n_Steps={i}")
ax[0].set_xlabel("x")

ax[1].set_xlabel("x")
ax[0].set_ylabel(r"$log_{10}(\epsilon_i)$")
ax[1].set_ylabel(r"$log_{10}(\Delta_i)$")
plt.legend()
plt.savefig("error.pdf")
plt.show()
