import numpy as np
import matplotlib.pyplot as plt

n = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]

avg_general_array = []
avg_special_array = []

for i in n:
    # Loading Data file
    general_data = np.loadtxt("Time_General_solution_" + str(i) + ".txt")
    special_data = np.loadtxt("Time_Special_solution_" + str(i) + ".txt")

    # average of all the runs in each file
    avg_general_array.append(np.mean(general_data))
    avg_special_array.append(np.mean(special_data))


plt.plot(n, avg_general_array, linestyle="--", marker="o", label=f"General Algorithm")
plt.plot(n, avg_special_array, linestyle="--", marker="o", label=f"Special Algorithm")

plt.xlabel("Number of steps")
plt.ylabel("Avg. times for algorithm executions (s)")

plt.legend()
plt.savefig("Problem10.pdf")
plt.show()
