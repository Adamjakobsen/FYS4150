import numpy as np
import matplotlib.pyplot as plt


def std_dev(
    data, avg
):  # Will give us a standard deviation for each set of measurements
    diffs = [(x - avg) ** 2 for x in data]
    s2 = (1 / (len(data) - 1)) * np.sum(diffs)  # sample variance
    return np.sqrt(s2)


n = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]  # steps

avg_general_array = []
avg_special_array = []
dev_general_array = []
dev_special_array = []

for i in n:
    # Loading Data file
    general_data = np.loadtxt("Time_General_solution_" + str(i) + ".txt")
    special_data = np.loadtxt("Time_Special_solution_" + str(i) + ".txt")

    avg_general = np.mean(general_data)
    avg_special = np.mean(special_data)

    dev_general = std_dev(general_data, avg_general)
    dev_special = std_dev(special_data, avg_special)

    # average of all the runs in each file
    avg_general_array.append(avg_general)
    avg_special_array.append(avg_special)

    dev_general_array.append(dev_general)
    dev_special_array.append(dev_special)


plt.plot(
    n,
    avg_general_array,
    linestyle="--",
    marker="o",
    label=f"General Algorithm",
    markersize=5,
)
plt.errorbar(n, avg_general_array, dev_general_array, linestyle="None", color="black")

plt.plot(
    n,
    avg_special_array,
    linestyle="--",
    marker="o",
    label=f"Special Algorithm",
    markersize=5,
)
plt.errorbar(
    n,
    avg_special_array,
    dev_special_array,
    linestyle="None",
    color="black",
    label=f"Error Bar",
)

plt.xlabel("Number of steps")
plt.ylabel("Avg. times for algorithm executions (s)")

plt.legend()
plt.savefig("Problem10.pdf")
plt.show()
