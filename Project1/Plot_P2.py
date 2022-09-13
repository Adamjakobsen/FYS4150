import numpy as np
import matplotlib.pyplot as plt

# inspired by code for reading file https://stackoverflow.com/questions/38532298/how-can-you-plot-data-from-a-txt-file-using-matplotlib


data = np.loadtxt("data.txt")

x = np.zeros(len(data))
y = np.zeros(len(data))

for i in range(len(data)):
    x[i] = data[i][0]
    y[i] = data[i][1]


print(x)
print("\n")
print(y)

plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("u(x)")
plt.show()
