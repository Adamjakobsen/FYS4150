import numpy as np
import matplotlib.pyplot as plt
import os
print("Working dir:", os.getcwd())

# inspired by code for reading file https://stackoverflow.com/questions/38532298/how-can-you-plot-data-from-a-txt-file-using-matplotlib

data1 = np.loadtxt("NOperationsP5Dense.txt")
data2 = np.loadtxt("NOperationsP5Tridiag.txt")

x = np.zeros(len(data1))
y1 = np.zeros(len(data1))
y2 = np.zeros(len(data2))

for i in range(len(data1)):
    x[i] = data1[i][0]
    y1[i] = data1[i][1]
    y2[i] = data2[i][1]


plt.plot(x, y1)
plt.plot(x, y2)

plt.xlabel("x")
plt.ylabel("Number of operations")
plt.show()
