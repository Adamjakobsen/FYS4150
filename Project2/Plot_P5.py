import numpy as np
import matplotlib.pyplot as plt
import os


# inspired by code for reading file https://stackoverflow.com/questions/38532298/how-can-you-plot-data-from-a-txt-file-using-matplotlib

data1 = np.loadtxt('NOperationsP5Dense.txt ', float)
data2 = np.loadtxt('NOperationsP5Tridiag.txt ',float)

x = np.zeros(len(data1))
y1 = np.zeros(len(data1))
y2 = np.zeros(len(data2))

for i in range(len(data1)):
    x[i] = data1[i][0]
    y1[i] = data1[i][1]
    y2[i] = data2[i][1]


plt.plot(x, y1,label='Dense')
plt.plot(x, y2,label='Tridiagonal')

plt.xlabel("x")
plt.ylabel("Number of operations")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig('P5_fig.pdf')
plt.show()
