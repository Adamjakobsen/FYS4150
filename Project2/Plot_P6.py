import numpy as np
import matplotlib.pyplot as plt

#Get data from P6eigenboys.txt, first column is v1,second column is v2, third column is v3

data = np.loadtxt("P6eigenboys.txt",float)

v1 = data[:,0]
v2 = data[:,1]
v3 = data[:,2]
#add 0 to start of v1, v2, v3
v1 = np.insert(v1,0,0)
v2 = np.insert(v2,0,0)
v3 = np.insert(v3,0,0)
#add 0 to end of v1, v2, v3
v1 = np.append(v1,0)
v2 = np.append(v2,0)
v3 = np.append(v3,0)

x = np.linspace(0,1,len(v1))
#Plot the data
plt.plot(x,v1,label="v1")
plt.plot(x,v2,label="v2")
plt.plot(x,v3,label="v3")
plt.xlabel(r"$\^x$")
plt.ylabel("v")
plt.title("Solutions to the eigenvalue problem")
plt.legend()
plt.savefig('P6_fig'+f'_{len(x)}'+ '.pdf')
