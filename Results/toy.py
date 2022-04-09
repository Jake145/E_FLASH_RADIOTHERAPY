import numpy as np
import matplotlib.pyplot as plt

x,y,dy=np.loadtxt("ej212_measured_corrected.txt",unpack=True)

plt.figure(1)
plt.errorbar(x,y,yerr=dy)
plt.scatter(x,y)
plt.show()
