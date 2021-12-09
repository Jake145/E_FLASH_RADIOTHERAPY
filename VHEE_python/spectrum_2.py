import numpy as np
import os
import matplotlib.pyplot as plt
import csv
path_="C:/Users/pensa/Downloads"
filename="130MeV_Spectrum_xJake.txt"

path=os.path.join(path_,filename)

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:f}'.format})

energies,value=np.loadtxt(path,unpack=True)

max_energy=160
try:
    scaling_factor=max_energy - np.max(energies)
except:
    scaling_factor = 0
    print("scaling factor error")
energies=energies + scaling_factor
value=value/(np.sum(value))

command=np.array(["/gps/hist/point" for i in range(0,len(energies))])

plt.figure(1)
plt.bar(energies,value)
plt.grid()
plt.xlabel("Energy [MeV]")
plt.ylabel("Frequency normalized sum to 1")
plt.title("Spectrum 130 MeV")
plt.show()
path_2="C:/Users/pensa/Desktop"

ab = np.zeros(command.size, dtype=[('var1', 'U32'), ('var2', np.float32), ('var3', np.float32)])
ab['var1'] = command
ab['var2'] = energies
ab['var3'] = value

with open("C:/Users/pensa/Desktop/newfilePath.csv", "w+",newline="") as f:
    writer = csv.writer(f,delimiter="\t")
    for row in zip(ab['var1'],ab['var2'],ab['var3']):
        writer.writerow(row)
f.close()

a=np.column_stack((energies, value))
en=np.array([a[i][0] for i in range(0,len(a))])
val=np.array([a[i][1] for i in range(0,len(a))])
plt.figure(2)
plt.bar(en,val)
plt.grid()
plt.xlabel("Energy [MeV]")
plt.ylabel("Frequency normalized sum to 1")
plt.title("Spectrum 130 MeV")
plt.show()

