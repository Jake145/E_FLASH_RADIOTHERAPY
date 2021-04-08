import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob


CURR_DIR = "../PhotonsAndField/Field"
event_directories=next(os.walk(CURR_DIR))[1]
field_directories=glob.glob(os.path.join(CURR_DIR,"*"))
#cerenkov_directories=glob.glob(os.path.join(CURR_DIR,"Cerenkov_*"))
datapoints=[]
number_of_events=[]
Cerenkov=[]

def area(x):
    return np.pi*x**2
for index,paths in enumerate(field_directories):

    number_of_events.append(area(int(paths.split("/")[-1].split("cm")[0])))
    path=paths
    optic_data=glob.glob(os.path.join(path,"Optic_*.csv"))
    assert(len(optic_data)!=0)
    cerenkov=0
    scintillation=0
    for f in optic_data:
        df=pd.read_csv(f,names=["column"])
        assert(len(df["column"])!=0)
        for i,item in enumerate(df["column"]):
           

            
               
            cerenk=str(item).split("\t")[2]
            cerenk_num=int(cerenk.split(":")[1])
            cerenkov+=cerenk_num		    
            scint=str(item).split("\t")[1]
            scint_num=int(scint.split(":")[1])
            scintillation+=scint_num
    
    datapoints.append(scintillation/cerenkov)
    Cerenkov.append(cerenkov)     
	 
indexes=np.argsort(np.array(number_of_events))


plt.figure("Field",figsize=(25,10))
plt.title("Varying Field dimension at z= 2.5 cm")

plt.subplot(1,2,1)
plt.title("scintillation/cerenkov ratio")
plt.plot(np.array(number_of_events)[indexes],np.array(datapoints)[indexes],marker='o',linestyle='dashed',color='red',label='scint/cerenkov ratio')
plt.xlabel("Area Campo [cm^2]")
plt.ylabel("rapporto scintillazione/cerenkov")
#plt.legend()
plt.grid()

plt.subplot(1,2,2)
plt.title("Cerenkov created in scintillator+fiber")
plt.plot(np.array(number_of_events)[indexes],np.array(Cerenkov)[indexes],marker='o',linestyle='dashed',color='green',label='cerenkov number')
plt.xlabel("Area Campo [cm^2]")
plt.ylabel("numero di cerenkov")
#plt.legend()
plt.grid()

plt.show()

