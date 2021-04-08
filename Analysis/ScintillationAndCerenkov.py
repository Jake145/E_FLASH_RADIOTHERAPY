import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd

CURR_DIR = "../PhotonsAndField/Photon"
event_directories=next(os.walk(CURR_DIR))[1]
scintillation_directories=glob.glob(os.path.join(CURR_DIR,"*"))
#cerenkov_directories=glob.glob(os.path.join(CURR_DIR,"Cerenkov_*"))
datapoints=[]
number_of_events=[]

#print(scintillation_directories)

for index,paths in enumerate(scintillation_directories):

    number_of_events.append(int(paths.split("/")[-1])/1e4)
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
    #print("Scintillation: ",scintillation,"Cerenkov: ",cerenkov)        
    #print(scintillation/cerenkov)
    datapoints.append(scintillation/cerenkov)     
	 
indexes=np.argsort(np.array(number_of_events))


plt.figure("Optic Events")
plt.title("Produced Scintillation and Cerenkov Ratio at 90° and z= 2.5 cm")
plt.plot(np.array(number_of_events)[indexes],np.array(datapoints)[indexes],marker='o',linestyle='dashed',color='blue',label='scint/cerenkov ratio')
plt.xlabel("Numero di eventi [10^4]")
plt.ylabel("rapporto scintillazione/cerenkov")
plt.grid()
plt.show()
