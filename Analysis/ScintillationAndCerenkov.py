import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd

CURR_DIR = "C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Analysis"
datapaths=next(os.walk(CURR_DIR))[1][2:]
datapoints=[]
number_of_events=[]
for paths in datapaths:
    number_of_events.append(paths.split("_")[1])
    path=os.path.join(CURR_DIR,paths)
    optic_fiber_data=glob.glob(os.path.join(path,"Optic_fiber*.csv"))
    scintil_count=0
    cerenkov_count=0
    for f in optic_fiber_data:
        df=pd.read_csv(f,names=["column"])

        for i,item in enumerate(df["column"]):
            if item.split("\t")[1] == "Scintillation in core":
                scintil_count+=1
            elif item.split("\t")[1] == "Cerenkov in core":
                cerenkov_count+=1
            else:
                print("Unkwown result")
                pass
    print("incoming optical photons scintillation/cerenkov ratio:",scintil_count/cerenkov_count)

    optic_data=glob.glob(os.path.join(path,"Optic_*.csv"))[0:12]
    cerenkov_count_created=0
    for f in optic_data:
        df=pd.read_csv(f,names=["column"])
        try:
            cerenk=str(df["column"][1]).split("/t")[2]
            cerenk_num=int(cerenk.split(":")[1])
            cerenkov_count_created+=cerenk_num
        except:
            pass
    total_cerenkov=cerenkov_count_created+cerenkov_count
    print("total scintillation/cerenkov ratio",scintil_count/total_cerenkov)
    datapoints.append(scintil_count/total_cerenkov)


plt.figure("Optic Events")
plt.title("Scintillation and Cerenkov Ratio at 90°")
plt.plot(number_of_events[::-1],datapoints[::-1],marker='o',linestyle='dashed',color='blue',label='scint/cerenkov ratio')
plt.xlabel("Numero di eventi")
plt.ylabel("rapporto scintillazione/cerenkov")
plt.grid()
plt.show()
