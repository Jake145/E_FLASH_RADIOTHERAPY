
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import glob

CURR_DIR = "./7MeV"

event_directories=next(os.walk(CURR_DIR))[1]

datapoints=[]
depth=[]
incident_electrons=[]
fluxes=[]
backscatter=[]
#Scintillation and cerenkov

for index,paths in enumerate(event_directories):
    paths=os.path.join(CURR_DIR,paths)
    print(paths)
    depth.append(int(paths.split("/")[-1].split("mm")[0])) #You gotta put the os.path.join thing
    path=paths
    optic_data=glob.glob(os.path.join(path,"Optic_*.csv"))
    assert(len(optic_data)!=0)
    cerenkov=0
    scintillation=0
    for f in optic_data:
        try:
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
            
        except:
            pass 
        if cerenkov!=0:
            datapoints.append(scintillation/cerenkov)
        else:
            datapoints.append(0)   
    electron_datas=glob.glob(os.path.join(path,"Kinetic_*.csv"))
    assert(len(electron_datas)!=0)
    number_e=0
    energies=[]
    for f in electron_datas:
        df=pd.read_csv(f,names=["column"])
        try:
            assert(len(df["column"])!=0)
            for i,item in enumerate(df["column"]):
                name=str(item).split("\t")[0]
                print(name)
                if name=="Incident Electron":
                    number_e+=1
                elif name=="Incoming Energy":
                    energy=float(str(item).split("\t")[2])
                    print(energy)
                    energies.append(energy)
                else:
                    print("What is wrong with this file man?")
        except:
            pass
    try:
        assert(len(energies)!=0)
        plt.figure("profondità: %.2f mm"%int(paths.split("/")[-1].split("mm")[0]))
        plt.title("Spettro Energetico primari a %.2f mm"%int(paths.split("/")[-1].split("mm")[0]))
        bin_heights, bin_borders, _=plt.hist(energies,round(np.sqrt(len(energies))),facecolor='g',ec='black',alpha=0.5,label='histogram data',density=True)
        #plt.title("Spettro energetico primari incidenti allo scintillatore")
        plt.ylabel("Frequenza")
        plt.grid()
        plt.xlabel("Energie [Mev]")
        figname=("Spettro_%.2f.png"%int(paths.split("/")[-1].split("mm")[0]))
        plt.savefig(os.path.join(CURR_DIR,figname))
    except:
        pass
    fluxfile="flux.csv"
    df_ = pd.read_csv(os.path.join(path,fluxfile),names=['X','Y','Z','flux [Gy]','fluxsq [Gy^2]','entry'])
    flux_number=0
    try:
        assert(len(df_["entry"])!=0)
        for i,item in enumerate(df_["entry"]):
            flux_number+=item
    except:
        pass
    fluxes.append(flux_number)
    backscatt=(flux_number -number_e)/2
    backscatter.append(number_e/backscatt) 
    
    
    
indexes=np.argsort(np.array(depth))


plt.figure("Optic Events")
plt.title("Produced Scintillation and Cerenkov Ratio")
plt.plot(np.array(depth)[indexes],np.array(datapoints)[indexes],marker='o',linestyle='dashed',color='blue',label='scint/cerenkov ratio')
plt.xlabel("profondità [mm]")
plt.ylabel("rapporto scintillazione/cerenkov")
plt.grid()
plt.show()

plt.figure("Backscatter")
plt.title("Backscatter electrons")
plt.plot(np.array(depth)[indexes],np.array(backscatter)[indexes],
marker='o',linestyle='dashed',color='blue',label='backscatters')
plt.xlabel("profondità [mm]")
plt.ylabel("rapporto incidenti/backscatter")
plt.grid()
plt.show()
##Backscatter



