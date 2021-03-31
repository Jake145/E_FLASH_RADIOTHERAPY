import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def PDD_plotter(energy_det,x_det,dosepath,j,seed_number=70,evt_number=1000000,seed_or_evt=True):
    df = pd.read_csv(dosepath,names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])
    distance=np.linspace(0,100,len(np.array(df['dose [Gy]'][np.logical_and(df['Z']==0,df['Y']==0)])))
    distance_val,validation_dose=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/novac11PDD.txt',unpack=True)


    D=np.array(df['dose [Gy]'][np.logical_and(df['Z']==j,df['Y']==j)])

    ds=np.array(df['dosesq [Gy^2]'][np.logical_and(df['Z']==j,df['Y']==j)])

    indexes=np.where(np.in1d(np.round(distance),x_det ))[0]
    dose=D[indexes]

    difference=np.abs(((np.array(dose)/np.max(np.array(dose)))-(energy_det[::-1]/np.max(energy_det)))*100)
    if seed_or_evt:
        plt.figure(f"Seed:{seed_number}")
        #plt.xscale('log')

        #plt.yscale('log')
        plt.title(f'PDD rivelatore LYSO con seed: {seed_number}')
    elif seed_or_evt==False:
        plt.figure(f"Events:{evt_number}")
        #plt.xscale('log')

        #plt.yscale('log')
        plt.title(f'PDD rivelatore LYSO con numero di eventi: {evt_number}')
    else:
        print("What?")
    plt.ylabel('Dose relativa [%]')
    plt.xlabel('Profondità [mm]')
    plt.plot(x_det,(energy_det[::-1]/np.max(energy_det))*100,marker='o',linestyle='dashed',color='blue',label='Con rivelatore')

    plt.plot(x_det,(np.array(dose)/np.max(np.array(dose)))*100,marker='o',linestyle='dashed',color='red',label='Senza rivelatore')


    plt.plot(x_det,difference,marker='^',linestyle='dashed',color='green',label='differenza (in valore assoluto)')

    plt.plot(distance_val,validation_dose,label='PDD sperimentale in acqua',color='green',linestyle='dashed',alpha=0.6)

    plt.grid()
    plt.legend()
    plt.show()
    return (energy_det[::-1]/np.max(energy_det))*100,(np.array(dose)/np.max(np.array(dose)))*100,difference
if __name__ == "__main__":

    x=np.array([i*4 for i in range(1,26)])
    energy=np.array([5.1257,6.57224,9.56792,10.203,11.8436,9.05542,10.3728,13.0155,15.6298,15.013,14.6585,13.8702,18.6531,25.4665,70.4515,158.052,329.911,471.795,557.863,662.301,703.98,705.362,720.122,751.33,736.001])
    water_dose_path='C:/Users/pensa/Desktop/dose_.csv'
    index_=1

    en_,dose_,diff_=PDD_seed_first=PDD_plotter(energy,x,water_dose_path,index_,548235486)

    energy_seed70=np.array([6.88713,7.47544,8.38924,9.78076,8.41548,9.27428,11.4792,13.8505,11.1641,12.3578,13.8565,14.0372,16.2199,30.3783,76.0533,177.326,343.558,439.645,590.097,685.512,733.006,751.9,795.296,730.024,783.761])
    water_dose_path_70="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/PDDS_SEED/Seed70/dose_seed70_.csv"
    index_70=0

    en_70,dose_70,diff_70=PDD_seed_70=PDD_plotter(energy_seed70,x,water_dose_path,index_70,70)

    energy_seed42 = np.array([9.283347,10.044,8.87046,7.50841,8.50056,8.51872,11.4371,11.7872,11.8205,11.1162,13.5585,14.8237,18.8701,29.3505,76.5672,167.397,293.762,430.718,598.858,675.931,759.087,805.082,850.455,805.04,785.721])
    water_dose_path_42="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/PDDS_SEED/Seed42/dose_seed42.csv"
    index_42=0

    en_42,dose_42,diff_42=PDD_plotter(energy_seed42,x,water_dose_path_42,index_42,42)

    energy_seed70_10000=np.array([0,84.7644,36.4487,133.318,228.81,86.9574,0,115.54,35.8759,53.9828,106.37,526.958,0,47.5125,1820.62,286.211,2739.01,4384.2,1522.76,4638.43,7181.18,7098.72,12793.4,10394.5,8250.61])#kev
    water_dose_path_10000="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/PDDS_EVT/10000/dose.csv"

    energy_seed70_100000=np.array([0.921058,0.407367,0.339114,0.722062,1.11582,0.898083,1.00956,0.891968,1.01412,1.11796,1.22934,1.4381,1.15956,1.29224,5.54584,17.3882,28.548,35.8041,62.6658,64.55,71.9706,85.016,87.1357,82.0007,80.2609])#Mev
    water_dose_path_100000="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/PDDS_EVT/100000/dose.csv"
    energy_seed70_500000=np.array([4.94434,3.60706,3.27444,5.374565,4.74935,6.23723,5.36259,6.20096,5.87947,5.42365,7.12847,5.19247,7.01484,14.4978,42.7781,100.664,154.022,209.323,280.975,355.513,354.334,390.271,386.177,357.522,375.922])
    water_dose_path_500000="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/PDDS_EVT/500000/dose.csv"

    en_70_10000,dose_70_10000,diff_70_10000=PDD_seed_70=PDD_plotter(energy_seed70_10000,x,water_dose_path_10000,index_70,70,10000,False)

    en_70_100000,dose_70_100000,diff_70_100000=PDD_seed_70=PDD_plotter(energy_seed70_100000,x,water_dose_path_100000,index_70,70,100000,False)

    en_70_500000,dose_70_500000,diff_70_500000=PDD_seed_70=PDD_plotter(energy_seed70_500000,x,water_dose_path_500000,index_70,70,500000,False)
    distance_val,validation_dose=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/novac11PDD.txt',unpack=True)


    plt.figure("PDD LYSO al variare del seed",figsize=(25,10))
    plt.title("PDD LYSO al variare del seed")
    plt.subplot(2,3,1)
    plt.title(f'PDD rivelatore LYSO con seed: 548235486')
    plt.ylabel('Dose relativa [%]')
    plt.xlabel('Profondità [mm]')
    plt.plot(x,en_,marker='o',linestyle='dashed',color='blue',label='Con rivelatore')

    plt.plot(x,dose_,marker='o',linestyle='dashed',color='red',label='Senza rivelatore')


    plt.plot(x,diff_,marker='^',linestyle='dashed',color='green',label='differenza (in valore assoluto)')

    plt.plot(distance_val,validation_dose,label='PDD sperimentale in acqua',color='green',linestyle='dashed',alpha=0.6)
    plt.grid()
    plt.legend()
    plt.subplot(2,3,2)
    plt.title(f'PDD rivelatore LYSO con seed: 70')
    plt.ylabel('Dose relativa [%]')
    plt.xlabel('Profondità [mm]')
    plt.plot(x,en_70,marker='o',linestyle='dashed',color='blue',label='Con rivelatore')

    plt.plot(x,dose_70,marker='o',linestyle='dashed',color='red',label='Senza rivelatore')


    plt.plot(x,diff_70,marker='^',linestyle='dashed',color='green',label='differenza (in valore assoluto)')

    plt.plot(distance_val,validation_dose,label='PDD sperimentale in acqua',color='green',linestyle='dashed',alpha=0.6)
    plt.grid()
    plt.legend()
    plt.subplot(2,3,3)
    plt.title(f'PDD rivelatore LYSO con seed: 42')
    plt.ylabel('Dose relativa [%]')
    plt.xlabel('Profondità [mm]')
    plt.plot(x,en_42,marker='o',linestyle='dashed',color='blue',label='Con rivelatore')

    plt.plot(x,dose_42,marker='o',linestyle='dashed',color='red',label='Senza rivelatore')


    plt.plot(x,diff_42,marker='^',linestyle='dashed',color='green',label='differenza (in valore assoluto)')

    plt.plot(distance_val,validation_dose,label='PDD sperimentale in acqua',color='green',linestyle='dashed',alpha=0.6)
    plt.grid()
    plt.legend()

    plt.savefig("C:/Users/pensa/Desktop/AnalysisThesis/Images_of_analysis/PDD_seeds.png")





    plt.figure("PDD LYSO al variare del numero di eventi",figsize=(25,10))
    plt.title("PDD LYSO al variare del numero di eventi")
    plt.subplot(2,3,1)
    plt.title('PDD rivelatore LYSO con numero di eventi: 10000')
    plt.ylabel('Dose relativa [%]')
    plt.xlabel('Profondità [mm]')
    plt.plot(x,en_70_10000,marker='o',linestyle='dashed',color='blue',label='Con rivelatore')

    plt.plot(x,dose_70_10000,marker='o',linestyle='dashed',color='red',label='Senza rivelatore')


    plt.plot(x,diff_70_10000,marker='^',linestyle='dashed',color='green',label='differenza (in valore assoluto)')

    plt.plot(distance_val,validation_dose,label='PDD sperimentale in acqua',color='green',linestyle='dashed',alpha=0.6)
    plt.grid()
    plt.legend()
    plt.subplot(2,3,2)
    plt.title('PDD rivelatore LYSO con numero di eventi: 100000')
    plt.ylabel('Dose relativa [%]')
    plt.xlabel('Profondità [mm]')
    plt.plot(x,en_70_100000,marker='o',linestyle='dashed',color='blue',label='Con rivelatore')

    plt.plot(x,dose_70_100000,marker='o',linestyle='dashed',color='red',label='Senza rivelatore')


    plt.plot(x,diff_70_100000,marker='^',linestyle='dashed',color='green',label='differenza (in valore assoluto)')

    plt.plot(distance_val,validation_dose,label='PDD sperimentale in acqua',color='green',linestyle='dashed',alpha=0.6)
    plt.grid()
    plt.legend()
    plt.subplot(2,3,3)
    plt.title('PDD rivelatore LYSO con numero di eventi: 500000')
    plt.ylabel('Dose relativa [%]')
    plt.xlabel('Profondità [mm]')
    plt.plot(x,en_70_500000,marker='o',linestyle='dashed',color='blue',label='Con rivelatore')

    plt.plot(x,dose_70_500000,marker='o',linestyle='dashed',color='red',label='Senza rivelatore')


    plt.plot(x,diff_70_500000,marker='^',linestyle='dashed',color='green',label='differenza (in valore assoluto)')
    plt.plot(distance_val,validation_dose,label='PDD sperimentale in acqua',color='green',linestyle='dashed',alpha=0.6)
    plt.grid()
    plt.legend()

    plt.savefig("C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY7Analysis/Images_of_analysis/PDD_evts.png")

