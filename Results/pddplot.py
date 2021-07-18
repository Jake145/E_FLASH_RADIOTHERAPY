import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import glob
import os
import re

def PDD_plotter_out(dosepath,j,d): #j is the index of the column, d is the x lenght of the mesh in mm
    df = pd.read_csv(dosepath,names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])
    distance=np.linspace(0,d,len(np.array(df['dose [Gy]'][np.logical_and(df['Z']==0,df['Y']==0)])))


    D=np.array(df['dose [Gy]'][np.logical_and(df['Z']==j,df['Y']==j)])

    ds=np.array(df['dosesq [Gy^2]'][np.logical_and(df['Z']==j,df['Y']==j)])

    return D,ds,distance

if __name__ == "__main__":
    ##7 MeV Water
    simulated_energy_7=np.array([3.7348,3.8707,3.71136,3.43725,3.04878,2.38526,1.44811,0.865002,0.584597,0.28643,0.0731164,0.00361251])#GeV
    simulated_distance_7=np.array([2,4,10,12,15,18,22,25,27,30,35,40])#mm
    Run2_simulated_energy=np.array([3.82854,3.4693,1.82775,0.868089])##GeV
    Run2_simulated_events=np.array([2490,2801,2300,1632])
    Run2_distance=np.array([4,12.12,20.2,25])#mm
    validation_dose=np.array([50,90,100])
    validation_distance=np.array([27,18,12])
    path="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/7mevNOVAC_10000000_4cm_100bin_10000000.csv"
    dose_,sq,distance=PDD_plotter_out(path,0,80)

    print(f"R100 in water: {np.round(distance[15],2)} mm ")
    print(f"R90 in water: {np.round(distance[22],2)} mm ")
    print(f"R50 in water: {np.round(distance[31],2)} mm ")

    df = pd.read_csv(path,names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])

    indexes=[5,15,25,31]
    mass_lyso=7.4*1*0.2*0.2
    mass_lyso=mass_lyso *0.001
    coefficients=[]
    for i,item in enumerate(indexes):
        dose=df['dose [Gy]'][item]
        evts=1e7
        dose_exp=(Run2_simulated_energy[i] * 1.6e-10)/mass_lyso
        dose_norm=dose/evts
        dose_exp_norm=dose_exp/2e6
        coefficients.append(dose_norm/dose_exp_norm)
    coefficients=np.array(coefficients)
    print("Obtained coefficients: ",coefficients)
    ## Calculating spectra mean energy
    mainpath="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY"
    secondary_path="Coefficients"
    #outpath="C:/Users/pensa/Desktop/AnalysisThesis/Images_of_analysis"
    paths=glob.glob(os.path.join(mainpath,secondary_path))
    pattern="_[0-9]*"
    AvgEn=[]
    X=[]
    for p in paths:
        #print(p)
        datapaths=next(os.walk(p))[1]
        for datasetpath in datapaths:
            filenames=glob.glob(os.path.join(p,datasetpath,"Kinetic_*"))

            En=[]
            energies=[]
            for f in filenames:
                #print(f)
                df = pd.read_csv(f,names=['column'])
                for i,item in enumerate(df['column']):
                    if str(item).split("\t")[0] != 'Incoming Energy':
                        pass
                    else:
                    #print(item)
                        energies.append(float(str(item).split("\t")[2]))


            En=np.array(energies)
            depth=datasetpath
            X.append(depth)
            try:
                assert(len(En)!=0)
                AvgEn.append(En.mean())
                plt.figure(f"profondità: {depth}")
                plt.title(f"Spettro Energetico primari a {depth}")
                bin_heights, bin_borders, _=plt.hist(En,10,facecolor='g',ec='black',alpha=0.5,label='histogram data',density=True)
                #plt.title("Spettro energetico primari incidenti allo scintillatore")
                plt.ylabel("Frequenza")
                plt.grid()
                plt.xlabel("Energie [Mev]")

            except:
                pass
    plt.show()



    ##Plotting Coefficients in function of distance
    plt.figure("coefficients")
    plt.plot(Run2_distance,coefficients,marker='o',linestyle='dashed',color='blue')
    plt.xlabel("profondità [mm]")
    plt.ylabel("Dw/Dl")
    plt.title("coefficienti acqua-equivalenza")
    plt.grid()
    plt.show()

##Plotting Coefficients in function of distance
    Y_=[6.48,3.25,1.945,1.46]
    plt.figure("coefficients_en")
    plt.plot(Y_,coefficients,marker='o',linestyle='dashed',color='blue')
    plt.xlabel("Energia primari incidenti [MeV]")
    plt.ylabel("Dw/Dl")
    plt.title("coefficienti acqua-equivalenza")
    plt.grid()
    plt.show()

##Plotting PDD 7MeV
    x_r=np.array([4,12.12,25])
    y_r=np.array([simulated_energy_7[1]*coefficients[0],simulated_energy_7[3]*coefficients[1],simulated_energy_7[7]*coefficients[3]])
    y_r=y_r/np.max(y_r)*100
    plt.figure(1)
    plt.plot(distance,100*dose_/np.max(dose_),marker='.',linestyle='dashed',color='blue',label='acqua')
    plt.plot(simulated_distance_7,100*simulated_energy_7/np.max(simulated_energy_7),marker='.',linestyle='dashed',color='green',label='Lyso')

    plt.scatter(validation_distance,validation_dose,marker='.',color='red',label='R100, R90, R50 sperimentale acqua')

    plt.scatter(x_r,y_r,marker='.',color='magenta',label='Water Equivalence corrected Points')

    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("7 MeV Water Dose Distribution")

    plt.vlines([x_r[0]], 0, y_r[0], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[1]], 0, y_r[1], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[2]], 0, y_r[2], linestyles='dashed', colors='magenta',alpha=0.2)

    plt.hlines([y_r[0]], 0, x_r[0], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[1]], 0, x_r[1], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[2]], 0, x_r[2], linestyles='dashed', colors='magenta',alpha=0.2)


    plt.vlines([12], 0, 100, linestyles='dashed', colors='red',alpha=0.2)
    plt.vlines([18], 0, 90, linestyles='dashed', colors='red',alpha=0.2)
    plt.vlines([27], 0, 50, linestyles='dashed', colors='red',alpha=0.2)

    plt.hlines([100], 0, 12, linestyles='dashed', colors='red',alpha=0.2)
    plt.hlines([90], 0, 18, linestyles='dashed', colors='red',alpha=0.2)
    plt.hlines([50], 0, 27, linestyles='dashed', colors='red',alpha=0.2)

    #plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()
    ##9Mev and 7Mev PMMA
    path_9="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/9MevEF_PMMA_10mil_8cm_200bin.csv"
    dose_9,sq_9,distance_9=PDD_plotter_out(path_9,0,80)
    path_7="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/7MevEF_custom_200bin_8cm.csv"
    dose_7,sq_7,distance_7=PDD_plotter_out(path_7,0,80)
    #simulated_energy_7=np.array([3.80317,3.74554,3.77554,3.409,2.21552,0.805606])
    simulated_energy_7=np.array([3.74886,3.68112,3.13143,3.01489,2.13634,0.95427])

    simulated_energy_9=np.array([3.99749,3.96411,3.96835,3.83761,3.5549,2.49916])#gev lyso
    simulated_energy_9_V2=np.array([3.89034,3.90884,4.13303,4.08977,3.62473,2.91254])

    exp_distance_9=np.array([2,5,10,12,17,22])#mm lyso
    validation_distance_9,validation_dose_9=np.loadtxt("C:/Users/pensa/Desktop/EF_Validationdata.txt",unpack=True)
#seed 42
    #ej212_energy_sim=np.array([9.3235,9.55532,9.58284,9.43716,9.43479,9.20395,8.0759,6.54801,3.88037,2.31994,0.581116,0.0558538])#gev
    ej212_energy_sim_2=np.array([10.1766,10.5808,10.8589,10.5277,10.6088,10.6537,9.50409,7.40045,5.59706,2.67905,0.90506])#gev (starts from 6mm
    ej212_energy_sim_3=np.array([10.0883,10.331,10.5256,10.5953,10.4832,10.1652,9.11255,6.69084,4.84278,2.12886,0.436584])#gev (starts from 6mm
#seed 145
    ej_212_energy_sim_145=np.array([10.0633,10.458,10.4576,10.508,10.574,10.1423,8.96335,6.74539,4.97139,2.09613,0.439371])
#seed 111555875875
    ej_212_energy_sim_lol=np.array([10.1578,10.3613,10.3675,10.5164,10.4649,10.4815,8.97829,6.60971,4.80155,2.12145,0.42088])
#seed 8675309
    ej_212_energy_sim_jenny=np.array([10.0794,10.5218,10.4084,10.4092,10.4342,10.424,8.94853,6.60155,4.76539,2.12755,0.412595])
#seed 106
    ej_212_energy_sim_106=np.array([10.0259,10.2793,10.4171,10.5105,10.5316,10.2792,9.04665,6.72886,4.8282,2.0578,0.420623])

    ej212_distamce_sim_2=np.array([6,13,27])
    ej212_distance=np.array([4,6,9,11,13,15,17,22,27,30,35,40])+1#mm
    ej212_measured_charge=np.array([np.mean([7.220,7.210,7.150]),np.mean([7.360,7.310,7.260]),np.mean([7.510,7.460,7.420]),np.mean([7.580,7.530,7.480]),np.mean([7.610,7.550,7.520]),np.mean([7.490,7.460,7.400]),np.mean([7.300,7.240,7.210]),np.mean([6.360,6.370,6.330]),np.mean([4.450,4.460,4.480]),np.mean([2.930,2.960,2.980]),np.mean([0.580,0.620,0.580]),np.mean([0.480,0.480,0.460])])+0.57
    ej212_measured_monitor_units=np.array([18.43,18.29,18.24,18.20,18.10,18.07,18.10,18.09,18.15,18.06,18.16,18.21])
    ej212_measured_pulse=4
    ej212_charge_over_pulse=ej212_measured_charge/ej212_measured_pulse
    ej212_charge_over_mu=ej212_measured_charge/ej212_measured_monitor_units
##9Mev_ej212
    indexes_ej212=np.array([15,23,28,33,38,41,55,67,75,87,100])+2
    mass_ej212=1*1*2*0.2
    mass_ej212=mass_ej212 *0.001
    df_ej212 = pd.read_csv(path_9,names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])

    coefficients_ej212=[]
    for i,item in enumerate(indexes_ej212):
        dose_ej212=df_ej212['dose [Gy]'][item]
        evts_ej212=1e7
        dose_exp_ej212=(ej212_energy_sim_3[i] * 1.6e-10)/mass_ej212
        dose_norm_ej212=dose_ej212/evts_ej212
        dose_exp_norm_ej212=dose_exp_ej212/2e6
        coefficients_ej212.append(dose_norm_ej212/dose_exp_norm_ej212)
    coefficients_ej212=np.array(coefficients_ej212)
    print("Obtained coefficients: ",coefficients_ej212)


    x_r=ej212_distance[1:]
    y_r=ej212_energy_sim_3*coefficients_ej212
    y_r=y_r/np.max(y_r)*100
    c_corr_ej=ej212_charge_over_mu[1:]*coefficients_ej212
    c_corr_ej=c_corr_ej/(1/0.94-0.25*ej212_charge_over_mu[1:])

    plt.figure("ej212")
    #plt.plot(distance_9,100*dose_9/np.max(dose_9),marker='.',linestyle='dashed',color='blue',label='pmma simulato 9 MeV ')
    #plt.plot(ej212_distance+5,100*ej212_energy_sim/np.max(ej212_energy_sim),marker='.',linestyle='dashed',color='green',label='ej212 simulato 9 MeV ')
    plt.plot(ej212_distance,100*ej212_charge_over_mu/np.max(ej212_charge_over_mu),marker='.',linestyle='dashed',color='black',label='ej212 misurato 9 MeV ')
    #plt.plot(ej212_distance[1:]-1,ej212_energy_sim_2*100/np.max(ej212_energy_sim_2),marker='.',linestyle='dashed',color='green',label='ej212 simulato 9 MeV')
    #plt.scatter(x_r,y_r,marker='.',color='magenta',label='pmma Equivalence corrected Points')

    plt.plot(validation_distance_9,validation_dose_9,marker='.',linestyle='dashed',color='red',label='pmma misurato')

    plt.plot(ej212_distance[1:],100*c_corr_ej/np.max(c_corr_ej),marker='.',linestyle='dashed',color='green',label='ej212 corretto saturazione e spettro ')


    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash")



    #plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()
## 9mev lyso
    indexes=[5,13,25,30,42,55]
    mass_lyso=7.25*1*0.2*0.2
    mass_lyso=mass_lyso *0.001
    df = pd.read_csv(path_9,names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])

    coefficients=[]
    for i,item in enumerate(indexes):
        dose=df['dose [Gy]'][item]
        evts=1e7
        dose_exp=(simulated_energy_9_V2[i] * 1.6e-10)/mass_lyso
        dose_norm=dose/evts
        dose_exp_norm=dose_exp/2e6
        coefficients.append(dose_norm/dose_exp_norm)
    coefficients=np.array(coefficients)
    print("Obtained coefficients: ",coefficients)


    x_r=exp_distance_9
    y_r=simulated_energy_9_V2*coefficients
    y_r=y_r/np.max(y_r)*100

    measured_charge=np.array([4.854,4.834,4.9,4.556,4.214,3.12,1.46])+2.4
    measured_distance=np.array([0,2,5,10,12,17,22])
    measured_monitor_units=np.array([426.849,427.842,429.340,426.420,425.680,426.660,423.140])
    measured_charge_per_pulse=np.array([20.6,20.7,20.4,22,23.735,32.1,68.5])
##Saturation Correction
    c_mu=measured_charge/measured_monitor_units
    c_corr=c_mu[1:]*coefficients
    #c_corr=(measured_charge_per_pulse[1:]/(1/0.56-0.4*measured_charge_per_pulse[1:]))*c_corr
    c_corr=c_corr/np.max(c_corr)
    dis_corr=measured_distance[1:]
##

    plt.figure(2)
    plt.plot(distance_9,100*dose_9/np.max(dose_9),marker='.',linestyle='dashed',color='blue',label='pmma simulato 9 MeV ')
    plt.plot(exp_distance_9,100*simulated_energy_9_V2/np.max(simulated_energy_9_V2),marker='.',linestyle='dashed',color='green',label='lyso simulato 9 MeV ')

    plt.scatter(x_r,y_r,marker='.',color='magenta',label='Water Equivalence corrected Points')

    plt.plot(validation_distance_9,validation_dose_9,marker='.',linestyle='dashed',color='red',label='pmma misurato')

    plt.vlines([x_r[0]], 0, y_r[0], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[1]], 0, y_r[1], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[2]], 0, y_r[2], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[3]], 0, y_r[3], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[4]], 0, y_r[4], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[5]], 0, y_r[5], linestyles='dashed', colors='magenta',alpha=0.2)

    plt.hlines([y_r[0]], 0, x_r[0], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[1]], 0, x_r[1], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[2]], 0, x_r[2], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[3]], 0, x_r[3], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[4]], 0, x_r[4], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[5]], 0, x_r[5], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash")



    #plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()
##Corrected_exp
    plt.figure(3)
    plt.plot(exp_distance_9,100*simulated_energy_9_V2/np.max(simulated_energy_9_V2),marker='.',linestyle='dashed',color='yellow',label='lyso simulato 9 MeV ')

    plt.plot(validation_distance_9,validation_dose_9,marker='.',linestyle='dashed',color='red',label='pmma misurato')
    plt.plot(measured_distance,100*c_mu/np.max(c_mu),marker='.',linestyle='dashed',color='blue',label='lyso misurato 9 MeV ')
    plt.plot(dis_corr,100*c_corr,marker='.',linestyle='dashed',color='green',label='lyso corretto 9 MeV ')
    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash")



    #plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()

## 7Mev
    indexes=[5,13,25,30,42,55]
    mass_lyso=7.25*1*0.2*0.2
    mass_lyso=mass_lyso *0.001
    df = pd.read_csv(path_7,names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])

    coefficients=[]
    for i,item in enumerate(indexes):
        dose=df['dose [Gy]'][item]
        evts=1e7
        dose_exp=(simulated_energy_7[i] * 1.6e-10)/mass_lyso
        dose_norm=dose/evts
        dose_exp_norm=dose_exp/2e6
        coefficients.append(dose_norm/dose_exp_norm)
    coefficients=np.array(coefficients)
    print("Obtained coefficients: ",coefficients)
##
    x_r=exp_distance_9
    y_r=simulated_energy_7*coefficients
    y_r=y_r/np.max(y_r)*100
    validation_distance_7,_,validation_dose_7=np.loadtxt("C:/Users/pensa/Desktop/EF_val_7mev.txt",unpack=True)

    plt.figure(5)
    plt.plot(distance_7,100*dose_7/np.max(dose_7),marker='.',linestyle='dashed',color='blue',label='pmma simulato 7 MeV ')
    plt.plot(exp_distance_9,100*simulated_energy_7/np.max(simulated_energy_7),marker='.',linestyle='dashed',color='green',label='lyso simulato 7 MeV ')

    plt.scatter(x_r,y_r,marker='.',color='magenta',label='Water Equivalence corrected Points')

    plt.plot(validation_distance_7,validation_dose_7,marker='.',linestyle='dashed',color='red',label='pmma misurato')

    plt.vlines([x_r[0]], 0, y_r[0], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[1]], 0, y_r[1], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[2]], 0, y_r[2], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[3]], 0, y_r[3], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[4]], 0, y_r[4], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.vlines([x_r[5]], 0, y_r[5], linestyles='dashed', colors='magenta',alpha=0.2)

    plt.hlines([y_r[0]], 0, x_r[0], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[1]], 0, x_r[1], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[2]], 0, x_r[2], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[3]], 0, x_r[3], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[4]], 0, x_r[4], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.hlines([y_r[5]], 0, x_r[5], linestyles='dashed', colors='magenta',alpha=0.2)
    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash")



    #plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()
##
    x_r=exp_distance_9
    y_r=simulated_energy_9*coefficients
    y_r=y_r/np.max(y_r)*100
    validation_distance_7,_,validation_dose_7=np.loadtxt("C:/Users/pensa/Desktop/EF_val_7mev.txt",unpack=True)

    measured_charge=np.array([4.854,4.834,4.9,4.556,4.214,3.12,1.46])
    measured_charge_per_pulse=np.array([20.6,20.7,20.4,22,23.735,32.1,68.5])
    measured_distance=np.array([0,2,5,10,12,17,22])
    measured_monitor_units=np.array([426.849,427.842,429.340,426.420,425.680,426.660,423.140])
    c_mu=measured_charge/measured_monitor_units
    c_corr=c_mu[1:]*coefficients
    #c_corr=(measured_charge_per_pulse[1:]/(1/0.56-0.4*measured_charge_per_pulse[1:]))
    c_corr=c_corr/np.max(c_corr)
    dis_corr=measured_distance[1:]
    plt.figure(3)
    plt.plot(exp_distance_9,100*simulated_energy_7/np.max(simulated_energy_7),marker='.',linestyle='dashed',color='yellow',label='lyso simulato 7 MeV ')
    plt.plot(distance_7,100*dose_7/np.max(dose_7),marker='.',linestyle='dashed',color='blue',label='pmma simulato 7 MeV ')

    plt.plot(validation_distance_7,validation_dose_7,marker='.',linestyle='dashed',color='red',label='pmma misurato')
    plt.plot(measured_distance,100*c_mu/np.max(c_mu),marker='.',linestyle='dashed',color='green',label='lyso misurato 7 MeV ')
    plt.scatter(dis_corr,100*c_corr,marker='.',color='magenta',label='lyso corretto 7 MeV ')
    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash")



    #plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()