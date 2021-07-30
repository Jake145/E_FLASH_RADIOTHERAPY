"""This script plots the various MC test performed"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import glob
import os
import re
from uncertainties import ufloat,unumpy
from uncertainties.umath import *

def PDD_plotter_out(dosepath,d,j=0):
    """This function returns the dose, the dose squared and the depth arrays of the water/PMMA
       simulated PDD

    :type dosepath: str
    :param dosepath: path of the .csv file with the simulation data

    :type d: float
    :param d: lenght in depth direction of scoring mesh

    :type j: int
    :param j: bin index of beam axis

    :returns: returns array of dose, dose squared and depth corresponding to each voxel
    :rtype: array
    """

    df = pd.read_csv(dosepath,names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])

    bins=len(np.array(df['dose [Gy]'][np.logical_and(df['Z']==0,df['Y']==0)]))

    lenght= d/bins

    distance=np.array([lenght*x for x in range(1,bins+1)])

    D=np.array(df['dose [Gy]'][np.logical_and(df['Z']==j,df['Y']==j)])

    ds=np.array(df['dosesq [Gy^2]'][np.logical_and(df['Z']==j,df['Y']==j)])

    return D,ds,distance

def mean_array_calculator(*args):
    """This function calculates the mean and std dev of the measurements
       to create the array with the means and uncertainties

    :type *args: array
    :param *args: arrays to stack and calculate the mean

    :returns: array of means with uncertainties
    :rtype: unumpy array
    """
    x = np.matrix(np.vstack(args))
    return unumpy.uarray(np.array(x.mean(0))[0],np.array(x.std(0))[0])

if __name__ == "__main__":

    path_9="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/9MevEF_PMMA_10mil_8cm_200bin.csv"
    dose_9_,sq_9,distance_9=PDD_plotter_out(path_9,80)
    path_7="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/7MevEF_custom_200bin_8cm.csv"
    dose_7,sq_7,distance_7=PDD_plotter_out(path_7,80)#seed 42
    path_106="C:/Users/pensa/Desktop/Send/phantomstats/dose_106.csv" #seed 106
    path_145="C:/Users/pensa/Desktop/Send/phantomstats/dose_145.csv" #seed 145

    dose_106,sq_106,distance_106=PDD_plotter_out(path_106,80)
    dose_145,sq_145,distance_145=PDD_plotter_out(path_145,80)

    assert(len(distance_9)==len(distance_106))
    assert(len(distance_9)==len(distance_145))

    dose_9=mean_array_calculator(dose_9_,dose_106,dose_145)

    dose_9_means=unumpy.nominal_values(dose_9)
    dose_9_std=unumpy.std_devs(dose_9)

    path_9_pen="C:/Users/pensa/Desktop/Send/physicslistsim/water_penelope_10mil.csv"
    dose_9_pen,sq_9_pen,distance_9_pen=PDD_plotter_out(path_9_pen,80)

    path_9_1mm="C:/Users/pensa/Desktop/Send/physicslists/dose_pen_step_1mm.csv"
    dose_9_1mm,sq_9_pen,distance_9_pen=PDD_plotter_out(path_9_1mm,80)

    path_9_liv="C:/Users/pensa/Desktop/Send/physicslistsim/livermore_10_mil.csv"
    dose_9_liv,sq_9_liv,distance_9_liv=PDD_plotter_out(path_9_liv,80)

    path_9_std="C:/Users/pensa/Desktop/Send/physicslistsim/std_opt_10mil_pdd.csv"
    dose_9_std,sq_9_std,distance_9_std=PDD_plotter_out(path_9_std,80)

    validation_distance_9,validation_dose_9=np.loadtxt("C:/Users/pensa/Desktop/Send/EF_9Mev_Water.txt",unpack=True)


    plt.figure("test")

    plt.plot(distance_9_pen,100*dose_9_pen/np.max(dose_9_pen),linestyle='dashed',color='red',label='Simulated Water 9 MeV Penelope')

    #plt.plot(distance_9_pen,100*dose_9_1mm/np.max(dose_9_1mm),linestyle='dashed',color='blue',label=' penelope maxstep 1 mm ')

    plt.plot(distance_9_pen,100*dose_9_liv/np.max(dose_9_liv),linestyle='dashed',color='blue',label='Simulated Water 9 MeV Livermore ')

    plt.plot(distance_9_pen,100*dose_9_std/np.max(dose_9_std),linestyle='dashed',color='orange',label='Simulated Water 9 MeV Std opt 4 ')

   # plt.plot(validation_distance_9,validation_dose_9,linestyle='dashed',color='green',label='Validation 9 MeV')



    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution Physics Lists")



    #plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()


    path_9="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/9MevEF_PMMA_10mil_8cm_200bin.csv"
    dose_9_,sq_9,distance_9=PDD_plotter_out(path_9,80)
    path_7="C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/7MevEF_custom_200bin_8cm.csv"
    dose_7,sq_7,distance_7=PDD_plotter_out(path_7,80)#seed 42
    path_106="C:/Users/pensa/Desktop/Send/phantomstats/dose_106.csv" #seed 106
    path_145="C:/Users/pensa/Desktop/Send/phantomstats/dose_145.csv" #seed 145

    dose_106,sq_106,distance_106=PDD_plotter_out(path_106,80)
    dose_145,sq_145,distance_145=PDD_plotter_out(path_145,80)

    dose_9=mean_array_calculator(dose_9_,dose_106,dose_145)

    dose_9_means=unumpy.nominal_values(dose_9)
    dose_9_std=unumpy.std_devs(dose_9)

    plt.figure("seed")

    plt.plot(distance_9,100*dose_9_/np.max(dose_9_),linestyle='dashed',color='red',label=' Seed = 42')


    plt.plot(distance_106,100*dose_106/np.max(dose_106),linestyle='dashed',color='blue',label=' Seed = 106')

    plt.plot(distance_145,100*dose_145/np.max(dose_145),linestyle='dashed',color='green',label='Seed = 145 ')



    textstr = (r' $\frac{\Delta dose}{dose}|_{buildup}=%.2f $ ' % (100*dose_9_std[np.where(dose_9_means==np.max(dose_9_means))]/np.max(dose_9_means))+"%")

    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

    plt.text( 45,70, textstr,  fontsize=14,
            verticalalignment='top', bbox=props)

    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution 9 MeV PMMA at various initial seeds")



    plt.legend()

    plt.grid()
    plt.show()


    distance_val,validation_dose=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/novac11PDD.txt',unpack=True)
    distancer50,doser50=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/r50.txt',unpack=True)
    distancer80,doser80=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/r80.txt',unpack=True)
    distancer100,doser100=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/r100.txt',unpack=True)
    path_novac11_1="C:/Users/pensa/Desktop/Send/simulations/novac_11_11mev_PDD.csv" #seed 106
    path_novac11_2="C:/Users/pensa/Desktop/Send/simulations/novac11_pdd_seed42_11mev_200bin.csv" #seed 42

    dose_novac11_1,sq_novac11_1,distance_novac11_1=PDD_plotter_out(path_novac11_1,80)
    dose_novac11_2,sq_novac11_2,distance_novac11_2=PDD_plotter_out(path_novac11_2,80)

    dose_novac11=mean_array_calculator(dose_novac11_1,dose_novac11_2)

    dose_novac11_means=unumpy.nominal_values(dose_novac11)
    dose_novac11_std=unumpy.std_devs(dose_novac11)




###
    path_1mm="C:/Users/pensa/Desktop/Send/prodcuts/water9Mevcut1mm.csv"
    path_1cm="C:/Users/pensa/Desktop/Send/prodcuts/1cmcuts.csv"
    dose_1mm,_,distance_=PDD_plotter_out(path_1mm,80)
    dose_1cm,_,_=PDD_plotter_out(path_1cm,80)


    plt.figure("cuts")


    plt.plot(distance_,100*dose_1mm/np.max(dose_1mm),linestyle='dashed',color='blue',label='cut = 1 mm')




    plt.plot(distance_,100*dose_1cm/np.max(dose_1cm),linestyle='dashed',color='red',label='cut = 1 cm')
    plt.plot(distance_,100*dose_9_pen/np.max(dose_9_pen),linestyle='solid',color='green',label='cut = 0.1 mm')




    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution Production Cuts 9 MeV water")



    plt.legend()

    plt.grid()
    plt.show()

    distance=np.array([7,10,28,31])

    ej212_energy_sim_3=np.array([10.0883,10.331,6.69084,4.84278])#good cuts
    ej212_01mm=np.array([10.1261,10.3468,6.67481,4.90235])
    ej212_1mm=np.array([10.163,10.3414,6.80742,4.8673])
    ej212_1cm=np.array([10.753,10.8904,6.92844,4.95998])

###

    plt.figure("cuts_scint")


    plt.plot(distance,100*ej212_energy_sim_3/np.max(ej212_energy_sim_3),linestyle='solid',color='red',marker='.',label='selected cuts')




    plt.plot(distance,100*ej212_01mm/np.max(ej212_01mm),linestyle='dashed',color='green',label='0.1 mm cut')
    plt.plot(distance,100*ej212_1mm/np.max(ej212_1mm),linestyle='dashed',color='blue',label='1 mm cut')
    plt.plot(distance,100*ej212_1cm/np.max(ej212_1cm),linestyle='dashed',color='orange',label='1 cm cut')

    plt.vlines([distance[0]], 0, (100*ej212_01mm/np.max(ej212_01mm))[0], linestyles='dashed', colors='magenta',alpha=0.5)
    plt.vlines([distance[1]], 0, (100*ej212_01mm/np.max(ej212_01mm))[1], linestyles='dashed', colors='magenta',alpha=0.5)
    plt.vlines([distance[2]], 0, (100*ej212_01mm/np.max(ej212_01mm))[2], linestyles='dashed', colors='magenta',alpha=0.5)
    plt.vlines([distance[3]], 0, (100*ej212_01mm/np.max(ej212_01mm))[3], linestyles='dashed', colors='magenta',alpha=0.5)






    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution EJ212 Production Cuts 9 MeV water")



    plt.legend()

    plt.grid()
    plt.show()























