"""This script plots the validations for Novac 11 in water"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy
from scipy.signal import savgol_filter
from tsmoothie.smoother import *

if __name__ == "__main__":


    distance_val,validation_dose=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/novac11PDD.txt',unpack=True)
    distancer50,doser50=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/r50.txt',unpack=True)
    distancer80,doser80=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/r80.txt',unpack=True)
    distancer100,doser100=np.loadtxt('C:/Users/pensa/Desktop/E_FLASH_RADIOTHERAPY/Flash_ex_novo/VALIDATION/r100.txt',unpack=True)

    df = pd.read_csv('C:/Users/pensa/Desktop/Send/validation_data/dose_novac11_profiles_10mil.csv',names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])




    distanceprofile=np.linspace(-80,80,len(np.array(df['dose [Gy]'][np.logical_and(df['X']==0,df['Y']==0)])))


    j=50
    bin=19
    plt.figure('profilerhor',figsize=(8,6))

    plt.subplot(3,1,1)
    plt.title('Novac 11 Profiles at 10 MeV at R100, R80 and R50')
    a=0
    b=5
    k=j


    # operate smoothing
    smoother = ConvolutionSmoother(window_len=5, window_type='ones')
    smoother.smooth(np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.array(df['dose [Gy]'][np.logical_and(np.logical_and(df['X']==bin,df['Y']==k),df['Z']==50)])*100)

    # generate intervals
    low, up = smoother.get_intervals('sigma_interval', n_sigma=5)

    plt.xlabel('distance [mm]')
    plt.ylabel('dose [$\%$]')
    plt.grid()
    plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.array(df['dose [Gy]'][np.logical_and(np.logical_and(df['X']==bin,df['Y']==k),df['Z']==50)])*100,marker='.',linestyle='dashed',color='blue',label='simulation')
    plt.plot(distancer100,doser100/np.max(doser100)*100,label='R100 profile',color='green')
    plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed', color='red',label='smoothed curve')
    plt.legend()

    bin=bin+16
    plt.subplot(3,1,2)
    smoother.smooth(np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.array(df['dose [Gy]'][np.logical_and(np.logical_and(df['X']==bin,df['Y']==k),df['Z']==50)])*100)

    # generate intervals
    low, up = smoother.get_intervals('sigma_interval', n_sigma=5)
    plt.xlabel('distance [mm]')
    plt.ylabel('dose [$\%$]')
    plt.grid()
    plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.array(df['dose [Gy]'][np.logical_and(np.logical_and(df['X']==bin,df['Y']==k),df['Z']==50)])*100,marker='.',linestyle='dashed',color='gray',label='simulation')
    plt.plot(distancer80,doser80/np.max(doser80)*100,label='R80 profile',color='orange')
    plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed', color='magenta',label='smoothed curve')
    plt.legend()

    bin=bin+11
    plt.subplot(3,1,3)
    smoother.smooth(np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.array(df['dose [Gy]'][np.logical_and(np.logical_and(df['X']==bin,df['Y']==k),df['Z']==50)])*100)

    # generate intervals
    low, up = smoother.get_intervals('sigma_interval', n_sigma=5)
    plt.xlabel('distance [mm]')
    plt.ylabel('dose [$\%$]')
    plt.grid()
    plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.array(df['dose [Gy]'][np.logical_and(np.logical_and(df['X']==bin,df['Y']==k),df['Z']==50)])*100,marker='.',linestyle='dashed',color='olive',label='simulation')
    plt.plot(distancer50,doser50/np.max(doser50)*100,label='R50 profile',color='cyan')
    plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed', color='orange',label='smoothed curve')
    plt.legend()
    plt.show()


    plt.figure("novac11")

    plt.plot(distance_novac11_1,100*dose_novac11_means/np.max(dose_novac11_means),linestyle='dashed',color='blue',marker='.',label='Simulation')




    plt.plot(distance_val,100*validation_dose/np.max(validation_dose),linestyle='dashed',color='red',marker='.',label='Validation')




    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution Novac11 10 MeV Water Validation")



    plt.legend()

    plt.grid()
    plt.show()
