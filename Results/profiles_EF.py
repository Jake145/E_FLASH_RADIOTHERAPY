import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy
from scipy.signal import savgol_filter
from tsmoothie.smoother import *

#so far the best one is dose_newcuts with binning 40 40 40 and size 4 8 8 cm

#df = pd.read_csv('dose_newcuts_5deg_newlist.csv',names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])
df = pd.read_csv('C:/Users/pensa/Desktop/Results/profile_9_mev_ef_more_stats_.csv',names=['X','Y','Z','dose [Gy]','dosesq [Gy^2]','entry'])

#distance_val,validation_dose=np.loadtxt('novac11PDD.txt',unpack=True)
distancer_hor,doser_hor=np.loadtxt('C:/Users/pensa/Desktop/Results/efProfiles9Mev.txt',unpack=True)
#distancer_vert,doser_vert=np.loadtxt('r80.txt',unpack=True)
#distancer100,doser100=np.loadtxt('r100.txt',unpack=True)

#distance=np.linspace(0,80,len(np.array(df['dose [Gy]'][np.logical_and(df['Z']==0,df['Y']==0)])))

distanceprofile=np.linspace(-70,70,len(np.array(df['dose [Gy]'][np.logical_and(df['X']==0,df['Y']==0)])))


j=75
bin=18
plt.figure('profiler50')
plt.title('profile r50 simulation vs validation data')

a=0
b=5

for i in range(a,b):
    k=j+i
    plt.subplot(2,3,i+1)
    plt.title('%2.f'%i)

    #yhat = savgol_filter(np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)]))*100, 39, 7) # window size 39, polynomial order 7

# operate smoothing
    #smoother = ConvolutionSmoother(window_len=5, window_type='ones')
    #smoother.smooth(np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)]))*100)

# generate intervals
    #low, up = smoother.get_intervals('sigma_interval', n_sigma=5)

    plt.xlabel('distance [mm]')
    plt.ylabel('dose [$/%$]')
    plt.grid()
    plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==bin,df['Y']==k)])/np.array(df['dose [Gy]'][np.logical_and(np.logical_and(df['X']==bin,df['Y']==k),df['Z']==70)])*100,marker='.',linestyle='dashed',color='blue',label='simulation')
    plt.plot(distancer_hor,doser_hor/np.max(doser_hor)*100,label='validation',color='green')
    #plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed', color='red',label='smoothed curve')
    plt.legend()
plt.show()
##plt.savefig('r50.png')


plt.figure('profile50 to save')
plt.title('R50 profile')
plt.xlabel('distance [mm]')
plt.ylabel('dose [$/%$]')
plt.grid()
plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==17,df['Y']==j)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==17,df['Y']==j)]))*100,marker='o',linestyle='dashed',color='blue',label='simulation')
plt.plot(distancer_hor,doser_hor/np.max(doser_hor)*100,label='validation',color='green',linestyle='dashed')
#plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed', color='red',label='smoothed curve')
plt.legend()
#plt.savefig('r50.png')
plt.show()


plt.figure('profiler80')
plt.title('profile r80 simulation vs validation data')



for i in range(a,b):
    k=j+i
    plt.subplot(2,3,i+1)
    yhat = savgol_filter(np.array(df['dose [Gy]'][np.logical_and(df['X']==12,df['Y']==k)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==12,df['Y']==k)]))*100, 39, 7) # window size 39, polynomial order 7

# operate smoothing
    smoother = ConvolutionSmoother(window_len=5, window_type='ones')
    smoother.smooth(np.array(df['dose [Gy]'][np.logical_and(df['X']==12,df['Y']==k)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==12,df['Y']==k)]))*100)

# generate intervals
    low, up = smoother.get_intervals('sigma_interval', n_sigma=5)
    plt.title('%2.f'%i)

    plt.xlabel('distance [mm]')
    plt.ylabel('dose [$/%$]')
    plt.grid()
    plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==12,df['Y']==k)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==12,df['Y']==k)]))*100,marker='o',linestyle='dashed',color='blue',label='simulation')
    plt.plot(distancer80,doser80/np.max(doser80)*100,label='validation',color='green')
    plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed',           color='red',label='smoothed curve')
    plt.legend()
    ##plt.savefig('r80.png')
plt.show()

plt.figure('profile80 to save')
plt.title('R80 profile')

plt.xlabel('distance [mm]')
plt.ylabel('dose [$/%$]')
plt.grid()
plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==12,df['Y']==j)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==12,df['Y']==j)]))*100,marker='o',linestyle='dashed',color='blue',label='simulation')
plt.plot(distancer80,doser80/np.max(doser80)*100,label='validation',color='green',linestyle='dashed')
#plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed', color='red',label='smoothed curve')
plt.legend()
#plt.savefig('r80.png')
plt.show()

plt.figure('profiler100')
plt.title('profile r100 simulation vs validation data')



for i in range(a,b):
    k=j+i
    plt.subplot(2,3,i+1)
    yhat = savgol_filter(np.array(df['dose [Gy]'][np.logical_and(df['X']==6,df['Y']==k)])/          np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==6,df['Y']==k)]))*100, 39, 7) # window size 39, polynomial order 7

# operate smoothing
    smoother = ConvolutionSmoother(window_len=5, window_type='ones')
    smoother.smooth(np.array(df['dose [Gy]'][np.logical_and(df['X']==6,df['Y']==k)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==6,df['Y']==k)]))*100)

# generate intervals
    low, up = smoother.get_intervals('sigma_interval', n_sigma=5)
    plt.title('%2.f'%i)

    plt.xlabel('distance [mm]')
    plt.ylabel('dose [$/%$]')
    plt.grid()
    plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==6,df['Y']==k)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==6,df['Y']==k)]))*100,marker='o',linestyle='dashed',color='blue',label='simulation')
    plt.plot(distancer100,doser100/np.max(doser100)*100,label='validation',color='green')
    plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed', color='red',label='smoothed curve')

    plt.legend()
##plt.savefig('r100.png')
plt.show()

plt.figure('profile100 to save')
plt.title('R100 profile')

plt.xlabel('distance [mm]')
plt.ylabel('dose [$/%$]')
plt.grid()
plt.plot(distanceprofile,np.array(df['dose [Gy]'][np.logical_and(df['X']==6,df['Y']==j)])/np.max(np.array(df['dose [Gy]'][np.logical_and(df['X']==6,df['Y']==j)]))*100,marker='o',linestyle='dashed',color='blue',label='simulation')
plt.plot(distancer100,doser100/np.max(doser100)*100,label='validation',color='green',linestyle='dashed')
#plt.plot(distanceprofile,smoother.smooth_data[0], linestyle='dashed', color='red',label='smoothed curve')
plt.legend()
#plt.savefig('r100.png')
plt.show()
