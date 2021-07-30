"""This script fits the experimental data and creates the saturation plots"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import glob
import os
import re
import scipy.stats as stats
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from labellines import labelLine, labelLines
from uncertainties import ufloat,unumpy
from uncertainties.umath import *

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

def model(x,alpha,k):
    """This function creates the model that relates the DPP to the CPP

    :type x: array
    :param x: CPP

    :type alpha: float
    :param alpha: float

    :type k: float
    :param k: float

    :returns: values calculated for the given alpha and k
    :rtype: unumpy array
    """
    return x/((1/k-alpha*x))

def model_inv(x,alpha,k):
    """This function creates the model that relates the CPP to the DPP

    :type x: array
    :param x: DPP

    :type alpha: float
    :param alpha: float

    :type k: float
    :param k: float

    :returns: values calculated for the given alpha and k
    :rtype: unumpy array
    """
    return k*x/((1+alpha*x))

def get_index(x,y):
    """This function finds the indexes of the first array that is the same as the
       elements of the second one

    :type x: array
    :param x: first array

    :type y: array
    :param y: second array

    :returns: indexes
    :rtype: list
    """
    indexes=[]
    for i,item in enumerate(x):
        for j,jtem in enumerate(y):
            try:
                assert(item==jtem)
                indexes.append(j)
            except:
                pass


    return indexes



def plotter(cpp_mean,cpp_sigma,p0_=[0.5,2],detector_name="EJ212",cord_1=10,cord_2=2):

    """This function does the curve fitting of CPP vs DPP and plots the resulting fit

    :type cpp_mean: array
    :param cpp_mean: mean values of the CPP

    :type cpp_sigma: array
    :param cpp_sigma: uncertainty values of the CPP

    :type p0_: list
    :param p0_: initial values of alpha and K

    :type detector_name: str
    :param detector_name: name of the detector

    :type cord_1: float
    :param cord_1: x coordinate for the results text box

    :type cord_2: float
    :param cord_2: y coordinate for the results text box


    :returns: float
    :rtype: returns alpha, k and their uncertainties
    """

    tmp,sigma_tmp= cpp_mean,cpp_sigma

    popt, pcov = curve_fit(model_inv,ddp_gaf_pl,tmp,p0_,sigma=sigma_tmp, absolute_sigma=False)

    sigma_alpha,sigma_k=np.diagonal(pcov)

    correlation_alpha_k,correlation_k_alpha=np.diagonal(np.fliplr(pcov))

    alpha_,k_=popt

    stats_,pvalue=stats.chisquare(f_obs=ddp_gaf_pl, f_exp=model_inv(tmp,alpha_,k_))

    r2=r2_score(tmp, model_inv(ddp_gaf_pl,alpha_,k_))

    plt.figure(f"pulse lenght: {pulse_lenght}")

    plt.subplot(2,1,1)

    plt.xlabel("Dose per pulse [Gy/p]")

    plt.ylabel("Charge per pulse [nC/p]")


    plt.xlim(0,14)

    plt.title(f"Saturation for {detector_name} at pulse lenght {pulse_lenght} $\mu s$")

    plt.scatter(ddp_gaf_pl,tmp,marker='.',color='red',label="experimental data")

    plt.errorbar(ddp_gaf_pl,tmp,yerr=sigma_tmp, xerr=gaf_uncertainty[gaf_points], ls='')

    plt.plot(np.linspace(0,20,100),model_inv(np.linspace(0,20,100),alpha_,k_),color='magenta',linestyle='--')


    textstr = '\n'.join((

        r'$y = k\cdot\frac{x}{1 + \alpha x}$',

        r'$\alpha=%.4f \pm %.4f $' % (alpha_,sigma_alpha ),

        r'$k=%.4f \pm %.4f $' % (k_,sigma_k ),


        r'$R^2$ score = %.3f ' % (r2)))

    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

    plt.text( cord_1,cord_2, textstr,  fontsize=14,
            verticalalignment='top', bbox=props)

    plt.grid()

    plt.subplot(2,1,2)


    plt.xlabel("Dose per pulse [Gy/p]")
    plt.ylabel(r"$\frac{q-f(x)}{\sigma}$")
    res = (tmp - model_inv(ddp_gaf_pl,alpha_,k_))/sigma_tmp
    plt.scatter(ddp_gaf_pl,res,marker='o', color='green',label='residuals')
    plt.legend()
    plt.grid()


    plt.figure(f"pulse lenght: {pulse_lenght} residuals")
    plt.title(f"Residuals for {detector_name} at pulse lenght {pulse_lenght} $\mu s$")
    plt.xlabel("Dose per pulse [Gy/p]")
    plt.ylabel(r"$\frac{q-f(x)}{\sigma}$")
    res = (tmp - model_inv(ddp_gaf_pl,alpha_,k_))/sigma_tmp
    plt.scatter(ddp_gaf_pl,res,marker='o', color='green',label='residuals')
    plt.legend()
    plt.grid()


    plt.show()


    return alpha_,k_,sigma_alpha,sigma_k



if __name__ == "__main__":

## EJ212

    ddp_gaf=np.array([12.5,5.8,4.5,3.5,3.1,2.9,2.7,2,1.2,0.7,0.5])

    ssd_gaf=np.array([73,90,106,105,114,122,123,134,162,190,217])

    ddp_gaf=ddp_gaf[np.argsort(ssd_gaf)]

    ssd_gaf=ssd_gaf[np.argsort(ssd_gaf)]

    gaf_uncertainty=0.03*ddp_gaf

    ssd=np.array([73,90,105,106,114,122,134,162,190,217])

    gaf_points=get_index(ssd,ssd_gaf)

    ssd_gaf_pl=ssd_gaf[gaf_points]

    ddp_gaf_pl=ddp_gaf[gaf_points]

    detector_name="EJ212"

    pulses=np.array([1,1,1,4,4,4,4,4,4,4])




##pulse lenght=4 mus
    pulse_lenght=4

    measured_charge_first=np.array([3.56,2.600,1.560,8.370,7.360,6.040,4.910,2.350,0.910,0.450])

    measured_charge_second=np.array([3.59,2.610,1.570,8.330,7.340,6.030,4.930,2.390,0.990,0.450])

    measured_charge_third=np.array([3.59,2.610,1.570,8.08,7.310,6.010,4.910,2.380,0.990,0.440])

    measured_charge_no_noise_4p=mean_array_calculator(measured_charge_first,measured_charge_second,measured_charge_third)

    measured_charge_mean_no_noise_4p=unumpy.nominal_values(measured_charge_no_noise_4p)

    measured_charge_sigma_no_noise_4p=unumpy.std_devs(measured_charge_no_noise_4p)

    dark_noise=unumpy.uarray([-0.19,-0.19,-0.19,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44],[0.03,0.03,0.03,0.01,0.01,0.01,0.01,0.01,0.01,0.01])
    measured_charge_4p = measured_charge_no_noise_4p - dark_noise

    measured_charge_per_pulse_4p = measured_charge_4p/pulses

    measured_charge_per_pulse_mean_4p=unumpy.nominal_values(measured_charge_per_pulse_4p)

    measured_charge_per_pulse_sigma_4p=unumpy.std_devs(measured_charge_per_pulse_4p)

    alfa_4p,kappa_4p,sig_alfa_4p,sig_kappa_4p=plotter(measured_charge_per_pulse_mean_4p,measured_charge_per_pulse_sigma_4p,[0.01,0.8],cord_2=3)


##pulse lenght=3 mus
    pulse_lenght=3

    measured_charge_first=np.array([3.21,2.06,1.22,6.49,5.720,4.68,3.820,1.84,0.750,0.33])

    measured_charge_second=np.array([3.21,2.05,1.23,6.47,5.710,4.67,3.810,1.85,0.760,0.34])

    measured_charge_third=np.array([3.20,2.05,1.22,6.45,5.700,4.66,3.820,1.84,0.750,0.33])

    measured_charge_no_noise_3p=mean_array_calculator(measured_charge_first,measured_charge_second,measured_charge_third)

    measured_charge_mean_no_noise_3p=unumpy.nominal_values(measured_charge_no_noise_3p)

    measured_charge_sigma_no_noise_3p=unumpy.std_devs(measured_charge_no_noise_3p)

    dark_noise=unumpy.uarray([-0.14,-0.14,-0.14,-0.33,-0.33,-0.33,-0.33,-0.33,-0.33,-0.33],[0,0,0,0.01,0.01,0.01,0.01,0.01,0.01,0.01])
    measured_charge_3p = measured_charge_no_noise_3p - dark_noise

    measured_charge_per_pulse_3p = measured_charge_3p/pulses

    measured_charge_per_pulse_mean_3p=unumpy.nominal_values(measured_charge_per_pulse_3p)

    measured_charge_per_pulse_sigma_3p=unumpy.std_devs(measured_charge_per_pulse_3p)


    alfa_3p,kappa_3p,sig_alfa_3p,sig_kappa_3p=plotter(measured_charge_per_pulse_mean_3p,measured_charge_per_pulse_sigma_3p,p0_=[0.16,1])

##pulse lenght=2 mus
    pulse_lenght=2

    measured_charge_first=np.array([2.400,1.460,0.86,4.56,4.00,3.28,2.75,1.27,0.49,0.21])

    measured_charge_second=np.array([2.410,1.460,0.86,4.55,4.01,3.29,2.75,1.26,0.50,0.22])

    measured_charge_third=np.array([2.400,1.450,0.86,4.54,4.01,3.29,2.74,1.27,0.50,0.205])

    measured_charge_no_noise_2p=mean_array_calculator(measured_charge_first,measured_charge_second,measured_charge_third)

    measured_charge_mean_no_noise_2p=unumpy.nominal_values(measured_charge_no_noise_2p)

    measured_charge_sigma_no_noise_2p=unumpy.std_devs(measured_charge_no_noise_2p)

    dark_noise=unumpy.uarray([-0.11,-0.11,-0.11,-0.23,-0.23,-0.23,-0.23,-0.23,-0.23,-0.23],[0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])

    measured_charge_2p = measured_charge_no_noise_2p - dark_noise

    measured_charge_per_pulse_2p = measured_charge_2p/pulses

    measured_charge_per_pulse_mean_2p=unumpy.nominal_values(measured_charge_per_pulse_2p)

    measured_charge_per_pulse_sigma_2p=unumpy.std_devs(measured_charge_per_pulse_2p)

    alfa_2p,kappa_2p,sig_alfa_2p,sig_kappa_2p=plotter(measured_charge_per_pulse_mean_2p,measured_charge_per_pulse_sigma_2p,cord_2=1.7)

##pulse lenght=1 mus
    pulse_lenght=1


    measured_charge_first=np.array([1.39,0.83,0.49,2.58,2.27,1.85,1.52,0.69,0.25,0.09])

    measured_charge_second=np.array([1.39,0.83,0.49,2.57,2.26,1.84,1.53,0.68,0.25,0.10])

    measured_charge_third=np.array([1.38,0.82,0.48,2.58,2.27,1.86,1.53,0.68,0.26,0.11])

    measured_charge_no_noise_1p=mean_array_calculator(measured_charge_first,measured_charge_second,measured_charge_third)

    measured_charge_mean_no_noise_1p=unumpy.nominal_values(measured_charge_no_noise_1p)

    measured_charge_sigma_no_noise_1p=unumpy.std_devs(measured_charge_no_noise_1p)

    dark_noise=unumpy.uarray([-0.08,-0.08,-0.08,-0.13,-0.13,-0.13,-0.13,-0.13,-0.13,-0.13],[0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])


    measured_charge_1p = measured_charge_no_noise_1p - dark_noise

    measured_charge_per_pulse_1p = measured_charge_1p/pulses

    measured_charge_per_pulse_mean_1p=unumpy.nominal_values(measured_charge_per_pulse_1p)

    measured_charge_per_pulse_sigma_1p=unumpy.std_devs(measured_charge_per_pulse_1p)



    alfa_1p,kappa_1p,sig_alfa_1p,sig_kappa_1p=plotter(measured_charge_per_pulse_mean_1p,measured_charge_per_pulse_sigma_1p,p0_=[0.1,1],cord_1=9,cord_2=1)

## fixed dpp, charge vs pulse charge

    pulses_plastic= 4

    pulse_lenght_plastic=np.array([1,2,3,4]).reshape((-1, 1))

    charge_plastic=np.array([measured_charge_per_pulse_mean_1p,measured_charge_per_pulse_mean_2p,measured_charge_per_pulse_mean_3p,measured_charge_per_pulse_mean_4p])

    charge_sigma_plastic=np.array([measured_charge_per_pulse_sigma_1p,measured_charge_per_pulse_sigma_2p,measured_charge_per_pulse_sigma_3p,measured_charge_per_pulse_sigma_4p])


    plt.figure("bigplot")
    plt.title("EJ212 Charge vs Pulse Lenght")
    for i in range(0,len(ssd_gaf_pl)):
        plt.subplot(5,2,i+1,)

        charge=np.array([measured_charge_per_pulse_mean_1p[i],measured_charge_per_pulse_mean_2p[i],measured_charge_per_pulse_mean_3p[i],measured_charge_per_pulse_mean_4p[i]])
        model = LinearRegression().fit(pulse_lenght_plastic, charge)
        r_sq = model.score(pulse_lenght_plastic, charge)
        q_=model.intercept_
        m_=model.coef_[0]
        x=np.linspace(0,max(pulse_lenght_plastic)+1,100)

        plt.xlabel(r"pulse lenght [$\mu$s]")
        plt.ylabel("CPP [nC/p]")
        plt.scatter(pulse_lenght_plastic, charge, marker='o',color='red',label=f"EJ212 at DPP = {ddp_gaf_pl[i]} Gy/p")
        plt.plot(x,m_*x+q_,linestyle='--',color='magenta',label="linear model")
        plt.legend()
        plt.grid()
    plt.show()

##Plot with all the points
    aplhas=np.array([alfa_1p,alfa_2p,alfa_3p,alfa_4p])
    kappas=np.array([kappa_1p,kappa_2p,kappa_3p,kappa_4p])
    colors=['red','blue','green','magenta']
    plt.figure("Ej212 Total")
    for i,item in enumerate(charge_plastic):
        y=item
        sigma_y=charge_sigma_plastic[i]
        alpha=aplhas[i]
        kappa=kappas[i]
        pl=pulse_lenght_plastic[i][0]
        plt.xlabel("Dose per pulse [Gy/p]")

        plt.ylabel("Charge per pulse [nC/p]")



        plt.title(f"Saturation for {detector_name} at various pulse lenghts")


        plt.plot(np.linspace(0,20,100),model_inv(np.linspace(0,20,100),alpha,kappa),linestyle='-',color=colors[i],label=str(pl)+r'$\mu$s')



        plt.scatter(ddp_gaf_pl,y,marker='.',color=colors[i])

        plt.errorbar(ddp_gaf_pl,y,yerr=sigma_y, xerr=gaf_uncertainty[gaf_points], ls='',color=colors[i])

    labelLines(plt.gca().get_lines(),zorder=2.5)
    plt.grid()
    plt.show()


##LYSO
    ssd=np.array([73,162,190,217,114,122])

    ssd.sort()

    gaf_points=get_index(ssd,ssd_gaf)

    ssd_gaf_pl=ssd_gaf[gaf_points]

    ddp_gaf_pl=ddp_gaf[gaf_points]

    detector_name_="LYSO"

    pulses=np.array([1,2,2,4,4,4])

    dark_noise_means = np.array([0.03,-0.03,-0.03,-0.19,-0.19,-0.18])

    dark_noise_sigmas = np.array([0.001,0.01,0.01,0.01,0.01,0.03])

##Pulse lenght = 4
    pulse_lenght=4

    measured_charge_mean_4l=np.array([np.mean([3.74,3.74,3.76,3.75,3.75]),np.mean([5.27,5.26,5.24]),np.mean([4.55,4.55,4.56]),np.mean([4.39,4.39,4.39,4.37,4.38]),np.mean([2.66,2.64,2.64]),

        np.mean([1.37,1.38,1.38])])

    measured_charge_sigma_4l=np.array([np.std([3.74,3.74,3.76,3.75,3.75]),np.std([5.27,5.26,5.24]),np.std([4.55,4.55,4.56]),np.std([4.39,4.39,4.39,4.37,4.38]),np.std([2.66,2.64,2.64]),np.std([1.37,1.38,1.38])])



    dark_noise=unumpy.uarray(dark_noise_means,dark_noise_sigmas)

    measured_charge_4l = unumpy.uarray(measured_charge_mean_4l,measured_charge_sigma_4l) - dark_noise

    measured_charge_per_pulse_4l=measured_charge_4l/pulses

    measured_charge_per_pulse_mean_4l=unumpy.nominal_values(measured_charge_per_pulse_4l)

    measured_charge_per_pulse_sigma_4l=unumpy.std_devs(measured_charge_per_pulse_4l)


    alfa_4l,kappa_4l,sig_alfa_4l,sig_kappa_4l=plotter(measured_charge_per_pulse_mean_4l,measured_charge_per_pulse_sigma_4l,detector_name=detector_name_)

##Pulse lenght = 3

    pulse_lenght=3

    measured_charge_mean_3l=np.array([np.mean([3.57,3.58,3.58,3.59,3.59]),np.mean([4.13,4.13,4.14]),np.mean([3.57,3.56,3.56]),np.mean([3.44,3.45,3.43]),np.mean([2.07,2.06,2.07]),

        np.mean([1.08,1.09,1.08])])

    measured_charge_sigma_3l=np.array([np.std([3.57,3.58,3.58,3.59,3.59]),np.std([4.13,4.13,4.14]),np.std([3.57,3.56,3.56]),np.std([3.44,3.45,3.43]),np.std([2.06,2.07,2.07]),np.std([1.08,1.09,1.08])])



    measured_charge_3l = unumpy.uarray(measured_charge_mean_3l,measured_charge_sigma_3l) - dark_noise

    measured_charge_per_pulse_3l=measured_charge_3l/pulses

    measured_charge_per_pulse_mean_3l=unumpy.nominal_values(measured_charge_per_pulse_3l)

    measured_charge_per_pulse_sigma_3l=unumpy.std_devs(measured_charge_per_pulse_3l)


    alfa_3l,kappa_3l,sig_alfa_3l,sig_kappa_3l=plotter(measured_charge_per_pulse_mean_3l,measured_charge_per_pulse_sigma_3l,detector_name=detector_name_)


##Pulse lenght = 2

    pulse_lenght=2

    measured_charge_mean_2l=np.array([np.mean([3.36,3.36,3.36,3.37,3.37]),np.mean([2.91,2.92,2.91]),np.mean([2.5,2.51,2.51]),np.mean([2.4,2.41,2.43]),np.mean([1.44,1.45,1.45]),

        np.mean([0.76,0.76,0.77])])

    measured_charge_sigma_2l=np.array([4*np.std([3.36,3.36,3.36,3.37,3.37]),4*np.std([2.91,2.92,2.91]),4*np.std([2.5,2.51,2.51]),np.std([2.4,2.41,2.43]),np.std([1.44,1.45,1.45]),np.std([0.76,0.76,0.77])])

    measured_charge_2l = unumpy.uarray(measured_charge_mean_2l,measured_charge_sigma_2l) - dark_noise

    measured_charge_per_pulse_2l=measured_charge_2l/pulses

    measured_charge_per_pulse_mean_2l=unumpy.nominal_values(measured_charge_per_pulse_2l)

    measured_charge_per_pulse_sigma_2l=unumpy.std_devs(measured_charge_per_pulse_2l)


    alfa_2l,kappa_2l,sig_alfa_2l,sig_kappa_2l=plotter(measured_charge_per_pulse_mean_2l,measured_charge_per_pulse_sigma_2l,detector_name=detector_name_)


##Pulse lenght = 1

    pulse_lenght=1

    measured_charge_mean_1l=np.array([np.mean([2.71,2.72,2.71,2.71,2.71]),np.mean([1.64,1.65,1.65]),np.mean([1.42,1.41,1.42]),np.mean([1.37,1.36,1.36]),np.mean([0.81,0.82,0.81]),

        np.mean([0.43,0.43,0.44])])

    measured_charge_sigma_1l=np.array([4*np.std([2.71,2.72,2.71,2.71,2.71]),4*np.std([1.64,1.65,1.65]),4*np.std([1.42,1.41,1.42]),np.std([1.37,1.36,1.36]),np.std([0.81,0.82,0.81]),np.std([0.43,0.43,0.44])])

    measured_charge_1l = unumpy.uarray(measured_charge_mean_1l,measured_charge_sigma_1l) - dark_noise

    measured_charge_per_pulse_1l=measured_charge_1l/pulses

    measured_charge_per_pulse_mean_1l=unumpy.nominal_values(measured_charge_per_pulse_1l)

    measured_charge_per_pulse_sigma_1l=unumpy.std_devs(measured_charge_per_pulse_1l)


    alfa_1l,kappa_1l,sig_alfa_1l,sig_kappa_1l=plotter(measured_charge_per_pulse_mean_1l,measured_charge_per_pulse_sigma_1l,detector_name=detector_name_)



## fixed dpp, charge vs pulse charge LYSO


    pulse_lenght_lyso=np.array([1,2,3,4]).reshape((-1, 1))

    charge_lyso=np.array([measured_charge_per_pulse_mean_1l,measured_charge_per_pulse_mean_2l,measured_charge_per_pulse_mean_3l,measured_charge_per_pulse_mean_4l])

    charge_sigma_lyso=np.array([measured_charge_per_pulse_sigma_1l,measured_charge_per_pulse_sigma_2l,measured_charge_per_pulse_sigma_3l,measured_charge_per_pulse_sigma_4l])


    plt.figure("bigplot")
    plt.title("LYSO Charge vs Pulse Lenght")
    for i in range(0,len(ssd_gaf_pl)):
        plt.subplot(5,2,i+1,)

        charge=np.array([measured_charge_per_pulse_mean_1l[i],measured_charge_per_pulse_mean_2l[i],measured_charge_per_pulse_mean_3l[i],measured_charge_per_pulse_mean_4l[i]])
        model = LinearRegression().fit(pulse_lenght_lyso, charge)
        r_sq = model.score(pulse_lenght_lyso, charge)
        q_=model.intercept_
        m_=model.coef_[0]
        x=np.linspace(0,max(pulse_lenght_lyso)+1,100)
        plt.xlabel(r"pulse lenght [$\mu$s]")
        plt.ylabel("CPP [nC/p]")
        plt.scatter(pulse_lenght_lyso, charge, marker='o',color='red',label=f"LYSO at DPP = {ddp_gaf_pl[i]} Gy/p")
        plt.plot(x,m_*x+q_,linestyle='--',color='magenta',label="linear model")
        plt.legend()
        plt.grid()
    plt.show()


##Plot with all the points
    aplhas=np.array([alfa_1l,alfa_2l,alfa_3l,alfa_4l])
    kappas=np.array([kappa_1l,kappa_2l,kappa_3l,kappa_4l])
    colors=['red','blue','green','magenta']
    plt.figure("LYSO Total")
    for i,item in enumerate(charge_lyso):
        y=item
        sigma_y=(charge_sigma_lyso[i])
        alpha=aplhas[i]
        kappa=kappas[i]
        pl=pulse_lenght_lyso[i][0]
        plt.xlabel("Dose per pulse [Gy/p]")

        plt.ylabel("Charge per pulse [nC/p]")


        plt.title(f"Saturation for LYSO at various pulse lenghts")


        plt.plot(np.linspace(0,20,100),model_inv(np.linspace(0,20,100),alpha,kappa),linestyle='-',color=colors[i],label=str(pl)+r'$\mu$s')



        plt.scatter(ddp_gaf_pl,y,marker='.',color=colors[i])

        plt.errorbar(ddp_gaf_pl,y,yerr=sigma_y, xerr=gaf_uncertainty[gaf_points], ls='',color=colors[i])

    labelLines(plt.gca().get_lines(),zorder=2.5)
    plt.grid()
    plt.show()

