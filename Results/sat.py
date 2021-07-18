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
def model(x,alpha,k):
    return x/((1/k-alpha*x))

def model_inv(x,alpha,k):
    return k*x/((1+alpha*x))

def get_index(x,y):
    indexes=[]
    for i,item in enumerate(x):
        for j,jtem in enumerate(y):
            try:
                assert(item==jtem)
                indexes.append(j)
            except:
                pass


    return indexes

def tmp_maker(x,y,z,w):
    tmp=(x - y)/z

    sigma_tmp=w/z

    return tmp,sigma_tmp

def plotter(measured_charge_mean,pulses,dark_noise,measured_charge_sigma,p0_=[0.5,2],detector_name="EJ212",cord_1=10,cord_2=2):

    tmp,sigma_tmp= tmp_maker(measured_charge_mean,dark_noise,pulses,measured_charge_sigma)

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

    #plt.ylim(0,200)

    #plt.xlim(-0.5,5)

    plt.title(f"Saturation for {detector_name} at pulse lenght {pulse_lenght} $\mu s$")

    plt.scatter(ddp_gaf_pl,tmp,marker='.',color='red',label="experimental data")

    plt.errorbar(ddp_gaf_pl,tmp,yerr=sigma_tmp, xerr=gaf_uncertainty[gaf_points], ls='')

    plt.plot(np.linspace(0,20,100),model_inv(np.linspace(0,20,100),alpha_,k_),color='magenta',linestyle='--')



    textstr = '\n'.join((

        r'$y = k\cdot \frac{x}{1 + \alpha x}$',

        r'$\alpha=%.3f \pm %.3f $' % (alpha_,sigma_alpha ),

        r'$k=%.3f \pm %.3f $' % (k_,sigma_k ),

        #r'$\chi ^2 redux$ = %.3f , pvalue = %.3f ' %(stats_/(len(tmp)-2),pvalue)))

        r'$R^2$ score = %.3f ' % (r2)))

    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

    # place a text box in upper left in axes coords
    plt.text( cord_1,cord_2, textstr,  fontsize=14,
            verticalalignment='top', bbox=props)

    plt.grid()

    plt.subplot(2,1,2)

    #plt.figure(f"pulse lenght: {pulse_lenght} residuals")
    #plt.title(f"Residuals for {detector_name} at pulse lenght {pulse_lenght} $\mu s$")
    plt.xlabel("Dose per pulse [Gy/p]")
    plt.ylabel(r"$\frac{q-f(x)}{\sigma}$")
    res = (tmp - model_inv(ddp_gaf_pl,alpha_,k_))/sigma_tmp
    plt.scatter(ddp_gaf_pl,res,marker='o', color='green',label='residuals')
    plt.legend()
    plt.grid()

    #plt.savefig(f"C:/Users/pensa/Desktop/Results/{detector_name}_{pulse_lenght}.png")

    plt.figure(f"pulse lenght: {pulse_lenght} residuals")
    plt.title(f"Residuals for {detector_name} at pulse lenght {pulse_lenght} $\mu s$")
    plt.xlabel("Dose per pulse [Gy/p]")
    plt.ylabel(r"$\frac{q-f(x)}{\sigma}$")
    res = (tmp - model_inv(ddp_gaf_pl,alpha_,k_))/sigma_tmp
    plt.scatter(ddp_gaf_pl,res,marker='o', color='green',label='residuals')
    plt.legend()
    plt.grid()


    #plt.show()


    return alpha_,k_,sigma_alpha,sigma_k



if __name__ == "__main__":

    ddp_gaf=np.unique(np.array([3.1,2.9,5.8,3.5,12.5,5.8,4.5,3.5,3.1,2.9,2.7,2,1.2,0.7,0.5]))

    ssd_gaf=np.unique(np.array([114,122,90,105,73,90,106,105,114,122,123,134,162,190,217]))

    ddp_gaf=ddp_gaf[np.argsort(ssd_gaf)]

    ssd_gaf=ssd_gaf[np.argsort(ssd_gaf)][::-1]

    gaf_uncertainty=0.03*ddp_gaf

    ssd=np.array([73,90,105,106,114,122,134,162,190,217])

    gaf_points=get_index(ssd,ssd_gaf)

    ssd_gaf_pl=ssd_gaf[gaf_points]

    ddp_gaf_pl=ddp_gaf[gaf_points]

    detector_name="EJ212"



##pulse lenght=4 mus
    pulse_lenght=4

    measured_charge_mean_4p=np.array([4*np.mean([3.56,3.59,3.59]),4*np.mean([2.600,2.610,2.610]),4*np.mean([1.560,1.570,1.570]),np.mean([8.370,8.330,8.270,8.08]),np.mean([7.360,7.340,7.310]),

        np.mean([6.040,6.030,6.010]),np.mean([4.910,4.930,4.910]),np.mean([2.350,2.390,2.380]),

        np.mean([0.910,0.990,0.990]),np.mean([0.450,0.450,0.440])])

    measured_charge_sigma_4p=np.array([4*np.std([3.56,3.59,3.59]),4*np.std([2.600,2.610,2.610]),4*np.std([1.560,1.570,1.570]),np.std([8.370,8.330,8.270,8.08]),np.std([7.360,7.340,7.310]),np.std([6.040,6.030,6.010]),

        np.std([4.910,4.930,4.910]),np.std([2.350,2.390,2.380]),np.std([0.910,0.990,0.990]),np.std([0.450,0.450,0.440])])

    pulses=4

    dark_noise_4p=np.array([-0.19*4,-0.19*4,-0.19*4,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44])

    alfa_4p,kappa_4p,sig_alfa_4p,sig_kappa_4p=plotter(measured_charge_mean_4p,pulses,dark_noise_4p,measured_charge_sigma_4p,[0.3,1])
##pulse lenght=3 mus
    pulse_lenght=3

    measured_charge_mean_3p=np.array([4*np.mean([3.21,3.21,3.20]),4*np.mean([2.06,2.05,2.05]),4*np.mean([1.22,1.23,1.22]),np.mean([6.49,6.47,6.45]),np.mean([5.720,5.710,5.700]),

        np.mean([4.68,4.67,4.66]),np.mean([3.820,3.810,3.820]),np.mean([1.84,1.85,1.84]),

        np.mean([0.750,0.760,0.750]),np.mean([0.33,0.34,0.33])])

    measured_charge_sigma_3p=np.array([4*np.std([3.21,3.21,3.20]),4*np.std([2.06,2.05,2.05]),4*np.std([1.22,1.23,1.22]),np.std([6.49,6.47,6.45]),np.std([5.720,5.710,5.700]),np.std([4.68,4.67,4.66]),

        np.std([3.820,3.810,3.820]),np.std([1.84,1.85,1.84]),np.std([0.750,0.760,0.750]),np.std([0.33,0.34,0.33])])

    pulses=4

    dark_noise_3p=np.array([-0.14*4,-0.14*4,-0.14*4,-0.33,-0.33,-0.33,-0.33,-0.33,-0.33,-0.33])

    alfa_3p,kappa_3p,sig_alfa_3p,sig_kappa_3p=plotter(measured_charge_mean_3p,pulses,dark_noise_3p,measured_charge_sigma_3p,[0.16,1])

##pulse lenght=2 mus
    pulse_lenght=2

    measured_charge_mean_2p=np.array([4*np.mean([2.400,2.410,2.400]),4*np.mean([1.460,1.460,1.450]),4*np.mean([0.87,0.86,0.86]),np.mean([4.56,4.55,4.54]),np.mean([4.00,4.01,4.01]),

        np.mean([3.28,3.29,3.29]),np.mean([2.75,2.75,2.74]),np.mean([1.27,1.26,1.27]),

        np.mean([0.49,0.50,0.50]),np.mean([0.21,0.22,0.205])])

    measured_charge_sigma_2p=np.array([4*np.std([2.400,2.410,2.400]),4*np.std([1.460,1.460,1.450]),4*np.std([0.87,0.86,0.86]),np.std([4.56,4.55,4.54]),np.std([4.00,4.01,4.01]),np.std([3.28,3.29,3.29]),

        np.std([2.75,2.75,2.74]),np.std([1.27,1.26,1.27]),np.std([0.49,0.50,0.50]),np.std([0.21,0.22,0.205])])

    pulses=4

    dark_noise_2p=np.array([-0.11*4,-0.11*4,-0.11*4,-0.23,-0.23,-0.23,-0.23,-0.23,-0.23,-0.23])

    alfa_2p,kappa_2p,sig_alfa_2p,sig_kappa_2p=plotter(measured_charge_mean_2p,pulses,dark_noise_2p,measured_charge_sigma_2p)

##pulse lenght=1 mus
    pulse_lenght=1

    measured_charge_mean_1p=np.array([4*np.mean([1.39,1.39,1.38]),4*np.mean([0.83,0.83,0.82]),4*np.mean([0.49,0.49,0.48]),np.mean([2.58,2.57,2.58]),np.mean([2.27,2.26,2.27]),

        np.mean([1.85,1.84,1.86]),np.mean([1.52,1.53,1.53]),np.mean([0.69,0.68,0.68]),

        np.mean([0.25,0.25,0.26]),np.mean([0.09,0.10,0.11])])

    measured_charge_sigma_1p=np.array([4*np.std([1.39,1.39,1.38]),4*np.std([0.83,0.83,0.82]),4*np.std([0.49,0.49,0.48]),np.std([2.58,2.58,2.57]),np.std([2.27,2.26,2.27]),np.std([1.85,1.84,1.86]),

        np.std([1.53,1.53,1.52]),np.std([0.69,0.68,0.68]),np.std([0.25,0.25,0.26]),np.std([0.09,0.10,0.11])])

    pulses=4

    dark_noise_1p=np.array([-0.08*4,-0.08*4,-0.08*4,-0.13,-0.13,-0.13,-0.13,-0.13,-0.13,-0.13])

    alfa_1p,kappa_1p,sig_alfa_1p,sig_kappa_1p=plotter(measured_charge_mean_1p,pulses,dark_noise_1p,measured_charge_sigma_1p,[0.1,1],cord_1=9,cord_2=1)

## fixed dpp, charge vs pulse charge

    pulses_plastic= 4

    pulse_lenght_plastic=np.array([1,2,3,4]).reshape((-1, 1))

    charge_plastic=np.array([measured_charge_mean_1p,measured_charge_mean_2p,measured_charge_mean_3p,measured_charge_mean_4p])

    charge_sigma_plastic=np.array([measured_charge_sigma_1p,measured_charge_sigma_2p,measured_charge_sigma_3p,measured_charge_sigma_4p])

    dark_plastic=np.array([dark_noise_1p,dark_noise_2p,dark_noise_3p,dark_noise_4p])

    plt.figure("bigplot")
    plt.title("EJ212 Charge vs Pulse Lenght")
    for i in range(0,len(ssd_gaf_pl)):
        plt.subplot(5,2,i+1,)

        charge=np.array([measured_charge_mean_1p[i],measured_charge_mean_2p[i],measured_charge_mean_3p[i],measured_charge_mean_4p[i]])
        model = LinearRegression().fit(pulse_lenght_plastic, charge)
        r_sq = model.score(pulse_lenght_plastic, charge)
        q_=model.intercept_
        m_=model.coef_[0]
        x=np.linspace(0,max(pulse_lenght_plastic)+1,100)
        #plt.figure(f"Ej212_DPP={ddp_gaf_pl[i]}")
        #plt.title(f"DPP = {ddp_gaf_pl[i]} Gy/p")
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
        y=(item - dark_plastic[i])/pulses_plastic
        sigma_y=(charge_sigma_plastic[i])/pulses_plastic
        alpha=aplhas[i]
        kappa=kappas[i]
        pl=pulse_lenght_plastic[i][0]
        plt.xlabel("Dose per pulse [Gy/p]")

        plt.ylabel("Charge per pulse [nC/p]")

        #plt.ylim(0,200)

        #plt.xlim(-0.5,5)

        plt.title(f"Saturation for {detector_name} at various pulse lenghts")


        plt.plot(np.linspace(0,20,100),model_inv(np.linspace(0,20,100),alpha,kappa),linestyle='-',color=colors[i],label=str(pl)+r'$\mu$s')



        plt.scatter(ddp_gaf_pl,y,marker='.',color=colors[i])

        plt.errorbar(ddp_gaf_pl,y,yerr=sigma_y, xerr=gaf_uncertainty[gaf_points], ls='',color=colors[i])

    #plt.legend()
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
##Pulse lenght = 4

    pulse_lenght=4

    measured_charge_mean_4l=np.array([4*np.mean([3.74,3.74,3.76,3.75,3.75]),2*np.mean([5.27,5.26,5.24]),2*np.mean([4.55,4.55,4.56]),np.mean([4.39,4.39,4.39,4.37,4.38]),np.mean([2.66,2.64,2.64]),

        np.mean([1.37,1.38,1.38])])

    measured_charge_sigma_4l=np.array([4*np.std([3.74,3.74,3.76,3.75,3.75]),4*np.std([5.27,5.26,5.24]),4*np.std([4.55,4.55,4.56]),np.std([4.39,4.39,4.39,4.37,4.38]),np.std([2.66,2.64,2.64]),np.std([1.37,1.38,1.38])])

    pulses=4

    dark_noise_4l=np.array([-0.03,-0.03,-0.19,-0.19,-0.19,-0.19])

    alfa_4l,kappa_4l,sig_alfa_4l,sig_kappa_4l=plotter(measured_charge_mean_4l,pulses,dark_noise_4l,measured_charge_sigma_4l,detector_name=detector_name_)

##Pulse lenght = 3

    pulse_lenght=3

    measured_charge_mean_3l=np.array([4*np.mean([3.57,3.58,3.58,3.59,3.59]),2*np.mean([4.13,4.13,4.14]),2*np.mean([3.57,3.56,3.56]),np.mean([3.44,3.45,3.43]),np.mean([2.07,2.06,2.07]),

        np.mean([1.08,1.09,1.08])])

    measured_charge_sigma_3l=np.array([4*np.std([3.57,3.58,3.58,3.59,3.59]),4*np.std([4.13,4.13,4.14]),4*np.std([3.57,3.56,3.56]),np.std([3.44,3.45,3.43]),np.std([2.06,2.07,2.07]),np.std([1.08,1.09,1.08])])

    pulses=4

    dark_noise_3l=np.array([-0.03,-0.03,-0.19,-0.19,-0.19,-0.19])

    alfa_3l,kappa_3l,sig_alfa_3l,sig_kappa_3l=plotter(measured_charge_mean_3l,pulses,dark_noise_3l,measured_charge_sigma_3l,detector_name=detector_name_)


##Pulse lenght = 2

    pulse_lenght=2

    measured_charge_mean_2l=np.array([4*np.mean([3.36,3.36,3.36,3.37,3.37]),2*np.mean([2.91,2.92,2.91]),2*np.mean([2.5,2.51,2.51]),np.mean([2.4,2.41,2.43]),np.mean([1.44,1.45,1.45]),

        np.mean([0.76,0.76,0.77])])

    measured_charge_sigma_2l=np.array([4*np.std([3.36,3.36,3.36,3.37,3.37]),4*np.std([2.91,2.92,2.91]),4*np.std([2.5,2.51,2.51]),np.std([2.4,2.41,2.43]),np.std([1.44,1.45,1.45]),np.std([0.76,0.76,0.77])])

    pulses=4

    dark_noise_2l=np.array([-0.03,-0.03,-0.19,-0.19,-0.19,-0.19])

    alfa_2l,kappa_2l,sig_alfa_2l,sig_kappa_2l=plotter(measured_charge_mean_2l,pulses,dark_noise_2l,measured_charge_sigma_2l,detector_name=detector_name_)


##Pulse lenght = 1

    pulse_lenght=1

    measured_charge_mean_1l=np.array([4*np.mean([2.71,2.72,2.71,2.71,2.71]),2*np.mean([1.64,1.65,1.65]),2*np.mean([1.42,1.41,1.42]),np.mean([1.37,1.36,1.36]),np.mean([0.81,0.82,0.81]),

        np.mean([0.43,0.43,0.44])])

    measured_charge_sigma_1l=np.array([4*np.std([2.71,2.72,2.71,2.71,2.71]),4*np.std([1.64,1.65,1.65]),4*np.std([1.42,1.41,1.42]),np.std([1.37,1.36,1.36]),np.std([0.81,0.82,0.81]),np.std([0.43,0.43,0.44])])

    pulses=4

    dark_noise_1l=np.array([-0.03,-0.03,-0.19,-0.19,-0.19,-0.19])

    alfa_1l,kappa_1l,sig_alfa_1l,sig_kappa_1l=plotter(measured_charge_mean_1l,pulses,dark_noise_1l,measured_charge_sigma_1l,detector_name=detector_name_)



## fixed dpp, charge vs pulse charge LYSO

    pulses_lyso= 4

    pulse_lenght_lyso=np.array([1,2,3,4]).reshape((-1, 1))

    charge_lyso=np.array([measured_charge_mean_1l,measured_charge_mean_2l,measured_charge_mean_3l,measured_charge_mean_4l])

    charge_sigma_lyso=np.array([measured_charge_sigma_1l,measured_charge_sigma_2l,measured_charge_sigma_3l,measured_charge_sigma_4l])

    dark_lyso=np.array([dark_noise_1l,dark_noise_2l,dark_noise_3l,dark_noise_4l])

    plt.figure("bigplot")
    plt.title("LYSO Charge vs Pulse Lenght")
    for i in range(0,len(ssd_gaf_pl)):
        plt.subplot(5,2,i+1,)

        charge=np.array([measured_charge_mean_1l[i],measured_charge_mean_2l[i],measured_charge_mean_3l[i],measured_charge_mean_4l[i]])
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
        y=(item - dark_lyso[i])/pulses_lyso
        sigma_y=(charge_sigma_lyso[i])/pulses_lyso
        alpha=aplhas[i]
        kappa=kappas[i]
        pl=pulse_lenght_lyso[i][0]
        plt.xlabel("Dose per pulse [Gy/p]")

        plt.ylabel("Charge per pulse [nC/p]")

        #plt.ylim(0,200)

        #plt.xlim(-0.5,5)

        plt.title(f"Saturation for LYSO at various pulse lenghts")


        plt.plot(np.linspace(0,20,100),model_inv(np.linspace(0,20,100),alpha,kappa),linestyle='-',color=colors[i],label=str(pl)+r'$\mu$s')



        plt.scatter(ddp_gaf_pl,y,marker='.',color=colors[i])

        plt.errorbar(ddp_gaf_pl,y,yerr=sigma_y, xerr=gaf_uncertainty[gaf_points], ls='',color=colors[i])

    #plt.legend()
    labelLines(plt.gca().get_lines(),zorder=2.5)
    plt.grid()
    plt.show()

