import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from labellines import labelLine, labelLines
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from uncertainties import ufloat, unumpy
from uncertainties.umath import *


def PDD_plotter_out(dosepath, d, j=0):
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

    df = pd.read_csv(
        dosepath, names=["X", "Y", "Z", "dose [Gy]", "dosesq [Gy^2]", "entry"]
    )

    bins = len(np.array(df["dose [Gy]"][np.logical_and(df["Z"] == 0, df["Y"] == 0)]))

    lenght = d / bins

    distance = np.array([lenght * x for x in range(1, bins + 1)])

    D = np.array(df["dose [Gy]"][np.logical_and(df["Z"] == j, df["Y"] == j)])

    ds = np.array(df["dosesq [Gy^2]"][np.logical_and(df["Z"] == j, df["Y"] == j)])

    return D, ds, distance


def mean_array_calculator(*args):
    """This function calculates the mean and std dev of the measurements
       to create the array with the means and uncertainties

    :type *args: array
    :param *args: arrays to stack and calculate the mean

    :returns: array of means with uncertainties
    :rtype: unumpy array
    """
    x = np.vstack(args)
    return unumpy.uarray(np.array(x.mean(0)), np.array(x.std(0)))


def find_nearest(array, value):
    """This function finds the nearest element in the array to the selected value

    :type array: array
    :param array: selected array

    :type value: float
    :param value: value selected

    :returns: index of the nearest element
    :rtype: int
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def model(x, alpha, k):
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
    return x / ((1 / k - alpha * x))


def model_inv(x, alpha, k):
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
    return k * x / ((1 + alpha * x))


def get_index(x, y):
    """This function finds the indexes of the second array that is the same as the
       elements of the first one. Every element must not be repeated

    :type x: array
    :param x: first array

    :type y: array
    :param y: second array

    :returns: indexes
    :rtype: list
    """
    indexes = []
    for i, item in enumerate(x):
        for j, jtem in enumerate(y):
            try:
                assert item == jtem
                indexes.append(j)
            except:
                pass

    return indexes


def plotter(
    cpp_mean,
    cpp_sigma,
    ddp_gaf_pl,
    pulse_lenght,
    p0_=[0.5, 2],
    detector_name="EJ212",
    cord_1=10,
    cord_2=2,
):

    """This function does the curve fitting of CPP vs DPP and plots the resulting fit

    :type cpp_mean: array
    :param cpp_mean: mean values of the CPP

    :type cpp_sigma: array
    :param cpp_sigma: uncertainty values of the CPP

    :type ddp_gaf_pl: array
    :param ddp_gaf_pl: Gafchromic DPP

    :type pulse_lenght: int
    :param pulse_lenght: pulse lenght

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

    tmp, sigma_tmp = cpp_mean, cpp_sigma

    gaf_uncertainty = 0.03 * ddp_gaf_pl

    popt, pcov = curve_fit(
        model_inv, ddp_gaf_pl, tmp, p0_, sigma=sigma_tmp, absolute_sigma=False
    )

    sigma_alpha, sigma_k = np.diagonal(pcov)

    correlation_alpha_k, correlation_k_alpha = np.diagonal(np.fliplr(pcov))

    alpha_, k_ = popt

    stats_, pvalue = stats.chisquare(f_obs=ddp_gaf_pl, f_exp=model_inv(tmp, alpha_, k_))

    r2 = r2_score(tmp, model_inv(ddp_gaf_pl, alpha_, k_))

    plt.figure(f"pulse lenght: {pulse_lenght}")

    plt.subplot(2, 1, 1)

    plt.xlabel("Dose per pulse [Gy/p]")

    plt.ylabel("Charge per pulse [nC/p]")

    plt.xlim(0, 14)

    plt.title(f"Saturation for {detector_name} at pulse lenght {pulse_lenght} $\mu s$")

    plt.scatter(ddp_gaf_pl, tmp, marker=".", color="red", label="experimental data")

    plt.errorbar(ddp_gaf_pl, tmp, yerr=sigma_tmp, xerr=gaf_uncertainty, ls="")

    plt.plot(
        np.linspace(0, 20, 100),
        model_inv(np.linspace(0, 20, 100), alpha_, k_),
        color="magenta",
        linestyle="--",
    )

    textstr = "\n".join(
        (
            r"$y = k\cdot\frac{x}{1 + \alpha x}$",
            r"$\alpha=%.4f \pm %.4f $" % (alpha_, sigma_alpha),
            r"$k=%.4f \pm %.4f $" % (k_, sigma_k),
            r"$R^2$ score = %.3f " % (r2),
        )
    )

    props = dict(boxstyle="round", facecolor="white", alpha=0.5)

    plt.text(cord_1, cord_2, textstr, fontsize=14, verticalalignment="top", bbox=props)

    plt.grid()

    plt.subplot(2, 1, 2)

    plt.xlabel("Dose per pulse [Gy/p]")
    plt.ylabel(r"$\frac{q-f(x)}{\sigma}$")
    res = (tmp - model_inv(ddp_gaf_pl, alpha_, k_)) / sigma_tmp
    plt.scatter(ddp_gaf_pl, res, marker="o", color="green", label="residuals")
    plt.legend()
    plt.grid()

    plt.show()

    return alpha_, k_, sigma_alpha, sigma_k
