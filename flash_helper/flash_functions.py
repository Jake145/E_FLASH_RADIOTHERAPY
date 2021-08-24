import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.odr import *
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


def model_inv(beta, x):
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
    return beta[0] * x / ((1 + beta[1] * x))


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
    p0_=[2, 0.5],
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


    :returns: alpha, k and their uncertainties
    :rtype: float
    """

    tmp, sigma_tmp = cpp_mean, cpp_sigma

    gaf_uncertainty = 0.05 * ddp_gaf_pl

    saturation = Model(model_inv)

    mydata = RealData(ddp_gaf_pl, tmp, sx=gaf_uncertainty, sy=sigma_tmp)

    myodr = ODR(mydata, saturation, p0_)

    myodr.set_job(fit_type=0)

    myoutput = myodr.run()

    k__ = myoutput.beta[0]

    alpha__ = myoutput.beta[1]

    sigma_k__ = myoutput.sd_beta[0]

    sigma_alpha__ = myoutput.sd_beta[1]

    r2 = r2_score(tmp, model_inv([k__, alpha__], ddp_gaf_pl))

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
        model_inv([k__, alpha__], np.linspace(0, 20, 100)),
        color="magenta",
        linestyle="--",
    )

    textstr = "\n".join(
        (
            r"$y = k\cdot\frac{x}{1 + \alpha x}$",
            r"$\alpha=%.2f \pm %.2f $" % (alpha__, sigma_alpha__),
            r"$k=%.2f \pm %.2f $" % (k__, sigma_k__),
            r"$R^2$ score = %.3f " % (r2),
        )
    )

    props = dict(boxstyle="round", facecolor="white", alpha=0.5)

    plt.text(cord_1, cord_2, textstr, fontsize=14, verticalalignment="top", bbox=props)

    plt.grid()

    plt.subplot(2, 1, 2)

    plt.xlabel("Dose per pulse [Gy/p]")
    plt.ylabel(r"$\frac{q-f(x)}{\sigma}$")
    res_x = myoutput.delta
    res_y = myoutput.eps
    plt.scatter(ddp_gaf_pl, res_x, marker="o", color="green", label="X residuals")
    plt.scatter(ddp_gaf_pl, res_y, marker="o", color="red", label="Y residuals")

    plt.legend()
    plt.grid()

    plt.show()

    return alpha__, k__, sigma_alpha__, sigma_k__


def find_r(x, y, R, uncertainties=False):
    """This function finds the range value R of the PDD

    :type x: array
    :param x: depth array

    :type y: array
    :param y: dose normalized by the maximum dose (ranges from 0 to 1)

    :type R: int
    :param R: Range value desired (ex. 50, 80, 100)

    :type pulse_lenght: int
    :param pulse_lenght: pulse lenght

    :type uncertainties: bool
    :param uncertainties: set True if using Uncertainties

    :returns: index of x for the selected range
    :rtype: float
    """
    if uncertainties == False:
        idx = find_nearest(y * 100, R)
        if R != 100 and idx != 0:
            idx2 = find_nearest(y, np.partition(np.abs(y - (R / 100)), 2)[1])
            m = 100 * (y[idx] - y[idx2]) / (x[idx] - x[idx2])
            q = 100 * y[idx] - m * x[idx]
            return (R - q) / m

        else:
            return x[idx]
    else:
        idx = find_nearest(unumpy.nominal_values(y) * 100, R)

        if R != 100 and idx != 0:
            idx2 = find_nearest(
                unumpy.nominal_values(y) * 100,
                np.partition(np.abs(100 * unumpy.nominal_values(y) - (R)), 2)[1],
            )
            m = 100 * (y[idx] - y[idx2]) / (x[idx] - x[idx2])
            q = 100 * y[idx] - m * x[idx]
            return (R - q) / m

        else:
            return x[idx]
