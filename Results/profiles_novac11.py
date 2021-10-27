"""This script plots the validations for Novac 11 in water"""
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("TkAgg")
import sys

import matplotlib.pyplot as plt
import scipy
from scipy.signal import savgol_filter
from tsmoothie.smoother import *
from uncertainties import ufloat, unumpy
from uncertainties.umath import *

sys.path.insert(0, "../")
from flash_helper.flash_functions import (PDD_plotter_out, find_nearest,
                                          find_r, mean_array_calculator)

if __name__ == "__main__":

    distance_val, validation_dose = np.loadtxt(
        "../Flash_ex_novo/VALIDATION/novac11PDD.txt", unpack=True
    )
    distancer50, doser50 = np.loadtxt(
        "../Flash_ex_novo/VALIDATION/r50.txt", unpack=True
    )
    distancer80, doser80 = np.loadtxt(
        "../Flash_ex_novo/VALIDATION/r80.txt", unpack=True
    )
    distancer100, doser100 = np.loadtxt(
        "../Flash_ex_novo/VALIDATION/r100.txt", unpack=True
    )

    df = pd.read_csv(
        "Send/validation_data/dose_novac11_profiles_10mil.csv",
        names=["X", "Y", "Z", "dose [Gy]", "dosesq [Gy^2]", "entry"],
    )

    distanceprofile = np.linspace(
        -80,
        80,
        len(np.array(df["dose [Gy]"][np.logical_and(df["X"] == 0, df["Y"] == 0)])),
    )

    j = 50
    bin = 19
    plt.figure("profilerhor", figsize=(8, 6))

    plt.subplot(3, 1, 1)
    plt.title("Novac 11 Profiles at 10 MeV at R100, R80 and R50")
    a = 0
    b = 5
    k = j

    # operate smoothing
    smoother = ConvolutionSmoother(window_len=5, window_type="ones")
    smoother.smooth(
        np.array(df["dose [Gy]"][np.logical_and(df["X"] == bin, df["Y"] == k)])
        / np.array(
            df["dose [Gy]"][
                np.logical_and(
                    np.logical_and(df["X"] == bin, df["Y"] == k), df["Z"] == 50
                )
            ]
        )
        * 100
    )

    # generate intervals
    low, up = smoother.get_intervals("sigma_interval", n_sigma=5)

    plt.xlabel("distance [mm]")
    plt.ylabel("dose [$\%$]")
    plt.grid()
    plt.plot(
        distanceprofile,
        np.array(df["dose [Gy]"][np.logical_and(df["X"] == bin, df["Y"] == k)])
        / np.array(
            df["dose [Gy]"][
                np.logical_and(
                    np.logical_and(df["X"] == bin, df["Y"] == k), df["Z"] == 50
                )
            ]
        )
        * 100,
        marker=".",
        linestyle="dashed",
        color="blue",
        label="simulation",
    )
    plt.plot(
        distancer100,
        doser100 / np.max(doser100) * 100,
        label="R100 profile",
        color="green",
    )
    plt.plot(
        distanceprofile,
        smoother.smooth_data[0],
        linestyle="dashed",
        color="red",
        label="smoothed curve",
    )
    plt.legend()

    bin = bin + 16
    plt.subplot(3, 1, 2)
    smoother.smooth(
        np.array(df["dose [Gy]"][np.logical_and(df["X"] == bin, df["Y"] == k)])
        / np.array(
            df["dose [Gy]"][
                np.logical_and(
                    np.logical_and(df["X"] == bin, df["Y"] == k), df["Z"] == 50
                )
            ]
        )
        * 100
    )

    # generate intervals
    low, up = smoother.get_intervals("sigma_interval", n_sigma=5)
    plt.xlabel("distance [mm]")
    plt.ylabel("dose [$\%$]")
    plt.grid()
    plt.plot(
        distanceprofile,
        np.array(df["dose [Gy]"][np.logical_and(df["X"] == bin, df["Y"] == k)])
        / np.array(
            df["dose [Gy]"][
                np.logical_and(
                    np.logical_and(df["X"] == bin, df["Y"] == k), df["Z"] == 50
                )
            ]
        )
        * 100,
        marker=".",
        linestyle="dashed",
        color="gray",
        label="simulation",
    )
    plt.plot(
        distancer80,
        doser80 / np.max(doser80) * 100,
        label="R80 profile",
        color="orange",
    )
    plt.plot(
        distanceprofile,
        smoother.smooth_data[0],
        linestyle="dashed",
        color="magenta",
        label="smoothed curve",
    )
    plt.legend()

    bin = bin + 11
    plt.subplot(3, 1, 3)
    smoother.smooth(
        np.array(df["dose [Gy]"][np.logical_and(df["X"] == bin, df["Y"] == k)])
        / np.array(
            df["dose [Gy]"][
                np.logical_and(
                    np.logical_and(df["X"] == bin, df["Y"] == k), df["Z"] == 50
                )
            ]
        )
        * 100
    )

    # generate intervals
    low, up = smoother.get_intervals("sigma_interval", n_sigma=5)
    plt.xlabel("distance [mm]")
    plt.ylabel("dose [$\%$]")
    plt.grid()
    plt.plot(
        distanceprofile,
        np.array(df["dose [Gy]"][np.logical_and(df["X"] == bin, df["Y"] == k)])
        / np.array(
            df["dose [Gy]"][
                np.logical_and(
                    np.logical_and(df["X"] == bin, df["Y"] == k), df["Z"] == 50
                )
            ]
        )
        * 100,
        marker=".",
        linestyle="dashed",
        color="olive",
        label="simulation",
    )
    plt.plot(
        distancer50, doser50 / np.max(doser50) * 100, label="R50 profile", color="cyan"
    )
    plt.plot(
        distanceprofile,
        smoother.smooth_data[0],
        linestyle="dashed",
        color="orange",
        label="smoothed curve",
    )
    plt.legend()
    plt.show()

    distance_val, validation_dose = np.loadtxt(
        "../Flash_ex_novo/VALIDATION/novac11PDD.txt", unpack=True
    )
    distancer50, doser50 = np.loadtxt(
        "../Flash_ex_novo/VALIDATION/r50.txt", unpack=True
    )
    distancer80, doser80 = np.loadtxt(
        "../Flash_ex_novo/VALIDATION/r80.txt", unpack=True
    )
    distancer100, doser100 = np.loadtxt(
        "../Flash_ex_novo/VALIDATION/r100.txt", unpack=True
    )
    path_novac11_1 = "Send/simulations/novac_11_11mev_PDD.csv"  # seed 106
    path_novac11_2 = "Send/simulations/novac11_pdd_seed42_11mev_200bin.csv"  # seed 42

    dose_novac11_1, sq_novac11_1, distance_novac11_1 = PDD_plotter_out(
        path_novac11_1, 80
    )
    dose_novac11_2, sq_novac11_2, distance_novac11_2 = PDD_plotter_out(
        path_novac11_2, 80
    )

    dose_novac11 = mean_array_calculator(dose_novac11_1, dose_novac11_2)

    dose_novac11_means = unumpy.nominal_values(dose_novac11)
    dose_novac11_std = unumpy.std_devs(dose_novac11)

    r100_val = find_r(distance_val, validation_dose / 100, 100)
    r90_val = find_r(
        distance_val[distance_val > r100_val],
        validation_dose[distance_val > r100_val] / 100,
        90,
    )
    r50_val = find_r(distance_val, validation_dose / 100, 50)

    r100_water = find_r(
        distance_novac11_1, dose_novac11 / max(dose_novac11_means), 100, True
    )
    r90_water = find_r(
        distance_novac11_1[distance_novac11_1 > r100_water],
        ((dose_novac11) / max(dose_novac11_means))[distance_novac11_1 > r100_water],
        90,
        True,
    )
    r50_water = find_r(
        distance_novac11_1[distance_novac11_1 > r100_water],
        ((dose_novac11) / max(dose_novac11_means))[distance_novac11_1 > r100_water],
        50,
        True,
    )
    print(f"Validation data: R100 = {r100_val}, R90={r90_val}, R50={r50_val}")
    print(f"water Sim data: R100 = {r100_water}, R90={r90_water}, R50={r50_water}")

    plt.figure("novac11")

    plt.plot(
        distance_novac11_1,
        100 * dose_novac11_means / np.max(dose_novac11_means),
        linestyle="dashed",
        color="blue",
        marker=".",
        label="Simulation",
    )

    plt.plot(
        distance_val,
        100 * validation_dose / np.max(validation_dose),
        linestyle="dashed",
        color="red",
        marker=".",
        label="Validation",
    )

    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution Novac11 10 MeV Water Validation")

    plt.legend()

    plt.grid()
    plt.show()
