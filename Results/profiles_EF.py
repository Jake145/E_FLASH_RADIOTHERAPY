"""This script plots the profiles and pdd for EF 9 MeV validations in water"""
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("TkAgg")
import sys

import matplotlib.pyplot as plt
import scipy
from scipy.signal import savgol_filter
from tsmoothie.smoother import *

sys.path.insert(0, "../")
from flash_helper.flash_functions import PDD_plotter_out, find_nearest, find_r

if __name__ == "__main__":

    df = pd.read_csv(
        "Send/Profile_ef_9_mev_10mil_100Bin.csv",
        names=["X", "Y", "Z", "dose [Gy]", "dosesq [Gy^2]", "entry"],
    )

    distancer_hor, doser_hor = np.loadtxt("efProfiles9Mev.txt", unpack=True)

    distancer_vert, doser_vert = np.loadtxt("efProfiles9Mev_vert.txt", unpack=True)

    distanceprofile = np.linspace(
        -70,
        70,
        len(np.array(df["dose [Gy]"][np.logical_and(df["X"] == 0, df["Y"] == 0)])),
    )

    j = 50
    bin = 18
    plt.figure("profilerhor", figsize=(8, 4))

    plt.subplot(2, 1, 1)
    plt.title("ElectronFlash Profiles at 9 MeV (18mm)")
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
        distancer_hor,
        doser_hor / np.max(doser_hor) * 100,
        label="horizontal profile",
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

    plt.subplot(2, 1, 2)
    smoother.smooth(
        np.array(df["dose [Gy]"][np.logical_and(df["X"] == bin, df["Z"] == k)])
        / np.array(
            df["dose [Gy]"][
                np.logical_and(
                    np.logical_and(df["X"] == bin, df["Z"] == k), df["Y"] == 50
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
        np.array(df["dose [Gy]"][np.logical_and(df["X"] == bin, df["Z"] == k)])
        / np.array(
            df["dose [Gy]"][
                np.logical_and(
                    np.logical_and(df["X"] == bin, df["Z"] == k), df["Y"] == 50
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
        distancer_hor,
        doser_hor / np.max(doser_hor) * 100,
        label="vertical profile",
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
    plt.show()
    ##pdd
    d = 80

    j = 50

    path_9_pen = "Send/physicslistsim/water_penelope_10mil.csv"

    dose_9_pen, sq_9_pen, distance_9_pen = PDD_plotter_out(path_9_pen, 80)

    bins = len(np.array(df["dose [Gy]"][np.logical_and(df["Z"] == 0, df["Y"] == 0)]))

    lenght = d / bins

    distance = np.array([lenght * x for x in range(1, bins + 1)])

    D = np.array(df["dose [Gy]"][np.logical_and(df["Z"] == j, df["Y"] == j)])

    ds = np.array(df["dosesq [Gy^2]"][np.logical_and(df["Z"] == j, df["Y"] == j)])

    validation_distance_9, validation_dose_9 = np.loadtxt(
        "Send/EF_9Mev_Water.txt", unpack=True
    )
    r100_val = find_r(validation_distance_9, validation_dose_9 / 100, 100)
    r90_val = find_r(
        validation_distance_9[validation_distance_9 > r100_val],
        validation_dose_9[validation_distance_9 > r100_val] / 100,
        90,
    )
    r50_val = find_r(validation_distance_9, validation_dose_9 / 100, 50)

    r100_water = find_r(distance_9_pen, dose_9_pen / max(dose_9_pen), 100)
    r90_water = find_r(
        distance_9_pen[distance_9_pen > r100_water],
        ((dose_9_pen) / max(dose_9_pen))[distance_9_pen > r100_water],
        90,
    )
    r50_water = find_r(
        distance_9_pen[distance_9_pen > r100_water],
        ((dose_9_pen) / max(dose_9_pen))[distance_9_pen > r100_water],
        50,
    )
    print(f"Validation data: R100 = {r100_val}, R90={r90_val}, R50={r50_val}")
    print(f"water Sim data: R100 = {r100_water}, R90={r90_water}, R50={r50_water}")

    plt.figure("pdd")
    plt.plot(
        validation_distance_9,
        validation_dose_9,
        linestyle="dashed",
        marker=".",
        color="green",
        label="Validation 9 MeV",
    )

    plt.plot(
        distance_9_pen,
        100 * dose_9_pen / np.max(dose_9_pen),
        linestyle="dashed",
        color="blue",
        marker=".",
        label="Simulated Water 9 MeV ",
    )

    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash 9 MeV water")

    plt.legend()

    plt.grid()
    plt.show()
