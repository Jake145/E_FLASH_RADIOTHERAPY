"""This scrip is used for medium-equivalence coefficient calculation and PDD plotting"""

import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from uncertainties import ufloat, unumpy
from uncertainties.umath import *

sys.path.insert(0, "../")
from flash_helper.flash_functions import (PDD_plotter_out, find_nearest,
                                          mean_array_calculator)

if __name__ == "__main__":
    print("go")

    ##7 MeV Water, This section calculates the coefficients for the Novac7
    simulated_energy_7 = np.array(
        [
            3.7348,
            3.8707,
            3.71136,
            3.43725,
            3.04878,
            2.38526,
            1.44811,
            0.865002,
            0.584597,
            0.28643,
            0.0731164,
            0.00361251,
        ]
    )  # GeV
    simulated_distance_7 = np.array(
        [2, 4, 10, 12, 15, 18, 22, 25, 27, 30, 35, 40]
    )  # mm
    Run2_simulated_energy = np.array([3.82854, 3.4693, 1.82775, 0.868089])  ##GeV
    Run2_simulated_events = np.array([2490, 2801, 2300, 1632])
    Run2_distance = np.array([4, 12.12, 20.2, 25])  # mm
    validation_dose = np.array([50, 90, 100])
    validation_distance = np.array([26.5, 18, 12])
    path = "../Flash_ex_novo/VALIDATION/7mevNOVAC_10000000_4cm_100bin_10000000.csv"  # seed 42
    dose_, sq, distance = PDD_plotter_out(path, 80)

    r100 = np.round(distance[find_nearest(dose_ / np.max(dose_), 1)], 2)
    distance_2 = distance[distance > r100]
    r90 = np.round(distance_2[find_nearest(dose_ / np.max(dose_), 0.9)], 2)
    r50 = np.round(distance[find_nearest(dose_ / np.max(dose_), 0.5)], 2)
    print(f"R100 in water: {r100} mm ")
    print(f"R90 in water: {r90} mm ")
    print(f"R50 in water: {r50} mm ")

    df = pd.read_csv(path, names=["X", "Y", "Z", "dose [Gy]", "dosesq [Gy^2]", "entry"])
    tic = time.perf_counter_ns()
    indexes = [5, 15, 25, 31]
    mass_lyso = 7.4 * 1 * 0.2 * 0.2
    mass_lyso = mass_lyso * 0.001
    tic = time.perf_counter_ns()
    dose = np.array(df["dose [Gy]"])[indexes]
    x = np.divide(dose / 1e7, (Run2_simulated_energy * 1.6e-10 / mass_lyso) / 2e6)
    coefficients = np.array(x)
    toc = time.perf_counter_ns()
    print(
        "Obtained coefficients: ",
        coefficients,
        "time elapsed: ",
        toc - tic,
        "nanoseconds",
    )

    ##Plotting PDD 7MeV Novac 7
    x_r = np.array([4, 12.12, 25])
    y_r = np.array(
        [
            simulated_energy_7[1] * coefficients[0],
            simulated_energy_7[3] * coefficients[1],
            simulated_energy_7[7] * coefficients[3],
        ]
    )
    y_r = y_r / np.max(y_r) * 100
    plt.figure(1)
    plt.plot(
        distance,
        100 * dose_ / np.max(dose_),
        marker=".",
        linestyle="dashed",
        color="blue",
        label="Simulated Water 7 MeV",
    )
    plt.plot(
        simulated_distance_7,
        100 * simulated_energy_7 / np.max(simulated_energy_7),
        marker=".",
        linestyle="dashed",
        color="green",
        label="Lyso",
    )

    plt.scatter(
        validation_distance,
        validation_dose,
        marker=".",
        color="red",
        label="R100, R90, R50 Water Validation",
    )

    plt.scatter(
        x_r,
        y_r,
        marker=".",
        color="magenta",
        label="Water Equivalence corrected Points",
    )

    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Novac7 7 MeV Lyso dosimeter in water phantom simulation")

    plt.vlines([x_r[0]], 0, y_r[0], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.vlines([x_r[1]], 0, y_r[1], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.vlines([x_r[2]], 0, y_r[2], linestyles="dashed", colors="magenta", alpha=0.2)

    plt.hlines([y_r[0]], 0, x_r[0], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.hlines([y_r[1]], 0, x_r[1], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.hlines([y_r[2]], 0, x_r[2], linestyles="dashed", colors="magenta", alpha=0.2)

    plt.vlines([12], 0, 100, linestyles="dashed", colors="red", alpha=0.2)
    plt.vlines([18], 0, 90, linestyles="dashed", colors="red", alpha=0.2)
    plt.vlines([27], 0, 50, linestyles="dashed", colors="red", alpha=0.2)

    plt.hlines([100], 0, 12, linestyles="dashed", colors="red", alpha=0.2)
    plt.hlines([90], 0, 18, linestyles="dashed", colors="red", alpha=0.2)
    plt.hlines([50], 0, 27, linestyles="dashed", colors="red", alpha=0.2)

    # plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()
    ##9Mev and Dark PMMA, this section calculates the PDDs in PMMA for 9Mev EF and Dark spectrum
    path_9 = "../Flash_ex_novo/VALIDATION/9MevEF_PMMA_10mil_8cm_200bin.csv"
    dose_9_, sq_9, distance_9 = PDD_plotter_out(path_9, 80)
    path_7 = "../Flash_ex_novo/VALIDATION/7MevEF_custom_200bin_8cm.csv"
    path_7_1000 = "Send/darklyso/pmmalysodarkpdd.csv"
    path_7_54545454 = "Send/darklyso/pmmapdddarklyso2.csv"
    dose_7_1, sq_7_1, distance_7_1 = PDD_plotter_out(path_7, 80)  # seed 42
    dose_7_2, sq_7_2, distance_7_2 = PDD_plotter_out(path_7_1000, 80)  # seed 1000
    dose_7_3, sq_7_3, distance_7_3 = PDD_plotter_out(path_7_54545454, 80)  # seed 42
    dose_7 = mean_array_calculator(dose_7_1, dose_7_2, dose_7_3)
    dose_7_means = unumpy.nominal_values(dose_7)
    dose_7_stds = unumpy.std_devs(dose_7)

    path_7_1000 = "Send/darklyso/pmmalysodarkpdd.csv"
    path_7_54545454 = "Send/darklyso/pmmalysodarkpdd2.csv"

    path_106 = "Send/phantomstats/dose_106.csv"  # seed 106
    path_145 = "Send/phantomstats/dose_145.csv"  # seed 145

    dose_106, sq_106, distance_106 = PDD_plotter_out(path_106, 80)
    dose_145, sq_145, distance_145 = PDD_plotter_out(path_145, 80)

    assert len(distance_9) == len(distance_106)
    assert len(distance_9) == len(distance_145)

    dose_9 = mean_array_calculator(dose_9_, dose_106, dose_145)

    dose_9_means = unumpy.nominal_values(dose_9)
    dose_9_std = unumpy.std_devs(dose_9)

    # simulated_energy_7=np.array([3.80317,3.74554,3.77554,3.409,2.21552,0.805606])
    simulated_energy_7 = np.array(
        [3.74886, 3.68112, 3.13143, 3.01489, 2.13634, 0.95427]
    )  # lyso seed 42
    simulated_energy_7_seed_106 = np.array(
        [4.00161, 3.7975, 3.21129, 3.25056, 1.9335, 0.892416]
    )
    simulated_energy_7_seed_145 = np.array(
        [3.68242, 3.77736, 3.38892, 3.05226, 1.98661, 0.908911]
    )

    simulated_energy_7 = mean_array_calculator(
        simulated_energy_7, simulated_energy_7_seed_106, simulated_energy_7_seed_145
    )
    simulated_energy_7_means = unumpy.nominal_values(simulated_energy_7)
    simulated_energy_7_stds = unumpy.std_devs(simulated_energy_7)
    simulated_energy_9 = np.array(
        [3.99749, 3.96411, 3.96835, 3.83761, 3.5549, 2.49916]
    )  # gev lyso
    simulated_energy_9_V2 = np.array(
        [3.89034, 3.90884, 4.13303, 4.08977, 3.62473, 2.91254]
    )

    exp_distance_9 = np.array([2, 5, 10, 12, 17, 22]) + 1.5  # mm lyso
    validation_distance_9, validation_dose_9 = np.loadtxt(
        "EF_Validationdata.txt", unpack=True
    )
    # seed 42
    # ej212_energy_sim=np.array([9.3235,9.55532,9.58284,9.43716,9.43479,9.20395,8.0759,6.54801,3.88037,2.31994,0.581116,0.0558538])#gev
    # ej212_energy_sim_2=np.array([10.1766,10.5808,10.8589,10.5277,10.6088,10.6537,9.50409,7.40045,5.59706,2.67905,0.90506])#gev (starts from 6mm
    ej212_energy_sim_3 = np.array(
        [
            9.84701,
            10.0883,
            10.331,
            10.5256,
            10.5953,
            10.4832,
            10.1652,
            9.11255,
            6.69084,
            4.84278,
            2.12886,
            0.436584,
        ]
    )  # gev (starts from 5mm
    # seed 145
    ej_212_energy_sim_145 = np.array(
        [
            9.95242,
            10.0633,
            10.458,
            10.4576,
            10.508,
            10.574,
            10.1423,
            8.96335,
            6.74539,
            4.97139,
            2.09613,
            0.439371,
        ]
    )
    # seed 111555875875
    ej_212_energy_sim_lol = np.array(
        [
            9.95562,
            10.1578,
            10.3613,
            10.3675,
            10.5164,
            10.4649,
            10.4815,
            8.97829,
            6.60971,
            4.80155,
            2.12145,
            0.42088,
        ]
    )
    # seed 8675309
    ej_212_energy_sim_jenny = np.array(
        [
            9.76461,
            10.0794,
            10.5218,
            10.4084,
            10.4092,
            10.4342,
            10.424,
            8.94853,
            6.60155,
            4.76539,
            2.12755,
            0.412595,
        ]
    )
    # seed 106
    ej_212_energy_sim_106 = np.array(
        [
            9.79503,
            10.0259,
            10.2793,
            10.4171,
            10.5105,
            10.5316,
            10.2792,
            9.04665,
            6.72886,
            4.8282,
            2.0578,
            0.420623,
        ]
    )

    ej212_energy_sim = mean_array_calculator(
        ej212_energy_sim_3,
        ej_212_energy_sim_145,
        ej_212_energy_sim_lol,
        ej_212_energy_sim_jenny,
        ej_212_energy_sim_106,
    )

    ej212_energy_sim_means = unumpy.nominal_values(ej212_energy_sim)
    ej212_energy_sim_std = unumpy.std_devs(ej212_energy_sim)

    # ej212_distamce_sim_2=np.array([6,13,27])
    ej212_distance = np.array([4, 6, 9, 11, 13, 15, 17, 22, 27, 30, 35, 40]) + 1  # mm

    ej212_noise = ufloat(np.mean([-0.57, -0.54, -0.53]), np.std([-0.57, -0.54, -0.53]))

    ej212_measured_charge_first = np.array(
        [
            7.220,
            7.360,
            7.510,
            7.580,
            7.610,
            7.490,
            7.300,
            6.360,
            4.450,
            2.930,
            0.580,
            0.480,
        ]
    )

    ej212_measured_charge_second = np.array(
        [
            7.210,
            7.310,
            7.460,
            7.530,
            7.550,
            7.460,
            7.240,
            6.370,
            4.460,
            2.960,
            0.620,
            0.480,
        ]
    )

    ej212_measured_charge_third = np.array(
        [
            7.150,
            7.260,
            7.420,
            7.480,
            7.520,
            7.400,
            7.210,
            6.330,
            4.480,
            2.980,
            0.580,
            0.460,
        ]
    )

    ej212_measured_charge = (
        mean_array_calculator(
            ej212_measured_charge_first,
            ej212_measured_charge_second,
            ej212_measured_charge_third,
        )
        - ej212_noise
    )

    ej212_measured_charge_means = unumpy.nominal_values(ej212_measured_charge)

    ej212_measured_charge_std = unumpy.std_devs(ej212_measured_charge)

    ej212_measured_monitor_units = np.array(
        [
            18.43,
            18.29,
            18.24,
            18.20,
            18.10,
            18.07,
            18.10,
            18.09,
            18.15,
            18.06,
            18.16,
            18.21,
        ]
    )
    ej212_measured_pulse = 4
    ej212_charge_over_pulse = ej212_measured_charge / ej212_measured_pulse
    ej212_charge_over_mu = ej212_measured_charge / ej212_measured_monitor_units

    ej212_charge_over_mu_means = unumpy.nominal_values(ej212_charge_over_mu)
    ej212_charge_over_pulse_means = unumpy.nominal_values(ej212_charge_over_pulse)

    ##9Mev_ej212 EF medium-equivalence coefficients
    indexes_ej212 = np.array([10, 15, 22, 27, 32, 37, 40, 55, 67, 75, 87, 100]) + 2
    mass_ej212 = 1 * 1 * 2 * 0.2
    mass_ej212 = mass_ej212 * 0.001

    tic = time.perf_counter_ns()
    dose_ej212 = np.array(dose_9)[indexes_ej212]
    x = np.divide(dose_ej212 / 1e7, (ej212_energy_sim * 1.6e-10 / mass_ej212) / 2e6)
    coefficients_ej212 = np.array(x)
    toc = time.perf_counter_ns()
    print(
        "Obtained coefficients: ",
        coefficients_ej212,
        "time elapsed: ",
        toc - tic,
        "nanoseconds",
    )

    coefficients_ej212_means = unumpy.nominal_values(coefficients_ej212)

    coefficients_ej212_std = unumpy.std_devs(coefficients_ej212)

    alpha = 0.036
    kappa = 0.655
    sigma_alpha = 0.002
    sigma_kappa = 0.005
    alpha = ufloat(alpha, sigma_alpha)
    kappa = ufloat(kappa, sigma_kappa)

    x_r = ej212_distance
    y_r = ej212_energy_sim * coefficients_ej212
    y_r = y_r / np.max(y_r) * 100
    c_corr_ej = ej212_charge_over_pulse * coefficients_ej212
    # c_corr_ej=c_corr_ej/(kappa-alpha*ej212_charge_over_pulse)
    c_corr_ej_means = unumpy.nominal_values(c_corr_ej)
    c_corr_ej_err = unumpy.std_devs(c_corr_ej)
    plt.figure("ej212")
    plt.plot(
        distance_9,
        100 * unumpy.nominal_values(dose_9) / max(unumpy.nominal_values(dose_9)),
        marker=".",
        linestyle="dashed",
        color="blue",
        label="Simulated PMMA 9 MeV ",
    )
    plt.plot(
        ej212_distance,
        100
        * unumpy.nominal_values(ej212_energy_sim)
        / max(unumpy.nominal_values(ej212_energy_sim)),
        marker=".",
        linestyle="dashed",
        color="magenta",
        label="Simulated EJ212 9 MeV ",
    )
    plt.errorbar(
        ej212_distance,
        unumpy.nominal_values(
            100 * (ej212_charge_over_pulse / max(ej212_charge_over_pulse_means))
        ),
        yerr=unumpy.std_devs(
            100 * (ej212_charge_over_pulse / max(ej212_charge_over_pulse_means))
        ),
        marker=".",
        linestyle="dashed",
        color="black",
        label="Measured EJ212 9 MeV ",
    )

    plt.plot(
        validation_distance_9,
        validation_dose_9,
        marker=".",
        linestyle="dashed",
        color="red",
        label="Measured PMMA",
    )

    plt.errorbar(
        ej212_distance,
        100 * unumpy.nominal_values((c_corr_ej / max(c_corr_ej_means))),
        yerr=100 * unumpy.std_devs((c_corr_ej / max(c_corr_ej_means))),
        marker=".",
        linestyle="dashed",
        color="green",
        label="Corrected EJ212 ",
    )

    # plt.scatter(x_r,unumpy.nominal_values(y_r),marker='.',color='blue',label='pmma Equivalence corrected Points')
    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash")

    plt.legend()

    plt.grid()
    plt.show()
    ## 9mev lyso EF coefficients
    indexes = [5, 13, 25, 30, 42, 55]
    mass_lyso = 7.25 * 1 * 0.2 * 0.2
    mass_lyso = mass_lyso * 0.001
    df = pd.read_csv(
        path_9, names=["X", "Y", "Z", "dose [Gy]", "dosesq [Gy^2]", "entry"]
    )

    tic = time.perf_counter_ns()
    dose = np.array(df["dose [Gy]"])[indexes]
    x = np.divide(dose / 1e7, (simulated_energy_9_V2 * 1.6e-10 / mass_lyso) / 2e6)
    coefficients = np.array(x)
    toc = time.perf_counter_ns()
    print("Obtained coefficients: ", coefficients)

    x_r = exp_distance_9
    y_r = simulated_energy_9_V2 * coefficients
    y_r = y_r / np.max(y_r) * 100

    measured_charge = np.array([4.854, 4.834, 4.9, 4.556, 4.214, 3.12, 1.46])
    measured_distance = np.array([0, 2, 5, 10, 12, 17, 22]) + 1.5
    measured_monitor_units = np.array(
        [426.849, 427.842, 429.340, 426.420, 425.680, 426.660, 423.140]
    )
    measured_charge_per_pulse = np.array([20.6, 20.7, 20.4, 22, 23.735, 32.1, 68.5])
    ## Correction
    c_mu = measured_charge / measured_monitor_units
    c_corr = c_mu[1:] * coefficients
    # c_corr=(measured_charge_per_pulse[1:]/(1/0.56-0.4*measured_charge_per_pulse[1:]))*c_corr
    c_corr = c_corr / np.max(c_corr)
    dis_corr = measured_distance[1:]
    ## 9 Mev LYSO plot

    plt.figure(2)
    plt.plot(
        distance_9,
        100 * unumpy.nominal_values(dose_9) / max(unumpy.nominal_values(dose_9)),
        marker=".",
        linestyle="dashed",
        color="blue",
        label="pmma simulato 9 MeV ",
    )
    plt.plot(
        exp_distance_9,
        100 * simulated_energy_9_V2 / np.max(simulated_energy_9_V2),
        marker=".",
        linestyle="dashed",
        color="green",
        label="lyso simulato 9 MeV ",
    )

    plt.scatter(
        unumpy.nominal_values(x_r),
        unumpy.nominal_values(y_r),
        marker=".",
        color="magenta",
        label="Water Equivalence corrected Points",
    )

    plt.plot(
        validation_distance_9,
        validation_dose_9,
        marker=".",
        linestyle="dashed",
        color="red",
        label="pmma misurato",
    )

    plt.vlines([x_r[0]], 0, y_r[0], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.vlines([x_r[1]], 0, y_r[1], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.vlines([x_r[2]], 0, y_r[2], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.vlines([x_r[3]], 0, y_r[3], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.vlines([x_r[4]], 0, y_r[4], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.vlines([x_r[5]], 0, y_r[5], linestyles="dashed", colors="magenta", alpha=0.2)

    plt.hlines([y_r[0]], 0, x_r[0], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.hlines([y_r[1]], 0, x_r[1], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.hlines([y_r[2]], 0, x_r[2], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.hlines([y_r[3]], 0, x_r[3], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.hlines([y_r[4]], 0, x_r[4], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.hlines([y_r[5]], 0, x_r[5], linestyles="dashed", colors="magenta", alpha=0.2)
    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash")

    # plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()

    ## Dark spectrum coefficients
    indexes_lyso = [5, 13, 25, 30, 42, 55]

    tic = time.perf_counter_ns()
    dose_lyso = np.array(dose_7)[indexes_lyso]
    x = np.divide(dose_lyso / 1e7, (simulated_energy_7 * 1.6e-10 / mass_lyso) / 2e6)
    coefficients_lyso = np.array(x)
    toc = time.perf_counter_ns()
    print(
        "Obtained coefficients: ",
        coefficients_lyso,
        "time elapsed: ",
        toc - tic,
        "nanoseconds",
    )

    coefficients_lyso_means = unumpy.nominal_values(coefficients_lyso)

    coefficients_lyso_std = unumpy.std_devs(coefficients_lyso)
    print("Obtained coefficients: ", coefficients_lyso)

    ## Dark spectrum measurements and PDD plot
    measured_charge_lyso_1 = np.array([4.85, 4.90, 5.02, 4.64, 4.30, 3.28, 1.47])
    measured_charge_lyso_2 = np.array([4.89, 4.89, 4.89, 4.57, 4.26, 3.20, 1.45])
    measured_charge_lyso_3 = np.array([4.84, 4.81, 4.87, 4.58, 4.20, 3.06, 1.47])
    measured_charge_lyso_4 = np.array([4.80, 4.76, 4.88, 4.51, 4.16, 3.04, 1.49])
    measured_charge_lyso_5 = np.array([4.89, 4.81, 4.84, 4.48, 4.15, 3.02, 1.42])
    measured_charge = mean_array_calculator(
        measured_charge_lyso_1,
        measured_charge_lyso_2,
        measured_charge_lyso_3,
        measured_charge_lyso_4,
        measured_charge_lyso_5,
    )
    noise = ufloat(
        np.mean([0.310, 0.310, 0.310, 0.310, 0.320]),
        np.std([0.310, 0.310, 0.310, 0.310, 0.320]),
    )
    measured_distance = np.array([0, 2, 5, 10, 12, 17, 22]) + 1.5
    charge_per_pulse = (measured_charge - noise) / 100

    x_r = measured_distance[1:]
    y_r = charge_per_pulse[1:] * coefficients_lyso
    y_r = y_r / np.max(unumpy.nominal_values(y_r)) * 100
    validation_distance_7, _, validation_dose_7 = np.loadtxt(
        "EF_val_7mev.txt", unpack=True
    )

    plt.figure(5)
    plt.plot(
        distance_9,
        100 * dose_7_means / np.max(dose_7_means),
        marker=".",
        linestyle="dashed",
        color="blue",
        label="Simulated PMMA  ",
    )
    plt.plot(
        validation_distance_7,
        validation_dose_7,
        marker=".",
        linestyle="dashed",
        color="red",
        label="Measured PMMA validation",
    )

    plt.errorbar(
        measured_distance[1:],
        100
        * unumpy.nominal_values(charge_per_pulse)[1:]
        / np.max(unumpy.nominal_values(charge_per_pulse)[1:]),
        yerr=unumpy.std_devs(100 * charge_per_pulse[1:] / max(charge_per_pulse[1:])),
        marker=".",
        linestyle="--",
        color="black",
        label="Measured LYSO",
    )

    plt.plot(
        exp_distance_9,
        100
        * unumpy.nominal_values(simulated_energy_7)
        / np.max(unumpy.nominal_values(simulated_energy_7)),
        marker=".",
        linestyle="dashed",
        color="green",
        label="Simulated LYSO  ",
    )
    plt.errorbar(
        x_r,
        unumpy.nominal_values(y_r),
        yerr=unumpy.std_devs(y_r),
        marker=".",
        linestyle="",
        color="magenta",
        label="PMMA Equivalence corrected Points",
    )

    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.title("Dose Distribution ElectronFlash Dark Spectrum")

    # plt.xlim(0,50)
    plt.legend()

    plt.grid()
    plt.show()
    ## Dark spectrum plotting
    x, y = np.loadtxt("Send/efdark.txt.txt", unpack=True)
    plt.figure("dark spectrum")
    plt.title("ElectronFlash Dark Spectrum")
    plt.xlabel("energy [MeV]")
    plt.ylabel("intensity")
    plt.bar(x, y, color="orange")
    plt.xlim(0, 9)
    plt.grid()
    plt.show()