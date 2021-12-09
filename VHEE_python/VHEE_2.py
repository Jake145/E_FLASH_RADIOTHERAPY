import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from tsmoothie.smoother import *
from uncertainties import ufloat, unumpy
from uncertainties.umath import *
from labellines import labelLine, labelLines
###
import sys
sys.path.insert(0, "E_FLASH_RADIOTHERAPY")
from flash_helper.flash_functions import PDD_plotter_out, mean_array_calculator,find_nearest,find_r

if __name__ == "__main__":
    path_130 = "C:/Users/JakeHarold/Downloads/pencil160_spectrum_ssd700.csv"
    depth=200 #mm
    center=25 #bin
    evt=1e6
    measurement=20 #mm
    dose_130, sq_130, distance_130 = PDD_plotter_out(path_130, depth,center)
    r100_130 = find_r(distance_130, (dose_130 / np.max(dose_130)) , 100)
    r90_130 = find_r(
        distance_130[distance_130 > r100_130],
        (dose_130 / np.max(dose_130))[distance_130 > r100_130] ,
        90,
    )
    r50_130 = find_r(distance_130, (dose_130 / np.max(dose_130)) , 50)
    evt=1e6

    print(f"Pencil Beam 160 MeV - R100 = {r100_130} mm , R90 = {r90_130} mm")


    plt.figure("pencil")

    plt.plot(
        distance_130,
        100 * (dose_130 / np.max(dose_130)),
        linestyle="dashed",marker='o',
        color="red",
        label="Simulated Water 160 MeV",
    )



    plt.plot(

        distance_130[find_nearest(distance_130,measurement)],(100*dose_130 / np.max(dose_130))[find_nearest(distance_130,measurement)],
        linestyle="",marker='o',
        color="blue",
        label="Measurement point",
    )

    plt.xlabel("distance [mm]")
    plt.ylabel("relative dose[%]")
    plt.plot([],[],label=f"Dose at measurement: {dose_130[find_nearest(distance_130,measurement)]/evt} Gy/electron")
    plt.title("Dose Distribution 160 MeV Pencil Beam")

    print("number at max:",np.max(dose_130)/evt)
    plt.legend()

    plt.grid()
    plt.show()
