import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CURR_DIR = "./7MeV"

event_directories = next(os.walk(CURR_DIR))[1]

datapoints_p = []
scintillations_p = []
cerenkovs_p = []
scintillations_s = []
cerenkovs_s = []
datapoints_s = []
depth = []
incident_electrons = []
fluxes = []
backscatter = []
# Scintillation and cerenkov

for index, paths in enumerate(event_directories):
    paths = os.path.join(CURR_DIR, paths)
    print(paths)
    depth.append(
        int(paths.split("/")[-1].split("mm")[0])
    )  # You gotta put the os.path.join thing
    path = paths
    optic_data = glob.glob(os.path.join(path, "Optic_*.csv"))
    assert len(optic_data) != 0
    cerenkov_p = 0
    cerenkov_s = 0
    scintillation_p = 0
    scintillation_s = 0
    for f in optic_data:
        try:
            df = pd.read_csv(f, names=["column"])
            assert len(df["column"]) != 0
            for i, item in enumerate(df["column"]):

                cerenk_p = str(item).split("\t")[2]
                cerenk_num_p = int(cerenk_p.split(":")[1])
                cerenkov_p += cerenk_num_p
                cerenk_s = str(item).split("\t")[4]
                cerenk_num_s = int(cerenk_s.split(":")[1])
                cerenkov_s += cerenk_num_s
                scint_p = str(item).split("\t")[1]
                scint_num_p = int(scint_p.split(":")[1])
                scintillation_p += scint_num_p
                scint_s = str(item).split("\t")[3]
                scint_num_s = int(scint_s.split(":")[1])
                scintillation_s += scint_num_s

        # print(scintillation/cerenkov)

        except:
            pass
        if cerenkov_p != 0:
            datapoints_p.append(scintillation_p / cerenkov_p)
        else:
            datapoints_p.append(0)
        if cerenkov_s != 0:
            datapoints_s.append(scintillation_s / cerenkov_s)
        else:
            datapoints_s.append(0)
        scintillations_p.append(scintillation_p / (scintillation_p + scintillation_s))
        scintillations_s.append(scintillation_s / (scintillation_p + scintillation_s))
        cerenkovs_p.append(cerenkov_p / (cerenkov_p + cerenkov_s))
        cerenkovs_s.append(cerenkov_s / (cerenkov_p + cerenkov_s))
    # print("Primary Scintillation: ",scintillation_p,"Primary Cerenkov: ",cerenkov_p,"Secondary Scintillation: ",scintillation_s,"Secondary Cerenkov: ",cerenkov_s)
    electron_datas = glob.glob(os.path.join(path, "Kinetic_*.csv"))
    assert len(electron_datas) != 0
    number_e = 0
    energies = []
    for f in electron_datas:
        df = pd.read_csv(f, names=["column"])
        try:
            assert len(df["column"]) != 0
            for i, item in enumerate(df["column"]):
                name = str(item).split("\t")[0]

                if name == "Incident Electron":
                    number_e += 1
                elif name == "Incoming Energy":
                    energy = float(str(item).split("\t")[2])

                    energies.append(energy)
                else:
                    print("What is wrong with this file man?")
        except:
            pass
    try:
        assert len(energies) != 0
        plt.figure("profondità: %.2f mm" % int(paths.split("/")[-1].split("mm")[0]))
        plt.title(
            "Spettro Energetico primari a %.2f mm"
            % int(paths.split("/")[-1].split("mm")[0])
        )
        bin_heights, bin_borders, _ = plt.hist(
            energies,
            round(np.sqrt(len(energies))),
            facecolor="g",
            ec="black",
            alpha=0.5,
            label="histogram data",
            density=True,
        )
        # plt.title("Spettro energetico primari incidenti allo scintillatore")
        plt.ylabel("Frequenza")
        plt.grid()
        plt.xlabel("Energie [Mev]")
        figname = "Spettro_%.2f.png" % int(paths.split("/")[-1].split("mm")[0])
        plt.savefig(os.path.join(CURR_DIR, figname))
    except:
        pass
    fluxfile = "flux.csv"
    df_ = pd.read_csv(
        os.path.join(path, fluxfile),
        names=["X", "Y", "Z", "flux [Gy]", "fluxsq [Gy^2]", "entry"],
    )
    flux_number = 0
    try:
        assert len(df_["entry"]) != 0
        for i, item in enumerate(df_["entry"]):
            flux_number += item
    except:
        pass
    fluxes.append(flux_number)
    backscatt = (flux_number - number_e) / 2
    backscatter.append(number_e / backscatt)


indexes = np.argsort(np.array(depth))


plt.figure("Optic Events Primary")
plt.title("Produced Scintillation and Cerenkov Ratio")
plt.plot(
    np.array(depth)[indexes],
    np.array(datapoints_p)[indexes],
    marker="o",
    linestyle="dashed",
    color="blue",
    label="primary scint/cerenkov ratio",
)
plt.plot(
    np.array(depth)[indexes],
    np.array(datapoints_s)[indexes],
    marker="o",
    linestyle="dashed",
    color="green",
    label="secondary scint/cerenkov ratio",
)
plt.xlabel("profondità [mm]")
plt.ylabel("rapporto scintillazione/cerenkov")
plt.legend()
plt.grid()
plt.show()

plt.figure("Total scintillation")
plt.title("Scintillations over total scintillation")
plt.plot(
    np.array(depth)[indexes],
    np.array(scintillations_p)[indexes],
    marker="o",
    linestyle="dashed",
    color="blue",
    label="primary scintillation over total",
)
plt.plot(
    np.array(depth)[indexes],
    np.array(scintillations_s)[indexes],
    marker="o",
    linestyle="dashed",
    color="green",
    label="secondary scintillation over total",
)
plt.xlabel("profondità [mm]")
plt.ylabel("numero di scintillazioni")
plt.legend()
plt.grid()
plt.show()

plt.figure("Total cerenkov")
plt.title("Cerenkov over total cerenkov")
plt.plot(
    np.array(depth)[indexes],
    np.array(cerenkovs_p)[indexes],
    marker="o",
    linestyle="dashed",
    color="blue",
    label="primary cerenkov over total",
)
plt.plot(
    np.array(depth)[indexes],
    np.array(cerenkovs_s)[indexes],
    marker="o",
    linestyle="dashed",
    color="green",
    label="secondary cerenkov over total",
)
plt.xlabel("profondità [mm]")
plt.ylabel("numero di cerenkov")
plt.legend()
plt.grid()
plt.show()

plt.figure("Backscatter")
plt.title("Backscatter electrons")
plt.plot(
    np.array(depth)[indexes],
    np.array(backscatter)[indexes],
    marker="o",
    linestyle="dashed",
    color="blue",
    label="backscatters",
)
plt.xlabel("profondità [mm]")
plt.ylabel("rapporto incidenti/backscatter")
plt.grid()
plt.show()
##Backscatter
