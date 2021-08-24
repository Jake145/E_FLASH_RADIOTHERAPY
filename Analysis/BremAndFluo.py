import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CURR_DIR = "../Secondaries/BremFluo"

event_directories = next(os.walk(CURR_DIR))[1]

BremAndFLuo = glob.glob(os.path.join(CURR_DIR, "*"))


datapoints = []
depth = []


for index, paths in enumerate(BremAndFLuo):

    depth.append(int(paths.split("/")[-1].split("cm")[0]))
    path = paths
    brem_data = glob.glob(os.path.join(path, "Optic_fiber*.csv"))
    Kinetic_data = glob.glob(os.path.join(path, "Kinetic*.csv"))
    optic_data = glob.glob(os.path.join(path, "Optic_*.csv"))[0:12]
    assert len(optic_data) != 0
    Brem = 0
    Fluo = 0
    scintillation = 0
    for f in brem_data:
        event_id = []
        track_id = []
        df = pd.read_csv(f, names=["column"])
        try:
            assert len(df["column"]) != 0
            for i, item in enumerate(df["column"]):
                event = str(item).split("\t")[0]
                name = str(item).split("\t")[1]
                track = str(item).split("\t")[2]
                if name[0] == "B":
                    if event in event_id and track in track_id:
                        pass
                    else:
                        Brem += 1
                        event_id.append(event)
                        track_id.append(track)
                elif name[0] == "F":
                    if event in event_id and track in track_id:
                        pass
                    else:
                        Fluo += 1
                        event_id.append(event)
                        track_id.append(track)
        except:
            pass
    for f in optic_data:
        try:
            df = pd.read_csv(f, names=["column"])
            assert len(df["column"]) != 0
            for i, item in enumerate(df["column"]):

                scint = str(item).split("\t")[1]
                scint_num = int(scint.split(":")[1])
                scintillation += scint_num
        except:
            pass

    try:
        data = scintillation / (Brem + Fluo)

    except:
        data = -1
    datapoints.append(data)
datapoints = np.array(datapoints)
depth = np.array(depth)
for i, item in enumerate(datapoints):
    if item == -1:
        depth[i] == -1
datapoints = datapoints[datapoints > -1]
depth = depth[depth > -1]
indexes = np.argsort(np.array(depth))
plt.figure(1)
plt.title("Rapporto Scintillazione/(Bremsstrhalung + Fluorescenza)")
plt.plot(
    depth[indexes],
    datapoints[indexes],
    marker="o",
    linestyle="dashed",
    color="green",
    label="cerenkov number",
)
plt.xlabel("profondità [cm]")
plt.ylabel("rapporto S/(B+F)")
# plt.legend()
plt.grid()

plt.show()
