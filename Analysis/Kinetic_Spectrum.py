import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mainpath = ".."
secondary_path = "PDDS_SEED"
outpath = "Images_of_analysis"
paths = glob.glob(os.path.join(mainpath, secondary_path))
print(paths)
pattern = "_[0-9]*"
for p in paths:
    print(p)
    datapaths = next(os.walk(p))[1]
    for datasetpath in datapaths:
        filenames = glob.glob(os.path.join(p, datasetpath, "Kinetic_*"))

        En = []
        energies = []
        for f in filenames:
            # print(f)
            df = pd.read_csv(f, names=["column"])
            for i, item in enumerate(df["column"]):
                # print(item)
                energies.append(float(str(item).split("\t")[1]))

        En = np.array(energies)
        depthstring = re.findall(pattern, datasetpath)[-1].split("\\")[1]
        depth = float((100 - int(depthstring)) / 10)

        try:
            assert len(En) != 0
            plt.figure("Depth: %.2f cm" % depth)
            plt.title("Primaries Energy Spectrum at %.2f cm" % depth)
            bin_heights, bin_borders, _ = plt.hist(
                En,
                round(np.sqrt(len(En))),
                facecolor="g",
                ec="black",
                alpha=0.5,
                label="histogram data",
                density=True,
            )
            # plt.title("Spettro energetico primari incidenti allo scintillatore")
            plt.ylabel("Frequency")
            plt.grid()
            plt.xlabel("Energy [Mev]")
            figname = "Spettro_%.2f.png" % depth
            #plt.savefig(os.path.join(outpath, figname))
        except:
            pass
plt.show()
