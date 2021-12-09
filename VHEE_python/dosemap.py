import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.insert(0, "E_FLASH_RADIOTHERAPY")
from flash_helper.flash_functions import PDD_plotter_out, mean_array_calculator,find_nearest,find_r


def dosemap(path,depth,edge,len_,measurement_point,evt,axis=0):
    df=pd.read_csv(
            path, names=["X", "Y", "Z", "dose [Gy]", "dosesq [Gy^2]", "entry"])
    bins = len(np.array(df["dose [Gy]"][np.logical_and(df["Z"] == 0, df["Y"] == 0)]))

    length = depth / bins

    distance = np.array([length * x for x in range(1, bins + 1)])
    x=np.linspace(- len_,len_,edge)
    y=x
    z=[]
    X=[]
    Y=[]
    x1,y1=np.meshgrid(x,y)
    if axis==0:
        for i in range(0,len(x)):
            for j in range(0,len(y)):
                X.append(i)
                Y.append(j)
                z.append(float(df["dose [Gy]"][np.logical_and(np.logical_and(df["Z"] == i, df["Y"] == j),df["X"]==find_nearest(distance,measurement_point))])/evt)

        z=np.array(z)
        z=z.reshape((edge,edge))
        return x1,y1,z
if __name__ == "__main__":
    path_160 = "C:/Users/JakeHarold/Downloads/pencil160_spectrum_good.csv"
    depth_=200 #mm
    edge_=50 #bins
    lenght_=25 #mm
    measurement_points= [0,20,60,100,150,200] #mm
    evt_=1e6

    rows = 2
    cols = 3
    axes=[]
    fig=plt.figure()

    for i in range(rows*cols):
        mesh_x,mesh_y,dose=dosemap(path_160,depth_,edge_,lenght_,measurement_points[i],evt_)
        axes.append( fig.add_subplot(rows, cols, i+1) )
        subplot_title=(f"Dose distribution xy plane at {measurement_points[i]/10} cm")
        axes[-1].set_title(subplot_title)
        plt.imshow(dose,interpolation='none',cmap="jet",extent=[-lenght_,lenght_,-lenght_,lenght_])
        plt.colorbar(orientation='horizontal',pad = 0.3).set_label("Gy/e")

        plt.xlabel("x [mm]")
        plt.ylabel("y [mm]")
    fig.tight_layout()
    plt.show()

