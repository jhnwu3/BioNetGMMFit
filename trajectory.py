import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
from textwrap import wrap
import sys

file = "trajectoryData.csv" # idk we'll change this later.
colors = ['b', 'g', 'g', 'r', 'c', 'm', 'y', 'k']
mpl.rcParams['figure.dpi'] = 500
df = pd.read_csv(file)
data = df.to_numpy()
plt.xlabel("Time (unit)", fontdict = {'fontsize' : 12})
plt.ylabel("Abundance/Concentration", fontdict = {'fontsize' : 12})
totalCols = data.shape[1]
finalGroundTruthCol = int( 1 + ((data.shape[1] - 1) / 2))
nGTruthCols = finalGroundTruthCol - 1
print(data)
print(data.shape)
for i in range(1,finalGroundTruthCol):
    plt.plot(data[:,0], data[:,i], colors[i - 1], label="S" + str(i))

for i in range(finalGroundTruthCol, totalCols):
    labelIndex = i - nGTruthCols
    colorIndex = i - finalGroundTruthCol
    plt.plot(data[:,0], data[:,i],colors[colorIndex] +  '--' , label="Estimated S" + str(labelIndex))

plt.savefig("TrajectoryOverlapped.png")