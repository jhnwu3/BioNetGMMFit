import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
from textwrap import wrap
import sys
import os 


def bins_list(min, max, step_size):
    binss= []
    current = min
    binss.append(current)
    while current < max:
        current += step_size
        binss.append(current)
    return binss

directory = "yeast/X"#"4_prot_CD3_CD8_CD28/1min_2min/Y" #"3_prot_linear_sim/X"
file = "X.csv"
path = os.path.join(directory, file)
print(path)
outPath = os.path.join(directory, file[:-4])
data = np.loadtxt(path, delimiter=',')
print(data.shape)
np.savetxt("covariance.csv", np.cov(data.T), delimiter="," )
binwidth = 2    
true = data
# fig = plt.figure(figsize=(8.0, 8.0))
fig, ax = plt.subplots(1,data.shape[1], figsize=(12,4))
# column wise graphing
for i in range(true.shape[1]):
    binwidth = (np.max(true[:,i]) - np.min(true[:,i]) )/ 50.0

    if binwidth == 0:
        binwidth = 1
    print(i)
    # print(true.shape) range(int(min(true[:,i])), int(max(true[:,i])) + 2*binwidth, binwidth)
    ax[i].hist(true[:,i],bins=bins_list(true[:,i].min(), true[:,i].max(), binwidth))
    # ax.legend(loc='upper right')
    ax[i].set_xlabel("Abundance")
    if i == 0:
        ax[i].set_ylabel("Number of Cells")
    
plt.savefig(outPath + '.png')