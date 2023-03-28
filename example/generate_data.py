import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
from textwrap import wrap
import sys
import os, shutil

# column wise
def logNormalMatrix(mu, sigma, nSamples):
    return np.exp(np.random.multivariate_normal(mu, sigma, size=nSamples))


directory = "yeast"
nSamples = 5000
sigma = np.loadtxt("covariance.csv", delimiter=',')
sigma = sigma[:5,:5]
sigma/=100 # scale things down a little bit.
# let's make the order, T, T2, D, Da, m, P  
print(sigma.max())
mu = np.array([100, 2, 0, 0, 0])
# print(logNormalMatrix(mu, sigma,10) - logNormalMatrix(mu, sigma,10))

# shutil.rmtree(directory)
# os.mkdir(directory)
# os.mkdir(os.path.join(directory, "X"))
# os.mkdir(os.path.join(directory, "Y"))

pathX = os.path.join(directory, "X", "X.csv")
pathY = os.path.join(directory, "Y", "Y.csv")
X = logNormalMatrix(mu, sigma, nSamples)
Y = logNormalMatrix(mu, sigma, nSamples)
np.savetxt(pathX, X)
np.savetxt(pathY, Y)


