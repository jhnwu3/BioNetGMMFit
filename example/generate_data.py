import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
from textwrap import wrap
import sys
import os, shutil

# column wise
# let mu be the means we want
# sigma, just the sigma value (related to standard deviation)
def logNormalMatrix(mean, sigma, nSamples):
    
    # m = e^(mu + std^2)
    # mu = ln m - std^2
    mu = np.log(mean) - 0.5* (sigma * sigma)
    print(mu)
    return np.random.lognormal(mu,sigma,(nSamples, mu.shape[0])) #np.exp(np.random.multivariate_normal(mu, sigma, size=nSamples))


directory = "yeast"
nSamples = 10000
sigma = np.loadtxt("covariance.csv", delimiter=',')
sigma = sigma[:5,:5]
sigma = sigma / 1000 # scale things down a little bit.
# print(sigma)
sigma = 1 * np.ones(5)
# let's make the order, T, T2, D, Da, m, P  
print(sigma.max())
mu = np.array([100, 2, 0.000001, 0.000001, 0.0000001])


# print(logNormalMatrix(mu, sigma,10) - logNormalMatrix(mu, sigma,10))

# shutil.rmtree(directory)
# os.mkdir(directory)
# os.mkdir(os.path.join(directory, "X"))
# os.mkdir(os.path.join(directory, "Y"))

pathX = os.path.join(directory, "X", "X.csv")
pathY = os.path.join(directory, "Y", "Y.csv")
X = logNormalMatrix(mu, sigma, nSamples)
Y = logNormalMatrix(mu, sigma, nSamples)

X[:,2:] = 0
Y[:,2:] = 0
np.savetxt(pathX, X, delimiter=',')
np.savetxt(pathY, Y, delimiter=',')


