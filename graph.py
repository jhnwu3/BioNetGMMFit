import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import sys

# def plot_confidence_interval(x, values, z=1.96, color='#2187bb', horizontal_line_width=0.25):
#     mean = np.mean(x)
#     stdev = np.std(x)
#     confidence_interval = z * stdev / np.sqrt(len(values))

#     left = x - horizontal_line_width / 2
#     top = mean - confidence_interval
#     right = x + horizontal_line_width / 2
#     bottom = mean + confidence_interval
#     plt.plot([x, x], [top, bottom], color=color)
#     plt.plot([left, right], [top, top], color=color)
#     plt.plot([left, right], [bottom, bottom], color=color)
#     plt.plot(x, mean, 'o', color='#f44336')

#     return mean, confidence_interval

def getFile(args):
    return args[args.index('-f') + 1]

def getName(args):
    return args[args.index('-n') + 1]

def getTime(args):
    return args[args.index('-t') + 1]

def getGraphType(args):
    return args[args.index('-g') + 1]

def getSimulatedRates(args):
    return args[args.index('-g') + 1]

class Graph:
    
    def plot_confidence_interval(x, values, z=1.96, color='#2187bb', horizontal_line_width=0.25):
        mean = np.mean(values)
        stdev = np.std(values)
        print(stdev)
        confidence_interval = z * stdev / np.sqrt(len(values))

        left = x - horizontal_line_width / 2
        top = mean - confidence_interval
        right = x + horizontal_line_width / 2
        bottom = mean + confidence_interval
        plt.plot([x, x], [top, bottom], color=color)
        plt.plot([left, right], [top, top], color=color)
        plt.plot([left, right], [bottom, bottom], color=color)
        plt.plot(x, mean, 'o', color='#f44336')

        return mean, confidence_interval
    
    # def __init__(self, dir):
    #     self.graphingDirectory = "/frontend/graph/" + dir +"/"
    #     self.fileName = dir
    #     # self.name = Graph_Writer.getFileName(argv)
    #     # self.t = Graph_Writer.getTime(argv)
    #     # self.estimates = np.genfromtxt(Graph_Writer.graphingDirectory + self.name + '_estimates.csv', delimiter=',')
    #     # print(self.estimates.shape)
    #     # self.moments = np.genfromtxt()
        
        
    def plotConfidenceIntervals(z, file, simulated = False):
        df = pd.read_csv(file)
        estimates = df.to_numpy()
        # estimates = np.genfromtxt(path, delimiter=',')
        nRates = estimates.shape[1] - 1
        categoriesNumeric = []
        for i in range(): 
            categoriesNumeric.append(i+1)
            
        plt.xticks(categoriesNumeric, df.columns)
        plt.title('Confidence Intervals For t=1, 6 cells, true thetas 0.24, 0.81')
        for i in range(nRates):
            Graph.plot_confidence_interval(i + 1, estimates[:,i], z)
        if simulated:
            plt.plot(1, 0.24, 'D', color='#013220')
            plt.plot(2, 0.81, 'D', color='#013220')
        plt.show()
        plt.savefig(file + '_estimates.png')
        
    def plotMomentsWithActualEvolvedMatrices(self, name, t): # this might make more sense overall actually. From here, we can get means, variances, and covariances.
        print("workinprogress")
        return 0
    
    # assume data format n moments X 2 columns (for X and Y) 
    def plotMoments(file): # get list of X, Y 
        df = pd.read_csv(file)
        moments = df.to_numpy()
        plt.title(file[:-4])
        plt.xlabel(df.columns[0])
        plt.ylabel(df.columns[1])
        plt.plot(moments[:,0],moments[:,1])
        plt.savefig(file[:-4] + '.png')
        
if "-f" not in sys.argv:
    print("Error Need to Specify Graph Files with -f")
    exit(0)
    
if "-g" not in sys.argv:
    print("Error Need to Specify Graph Files with -g")
    exit(0)

# graph.plotConfidenceIntervals(1.96)
if getGraphType(sys.argv) == 'CI':
    Graph.plotConfidenceIntervals(getFile(sys.argv))
else:
    Graph.plotMoments(getFile(sys.argv))