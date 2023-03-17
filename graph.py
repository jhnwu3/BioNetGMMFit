import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
from textwrap import wrap
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
def compute_moments(x):
    means = np.average(x,axis=0)
    variances = np.var(x, axis=0)
    # n means/variances, (n^2 -n)/2 covariances, 
    cov = []
    for i in range(x.shape[1]):
        for j in range(i+1, x.shape[1]):
           cov.append(np.cov(x[:,i],x[:,j])[0][1])
    return np.append(means, np.append(variances,cov))

def getFile(args, multi=False):
    if multi:
        return args[args.index('-f') + 1], args[args.index('-f') + 2]
    return args[args.index('-f') + 1]

def getName(args):
    if '-n' in args:
        return args[args.index('-n') + 1]
    return ""

def getMomentSubset(args):
    if '-m' in args:
        return args[args.index('-m') + 1]
    return ""

def getTime(args):
    return args[args.index('-t') + 1]

def getGraphType(args):
    return args[args.index('-g') + 1]

def getSimulatedRates(args):
    return args[args.index('-r') + 1]

def getSimulatedRates(args):
    return args[args.index('-r') + 1]

class Graph:
    
    def plot_confidence_interval(x, values, axes, z=1.96, color='#2187bb', horizontal_line_width=0.25):
        mean = np.mean(values)
        stdev = np.std(values)
        confidence_interval = z * stdev / np.sqrt(len(values))

        left = x - horizontal_line_width / 2
        top = mean - confidence_interval
        right = x + horizontal_line_width / 2
        bottom = mean + confidence_interval
        axes.plot([x, x], [top, bottom], color=color)
        axes.plot([left, right], [top, top], color=color)
        axes.plot([left, right], [bottom, bottom], color=color)
        axes.plot(x, mean, 'o', color='#f44336')

        return mean, confidence_interval
    
    # def __init__(self, dir):
    #     self.graphingDirectory = "/frontend/graph/" + dir +"/"
    #     self.fileName = dir
    #     # self.name = Graph_Writer.getFileName(argv)
    #     # self.t = Graph_Writer.getTime(argv)
    #     # self.estimates = np.genfromtxt(Graph_Writer.graphingDirectory + self.name + '_estimates.csv', delimiter=',')
    #     # print(self.estimates.shape)
    #     # self.moments = np.genfromtxt()
        
        
    def plotConfidenceIntervals(z, file, simulated = False, trueRatesFile='', title=""):
        df = pd.read_csv(file)
        estimates = df.to_numpy()
        # estimates = np.genfromtxt(path, delimiter=',')
        nRates = estimates.shape[1] - 1
        categoriesNumeric = []
        for i in range(estimates.shape[1]): 
            categoriesNumeric.append(i+1)
        fig, axes = plt.subplots(figsize=(5.0, 5.0))
        axes.spines.right.set_visible(False)
        axes.spines.top.set_visible(False)
        axes.set_xticks(categoriesNumeric, df.columns)
        axes.tick_params(axis='both', which='major', labelsize=12)
        axes.set_title(title,wrap=True, fontdict = {'fontsize' : 18})
        axes.set_ylabel("Estimate Value", fontdict = {'fontsize' : 12})
        groundTruth = []
        for i in range(nRates):
            Graph.plot_confidence_interval(i + 1, estimates[:,i], axes, z)
        if simulated:
            tRates = np.genfromtxt(trueRatesFile, delimiter=',')
            # print(tRates)
            for i in range(nRates):
                groundTruth.append(axes.plot(i+1,tRates[i], 'D', color='#013220'))
        # plt.style.use('ggplot')
        # plt.show()
        if simulated:
            axes.legend(groundTruth[0],["Ground Truth"])
        plt.savefig(file[:-4] + '_estimates.png')
        
    def plotMomentsWithActualEvolvedMatrices(xName, yName, gTitle=""): # this might make more sense overall actually. From here, we can get means, variances, and covariances.
        dfX = pd.read_csv(xName)
        dfY = pd.read_csv(yName)
        x = dfX.to_numpy()
        y = dfY.to_numpy()
        plt.title(xName[:-4] + " and " + yName[:-4])
        plt.xlabel("Predicted Moment")
        plt.ylabel("Observed Moment")
        xMoments = compute_moments(x)
        yMoments = compute_moments(y)
        plt.plot(np.unique(xMoments), np.poly1d(np.polyfit(xMoments, yMoments, 1))(np.unique(xMoments)))
        x123 = np.arange(0, np.max(xMoments))
        y123 = x123
        plt.plot(np.unique(x123), np.poly1d(np.polyfit(x123, y123, 1))(np.unique(x123)), color='#f44336')
        plt.scatter(xMoments, yMoments)
        plt.savefig(xName[:-7] + '.png')
    
    # assume data format n moments X 2 columns (for X and Y) 
    def plotMoments(file, title="", nSpecies=""): # get list of X, Y 
        df = pd.read_csv(file)
        moments = df.to_numpy()
        if title == "":
            title = file

        fig, axes = plt.subplots(figsize=(6.5, 6.0))
        axes.set_title(title, wrap=True,loc='center', fontdict = {'fontsize' : 20})    
        plt.xlabel("Predicted Moment", fontdict = {'fontsize' : 12})
        plt.ylabel("Observed Moment", fontdict = {'fontsize' : 12})
        axes.spines.right.set_visible(False)
        axes.spines.top.set_visible(False)
    
        x123 = np.arange(np.min(moments[:]), np.max(moments[:]))
        y123 = x123
        optimalLine, = axes.plot(np.unique(x123), np.poly1d(np.polyfit(x123, y123, 1))(np.unique(x123)),'--')
        a, b = np.polyfit(moments[:,0], moments[:,1], 1)
        bestFit, = axes.plot(np.unique(x123), np.poly1d(np.polyfit(moments[:,0], moments[:,1], 1))(np.unique(x123)), ':')
        
        totMoments = moments.shape[0]
        means = ""
        variances = ""
        covariances = ""
        if(nSpecies != ""):
            nSpecies = int(nSpecies)
            if(nSpecies*(nSpecies + 3) / 2 == totMoments):
                means = axes.scatter(moments[:nSpecies,0],moments[:nSpecies,1], s=80)
                variances = axes.scatter(moments[nSpecies:2*nSpecies,0], moments[nSpecies:2*nSpecies,1], s=80)
                covariances = axes.scatter(moments[2*nSpecies:totMoments,0], moments[2*nSpecies:totMoments,1], s=80)
                axes.legend([optimalLine, bestFit, means, variances, covariances], ["Perfect Fit",  "Best Fit of Observed vs. Predicted Line:" + " y= " + "{:.2f}".format(a) + "x + " + "{:.2f}".format(b), "Means", "Variances", "Covariances"])
            elif(2*nSpecies == totMoments):
                means = axes.scatter(moments[:nSpecies,0],moments[:nSpecies,1], s=80)
                variances = axes.scatter(moments[nSpecies:2*nSpecies,0], moments[nSpecies:2*nSpecies,1], s=80)
                axes.legend([optimalLine, bestFit, means, variances], ["Perfect Fit",  "Best Fit of Observed vs. Predicted Line:" + " y= " + "{:.2f}".format(a) + "x + " + "{:.2f}".format(b), "Means", "Variances"])      
            else:
                means = axes.scatter(moments[:,0],moments[:,1], s=80) 
                axes.legend([optimalLine, bestFit, means], ["Perfect Fit",  "Best Fit of Observed vs. Predicted Line:" + " y= " + "{:.2f}".format(a) + "x + " + "{:.2f}".format(b), "Means"])      
        else:
            axes.scatter(moments[:,0],moments[:,1] , s=80)   
            axes.legend([optimalLine, bestFit], ["Perfect Fit",  "Best Fit of Observed vs. Predicted Line:" + " y= " + "{:.2f}".format(a) + "x + " + "{:.2f}".format(b)])    
        
        plt.savefig(file[:-4] + '.png')
        
    def plotAllMoments(xFile, yFile, z=1.96, title=""): 
        # get data
        dfX = pd.read_csv(xFile)
       
        xMoments = dfX.to_numpy()
        yMoments = np.genfromtxt(yFile, delimiter=',')
        if title == "":
            title = xFile
        # compute confidence intervals 
        xAvgs = np.mean(xMoments, axis=0)
        xStds = np.std(xMoments, axis=0)
        nRuns = xMoments.shape[0]
        xErrorBars = z*xStds/ (np.sqrt(nRuns)) # 95% CI's
        print(xErrorBars)
        # plt.figure(frameon=False)
        fig, axes = plt.subplots(figsize=(6.5, 6.0))
        axes.set_title(title, wrap=True,loc='center', fontdict = {'fontsize' : 20})   
        plt.xlabel("Estimated Moment", fontdict = {'fontsize' : 12})
        plt.ylabel("Observed Moment", fontdict = {'fontsize' : 12})
        axes.spines.right.set_visible(False)
        axes.spines.top.set_visible(False)
        axes.scatter(xAvgs,yMoments)
       
        x123 = np.arange(np.min(xAvgs), np.max(xAvgs))
        y123 = x123
        optimalLine, = axes.plot(np.unique(x123), np.poly1d(np.polyfit(x123, y123, 1))(np.unique(x123)), color='red') # perfect fit
        bestFit, = axes.plot(np.unique(xAvgs), np.poly1d(np.polyfit(xAvgs, yMoments, 1))(np.unique(xAvgs)))
        axes.errorbar(xAvgs, yMoments, yerr = 10, fmt='-', elinewidth = 5, capsize=10)
        axes.legend([optimalLine, bestFit], ["Perfect Fit",  "Best Fit of Observed vs. Predicted Line"])
        plt.savefig(xFile[:-4] + 'intervalMomentsEstimated.png')
        
if __name__ == "__main__": 
    # print(sys.argv)
    if '-h' in sys.argv:
        print("Specify graphing type with -g <graph type> ")
        print("Specify files with -f <filename> <optional 2nd final name only if graph type is 'dMoments'>")
        exit(0)
            
    if "-f" not in sys.argv:
        print("Error Need to Specify Graph Files with -f")
        exit(0)
        
    if "-g" not in sys.argv:
        print("Error Need to Specify Graph Types with -g")
        exit(0)

    plt.rcParams['figure.constrained_layout.use'] = True
    plt.tight_layout()
    graphType = getGraphType(sys.argv)
    if graphType == 'CI':
        Graph.plotConfidenceIntervals(1.96, getFile(sys.argv), title=getName(sys.argv))
    elif graphType == 'CI_truth':
        Graph.plotConfidenceIntervals(1.96, getFile(sys.argv),simulated=True, trueRatesFile=getSimulatedRates(sys.argv), title=getName(sys.argv))
    elif graphType == 'Moments':
        Graph.plotMoments(getFile(sys.argv), title=getName(sys.argv), nSpecies= getMomentSubset(sys.argv))
    elif graphType == 'dMoments':
        dataX, dataY = getFile(sys.argv, multi=True)
        Graph.plotMomentsWithActualEvolvedMatrices(dataX,dataY)
    elif graphType =='allMoments':
        dataX, dataY = getFile(sys.argv, multi=True)
        Graph.plotAllMoments(dataX,dataY,title=getName(sys.argv))
    else:
        print("Error Invalid Graph Type Inputted:", graphType)
        print("")