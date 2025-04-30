import sys
import astropy

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.stats import sigma_clipped_stats
from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style

def setMatplotlibConf():
    rc_fonts = {
        "font.family": "serif",
        "font.size": 14,
        "font.weight" : "medium",
        # "text.usetex": True,  # laggs a little when generatin plots in my fedora36
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 8.0,
        "xtick.major.width": 2.8,
        "xtick.minor.size": 4.0,
        "xtick.minor.width": 2.5,
        "ytick.major.size": 8.0,
        "ytick.major.width": 1.8,
        "ytick.minor.size": 4.0,
        "ytick.minor.width": 1.8,
        "legend.handlelength": 3.0,
        "axes.linewidth" : 3.5,
        "xtick.major.pad" : 6,
        "ytick.major.pad" : 6,
        "legend.fancybox" : True,
        "mathtext.fontset" : "dejavuserif"
    }
    mpl.rcParams.update(rc_fonts)
    return(rc_fonts)


def configureAxis(ax, xlabel, ylabel, logScale=True):
    ax.xaxis.set_minor_locator(MultipleLocator(1000000))
    ax.yaxis.set_minor_locator(MultipleLocator(1000000))
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis='x', which='major', labelsize=25, pad=17)
    ax.tick_params(axis='y', which='major', labelsize=25, pad=17)
    ax.set_xlabel(xlabel, fontsize=30, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=30, labelpad=10)
    if(logScale): ax.set_yscale('log')

def readMetric(fileName):
    fileNames = []
    data = []

    with open(fileName, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    file  = parts[0]
                    value = float(parts[1])
                    fileNames.append(file)
                    data.append(value)
                except ValueError:
                    continue   
    return(np.array(fileNames), np.array(data))

def calculateFreedmanBins(data, initialValue = None):
    if (initialValue == None):
        bins = [np.nanmin(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= np.nanmax(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

setMatplotlibConf()

folderWithResiduals = sys.argv[1]
diagnosisFolder     = sys.argv[2]

fileWithMetric = folderWithResiduals + "/pixelsAdded.txt"
fileNames, values = readMetric(fileWithMetric)

# Detection of outliers (> 3sigma)
meanOfDistribution = np.nanmean(values)
stdOfDistribution  = np.nanstd(values)

outliers1 = np.array(np.where(values < (meanOfDistribution - 3*stdOfDistribution)))
outliers2 = np.array(np.where(values > (meanOfDistribution + 3*stdOfDistribution)))
outliersFileNames = np.concatenate((fileNames[outliers1[0]], fileNames[outliers2[0]]))
outliersValues    = np.concatenate((values[outliers1[0]], values[outliers2[0]]))

with open(diagnosisFolder + "/residualMetricOutliers.txt", 'w') as f:
    for i in range(len(outliersFileNames)):
        f.write(str(outliersFileNames[i]) + " " + str(outliersValues[i]) + "\n")

# Distribution plot
binsToUse = calculateFreedmanBins(values)
fig, ax = plt.subplots(1, 1, figsize=(15, 15))
configureAxis(ax, 'Metric (sum of px)', '', logScale=False)
n, bins, patches = ax.hist(values, bins=binsToUse, color="teal")
ax.hist(outliersValues, bins=binsToUse, color="crimson")
ax.vlines(x=meanOfDistribution + 3*stdOfDistribution, ymin=0, ymax=n.max(), color="black", lw=2, ls="--")
ax.vlines(x=meanOfDistribution - 3*stdOfDistribution, ymin=0, ymax=n.max(), color="black", lw=2, ls="--")

plt.savefig(diagnosisFolder + "/residualsMetricDistribution.png")
