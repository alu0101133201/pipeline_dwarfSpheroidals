import sys
import astropy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.io import fits
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

def calculateFreedmanBins(data, initialValue = None):
    if (initialValue == None):
        bins = [min(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= max(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

def extractFactorsFromFile(file):
    starsUsed = []
    try:
        with open(file, 'r') as f:
            for line in f:
                splittedLine = line.split()
                if (len(splittedLine) == 4): # If the file has been correctly created it has 4 fields
                    starsUsed.append(float(splittedLine[-1]))
    except (OSError, ValueError) as e:
        print(f"Error reading file {file}: {e}")
    return(np.array(starsUsed))

def saveHistogram(starsUsed, myBins, outputFileName):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_title("Number of stars used in each frame for calibration")
    configureAxis(ax, 'Number of Stars used', '', logScale=False)
    ax.hist(starsUsed, bins=myBins, color="teal")
    plt.savefig(outputFileName)


fileWithFactors = sys.argv[1]
outputFileName  = sys.argv[2]

setMatplotlibConf()

starsUsed = extractFactorsFromFile(fileWithFactors)
myBins = calculateFreedmanBins(starsUsed)
saveHistogram(starsUsed, myBins, outputFileName)