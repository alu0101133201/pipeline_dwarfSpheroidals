import os
import re
import sys
import glob
import math
import fnmatch
import astropy
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip
from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style

from datetime import datetime
import time

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

def computeMedianAndStd(values):
    mean, median, std = sigma_clipped_stats(values)
    return(median, std)

def calculateFreedmanBins(data, initialValue = None):
    if (initialValue == None):
        bins = [min(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= max(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

def saveHistogram(values, median, std, badValues, imageName, numOfStd, title, labelX):
    valuesToPlot = values[~np.isnan(values)]
    myBins = calculateFreedmanBins(valuesToPlot)
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.tight_layout(pad=7)
    configureAxis(ax, labelX, 'Number of frames', logScale=False)
    ax.set_title(title, fontsize=18, pad=17)

    counts, bins, patches = ax.hist(valuesToPlot, color="teal", bins=myBins)
    max_bin_height = counts.max() + 5
    ax.set_ylim(0, max_bin_height)

    ax.text(0.3755, 0.95, "Median: " + str(int(median)), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')
    ax.text(0.375, 0.9, "Std: " + str(int(std)), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')

    ax.hist(badValues, bins=myBins,color='darkred')

    ax.text(0.375, 0.85, "Rejected: " + str(len(badValues)), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')

    plt.savefig(imageName)
    return()


def obtainKeyWordFromFits(file, keyword):
    if os.path.exists(file):
        with fits.open(file) as hdul:
            header = hdul[HDU_TO_FIND_AIRMASS].header
            
            if keyword in header:
                keywordValue = header[keyword]
                return(keywordValue)
            else:
                raise Exception(f"Keyword '{keyword}' not found in the header.")
    else:
        raise Exception(f"File {file} does not exist.")


def obtainBackgroundStd(currentFile):
    backgroundValue = -1

    # First we read the background Value
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) != 1):
            raise Exception("File with the background estimation contains more that 1 line. Expected 1 line got " + str(len(lines)))
        
        splittedLine = lines[0].strip().split()
        numberOfFields = len(splittedLine)

        if (numberOfFields == 5):
            backgroundStd   = float(splittedLine[2])
        elif (numberOfFields == 1):
            return(float('nan'), float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 5 (constant estimation of the background), got " + str(numberOfFields))

    return(backgroundStd)


def identifyBadFrames(folderWithFrames, numberOfStdForRejecting):
    badFiles   = []
    allFiles   = []
    allStd     = []
    allBackgroundValues = []

    for currentFile in glob.glob(folderWithFrames + "/*.txt"):
        if fnmatch.fnmatch(currentFile, '*done*.txt'):
            continue
        currentStd = obtainBackgroundStd(currentFile)

        if (math.isnan(currentStd)):
            continue
        allFiles.append(currentFile)
        allStd.append(currentStd)

    allStd = np.array(allStd)
    std_mask = sigma_clip(allStd, sigma=numberOfStdForRejecting, cenfunc='median', stdfunc='std', maxiters=5, masked=True).mask

    allFiles = np.array(allFiles)
    badFiles = allFiles[std_mask]
    badStd   = allStd[std_mask]
    return(allFiles, allStd, badFiles, badStd)


HDU_TO_FIND_AIRMASS = 1

folderWithSkyEstimations      = sys.argv[1]
outputFolder                  = sys.argv[2]
outputFile                    = sys.argv[3]
numberOfStdForRejecting       = int(sys.argv[4])
nightNumber                   = int(sys.argv[5])

setMatplotlibConf()

allFiles, backgroundStds, badFiles, badStd = identifyBadFrames(folderWithSkyEstimations, numberOfStdForRejecting)

with open(outputFolder + "/" + outputFile, 'w') as file:
    for fileName in badFiles:
        file.write(fileName + '\n')

backgroundStdMedian, BackgroundStdStd = computeMedianAndStd(backgroundStds)
saveHistogram(backgroundStds, backgroundStdMedian, BackgroundStdStd, badStd,\
                outputFolder + f"/backgroundStdHist_beforeFlat_n{nightNumber}.png", numberOfStdForRejecting, "Excluding frames in flat construction", "Background STD (ADU)")