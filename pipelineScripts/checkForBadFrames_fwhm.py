import re
import sys
import glob
import math
import fnmatch
import astropy

import os
from astropy.io import fits

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from matplotlib.ticker import MultipleLocator

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



def retrieveFWHMValues(currentFile):
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) != 1):
            raise Exception("File " + currentFile + " with the FWHM estimation contains more that 1 line. Expected 1 line got " + str(len(lines)))
        
        splittedLine = lines[0].strip().split()
        numberOfFields = len(splittedLine)

        if (numberOfFields == 4):
            return(float(splittedLine[0]))
        elif (numberOfFields == 0):
            return(float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 4 (constant estimation of the background), got " + str(numberOfFields))

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

def saveHistogram(values, median, std, badValues, imageName, numOfStd, title):
    myBins = calculateFreedmanBins(values)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    configureAxis(ax, 'FWHM (px)', '', logScale=False)
    ax.set_title(title, fontsize=22, pad=17)
    counts, bins, patches = ax.hist(values, bins=myBins, color="teal")
    max_bin_height = counts.max() + 10
    ax.set_ylim(0, max_bin_height)

    ax.text(0.3755, 0.95, "Median: " + "{:.2f}".format(median), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')
    ax.text(0.375, 0.9, "Std: " + "{:.2f}".format(std), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')
    if len(badValues)!=0:
        counts_bad, bins_bad, patches_bad = ax.hist(badValues,bins=np.linspace(np.nanmin(badValues), np.nanmax(badValues), 4),color='red',label='Rejected FWHM')
        
        ax.legend()
    ax.text(0.3655,0.85,"Rej. Frames: "+str(len(badValues)),transform=ax.transAxes,
                fontsize=20,verticalalignment='top',horizontalalignment='left')
    plt.savefig(imageName)
    return()
    
def saveFWHMevol(allTable,badFiles,badFwhm,imageName):
    fig, ax = plt.subplots(1, 2, figsize=(10,10))
    configureAxis(ax[0], 'UTC', 'FWHM (pix)',logScale=False)
    configureAxis(ax[1], 'Airmass', 'FWHM (pix)',logScale=False)
    fig.suptitle('FWHM evolution',fontsize=22,pad=17)
    pattern=r"entirecamera_(\d+)"
    for row in range(len(allTable)):
        file=allTable.loc[row]['File']
        match=re.search(pattern,file)
        frame=float(match.group(1))
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,'DATE-OBS')
        air=obtainKeyWordFromFits(file,'AIRMASS')
        date_ok=datetime.fromisoformat(date)
        fwhm=allTable.loc[row]['FWHM']
        ax[0].scatter(date_ok,fwhm,marker='o',s=50,edgecolor='black',color='teal',zorder=5)
        ax[1].scatter(air,fwhm,marker='o',s=50,edgecolor='black',color='teal',zorder=5)

    
    if len(badFiles)!=0:
        for j in range(len(badFiles)):
            match=re.search(pattern,badFiles[j])
            frame=float(match.group(1))
            file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
            date=obtainKeyWordFromFits(file,'DATE-OBS')
            air=obtainKeyWordFromFits(file,'AIRMASS')
            date_ok=datetime.fromisoformat(date)
            ax[0].scatter(date_ok,badFwhm[j],marker='P',edgecolor='k',color='mediumorchid',s=80,zorder=6,label='Rejected FWHM')
            ax[1].scatter(air,badFwhm[j],marker='P',edgecolor='k',color='mediumorchid',s=80,zorder=6,label='Rejected FWHM')
            
        ax[0].legend()
    for label in ax[0].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')
    plt.tight_layout()
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
        raise Exception(f"File {fits_file_path} does not exist.")

HDU_TO_FIND_AIRMASS = 1

folderWithFWHM            = sys.argv[1]
outputFolder              = sys.argv[2]
outputFile                = sys.argv[3]
numberOfStdForRejecting    = float(sys.argv[4])
folderWithFramesWithAirmasses = sys.argv[5]

setMatplotlibConf()

# 1.- Obtain the FWHM values ------------------------
fwhmValues = np.array([])
for currentFile in glob.glob(folderWithFWHM + "/range_*.txt"):
    fwhmValue = retrieveFWHMValues(currentFile)
    if (not math.isnan(fwhmValue)):
        fwhmValues = np.concatenate((fwhmValues, [fwhmValue]))



def identifyBadFrames(folderWithFWHM, numberOfStdForRejecting):
    badFiles   = []
    allFiles   = []
    allFWHM     = []

    for currentFile in glob.glob(folderWithFWHM + "/range_*.txt"):
        if fnmatch.fnmatch(currentFile, '*done*.txt'):
            continue

        fwhmValue = retrieveFWHMValues(currentFile)
        if (math.isnan(fwhmValue)):
            continue
        allFiles.append(currentFile)
        allFWHM.append(fwhmValue)

    allFWHM = np.array(allFWHM)

    mask = sigma_clip(allFWHM, sigma=numberOfStdForRejecting, cenfunc='median', stdfunc='std', maxiters=5, masked=True).mask

    allFiles = np.array(allFiles)
    allTogether=pd.DataFrame({'File':allFiles,'FWHM':allFWHM})
    badFiles = allFiles[mask]
    badFWHM = allFWHM[mask]
    return(badFiles,badFWHM,allTogether)
    
# 2.- Identify what frames are outside the acceptance region -----------------------
badFiles,badFWHM,allData = identifyBadFrames(folderWithFWHM, numberOfStdForRejecting)
allData.to_csv(outputFolder+"/FileFWHMtable.csv")
saveFWHMevol(allData,badFiles,badFWHM,outputFolder+"/fwhmEvol.png")
pattern = r"entirecamera_\d+"
with open(outputFolder + "/" + outputFile, 'w') as file:
    for fileName in badFiles:
        match = re.search(pattern, fileName)
        result = match.group()
        file.write(result + '\n')

# 3.- Obtain the median and std and do teh histogram -------------------------------------
fwhmValueMean, fwhmValueStd = computeMedianAndStd(fwhmValues)
saveHistogram(fwhmValues, fwhmValueMean, fwhmValueStd, badFWHM, outputFolder + "/fwhmHist.png", numberOfStdForRejecting, "FWHM of frames")
