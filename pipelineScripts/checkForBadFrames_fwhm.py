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


def calculateFreedmanBins(data, initialValue = None):
    if (initialValue == None):
        bins = [min(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= max(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

def saveHistogram(values, fwhmRejectedIndices, astrometryRejectedIndices, imageName, title):
    myBins = calculateFreedmanBins(values)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    configureAxis(ax, 'FWHM (px)', '', logScale=False)
    ax.set_title(title, fontsize=22, pad=17)
    counts, bins, patches = ax.hist(values, bins=myBins, color="teal")
    max_bin_height = counts.max() + 20
    ax.set_ylim(0, max_bin_height)

    ax.text(0.6, 0.7,"Rejected Frames by fwhm: "+str(len(fwhmRejectedIndices)),transform=ax.transAxes,
                fontsize=20,verticalalignment='top',horizontalalignment='left')


    ax.hist(values[fwhmRejectedIndices], bins=myBins, color='mediumorchid', label='Rejected FWHM')
    ax.hist(values[astrometryRejectedIndices], bins=myBins, color='blue', label='Rejected astrometry')

    ax.legend(fontsize=15)
    plt.savefig(imageName)
    return()
    
def saveFWHMevol(allTable, fwhmRejectedIndices, astrometryRejectedIndices, fwhmThreshold, imageName):
    
    fig, ax = plt.subplots(2, 1, figsize=(20,10))
    configureAxis(ax[0], 'UTC', 'FWHM (pix)',logScale=False)
    configureAxis(ax[1], 'Airmass', 'FWHM (pix)',logScale=False)
    fig.suptitle('FWHM evolution',fontsize=22)
    pattern=r"(\d+).fits"

    alldate_ok = []
    allAir     = []
    for row in range(len(allTable)):
        file=allTable.loc[row]['File']
        match=re.search(pattern,file)
        frame=int(match.group(1))
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,'DATE-OBS')
        air=obtainKeyWordFromFits(file,'AIRMASS')
        date_ok=datetime.fromisoformat(date)
        fwhm=allTable.loc[row]['FWHM']
        ax[0].scatter(date_ok,fwhm,marker='o',s=50,edgecolor='black',color='teal',zorder=5)
        ax[1].scatter(air,fwhm,marker='o',s=50,edgecolor='black',color='teal',zorder=5)

        alldate_ok.append(date_ok)
        allAir.append(air)


    fwhmRejectedFiles  = [x for x in allTable['File'][fwhmRejectedIndices] ]
    fwhmRejectedValues = [x for x in allTable['FWHM'][fwhmRejectedIndices] ]

    astrometryRejectedFiles  = [x for x in allTable['File'][astrometryRejectedIndices] ]
    astrometryRejectedValues = [x for x in allTable['FWHM'][astrometryRejectedIndices] ]

    for j in range(len(fwhmRejectedFiles)):
        match=re.search(pattern, fwhmRejectedFiles[j])
        frame=int(match.group(1))
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,'DATE-OBS')
        air=obtainKeyWordFromFits(file,'AIRMASS')
        date_ok=datetime.fromisoformat(date)

        ax[0].scatter(date_ok, fwhmRejectedValues[j],marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=6,label='Rejected by FWHM' if (j==0) else "")
        ax[1].scatter(air, fwhmRejectedValues[j],marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=6,label='Rejected by FWHM'  if (j==0) else "")
        
    for j in range(len(astrometryRejectedFiles)):
        match=re.search(pattern, astrometryRejectedFiles[j])
        frame=int(match.group(1))
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,'DATE-OBS')
        air=obtainKeyWordFromFits(file,'AIRMASS')
        date_ok=datetime.fromisoformat(date)

        ax[0].scatter(date_ok, astrometryRejectedValues[j],facecolors='none', lw=1.5, edgecolor='blue',color='blue',s=350,zorder=6,label='Rejected by astrometry' if (j==0) else "")
        ax[1].scatter(air, astrometryRejectedValues[j],facecolors='none', lw=1.5, edgecolor='blue',color='blue',s=350,zorder=6,label='Rejected by astrometry'  if (j==0) else "")

    ax[0].hlines(fwhmThreshold, xmin=np.nanmin(alldate_ok), xmax=np.nanmax(alldate_ok), lw=1.5, ls="--", color="black", label="Rejection threshold")
    ax[1].hlines(fwhmThreshold, xmin=np.nanmin(allAir), xmax=np.nanmax(allAir), lw=1.5, ls="--", color="black")

    ax[0].legend(fontsize=15)
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
        raise Exception(f"File {file} does not exist.")

def getIndicesOfFiles(allData, filesNames):
    indices = []
    for i in filesNames:
        for j in range(len(allData["File"])):
            if (i == allData["File"][j]):
                indices.append(j)
    return(indices)


def identifyBadFrames(folderWithFWHM, fwhmThreshold):
    badFiles   = []
    allFiles   = []
    allFWHM     = []

    for currentFile in glob.glob(folderWithFWHM + "/range_*.txt"):
        if fnmatch.fnmatch(currentFile, '*done*.txt'):
            continue

        fwhmValue = retrieveFWHMValues(currentFile)
        if (math.isnan(fwhmValue)):
            continue
        
        allFiles.append(".".join(currentFile.split("/")[-1].split("_")[2].split(".")[:-1]))
        allFWHM.append(fwhmValue)

    allFWHM = np.array(allFWHM)

    badFiles = []
    badFWHM = []
    allTogether = []

    for i in range(len(allFWHM)):
        if (allFWHM[i] > fwhmThreshold):
            badFiles.append(allFiles[i])
            badFWHM.append(allFWHM[i])

    allTogether=pd.DataFrame({'File':allFiles,'FWHM':allFWHM})
    return(badFiles,badFWHM,allTogether)
   
def getRejectedFramesFromFile(file):
    badFrames = []

    with open(file, 'r') as f:
        for i in f:
            badFrames.append(i[:-1] + ".fits")
    return(np.array(badFrames))

HDU_TO_FIND_AIRMASS = 1

folderWithFWHM            = sys.argv[1]
outputFolder              = sys.argv[2]
outputFile                = sys.argv[3]
folderWithFramesWithAirmasses = sys.argv[4]
fwhmThreshold = float(sys.argv[5])

setMatplotlibConf()


# 0.- Get frames identified as candidates for rejection in previous step (astrometrisation, background value and background std)
rejectedAtrometryFile = outputFolder + "/identifiedBadFrames_astrometry.txt"
rejectedAstrometrisedFrames = getRejectedFramesFromFile(rejectedAtrometryFile)

# 1.- Obtain the FWHM values ------------------------
fwhmValues = np.array([])
for currentFile in glob.glob(folderWithFWHM + "/range_*.txt"):
    fwhmValue = retrieveFWHMValues(currentFile)
    if (not math.isnan(fwhmValue)):
        fwhmValues = np.concatenate((fwhmValues, [fwhmValue]))


 
# 2.- Identify what frames are outside the acceptance region -----------------------
badFilesFWHM, badFWHM, allData = identifyBadFrames(folderWithFWHM, fwhmThreshold)

with open(outputFolder + "/fwhmValues.dat", 'w') as f:
    for i in range(len(allData)):
        f.write(str(allData["File"][i]) + " " + str(allData["FWHM"][i]) + ("\n"))

# allData.to_csv(outputFolder+"/FileFWHMtable.csv")

fwhmRejectedIndices            = getIndicesOfFiles(allData, badFilesFWHM)
astrometryRejectedIndices      = getIndicesOfFiles(allData, rejectedAstrometrisedFrames)

saveFWHMevol(allData, fwhmRejectedIndices, astrometryRejectedIndices, fwhmThreshold, outputFolder+"/fwhmEvol.png")

pattern = r"\d+"
with open(outputFolder + "/" + outputFile, 'w') as file:
    for fileName in badFilesFWHM:
        match = re.search(pattern, fileName)
        result = match.group(0)
        file.write(result + '\n')


# 3.- Obtain the median and std and do teh histogram -------------------------------------
saveHistogram(allData["FWHM"], fwhmRejectedIndices, astrometryRejectedIndices, outputFolder + "/fwhmHist.png", "FWHM of frames")
