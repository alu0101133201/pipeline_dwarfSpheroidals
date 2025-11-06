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
from astropy.time import Time

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



def retrieveFWHMValues(currentFile,h,arcsecPerPix):
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) < h):
            raise Exception("File " + currentFile + " with the FWHM estimation contains less lines than number of extension. Expected at least" + str(h))
        
        splittedLine = lines[(h-1)].strip().split()
        numberOfFields = len(splittedLine)

        if (numberOfFields == 1):
            return(float(splittedLine[0])*arcsecPerPix)
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

def saveHistogram(allData, median, std, fwhmRejectedIndices, astrometryRejectedIndices,  imageName, numOfStd, title):
    values = allData["FWHM"]

    myBins = calculateFreedmanBins(values)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    configureAxis(ax, 'FWHM (arcsec)', '', logScale=False)
    ax.set_title(title, fontsize=22, pad=17)
    counts, bins, patches = ax.hist(values, bins=myBins, color="teal")
    max_bin_height = counts.max() + 20
    ax.set_ylim(0, max_bin_height)

    ax.text(0.6, 0.8, "Median: " + "{:.2f}".format(median), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')
    ax.text(0.6, 0.75, "Std: " + "{:.2f}".format(std), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')
    ax.text(0.6, 0.7,"Rejected Frames: "+str(len(fwhmRejectedIndices)),transform=ax.transAxes,
                fontsize=20,verticalalignment='top',horizontalalignment='left')


    ax.hist(values[fwhmRejectedIndices], bins=myBins, color='mediumorchid', label='Rejected FWHM')
    ax.hist(values[astrometryRejectedIndices], bins=myBins, color='blue', label='Rejected astrometry')

    ax.legend(fontsize=15)
    plt.savefig(imageName)
    return()
    
def saveFWHMevol(allTable, fwhmRejectedIndices, astrometryRejectedIndices,  imageName,airKey,dateKey,fwhmToReject):
    
    fig, ax = plt.subplots(2, 1, figsize=(20,10))
    configureAxis(ax[0], 'UTC', 'FWHM (arcsec)',logScale=False)
    configureAxis(ax[1], 'Airmass', 'FWHM (arcsec)',logScale=False)
    fig.suptitle('FWHM evolution',fontsize=22)
    pattern=r"(\d+).fits"

    for row in range(len(allTable)):
        file=allTable.loc[row]['File']
        match=re.search(pattern,file)
        frame=int(match.group(1))
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,dateKey)
        air=obtainKeyWordFromFits(file,airKey)
        if dateHeaderKey.startswith("DATE"):
            date_ok=datetime.fromisoformat(date)
        elif dateHeaderKey.startswith("MJD"):
            date_ok=Time(date,format='mjd').to_datetime()
        fwhm=allTable.loc[row]['FWHM']
        ax[0].scatter(date_ok,fwhm,marker='o',s=50,edgecolor='black',color='teal',zorder=5)
        ax[1].scatter(air,fwhm,marker='o',s=50,edgecolor='black',color='teal',zorder=5)


    fwhmRejectedFiles  = [x for x in allTable['File'][fwhmRejectedIndices] ]
    fwhmRejectedValues = [x for x in allTable['FWHM'][fwhmRejectedIndices] ]

    astrometryRejectedFiles  = [x for x in allTable['File'][astrometryRejectedIndices] ]
    astrometryRejectedValues = [x for x in allTable['FWHM'][astrometryRejectedIndices] ]

    
    for j in range(len(fwhmRejectedFiles)):
        match=re.search(pattern, fwhmRejectedFiles[j])
        frame=int(match.group(1))
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,dateKey)
        air=obtainKeyWordFromFits(file,airKey)
        if dateHeaderKey.startswith("DATE"):
            date_ok=datetime.fromisoformat(date)
        elif dateHeaderKey.startswith("MJD"):
            date_ok=Time(date,format='mjd').to_datetime()

        ax[0].scatter(date_ok, fwhmRejectedValues[j],marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=6,label='Rejected by FWHM' if (j==0) else "")
        ax[1].scatter(air, fwhmRejectedValues[j],marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=6,label='Rejected by FWHM'  if (j==0) else "")
        
    for j in range(len(astrometryRejectedFiles)):
        match=re.search(pattern, astrometryRejectedFiles[j])
        frame=int(match.group(1))
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,dateKey)
        air=obtainKeyWordFromFits(file,airKey)
        if dateHeaderKey.startswith("DATE"):
            date_ok=datetime.fromisoformat(date)
        elif dateHeaderKey.startswith("MJD"):
            date_ok=Time(date,format='mjd').to_datetime()

        ax[0].scatter(date_ok, astrometryRejectedValues[j],facecolors='none', lw=1.5, edgecolor='blue',color='blue',s=350,zorder=6,label='Rejected by astrometry' if (j==0) else "")
        ax[1].scatter(air, astrometryRejectedValues[j],facecolors='none', lw=1.5, edgecolor='blue',color='blue',s=350,zorder=6,label='Rejected by astrometry'  if (j==0) else "")
        
    ax[0].axhline(fwhmToReject,0,1,linestyle='--',linewidth=1,color='k',label='Rejection threshold')
    ax[1].axhline(fwhmToReject,0,1,linestyle='--',linewidth=1,color='k',label='Rejection threshold')
    ax[0].legend(fontsize=15, loc="upper right")
    for label in ax[0].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')
    plt.tight_layout()
    plt.savefig(imageName)
    return()

def obtainKeyWordFromFits(file, keyword):
    if os.path.exists(file):
        with fits.open(file) as hdul:
            header = hdul[HDU_TO_FIND_AIRMASS_DATE].header
            
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


def identifyBadFrames(folderWithFWHM, fwhmToReject):
    badFiles   = []
    allFiles   = []
    allFWHM     = []

    for currentFile in glob.glob(folderWithFWHM + "/fwhm_*.txt"):
        if fnmatch.fnmatch(currentFile, '*done*.txt'):
            continue

        fwhmValue = retrieveFWHMValues(currentFile,h,arcsecPerPix)
        if (math.isnan(fwhmValue)):
            continue
        
        allFiles.append(".".join(currentFile.split("/")[-1].split("_")[1].split(".")[:-1]))
        allFWHM.append(fwhmValue)

    allFWHM = np.array(allFWHM)

    mask = allFWHM>fwhmToReject
    allFiles = np.array(allFiles)
    allTogether=pd.DataFrame({'File':allFiles,'FWHM':allFWHM})
    badFiles = allFiles[mask]
    badFWHM = allFWHM[mask]
    return(badFiles,badFWHM,allTogether)
   
def getRejectedFramesFromFile(file):
    badFrames = []

    with open(file, 'r') as f:
        for i in f:
            badFrames.append(i[:-1] + ".fits")
    return(np.array(badFrames))

HDU_TO_FIND_AIRMASS_DATE = 0


folderWithFWHM            = sys.argv[1]
outputFolder              = sys.argv[2]
outputFile                = sys.argv[3]
fwhmToReject    = float(sys.argv[4])
folderWithFramesWithAirmasses = sys.argv[5]
h=int(sys.argv[6])
airMassKey=sys.argv[7]
dateHeaderKey=sys.argv[8]
arcsecPerPix=float(sys.argv[9])
setMatplotlibConf()
outputFolder_ccd=outputFolder+"/CCD"+str(h)

# 0.- Get frames identified as candidates for rejection in previous step (astrometrisation, background value and background std)
rejectedAtrometryFile = outputFolder+"/identifiedBadFrames_astrometry.txt"
rejectedAstrometrisedFrames = getRejectedFramesFromFile(rejectedAtrometryFile)



# 1.- Obtain the FWHM values ------------------------
fwhmValues = np.array([])
for currentFile in glob.glob(folderWithFWHM + "/fwhm_*.txt"):
    
    fwhmValue = retrieveFWHMValues(currentFile,h,arcsecPerPix)
    
    if (not math.isnan(fwhmValue)):
        fwhmValues = np.concatenate((fwhmValues, [fwhmValue]))

 
# 2.- Identify what frames are outside the acceptance region -----------------------
badFilesFWHM, badFWHM, allData = identifyBadFrames(folderWithFWHM, fwhmToReject)
allData.to_csv(outputFolder_ccd+"/FileFWHMtable.csv")

fwhmRejectedIndices            = getIndicesOfFiles(allData, badFilesFWHM)
astrometryRejectedIndices      = getIndicesOfFiles(allData, rejectedAstrometrisedFrames)

saveFWHMevol(allData, fwhmRejectedIndices, astrometryRejectedIndices,  outputFolder_ccd+"/fwhmEvol.png",airMassKey,dateHeaderKey,fwhmToReject)

pattern = r"\d+"
with open(outputFolder_ccd + "/" + outputFile, 'w') as file:
    for fileName in badFilesFWHM:
        match = re.search(pattern, fileName)
        result = match.group(0)
        file.write(result + '\n')

# 3.- Obtain the median and std and do teh histogram -------------------------------------
fwhmValueMean, fwhmValueStd = computeMedianAndStd(fwhmValues)
saveHistogram(allData, fwhmValueMean, fwhmValueStd, fwhmRejectedIndices, astrometryRejectedIndices,  outputFolder_ccd + "/fwhmHist.png", fwhmToReject, "FWHM of frames")
