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


def extractNumberFromName(filename):
    match = re.search(r"entirecamera_(\d+).txt", filename)
    if match:
        return int(match.group(1)) 
    else:
        raise Exception("Something went wrong with the fileNames")

def retrieveBackgroundValues(currentFile):
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) != 1):
            raise Exception("File with the background estimation contains more that 1 line. Expected 1 line got " + str(len(lines)))
        
        splittedLine = lines[0].strip().split()
        numberOfFields = len(splittedLine)

        if (numberOfFields == 3):
            return(float(splittedLine[1]))
        elif (numberOfFields == 1):
            return(float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 3 (constant estimation of the background), got " + str(numberOfFields))

def readAirMassesFromFile(fileName):
    values = []
    try:
        with open(fileName, 'r') as file:
            for line in file:
                print(line)
                values.append(float(line.strip()))
    except FileNotFoundError:
        print(f"Error: The file {fileName} was not found.")
    except ValueError:
        print("Error: Could not convert some lines to float. Please check the file content.")

    return values

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
    ax.set_title(title, fontsize=22, pad=17)

    counts, bins, patches = ax.hist(valuesToPlot, color="teal") #, bins=myBins)
    max_bin_height = counts.max() + 5
    ax.set_ylim(0, max_bin_height)

    ax.text(0.3755, 0.95, "Median: " + str(int(median)), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')
    ax.text(0.375, 0.9, "Std: " + str(int(std)), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')
    if len(badValues)!=0:
        if "STD" in labelX:
            type_rejection="STD"
        elif "counts" in labelX:
            type_rejection="Counts"
        counts_bad, bins_bad, patches_bad = ax.hist(badValues,bins=np.linspace(np.nanmin(badValues), np.nanmax(badValues), 4),color='red',label='Rejected'+type_rejection)

        
        ax.legend()
    ax.text(0.3655, 0.85, "Rejected: " + str(len(badValues)), transform=ax.transAxes, 
        fontsize=20, verticalalignment='top', horizontalalignment='left')
    # ax.text(0.1, 0.9, 'Mean: ' + "{:.0}".format(median), transform=ax.transAxes, 
    #     fontsize=18, verticalalignment='top', horizontalalignment='left')
    # ax.text(0.1, 0.85, 'Std: ' + "{:.0}".format(std), transform=ax.transAxes, 
    #     fontsize=18, verticalalignment='top', horizontalalignment='left')
    plt.savefig(imageName)
    return()


def obtainNumberFromFrame(currentFile):
    match = re.search(r'entirecamera_(\d+)\.txt', currentFile)
    if match:
        frameNumber = int(match.group(1))
        return(frameNumber)
    else:
        raise Exception("Something when wrong when trying to obtain the number frame (needed for accessing the airmass of it)")

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

def obtainAirmassFromFile(currentFile, airMassesFolder, airMassKeyWord):
    frameNumber = obtainNumberFromFrame(currentFile)

    fitsFileNamePatter = f"{frameNumber}.fits"
    fitsFilePath = os.path.join(airMassesFolder, fitsFileNamePatter)

    print("file for obtaining airmass: ", fitsFilePath)

    airMass = obtainKeyWordFromFits(fitsFilePath, airMassKeyWord)
    return(airMass)

def obtainNormalisedBackground(currentFile, folderWithAirMasses, airMassKeyWord):
    backgroundValue = -1

    # First we read the background Value
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) != 1):
            raise Exception("File with the background estimation contains more that 1 line. Expected 1 line got " + str(len(lines)))
        
        splittedLine = lines[0].strip().split()
        numberOfFields = len(splittedLine)

        if (numberOfFields == 3):
            backgroundValue = float(splittedLine[1])
            backgroundStd   = float(splittedLine[2])
        elif (numberOfFields == 1):
            return(float('nan'), float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 3 (constant estimation of the background), got " + str(numberOfFields))

    # Then we read the airmass
    try:
        airmass = obtainAirmassFromFile(currentFile, folderWithAirMasses, airMassKeyWord)
    except:
        print("something went wrong in obtaining the airmass, returning nans (file: " + str(currentFile) + ")")
        return(float('nan'), float('nan')) 
    return(backgroundValue / airmass, backgroundStd)

def saveBACKevol(allTable,badFiles,badBack,imageName):
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    configureAxis(ax, 'Frame Number', 'Background (ADU)',logScale=False)
    ax.set_title('Background evolution',fontsize=22,pad=17)
    pattern=r"entirecamera_(\d+)"
    for row in range(len(allTable)):
        file=allTable.loc[row]['File']
        match=re.search(pattern,file)
        frame=float(match.group(1))
        bck=allTable.loc[row]['Background']
        ax.scatter(frame,bck,marker='o',s=50,edgecolor='black',color='teal',zorder=5)
    
    if len(badFiles)!=0:
        for j in range(len(badFiles)):
            match=re.search(pattern,badFiles[j])
            frame=float(match.group(1))
            ax.scatter(frame,badBack[j],marker='X',edgecolor='k',color='darkred',s=80,zorder=6,label='Rejected back.')
        ax.legend()
    plt.tight_layout()
    plt.savefig(imageName)
    return()
def saveSTDevol(allTable,badFiles,badSTD,imageName):
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    configureAxis(ax, 'Frame Number', 'Background STD (ADU)',logScale=False)
    ax.set_title('STD evolution',fontsize=22,pad=17)
    pattern=r"entirecamera_(\d+)"
    for row in range(len(allTable)):
        file=allTable.loc[row]['File']
        match=re.search(pattern,file)
        frame=float(match.group(1))
        bck=allTable.loc[row]['STD']
        ax.scatter(frame,bck,marker='o',s=50,edgecolor='black',color='teal',zorder=5)
    
    if len(badFiles)!=0:
        for j in range(len(badFiles)):
            match=re.search(pattern,badFiles[j])
            frame=float(match.group(1))
            ax.scatter(frame,badSTD[j],marker='X',edgecolor='k',color='darkred',s=80,zorder=6,label='Rejected STD')
        ax.legend()
    plt.tight_layout()
    plt.savefig(imageName)
    return()
def identifyBadFrames(folderWithFrames, folderWithFramesWithAirmasses, airMassKeyWord, numberOfStdForRejecting):
    badFiles   = []
    allFiles   = []
    allStd     = []
    allBackgroundValues = []

    for currentFile in glob.glob(folderWithFrames + "/*.txt"):
        if fnmatch.fnmatch(currentFile, '*done*.txt'):
            continue
        currentValue, currentStd = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord)

        if (math.isnan(currentValue)):
            continue
        allFiles.append(currentFile)
        allStd.append(currentStd)
        allBackgroundValues.append(currentValue)

    allStd = np.array(allStd)
    allBackgroundValues = np.array(allBackgroundValues)

    std_mask = sigma_clip(allStd, sigma=numberOfStdForRejecting, cenfunc='median', stdfunc='std', maxiters=5, masked=True).mask
    values_mask = sigma_clip(allBackgroundValues, sigma=numberOfStdForRejecting, cenfunc='median', stdfunc='std', maxiters=5, masked=True).mask

    combined_mask = values_mask | std_mask
    allFiles = np.array(allFiles)
    allTogether=pd.DataFrame({'File':allFiles,'Background':allBackgroundValues,'STD':allStd})

    badFiles = allFiles[combined_mask]
    badValues=allBackgroundValues[values_mask]; badFilesBCK=allFiles[values_mask]
    badStd=allStd[std_mask]; badFilesSTD=allFiles[std_mask]
    return(badFiles,badValues,badStd,allTogether,badFilesBCK,badFilesSTD)


HDU_TO_FIND_AIRMASS = 1

folderWithSkyEstimations      = sys.argv[1]
folderWithFramesWithAirmasses = sys.argv[2] # The airmasses are in the header of the fits files that are in this folder
airMassKeyWord                = sys.argv[3]
outputFolder                  = sys.argv[4]
outputFile                    = sys.argv[5]
numberOfStdForRejecting       = float(sys.argv[6])

setMatplotlibConf()

# 1.- Obtain the normalised background values and std values ------------------
normalisedBackgroundValues = []
backgroundStds             = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    print("current file: ", currentFile)
    if fnmatch.fnmatch(currentFile, '*done*.txt'):
        continue
    currentValue, currentStd = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord)
    normalisedBackgroundValues.append(currentValue)
    backgroundStds.append(currentStd)
    
normalisedBackgroundValues = np.array(normalisedBackgroundValues)
backgroundStds = np.array(backgroundStds)

print("\n\n")


# 2.- Identify what frames are outside the acceptance region ------------------
badFiles,badValues,badStd,allData,badFilesBCK,badFilesSTD = identifyBadFrames(folderWithSkyEstimations,folderWithFramesWithAirmasses, airMassKeyWord, numberOfStdForRejecting)
saveBACKevol(allData,badFilesBCK,badValues,outputFolder+"/backgroundEvolution.png")
saveSTDevol(allData,badFilesSTD,badStd,outputFolder+"/stdEvolution.png")

with open(outputFolder + "/" + outputFile, 'w') as file:
    for fileName in badFiles:
        file.write(fileName + '\n')

# 3.- Obtain the median and std and do the histograms ------------------
backgroundValueMedian, backgroundValueStd = computeMedianAndStd(normalisedBackgroundValues)
saveHistogram(normalisedBackgroundValues, backgroundValueMedian, backgroundValueStd,badValues, \
                outputFolder + "/backgroundHist.png", numberOfStdForRejecting, "Background values normalised by the airmass", "Background counts (ADU)")

backgroundStdMedian, BackgroundStdStd = computeMedianAndStd(backgroundStds)
saveHistogram(backgroundStds, backgroundStdMedian, BackgroundStdStd, badStd,\
                outputFolder + "/backgroundStdHist.png", numberOfStdForRejecting, "Background std values", "Background STD (ADU)")
