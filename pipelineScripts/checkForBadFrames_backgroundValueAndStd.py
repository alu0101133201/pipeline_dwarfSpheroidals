import os
import re
import sys
import glob
import math
import fnmatch
import astropy

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.stats import sigma_clipped_stats


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

def saveHistogram(values, median, std, imageName, numOfStd, title):
    valuesToPlot = values[~np.isnan(values)]
    myBins = calculateFreedmanBins(valuesToPlot)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_title(title)
    counts, bins, patches = ax.hist(valuesToPlot, bins=myBins)
    max_bin_height = counts.max() + 10
    ax.set_ylim(0, max_bin_height)

    ax.vlines(median, ymin=0, ymax=max_bin_height, color="black", linestyle="--", linewidth=2.0, label="Median")
    ax.vlines(median + numOfStd * std, ymin=0, ymax=max_bin_height, color="red", linestyle="--", linewidth=2.0, label=str(numOfStd) + " std")
    ax.vlines(median - numOfStd * std, ymin=0, ymax=max_bin_height, color="red", linestyle="--", linewidth=2.0)

    ax.legend(fontsize=18)
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
            return(float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 3 (constant estimation of the background), got " + str(numberOfFields))

    # Then we read the airmass
    airmass = obtainAirmassFromFile(currentFile, folderWithAirMasses, airMassKeyWord)
    return(backgroundValue / airmass, backgroundStd)

HDU_TO_FIND_AIRMASS = 1

folderWithSkyEstimations      = sys.argv[1]
folderWithFramesWithAirmasses = sys.argv[2] # The airmasses are in the header of the fits files that are in this folder
airMassKeyWord                = sys.argv[3]
outputFolder                  = sys.argv[4]
outputFile                    = sys.argv[5]
numberOfStdForRejecting       = float(sys.argv[6])


# 1.- Obtain the normalised background values and std values ------------------
normalisedBackgroundValues = []
backgroundStds             = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if fnmatch.fnmatch(currentFile, '*done*.txt'):
        continue
    currentValue, currentStd = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord)
    normalisedBackgroundValues.append(currentValue)
    backgroundStds.append(currentStd)
    
normalisedBackgroundValues = np.array(normalisedBackgroundValues)
backgroundStds = np.array(backgroundStds)

# 2.- Obtain the median and std and do the histograms ------------------
backgroundValueMedian, backgroundValueStd = computeMedianAndStd(normalisedBackgroundValues)
saveHistogram(normalisedBackgroundValues, backgroundValueMedian, backgroundValueStd, \
                outputFolder + "/backgroundHist.png", numberOfStdForRejecting, "Background values normalised by the airmass")

backgroundStdMedian, BackgroundStdStd = computeMedianAndStd(backgroundStds)
saveHistogram(backgroundStds, backgroundStdMedian, BackgroundStdStd, \
                outputFolder + "/backgroundStdHist.png", numberOfStdForRejecting, "Background std values")


# 3.- Identify what frames are outside the acceptance region ------------------
badFiles = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if fnmatch.fnmatch(currentFile, '*done*.txt'):
        continue
    currentValue, currentStd = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord)

    if (math.isnan(currentValue)):
        continue
    if ((currentValue > (backgroundValueMedian + (numberOfStdForRejecting * backgroundValueStd))) or \
        (currentValue < (backgroundValueMedian - (numberOfStdForRejecting * backgroundValueStd))) or \
        (currentStd   > (backgroundStdMedian   + (numberOfStdForRejecting * BackgroundStdStd))) or \
        (currentStd   < (backgroundStdMedian   - (numberOfStdForRejecting * BackgroundStdStd)))):
        badFiles.append(currentFile.split('/')[-1])
        
with open(outputFolder + "/" + outputFile, 'w') as file:
    for fileName in badFiles:
        file.write(fileName + '\n')

