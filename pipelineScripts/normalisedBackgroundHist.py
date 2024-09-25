# This python file produces a histogram of the air-mass normalised background values in surface-brightness units

# This is something similar to the histogram produced in the detection of bad frames based on the background values (script checkForBadFrames_backgroundValue.py)
# But at that point we cannot generate this histogram because we need the data calibrated.
# The thing is a little bit tricky because the calibrated images are already background-subtracted. So if we want the background values in magnitudes/arcsec² we need
# To combine the background estimation done prior to its subtraction with the calibration factors

import re
import sys
import glob
import math
import fnmatch
import astropy

import numpy as np
import matplotlib.pyplot as plt

from astropy.stats import sigma_clipped_stats

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

def retrieveCalibrationFactors(currentFile):
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) != 1):
            raise Exception("File with calibration factor contains more that 1 line. Expected 1 line got " + str(len(lines)))
        
        splittedLine = lines[0].strip().split()
        numberOfFields = len(splittedLine)

        if (numberOfFields == 2):
            return(float(splittedLine[0]))
        elif (numberOfFields == 0):
            return(float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 3 (constant estimation of the background), got " + str(numberOfFields))


def countsToSurfaceBrightnessUnits(values, arcsecPerPx):
    magnitudes = -2.5 * np.log10(values) + 22.5 + 5*np.log10(arcsecPerPx)
    return(magnitudes)

def calculateFreedmanBins(data, initialValue = None):
    if (initialValue == None):
        bins = [min(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= max(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

def saveHistogram(values, imageName):
    myBins = calculateFreedmanBins(values)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    counts, bins, patches = ax.hist(values, bins=myBins)
    max_bin_height = counts.max()
    ax.set_ylim(0, max_bin_height)

    plt.savefig(imageName)
    return()

def filesMatch(file1, file2):
    calibrationPattern = r"_(\d+)_"
    match = re.search(calibrationPattern, file1)
    calibrationNumber = match.group(1)

    backgroundPattern = r"entirecamera_(\d+)"
    match = re.search(backgroundPattern, file2)
    backgroundNumber = match.group(1)

    if (calibrationNumber == backgroundNumber):
        return(True)
    return(False)

def applyCalibrationFactorsToBackgroundValues(backgroundValues, calibrationFactors):
    calibratedValues = []

    # We first need to match the values
    for i in calibrationFactors:
        calibrationFile = i[0]
        for j in backgroundValues:
            backgroundFile = j[0]
            if (filesMatch(calibrationFile, backgroundFile)):
                calibratedValues.append(i[1] * j[1])
    return(calibratedValues)

folderWithCalibrationFactors  = sys.argv[1]
folderWithSkyEstimations      = sys.argv[2]
arcsecPerPx                   = float(sys.argv[3])
destinationFolder             = sys.argv[4]





# This is going to change. I'm already normalising the flux somewhere else so this code should not change much but
# instead taking the background as I'm taking now I will take the normalised one





# Obtain the background values that are not calibrated
totalBackgroundValues = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if (fnmatch.fnmatch(currentFile, '*/done*')): continue
    backgroundValue = retrieveBackgroundValues(currentFile)

    if (not math.isnan(backgroundValue)):
        totalBackgroundValues.append([currentFile.split('/')[-1], backgroundValue])

# Obtain the calibration factors
totalCalibrationFactors = []
for currentFile in glob.glob(folderWithCalibrationFactors + "/alpha_*.txt"):
    calibrationFactor = retrieveCalibrationFactors(currentFile)
    if (not math.isnan(calibrationFactor)):
        totalCalibrationFactors.append([currentFile.split('/')[-1], calibrationFactor])
    

values = applyCalibrationFactorsToBackgroundValues(totalBackgroundValues, totalCalibrationFactors)
magnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(values, arcsecPerPx)
saveHistogram(magnitudesPerArcSecSq, destinationFolder + "/magnitudeHist.png")