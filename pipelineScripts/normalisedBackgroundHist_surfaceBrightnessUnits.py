# This python file produces a histogram of the air-mass normalised background values in surface-brightness units

# This is something similar to the histogram produced in the detection of bad frames based on the background values (script checkForBadFrames_backgroundValue.py)
# But at that point we cannot generate this histogram because we need the data calibrated.
# The thing is a little bit tricky because the calibrated images are already background-subtracted. So if we want the background values in magnitudes/arcsec² we need
# To combine the background estimation done prior to its subtraction with the calibration factors
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
        elif (numberOfFields == 1):
            return(float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 3 (constant estimation of the background), got " + str(numberOfFields))

    # Then we read the airmass
    airmass = obtainAirmassFromFile(currentFile, folderWithAirMasses, airMassKeyWord)
    return(backgroundValue / airmass)
    
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
    valuesToPlot = values[~np.isnan(values)]
    myBins = calculateFreedmanBins(valuesToPlot)
    myBins = np.linspace(np.nanmin(values), np.nanmax(values), 10)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    # ax.set_ylim(0, 60)
    counts, bins, patches = ax.hist(values, bins=myBins)
    max_bin_height = counts.max() + 10

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

    for j in backgroundValues:
        backgroundFile = j[0]
        found=False
        for i in calibrationFactors:
            calibrationFile = i[0]
            if (filesMatch(calibrationFile, backgroundFile)):
                calibratedValues.append(i[1] * j[1])
                found=True
                break
        if (not found):
            calibratedValues.append(np.nan)
                
    return(calibratedValues)

def removeBadFramesFromList(data, badFrames):
    noBadValues = []
    badFramesSet = set(badFrames)

    for i in data:
        match = re.search(r'(\d+)', i[0])
        if match:
            number = float(match.group(1))

            if number not in badFrames:
                noBadValues.append(i)
    return(noBadValues)

HDU_TO_FIND_AIRMASS = 1

folderWithSkyEstimations      = sys.argv[1]
folderWithFramesWithAirmasses = sys.argv[2] # The airmasses are in the header of the fits files that are in this folder
airMassKeyWord                = sys.argv[3]
folderWithCalibrationFactors  = sys.argv[4]
arcsecPerPx                   = float(sys.argv[5])
destinationFolder             = sys.argv[6]
rejectedFramesBackground      = sys.argv[7]
rejectedFramesFWHM            = sys.argv[8]


# 0.- Identify the files that have been identified as bad frames 
# This is needed because the data used comes from the noise-sky_it1 (since we need the background values) and the
# removal of bad frames is something done in posterior steps
badFrames = []
for currentFile in glob.glob(rejectedFramesBackground + "/*.fits"):
    match = re.search(r'entirecamera_(\d+)\.fits', currentFile)
    if match:
        number = int(match.group(1))  
        badFrames.append(number)
    else:
        raise Exception("Error identifying the number of the bad frames (background bad frames)")

for currentFile in glob.glob(rejectedFramesFWHM + "/*.fits"):
    match = re.search(r'entirecamera_(\d+)\.fits', currentFile)
    if match:
        number = int(match.group(1))  
        badFrames.append(number)
    else:
        raise Exception("Error identifying the number of the bad frames (FWHM bad frames)")


# 1.- Obtain the normalised background values ------------------
normalisedBackgroundValues = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if fnmatch.fnmatch(currentFile, '*done*.txt'):
        continue
    currentValue = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord)
    normalisedBackgroundValues.append([currentFile.split('/')[-1], currentValue])

# 2.- Obtain the calibration factors
totalCalibrationFactors = []
for currentFile in glob.glob(folderWithCalibrationFactors + "/alpha_*.txt"):
    calibrationFactor = retrieveCalibrationFactors(currentFile)
    if (not math.isnan(calibrationFactor)):
        totalCalibrationFactors.append([currentFile.split('/')[-1], calibrationFactor])
    
normalisedBackgroundValues = removeBadFramesFromList(normalisedBackgroundValues, badFrames)
totalCalibrationFactors    = removeBadFramesFromList(totalCalibrationFactors, badFrames)

values = applyCalibrationFactorsToBackgroundValues(normalisedBackgroundValues, totalCalibrationFactors)


magnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(values, arcsecPerPx)
saveHistogram(np.array(magnitudesPerArcSecSq), destinationFolder + "/magnitudeHist.png")

# Temporal code, just for checking the relation between background counts and magnitude
# This gives information of the calibration factors

x = []
x = [i[1] for i in normalisedBackgroundValues]

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.set_yscale('log')
ax.scatter(magnitudesPerArcSecSq, x)
plt.savefig(destinationFolder + "/countsVsMagnitudes.png")