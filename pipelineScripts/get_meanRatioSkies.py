import os
import re
import sys
import glob

import fnmatch


import numpy as np

import pandas as pd

from astropy.visualization import astropy_mpl_style

from astropy.io import fits
from astropy.stats import sigma_clipped_stats



HDU_TO_FIND_AIRMASS = 0
folderWithSkyEstimations = sys.argv[1]
commonCalibrationFactorFile = sys.argv[2]
outputFile = sys.argv[3]
ccd_ref=int(sys.argv[4])
num_ccd=int(sys.argv[5])
airMassKeyWord = sys.argv[6]
folderWithFramesWithAirmasses = sys.argv[7]
arcsecPerPx = float(sys.argv[8])
def obtainNormalisedBackground(currentFile, folderWithAirMasses, airMassKeyWord,h):
    backgroundValue = -1

    # First we read the background Value
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) < h):
            raise Exception("File with the background estimation contains more that 1 line. Expected 1 line got " + str(len(lines)))
        
        splittedLine = lines[(h-1)].strip().split()
        numberOfFields = len(splittedLine)

        if (numberOfFields == 5) or (numberOfFields == 3):
            backgroundValue = float(splittedLine[1])
        elif (numberOfFields == 1):
            return(float('nan'), float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 5 (constant estimation of the background), got " + str(numberOfFields))

    # Then we read the airmass
    try:
        airmass = obtainAirmassFromFile(currentFile, folderWithAirMasses, airMassKeyWord)
    except:
        print("Something went wrong in obtaining the airmass, returning nans (file " + str(currentFile) + ")")
        return(float('nan'), float('nan')) 
    return(backgroundValue, backgroundValue / airmass)
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
                if (keywordValue == "" or keywordValue == None):
                    keywordValue=np.nan
                return(keywordValue)
            else:
                raise Exception(f"Keyword '{keyword}' not found in the header.")
    else:
        raise Exception(f"File {file} does not exist.")

def obtainAirmassFromFile(currentFile, airMassesFolder, airMassKeyWord):
    frameNumber = obtainNumberFromFrame(currentFile)

    fitsFileNamePatter = f"{frameNumber}.fits"
    fitsFilePath = os.path.join(airMassesFolder, fitsFileNamePatter)

    airMass = obtainKeyWordFromFits(fitsFilePath, airMassKeyWord)
    return(airMass)
def filesMatch(file1, file2):
    calibrationPattern = r"_(\d+)."
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
            
            if ((not pd.isna(calibrationFile)) and (not pd.isna(backgroundFile))):
                if (filesMatch(calibrationFile, backgroundFile)):
                    cValues_row=[]
                    for h in range(1, num_ccd + 1):
                        cValues_row.append(i[1][h-1] * j[h])
                    calibratedValues.append(cValues_row)
                    found=True
                    break

        if (not found):
            calibratedValues=np.append(calibratedValues, [np.nan] * num_ccd)
                
    return(np.array(calibratedValues))
def countsToSurfaceBrightnessUnits(values, arcsecPerPx):
    magnitudes = -2.5 * np.log10(values) + 22.5 + 5*np.log10(arcsecPerPx)
    return(magnitudes)

pattern = os.path.join(folderWithFramesWithAirmasses, "*.fits")
totalNumberOfFrames = len(glob.glob(pattern))

originalBackgroundValues =  np.array([[np.nan] * (1+num_ccd) for _ in range(totalNumberOfFrames)], dtype=object)
normalisedBackgroundValues =  np.array([[np.nan] * (1+num_ccd) for _ in range(totalNumberOfFrames)], dtype=object)
files = np.full(totalNumberOfFrames, np.nan, dtype=object)
 
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if fnmatch.fnmatch(currentFile, '*done*.txt'):
        continue

    match = re.search(r'_(\d+)\.', currentFile)
    if match:
        number = int(match.group(1))  
    else:
        raise Exception(f"Number not found in the file name ({currentFile}). Something went wrong here")
    files[number-1] = currentFile
    normalisedBackgroundValues[number-1,0] = currentFile.split('/')[-1]
    originalBackgroundValues[number-1,0] = currentFile.split('/')[-1]
    for h in range(1, num_ccd + 1):
        originalBackground, normalisedBackground = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord,h)
        normalisedBackgroundValues[number-1,h] = normalisedBackground
        originalBackgroundValues[number-1,h] = originalBackground
commonCalibrationFactorValue=np.array([np.nan] * num_ccd, dtype=float)
with open(commonCalibrationFactorFile) as f:
   fileContent=f.read().split()
   for h in range(num_ccd):
       commonCalibrationFactorValue[h] = float(fileContent[2*h])

arrayWithCommonFactor = [[x[0], [cfac for cfac in commonCalibrationFactorValue]] for x in normalisedBackgroundValues]
valuesCalibratedOriginal = applyCalibrationFactorsToBackgroundValues(originalBackgroundValues, arrayWithCommonFactor)
valuesCalibratedNormalised = applyCalibrationFactorsToBackgroundValues(normalisedBackgroundValues, arrayWithCommonFactor)

magnitudesPerArcSecSqNormalised = countsToSurfaceBrightnessUnits(valuesCalibratedNormalised, arcsecPerPx)
diffs=np.zeros((totalNumberOfFrames,num_ccd), dtype=float)

for h in range(num_ccd):
    diffs[:,h]= magnitudesPerArcSecSqNormalised[:,ccd_ref-1] - magnitudesPerArcSecSqNormalised[:,h]

# Now we save the mean results
collapsed = np.mean(diffs, axis=0)
parameter=10**(-0.4*collapsed)

# Save to TXT
np.savetxt(outputFile, parameter, fmt="%.12f")