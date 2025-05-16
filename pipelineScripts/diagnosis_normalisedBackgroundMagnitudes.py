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
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

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

def obtainNormalisedBackground(currentFile, folderWithAirMasses, airMassKeyWord):
    backgroundValue = -1

    # First we read the background Value
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) != 1):
            raise Exception("File with the background estimation contains more that 1 line. Expected 1 line got " + str(len(lines)))
        
        splittedLine = lines[0].strip().split()
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
            return(float(np.nan)) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 3 (constant estimation of the background), got " + str(numberOfFields))


def countsToSurfaceBrightnessUnits(values, arcsecPerPx):
    magnitudes = -2.5 * np.log10(values) + 22.5 + 5*np.log10(arcsecPerPx)
    return(magnitudes)

def calculateFreedmanBins(data, initialValue = None):
    if (initialValue == None):
        bins = [np.nanmin(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= np.nanmax(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

def saveHistogram(values, rejectedAstrometryIndices, rejectedFWHMIndices, rejectedBackgroundIndices, title, imageName):
    myBins = calculateFreedmanBins(values[~np.isnan(values)])

    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_title(title, fontsize=22, pad=17)
    plt.tight_layout(pad=7.0)
    configureAxis(ax, 'Background (mag/arcsec^2)', '', logScale=False)
    counts, bins, patches = ax.hist(values, bins=myBins, color="teal")

    if (len(rejectedBackgroundIndices) > 0):
        ax.hist(values[rejectedBackgroundIndices - 1], bins=myBins, color="red", label="Rejected by background brightness")
    if (len(rejectedFWHMIndices)):
        ax.hist(values[rejectedFWHMIndices - 1], bins=myBins, color="mediumorchid", label="Rejected by fwhm")
    if (len(rejectedAstrometryIndices)):
        ax.hist(values[rejectedAstrometryIndices - 1], bins=myBins, color="blue", label="Rejected by astrometry")

    max_bin_height = counts.max() + 5
    ax.set_ylim(0, max_bin_height)
    ax.legend(fontsize=20, loc="upper right")
    plt.savefig(imageName)
    return()

def saveScatterFactors(factors, rejectedAstrometryIndices, rejectedFWHMIndices, rejectedBackgroundIndices, title, imageName, folderWithFramesWithAirmasses, destinationFolder):
    airMass  = []
    time     = []
    cfactors = []


    for i in factors:
        if (not pd.isna(i[0])): 
            match=re.search(r"_(\d+).",i[0])
            frame = match.group(1)
            file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
            date=obtainKeyWordFromFits(file,'DATE-OBS')
            air=obtainKeyWordFromFits(file,'AIRMASS')
            date_ok=datetime.fromisoformat(date)
            cfactors.append(i[1])
            airMass.append(air)
            time.append(date_ok)
        else:
            airMass.append(np.nan)
            time.append(np.nan)
            cfactors.append(np.nan)


    airMass=np.array(airMass,dtype='float')
    cfactors=np.array(cfactors,dtype='float')
    time = np.array(time)


    with open(destinationFolder + "/cfactors.txt", "w") as file:
        for a, b, c in zip(time, airMass, cfactors):
            file.write(f"{a}\t{b}\t{c}\n")  # Tab-separated columns




    fig, ax = plt.subplots(2, 1, figsize=(20,10))
    configureAxis(ax[0], 'UTC', 'Calibration Factor',logScale=False)
    configureAxis(ax[1], 'Airmass', 'Calibration Factor',logScale=False)
    fig.suptitle('Calibration factor evolution',fontsize=22)

    pattern=r"(\d+).fits"

    timeMask = ~pd.isna(time) & ~pd.isna(cfactors)
    ax[0].scatter(time[timeMask],cfactors[timeMask],marker='o',s=50,edgecolor='black',color='teal',zorder=0)
    ax[1].scatter(airMass,cfactors,marker='o',s=50,edgecolor='black',color='teal',zorder=0)

    if (len(rejectedAstrometryIndices) > 0):
        timeMask = ~pd.isna(time[rejectedAstrometryIndices - 1]) & ~pd.isna(cfactors[rejectedAstrometryIndices - 1])
        ax[0].scatter(time[rejectedAstrometryIndices - 1][timeMask], cfactors[rejectedAstrometryIndices - 1][timeMask], facecolors='none', lw=1.5, edgecolor='blue', s=350,zorder=1, label="Rejected by astrometry")
        ax[1].scatter(airMass[rejectedAstrometryIndices - 1], cfactors[rejectedAstrometryIndices - 1], facecolors='none', lw=1.5, edgecolor='blue', s=350,zorder=1, label="Rejected by astrometry")
    
    if (len(rejectedFWHMIndices) > 0):
        timeMask = ~pd.isna(time[rejectedFWHMIndices - 1]) & ~pd.isna(cfactors[rejectedFWHMIndices - 1])
        ax[0].scatter(time[rejectedFWHMIndices - 1][timeMask], cfactors[rejectedFWHMIndices - 1][timeMask], marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=1, label="Rejected by FWHM")
        ax[1].scatter(airMass[rejectedFWHMIndices - 1], cfactors[rejectedFWHMIndices - 1], marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=1, label="Rejected by FWHM")
   
    if (len(rejectedBackgroundIndices) > 0):
        timeMask = ~pd.isna(time[rejectedBackgroundIndices - 1]) & ~pd.isna(cfactors[rejectedBackgroundIndices - 1])
        ax[0].scatter(time[rejectedBackgroundIndices - 1][timeMask], cfactors[rejectedBackgroundIndices - 1][timeMask], marker='D',edgecolor='k',color='red',s=120,zorder=1, label="Rejected by background brightness")
        ax[1].scatter(airMass[rejectedBackgroundIndices - 1], cfactors[rejectedBackgroundIndices - 1], marker='D',edgecolor='k',color='red',s=120,zorder=1, label="Rejected by background brightness")



    ax[0].legend(loc="upper right", fontsize=20)
    
    ax[0].legend(loc="upper right", fontsize=20)
    fig.suptitle(title,fontsize=22)
    for label in ax[0].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')
    plt.tight_layout()
    
    plt.savefig(imageName)
    return()

def saveBackEvolution(magnitudesPerArcSecSq, rejectedAstrometryIndices, rejectedFWHMIndices, rejectedBackgroundIndices, title, imageName, folderWithFramesWithAirmasses, maxmimumBackgroundBrightness):
    airMass  = []
    time     = []
    

    for i in normalisedBackgroundValues:
        if (not pd.isna(i[0])): 
            match=re.search(r"_(\d+).",i[0])
            frame = match.group(1)
            file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
            date=obtainKeyWordFromFits(file,'DATE-OBS')
            air=obtainKeyWordFromFits(file,'AIRMASS')
            date_ok=datetime.fromisoformat(date)
            
            airMass.append(air)
            time.append(date_ok)
        else:
            airMass.append(np.nan)
            time.append(np.nan)

    airMass=np.array(airMass,dtype='float')
    time = np.array(time)

    fig, ax = plt.subplots(2, 1, figsize=(20,10))
    configureAxis(ax[0], 'UTC', '',logScale=False)
    configureAxis(ax[1], 'Airmass', '',logScale=False)
    fig.suptitle(title,fontsize=22)
    pattern=r"(\d+).fits"

    timeMask = ~pd.isna(time) & ~pd.isna(magnitudesPerArcSecSq)
    ax[0].scatter(time[timeMask],magnitudesPerArcSecSq[timeMask],marker='o',s=50,edgecolor='black',color='teal',zorder=0)
    ax[1].scatter(airMass,magnitudesPerArcSecSq,marker='o',s=50,edgecolor='black',color='teal',zorder=0)

    if (len(rejectedAstrometryIndices) > 0):
        timeMask = ~pd.isna(time[rejectedAstrometryIndices - 1]) & ~pd.isna(magnitudesPerArcSecSq[rejectedAstrometryIndices - 1])
        ax[0].scatter(time[rejectedAstrometryIndices - 1][timeMask], magnitudesPerArcSecSq[rejectedAstrometryIndices - 1][timeMask], facecolors='none', lw=1.5, edgecolor='blue', s=350,zorder=1, label="Rejected by astrometry")
        ax[1].scatter(airMass[rejectedAstrometryIndices - 1], magnitudesPerArcSecSq[rejectedAstrometryIndices - 1], facecolors='none', lw=1.5, edgecolor='blue', s=350,zorder=1, label="Rejected by astrometry")
    
    if (len(rejectedFWHMIndices) > 0):
        timeMask = ~pd.isna(time[rejectedFWHMIndices - 1]) & ~pd.isna(magnitudesPerArcSecSq[rejectedFWHMIndices - 1])
        ax[0].scatter(time[rejectedFWHMIndices - 1][timeMask], magnitudesPerArcSecSq[rejectedFWHMIndices - 1][timeMask], marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=1, label="Rejected by FWHM")
        ax[1].scatter(airMass[rejectedFWHMIndices - 1], magnitudesPerArcSecSq[rejectedFWHMIndices - 1], marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=1, label="Rejected by FWHM")
    
    if (len(rejectedBackgroundIndices) > 0):
        timeMask = ~pd.isna(time[rejectedBackgroundIndices - 1]) & ~pd.isna(magnitudesPerArcSecSq[rejectedBackgroundIndices - 1])
        ax[0].scatter(time[rejectedBackgroundIndices - 1][timeMask], magnitudesPerArcSecSq[rejectedBackgroundIndices - 1][timeMask], marker='D',edgecolor='k',color='red',s=120,zorder=1, label="Rejected by background brightness")
        ax[1].scatter(airMass[rejectedBackgroundIndices - 1], magnitudesPerArcSecSq[rejectedBackgroundIndices - 1], marker='D',edgecolor='k',color='red',s=120,zorder=1, label="Rejected by background brightness")


    time = time[~pd.isna(time)]
    ax[0].hlines(maxmimumBackgroundBrightness, xmin=np.nanmin(time), xmax=np.nanmax(time), color="black", ls="--", lw=2)
    ax[1].hlines(maxmimumBackgroundBrightness, xmin=np.nanmin(airMass), xmax=np.nanmax(airMass), color="black", ls="--", lw=2)


    ax[0].legend(loc="upper right", fontsize=20)
    
    for label in ax[0].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')
    fig.supylabel(r'Back. [mag arcsec$^{-2}$]',fontsize=30)
    plt.tight_layout()
    plt.savefig(imageName)
    return()

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
                    calibratedValues.append(i[1] * j[1])
                    found=True
                    break

        if (not found):
            calibratedValues.append(np.nan)
                
    return(calibratedValues)

def removeBadFramesFromList(data, badFrames_bck, badFrames_fwhm):
    noBadValues = []
    badFramesSet_bck = set(badFrames_bck)
    badFramesSet_fwhm = set(badFrames_fwhm)
    BadValues_bck = []
    BadValues_fwhm = []
    for i in data:
        
        match = re.search(r'_(\d+).', i[0])
        if match:
            number = float(match.group(1))

            if number in badFrames_bck:
                BadValues_bck.append(i)
            elif number in badFrames_fwhm:
                BadValues_fwhm.append(i)
            else:
                noBadValues.append(i)
    return(noBadValues,BadValues_bck,BadValues_fwhm)

def scatterPlotCountsVsMagnitudes(backgroundCounts, magnitudesPerArcSecSq, rejectedAstrometryIndices, rejectedFWHMIndices, fileName):
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_title("If calibration is correct, this should be a straight line", fontsize=20, pad=17)
    configureAxis(ax, 'Background (mag/arcsec^2)', 'Log(Background) (ADU)', logScale=False)

    plt.tight_layout(pad=6.0)
    ax.scatter(magnitudesPerArcSecSq, np.log10(backgroundCounts), s=40, color="teal")

    ax.scatter(magnitudesPerArcSecSq[rejectedAstrometryIndices], np.log10(backgroundCounts)[rejectedAstrometryIndices], facecolors='none', lw=1.5, edgecolor='blue',s=120, label="Rejected astrometry")
    ax.scatter(magnitudesPerArcSecSq[rejectedFWHMIndices], np.log10(backgroundCounts)[rejectedFWHMIndices], s=40, color="mediumorchid", label="Rejected FWHM")

    ax.legend(fontsize=18, loc="upper right")
    plt.savefig(fileName)

def getIndicesOfRejectedFrames(normalisedBackgroundValuesArray, rejectedFrames):
    rejectedFrames = []
    pattern = r"\d+"

    for j in rejectedFrames:
        rejectedFrames
    return(indices)

def identifyBadFrames(folderWithFrames,folderWithAirmasses,airMassKeyWord,folderWithCalibrationFactors):
    
    allFiles=[]
    allBackValues=[]
    for currentFile in glob.glob(folderWithFrames + "/*.txt"):
        if fnmatch.fnmatch(currentFile,'*done*.txt'):
            continue
        currentValue = obtainNormalisedBackground(currentFile, folderWithAirmasses, airMassKeyWord)
        matchID=re.search(r'_(\d+).',currentFile).group(1)
        calFactorFile=glob.glob(f"{folderWithCalibrationFactors}/alpha*Decals*{matchID}.txt")[0]
        calFactor=retrieveCalibrationFactors(calFactorFile)
        valueCalibrated = currentValue*calFactor

        magnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(valueCalibrated, arcsecPerPx)
        allFiles.append(currentFile)
        allBackValues.append(magnitudesPerArcSecSq)
    allFiles=np.array(allFiles)
    allBackValues=np.array(allBackValues)
    values_mask=allBackValues<20.5
    badFiles=allFiles[values_mask]
    badValues=allBackValues[values_mask]
    return(badFiles,badValues)

def identifyBadFrames(files, data, threshold):
    badFiles   = []
    badBackground = []
    for i in range(len(files)):
        if (data[i] < threshold):
            badFiles.append((files[i].split("/")[-1].split("_")[1].split(".")[:-1][0]))
            badBackground.append(data[i])
    badFiles = np.array(badFiles)
    badBackground = np.array(badBackground)

    allTogether=pd.DataFrame({'File':files,'FWHM':data})
    return(badFiles,badBackground,allTogether)


HDU_TO_FIND_AIRMASS = 1

folderWithSkyEstimations      = sys.argv[1]
folderWithFramesWithAirmasses = sys.argv[2] # The airmasses are in the header of the fits files that are in this folder
airMassKeyWord                = sys.argv[3]
folderWithCalibrationFactors  = sys.argv[4]
arcsecPerPx                   = float(sys.argv[5])
destinationFolder             = sys.argv[6]
maxmimumBackgroundBrightness  = float(sys.argv[7])
outputFile                    = sys.argv[8]

setMatplotlibConf()



pattern = os.path.join(folderWithFramesWithAirmasses, "*.fits")
totalNumberOfFrames = len(glob.glob(pattern))

# 0.- Identify the files that have been identified as bad frames 
# We have bad frames due to:
#   Astrometry
#   FWHM

rejectedFrames_astrometry = []
with open(destinationFolder + "/identifiedBadFrames_astrometry.txt", 'r') as f:
    for line in f:
        rejectedFrames_astrometry.append(int(line.strip()))
rejectedFrames_astrometry = np.array(rejectedFrames_astrometry)


rejectedFrames_FWHM = []
with open(destinationFolder + "/identifiedBadFrames_fwhm.txt", 'r') as f:
    for line in f:
        rejectedFrames_FWHM.append(int(line.strip()))
rejectedFrames_FWHM = np.array(rejectedFrames_FWHM)

# 1.- Obtain the normalised background values ------------------
originalBackgroundValues =  np.array([[np.nan, np.nan] for _ in range(totalNumberOfFrames)], dtype=object)
normalisedBackgroundValues =  np.array([[np.nan, np.nan] for _ in range(totalNumberOfFrames)], dtype=object)
files = np.full(totalNumberOfFrames, np.nan, dtype=object)
 
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if fnmatch.fnmatch(currentFile, '*done*.txt'):
        continue

    match = re.search(r'_(\d+)\.', currentFile)
    if match:
        number = int(match.group(1))  # or leave as string if needed
    else:
        raise Exception(f"Number not found in the file name ({currentFile}). Something went wrong here")

    originalBackground, normalisedBackground = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord)

    files[number-1] = currentFile
    normalisedBackgroundValues[number-1] = [currentFile.split('/')[-1], normalisedBackground]
    originalBackgroundValues[number-1] = [currentFile.split('/')[-1], originalBackground]

files = np.array(files)


# 2.- Obtain the calibration factors
totalCalibrationFactors = np.array([[np.nan, np.nan] for _ in range(totalNumberOfFrames)], dtype=object)

for currentFile in glob.glob(folderWithCalibrationFactors + "/alpha_*Decals*.txt"):
    match = re.search(r'_(\d+)\.', currentFile)
    if match:
        number = int(match.group(1))  # or leave as string if needed
    else:
        raise Exception(f"Number not found in the file name ({currentFile}). Something went wrong here")

    calibrationFactor = retrieveCalibrationFactors(currentFile)
    totalCalibrationFactors[number-1] = [currentFile.split('/')[-1], calibrationFactor]


valuesCalibrated = applyCalibrationFactorsToBackgroundValues(normalisedBackgroundValues, totalCalibrationFactors)

magnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(valuesCalibrated, arcsecPerPx)

with open(destinationFolder + "/backgroundMagnitudes.dat", 'w') as f:
    for i in range(len(magnitudesPerArcSecSq)):
        f.write(str(files[i]) + " " + str(magnitudesPerArcSecSq[i]) + ("\n"))

badFilesBackground, badBackground, allData = identifyBadFrames(files, magnitudesPerArcSecSq, maxmimumBackgroundBrightness)

pattern = r"\d+"
with open(destinationFolder + "/" + outputFile, 'w') as file:
    for fileName in badFilesBackground:
        match = re.search(pattern, fileName)
        result = match.group(0)
        file.write(result + '\n')


rejectedFrames_Background = np.array(badFilesBackground).astype(int)


saveHistogram(np.array(magnitudesPerArcSecSq), rejectedFrames_astrometry, rejectedFrames_FWHM, rejectedFrames_Background, \
                "Distribution of NORMALISED background magnitudes", destinationFolder + "/magnitudeHist.png")


saveScatterFactors(totalCalibrationFactors, rejectedFrames_astrometry, rejectedFrames_FWHM, rejectedFrames_Background, \
                "Evolution of calibration factors",destinationFolder + "/calibrationFactorEvolution.png", folderWithFramesWithAirmasses, destinationFolder)



originalValuesCalibrated = applyCalibrationFactorsToBackgroundValues(originalBackgroundValues, totalCalibrationFactors)
valuesCalibrated = applyCalibrationFactorsToBackgroundValues(normalisedBackgroundValues, totalCalibrationFactors)

originalMagnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(originalValuesCalibrated, arcsecPerPx)
normMagnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(valuesCalibrated, arcsecPerPx)


saveBackEvolution(originalMagnitudesPerArcSecSq,rejectedFrames_astrometry,rejectedFrames_FWHM, rejectedFrames_Background,\
    "Evolution of background magnitudes",destinationFolder+"/backgroundEvolution_originalMag",folderWithFramesWithAirmasses, maxmimumBackgroundBrightness)

saveBackEvolution(normMagnitudesPerArcSecSq,rejectedFrames_astrometry,rejectedFrames_FWHM, rejectedFrames_Background,\
    "Evolution of NORMALISED background magnitudes",destinationFolder+"/backgroundEvolution_normMag",folderWithFramesWithAirmasses, maxmimumBackgroundBrightness)

# x = [float(i) for i in valuesCalibrated]
# scatterPlotCountsVsMagnitudes(x, normMagnitudesPerArcSecSq, rejectedFrames_astrometry, rejectedFrames_FWHM, \
#                             destinationFolder + "/countsVsMagnitudes.png")
