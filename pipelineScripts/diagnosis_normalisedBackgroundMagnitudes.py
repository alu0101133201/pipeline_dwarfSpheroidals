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

        if (numberOfFields == 5):
            backgroundValue = float(splittedLine[1])
        elif (numberOfFields == 1):
            return(float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 5 (constant estimation of the background), got " + str(numberOfFields))

    # Then we read the airmass
    try:
        airmass = obtainAirmassFromFile(currentFile, folderWithAirMasses, airMassKeyWord)
    except:
        print("Something went wrong in obtaining the airmass, returning nans (file " + str(currentFile) + ")")
        return(float('nan'), float('nan')) 
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
        bins = [np.nanmin(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= np.nanmax(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

def saveHistogram(values, rejectedAstrometryIndices,  rejectedBackgroundValueIndices, rejectedBackgroundStdIndices, rejectedFWHMIndices, title, imageName):
    myBins = calculateFreedmanBins(values[~np.isnan(values)])

    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_title(title, fontsize=22, pad=17)
    plt.tight_layout(pad=7.0)
    configureAxis(ax, 'Background (mag/arcsec^2)', '', logScale=False)
    counts, bins, patches = ax.hist(values, bins=myBins, color="teal")

    ax.hist(values[rejectedBackgroundValueIndices], bins=myBins, color="darkred", label="Rejected by background value")
    ax.hist(values[rejectedBackgroundStdIndices], bins=myBins, color="gold", label="Rejected by background std")
    ax.hist(values[rejectedFWHMIndices], bins=myBins, color="mediumorchid", label="Rejected by fwhm")
    ax.hist(values[rejectedAstrometryIndices], bins=myBins, color="blue", label="Rejected by astrometry")

    max_bin_height = counts.max() + 5
    ax.set_ylim(0, max_bin_height)
    ax.legend(fontsize=20, loc="upper right")
    plt.savefig(imageName)
    return()

def saveScatterFactors(factors, rejectedAstrometryIndices, rejectedBackgroundValueIndices, rejectedBackgroundStdIndices, rejectedFWHMIndices, title, imageName, folderWithFramesWithAirmasses):
    airMass  = []
    time     = []
    cfactors = []

    for i in factors:
        match=re.search(r"_(\d+).",i[0])
        frame = match.group(1)
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,'DATE-OBS')
        air=obtainKeyWordFromFits(file,'AIRMASS')
        date_ok=datetime.fromisoformat(date)
        cfactors.append(i[1])
        airMass.append(air)
        time.append(date_ok)

    cfactors=np.array(cfactors,dtype='float')
    airMass=np.array(airMass,dtype='float')
    time = np.array(time)

    fig, ax = plt.subplots(2, 1, figsize=(20,10))
    configureAxis(ax[0], 'UTC', 'Calibration Factor',logScale=False)
    configureAxis(ax[1], 'Airmass', 'Calibration Factor',logScale=False)
    fig.suptitle('Calibration factor evolution',fontsize=22)
    pattern=r"(\d+).fits"

    ax[0].scatter(time,cfactors,marker='o',s=50,edgecolor='black',color='teal',zorder=0)
    ax[1].scatter(airMass,cfactors,marker='o',s=50,edgecolor='black',color='teal',zorder=0)

    ax[0].scatter(time[rejectedAstrometryIndices], cfactors[rejectedAstrometryIndices], facecolors='none', lw=1.5, edgecolor='blue', s=350,zorder=1, label="Rejected by astrometry")
    ax[0].scatter(time[rejectedBackgroundValueIndices], cfactors[rejectedBackgroundValueIndices], marker='X',edgecolor='k',color='darkred',s=120,zorder=1, label="Rejected by background")
    ax[0].scatter(time[rejectedBackgroundStdIndices], cfactors[rejectedBackgroundStdIndices], marker='D',edgecolor='k',color='gold',s=120,zorder=1, label="Rejected by std")
    ax[0].scatter(time[rejectedFWHMIndices], cfactors[rejectedFWHMIndices], marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=1, label="Rejected by FWHM")

    ax[1].scatter(airMass[rejectedAstrometryIndices], cfactors[rejectedAstrometryIndices], facecolors='none', lw=1.5, edgecolor='blue', s=350,zorder=1, label="Rejected by astrometry")
    ax[1].scatter(airMass[rejectedBackgroundValueIndices], cfactors[rejectedBackgroundValueIndices], marker='X',edgecolor='k',color='darkred',s=120,zorder=1, label="Rejected by background")
    ax[1].scatter(airMass[rejectedBackgroundStdIndices], cfactors[rejectedBackgroundStdIndices], marker='D',edgecolor='k',color='gold',s=120,zorder=1, label="Rejected by std")
    ax[1].scatter(airMass[rejectedFWHMIndices], cfactors[rejectedFWHMIndices], marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=1, label="Rejected by FWHM")

    ax[0].legend(loc="upper right", fontsize=20)
    
    for label in ax[0].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')
    plt.tight_layout()
    plt.savefig(imageName)
    return()

def saveBackEvolution(normalisedBackgroundValues,factors, rejectedAstrometryIndices, rejectedBackgroundValueIndices, rejectedBackgroundStdIndices, rejectedFWHMIndices, title, imageName, folderWithFramesWithAirmasses):
    airMass  = []
    time     = []
    

    for i in normalisedBackgroundValues:
        match=re.search(r"_(\d+).",i[0])
        frame = match.group(1)
        file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
        date=obtainKeyWordFromFits(file,'DATE-OBS')
        air=obtainKeyWordFromFits(file,'AIRMASS')
        date_ok=datetime.fromisoformat(date)
        
        airMass.append(air)
        time.append(date_ok)

    
    airMass=np.array(airMass,dtype='float')
    time = np.array(time)
    valuesCalibrated = applyCalibrationFactorsToBackgroundValues(normalisedBackgroundValues, factors)

    magnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(valuesCalibrated, arcsecPerPx)
    fig, ax = plt.subplots(2, 1, figsize=(20,10))
    configureAxis(ax[0], 'UTC', '',logScale=False)
    configureAxis(ax[1], 'Airmass', '',logScale=False)
    fig.suptitle(title,fontsize=22)
    pattern=r"(\d+).fits"

    ax[0].scatter(time,magnitudesPerArcSecSq,marker='o',s=50,edgecolor='black',color='teal',zorder=0)
    ax[1].scatter(airMass,magnitudesPerArcSecSq,marker='o',s=50,edgecolor='black',color='teal',zorder=0)

    ax[0].scatter(time[rejectedAstrometryIndices], magnitudesPerArcSecSq[rejectedAstrometryIndices], facecolors='none', lw=1.5, edgecolor='blue', s=350,zorder=1, label="Rejected by astrometry")
    ax[0].scatter(time[rejectedBackgroundValueIndices], magnitudesPerArcSecSq[rejectedBackgroundValueIndices], marker='X',edgecolor='k',color='darkred',s=120,zorder=1, label="Rejected by background")
    ax[0].scatter(time[rejectedBackgroundStdIndices], magnitudesPerArcSecSq[rejectedBackgroundStdIndices], marker='D',edgecolor='k',color='gold',s=120,zorder=1, label="Rejected by std")
    ax[0].scatter(time[rejectedFWHMIndices], magnitudesPerArcSecSq[rejectedFWHMIndices], marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=1, label="Rejected by FWHM")

    ax[1].scatter(airMass[rejectedAstrometryIndices], magnitudesPerArcSecSq[rejectedAstrometryIndices], facecolors='none', lw=1.5, edgecolor='blue', s=350,zorder=1, label="Rejected by astrometry")
    ax[1].scatter(airMass[rejectedBackgroundValueIndices], magnitudesPerArcSecSq[rejectedBackgroundValueIndices], marker='X',edgecolor='k',color='darkred',s=120,zorder=1, label="Rejected by background")
    ax[1].scatter(airMass[rejectedBackgroundStdIndices], magnitudesPerArcSecSq[rejectedBackgroundStdIndices], marker='D',edgecolor='k',color='gold',s=120,zorder=1, label="Rejected by std")
    ax[1].scatter(airMass[rejectedFWHMIndices], magnitudesPerArcSecSq[rejectedFWHMIndices], marker='P',edgecolor='k',color='mediumorchid',s=120,zorder=1, label="Rejected by FWHM")

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

def scatterPlotCountsVsMagnitudes(backgroundCounts, magnitudesPerArcSecSq, rejectedAstrometryIndices, rejectedBackgroundValueIndices, rejectedBackgroundStdIndices, rejectedFWHMIndices, fileName):
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_title("If calibration is correct, this should be a straight line", fontsize=20, pad=17)
    configureAxis(ax, 'Background (mag/arcsec^2)', 'Log(Background) (ADU)', logScale=False)

    plt.tight_layout(pad=6.0)
    ax.scatter(magnitudesPerArcSecSq, np.log10(backgroundCounts), s=40, color="teal")

    ax.scatter(magnitudesPerArcSecSq[rejectedAstrometryIndices], np.log10(backgroundCounts)[rejectedAstrometryIndices], facecolors='none', lw=1.5, edgecolor='blue',s=120, label="Rejected astrometry")
    ax.scatter(magnitudesPerArcSecSq[rejectedBackgroundValueIndices], np.log10(backgroundCounts)[rejectedBackgroundValueIndices], s=40, color="darkred", label="Rejected background value")
    ax.scatter(magnitudesPerArcSecSq[rejectedBackgroundStdIndices], np.log10(backgroundCounts)[rejectedBackgroundStdIndices], s=40, color="gold", label="Rejected background std")
    ax.scatter(magnitudesPerArcSecSq[rejectedFWHMIndices], np.log10(backgroundCounts)[rejectedFWHMIndices], s=40, color="mediumorchid", label="Rejected FWHM")

    ax.legend(fontsize=18, loc="upper right")
    plt.savefig(fileName)

def getIndicesOfRejectedFrames(normalisedBackgroundValuesArray, rejectedFrames):
    indices = []
    pattern = r"\d+"

    for j in rejectedFrames:
        for i in range(len(normalisedBackgroundValuesArray)):
            currentFrameName = normalisedBackgroundValuesArray[i][0]
            match = re.search(pattern, currentFrameName)
            result = match.group(0)
            if (result == j):
                indices.append(i)
    return(indices)


HDU_TO_FIND_AIRMASS = 1

folderWithSkyEstimations      = sys.argv[1]
folderWithFramesWithAirmasses = sys.argv[2] # The airmasses are in the header of the fits files that are in this folder
airMassKeyWord                = sys.argv[3]
folderWithCalibrationFactors  = sys.argv[4]
arcsecPerPx                   = float(sys.argv[5])
destinationFolder             = sys.argv[6]


setMatplotlibConf()


# 0.- Identify the files that have been identified as bad frames 
# We have bad frames due to:
#   Astrometry
#   Background value
#   Background Std
#   FWHM

rejectedFrames_astrometry = []
with open(destinationFolder + "/identifiedBadFrames_astrometry.txt", 'r') as f:
    for line in f:
        rejectedFrames_astrometry.append(line.strip())

rejectedFrames_backgroundValue = []
with open(destinationFolder + "/identifiedBadFrames_backgroundValue.txt", 'r') as f:
    for line in f:
        rejectedFrames_backgroundValue.append(line.strip())
 
rejectedFrames_backgroundStd = []
with open(destinationFolder + "/identifiedBadFrames_backgroundStd.txt", 'r') as f:
    for line in f:
        rejectedFrames_backgroundStd.append(line.strip())

rejectedFrames_FWHM = []
with open(destinationFolder + "/identifiedBadFrames_fwhm.txt", 'r') as f:
    for line in f:
        rejectedFrames_FWHM.append(line.strip())


# 1.- Obtain the normalised background values ------------------
normalisedBackgroundValues = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if fnmatch.fnmatch(currentFile, '*done*.txt'):
        continue
    currentValue = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord)
    normalisedBackgroundValues.append([currentFile.split('/')[-1], currentValue])


# 2.- Obtain the calibration factors
totalCalibrationFactors = []
for currentFile in glob.glob(folderWithCalibrationFactors + "/alpha_*Decals*.txt"):
    calibrationFactor = retrieveCalibrationFactors(currentFile)
    totalCalibrationFactors.append([currentFile.split('/')[-1], calibrationFactor])


rejectedAstrometryIndices      = getIndicesOfRejectedFrames(normalisedBackgroundValues, rejectedFrames_astrometry)
rejectedBackgroundValueIndices = getIndicesOfRejectedFrames(normalisedBackgroundValues, rejectedFrames_backgroundValue)
rejectedBackgroundStdIndices   = getIndicesOfRejectedFrames(normalisedBackgroundValues, rejectedFrames_backgroundStd)
rejectedFWHMIndices            = getIndicesOfRejectedFrames(normalisedBackgroundValues, rejectedFrames_FWHM)


valuesCalibrated = applyCalibrationFactorsToBackgroundValues(normalisedBackgroundValues, totalCalibrationFactors)

magnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(valuesCalibrated, arcsecPerPx)

saveHistogram(np.array(magnitudesPerArcSecSq), rejectedAstrometryIndices, rejectedBackgroundValueIndices, rejectedBackgroundStdIndices, rejectedFWHMIndices, \
                "Distribution of NORMALISED background magnitudes", destinationFolder + "/magnitudeHist.png")

saveScatterFactors(totalCalibrationFactors, rejectedAstrometryIndices, rejectedBackgroundValueIndices, rejectedBackgroundStdIndices, rejectedFWHMIndices, \
                "Evolution of calibration factors",destinationFolder + "/calibrationFactorEvolution.png", folderWithFramesWithAirmasses)

saveBackEvolution(normalisedBackgroundValues,totalCalibrationFactors,rejectedAstrometryIndices,rejectedBackgroundValueIndices, rejectedBackgroundStdIndices,rejectedFWHMIndices, \
    "Evolution of NORMALISED background magnitudes",destinationFolder+"/backgroundEvolution_magnitudes",folderWithFramesWithAirmasses)

x = [float(i) for i in valuesCalibrated]
scatterPlotCountsVsMagnitudes(x, magnitudesPerArcSecSq, rejectedAstrometryIndices, rejectedBackgroundValueIndices, rejectedBackgroundStdIndices, rejectedFWHMIndices, \
                            destinationFolder + "/countsVsMagnitudes.png")
