"""
This python script is for analysing the different calibrated background between CCDs on a multi-detector array.
The script will:
- Plot the background evolution in magnitudes per arsec^2 for each CCD
- Given a refference CCD, will compute the ratio to multiply each other CCD to mathc the refference sky
- Save this ratio
- Plot the background evolution after applying this ratio

"""
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
from astropy.time import Time
import random

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

def saveScatterPlot(data,parameter,ccd_ref, title, outputFile):
    time = []
    
    for i in normalisedBackgroundValues:
        
        if (not pd.isna(i[0])):
            match=re.search(r"_(\d+).",i[0])
            frame = match.group(1)
            file=folderWithFramesWithAirmasses+'/'+str(frame)+'.fits'
            date=obtainKeyWordFromFits(file,dateKey)
            if dateKey.startswith("DATE"):
                date_ok=datetime.fromisoformat(date)
            elif dateKey.startswith("MJD"):
                date_ok=Time(date,format='mjd').to_datetime() 
            
            time.append(date_ok)
            
        else:
            
            time.append(np.nan)
    time = np.array(time)
    
    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'H', 'X', 'd']
    fig , ax = plt.subplots(2,1,figsize=(20,10))
    configureAxis(ax[0],'','',logScale=False)
    configureAxis(ax[1],'UTC','',logScale=False)
    fig.suptitle(title, fontsize=22)
    pattern=r"(\d+).fits"
    for h in range(1, num_ccd + 1):
        timeMask = ~pd.isna(time) & ~pd.isna(data[:,h-1])
        marker=random.choice(markers)
        ax[0].scatter(time[timeMask], data[:,h-1][timeMask], marker=marker,s=50, alpha=0.7)
        label_legend = f"CCD {h} (ref.)" if h == ccd_ref else f"CCD {h}"
        ax[1].scatter(time[timeMask], data[:,h-1][timeMask]-2.5*np.log10(parameter[h-1]),label=label_legend, marker=marker,s=50, alpha=0.7)
    ax[1].legend(loc="best", fontsize=20)
    for label in ax[0].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')
    for label in ax[1].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')
    ax[0].set_title('Before gain correction', fontsize=20)
    ax[1].set_title('After gain correction', fontsize=20)
    fig.supylabel(r'Back. [mag arcsec$^{-2}$]',fontsize=30)
    plt.tight_layout()
    plt.savefig(outputFile)
    return()

HDU_TO_FIND_AIRMASS = 0

folderWithSkyEstimations = sys.argv[1]
destinationFolder = sys.argv[2]
commonCalibrationFactorFile = sys.argv[3]
arcsecPerPx = float(sys.argv[4])
dateKey = sys.argv[5]
airMassKeyWord = sys.argv[6]
folderWithFramesWithAirmasses = sys.argv[7]
num_ccd = int(sys.argv[8])
ccd_ref=int(sys.argv[9])  
outputRatioFile=sys.argv[10]
setMatplotlibConf()
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
   
 # 2.- Retrieve the calibration factors
commonCalibrationFactorValue=np.array([np.nan] * num_ccd, dtype=float)
with open(commonCalibrationFactorFile) as f:
   fileContent=f.read().split()
   for h in range(num_ccd):
       commonCalibrationFactorValue[h] = float(fileContent[2*h])

arrayWithCommonFactor = [[x[0], [cfac for cfac in commonCalibrationFactorValue]] for x in normalisedBackgroundValues]
valuesCalibratedOriginal = applyCalibrationFactorsToBackgroundValues(originalBackgroundValues, arrayWithCommonFactor)
valuesCalibratedNormalised = applyCalibrationFactorsToBackgroundValues(normalisedBackgroundValues, arrayWithCommonFactor)

magnitudesPerArcSecSqOriginal = countsToSurfaceBrightnessUnits(valuesCalibratedOriginal, arcsecPerPx)
magnitudesPerArcSecSqNormalised = countsToSurfaceBrightnessUnits(valuesCalibratedNormalised, arcsecPerPx)

diffs=np.zeros((totalNumberOfFrames,num_ccd), dtype=float)

for h in range(num_ccd):
    diffs[:,h]= magnitudesPerArcSecSqNormalised[:,ccd_ref-1] - magnitudesPerArcSecSqNormalised[:,h]

# Now we save the mean results
collapsed = np.mean(diffs, axis=0)
parameter=10**(-0.4*collapsed)
np.savetxt(outputRatioFile, parameter, fmt="%.12f")
saveScatterPlot(magnitudesPerArcSecSqNormalised,parameter,ccd_ref,"Evolution of Normalised Background magnitude per ccd",destinationFolder+"/backgroundCCDcomparison.png")