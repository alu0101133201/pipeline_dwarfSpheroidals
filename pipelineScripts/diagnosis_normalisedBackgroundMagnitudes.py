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
        bins = [min(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= max(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

def saveHistogram(values,rvalues_bck,rvalues_fwhm, title, imageName):
    valuesToPlot = values[~np.isnan(values)]
    
    rvaluesToPlot_fwhm=rvalues_fwhm[~np.isnan(rvalues_fwhm)]
    myBins = calculateFreedmanBins(valuesToPlot)
    myBins = np.linspace(np.nanmin(values), np.nanmax(values), 10)
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_title(title, fontsize=22, pad=17)
    plt.tight_layout(pad=7.0)
    configureAxis(ax, 'Background (mag/arcsec^2)', '', logScale=False)
    counts, bins, patches = ax.hist(values, bins=myBins, color="teal",label='Used frames')

    ##Check and hist rejected values
    if len(rvalues_bck)!=0:
        rvaluesToPlot_bck=rvalues_bck[~np.isnan(rvalues_bck)]
       #myBins_bck = calculateFreedmanBins(rvaluesToPlot_bck)
        myBins_bck = np.linspace(np.nanmin(rvalues_bck), np.nanmax(rvalues_bck), 10)
        counts_bck,bins_bck,patches_bck = ax.hist(rvalues_bck,bins=myBins_bck,color='darkred',label='Rejected Back.')
    
    if len(rvalues_fwhm)!=0:
        rvaluesToPlot_fwhm=rvalues_fwhm[~np.isnan(rvalues_fwhm)]
        #myBins_fwhm = calculateFreedmanBins(rvaluesToPlot_fwhm)
        myBins_fwhm = np.linspace(np.nanmin(rvalues_fwhm), np.nanmax(rvalues_fwhm), 10)
        counts_bck,bins_bck,patches_bck = ax.hist(rvalues_fwhm,bins=myBins_fwhm,color='mediumorchid',label='Rejected FWHM')
    max_bin_height = counts.max() + 10
    ax.legend()
    plt.savefig(imageName)
    return()

def saveScatterFactors(factors,rfactors_bck,rfactors_fwhm,title,imageName):
    frames=[]
    cfactors=[]
    for i in factors:
        match=re.search(r"_(\d+).",i[0])
        frames.append(match.group(1))
        cfactors.append(i[1])
    frames=np.array(frames,dtype='float')
    cfactors=np.array(cfactors,dtype='float')
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    ax.set_title(title,fontsize=22,pad=17)
    
    configureAxis(ax,'Number of Frame','Calibration Factor',logScale=False)
    ax.scatter(frames,cfactors,marker='o',edgecolor='k',color='teal',s=60,zorder=0)
    if len(rfactors_bck)!=0:
        frames_bck=[]; cfactors_bck=[]
        for i in rfactors_bck:
            match=re.search(r"_(\d+).",i[0])
            frames_bck.append(match.group(1))
            cfactors_bck.append(i[1])
        frames_bck=np.array(frames_bck,dtype='float')
        cfactors_bck=np.array(cfactors_bck,dtype='float')
        ax.scatter(frames_bck,cfactors_bck,marker='X',edgecolor='k',color='darkred',s=80,zorder=1,label='Rejected back.')
    if len(rfactors_fwhm)!=0:
        frames_fwhm=[]; cfactors_fwhm=[]
        for i in rfactors_fwhm:
            match=re.search(r"_(\d+).",i[0])
            frames_fwhm.append(match.group(1))
            cfactors_fwhm.append(i[1])
        frames_fwhm=np.array(frames_fwhm,dtype='float')
        cfactors_fwhm=np.array(cfactors_fwhm,dtype='float')
        ax.scatter(frames_fwhm,cfactors_fwhm,marker='P',edgecolor='k',color='mediumorchid',s=80,zorder=1,label='Rejected FWHM')

    ax.legend()
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

def scatterPlotCountsVsMagnitudes(backgroundCounts, magnitudesPerArcSecSq, fileName):
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    ax.set_title("If calibration is well-done it should be a straight line", fontsize=20, pad=17)
    configureAxis(ax, 'Background (mag/arcsec^2)', 'Log(Background) (ADU)', logScale=False)

    plt.tight_layout(pad=8.0)
    ax.scatter(magnitudesPerArcSecSq, np.log10(backgroundCounts), s=25, color="teal")
    plt.savefig(fileName)

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
# In order to trace reasons concerning the removal, we're gonna store it on different lists
badFrames_background = []
for currentFile in glob.glob(rejectedFramesBackground + "/*.fits"):
    match = re.search(r'entirecamera_(\d+)\.fits', currentFile)
    if match:
        number = int(match.group(1))  
        badFrames_background.append(number)
    else:
        raise Exception("Error identifying the number of the bad frames (background bad frames)")
badFrames_fwhm = []
for currentFile in glob.glob(rejectedFramesFWHM + "/*.fits"):
    match = re.search(r'entirecamera_(\d+)\.fits', currentFile)
    if match:
        number = int(match.group(1))  
        badFrames_fwhm.append(number)
    else:
        raise Exception("Error identifying the number of the bad frames (FWHM bad frames)")



setMatplotlibConf()

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
    if (not math.isnan(calibrationFactor)):
        totalCalibrationFactors.append([currentFile.split('/')[-1], calibrationFactor])

values,rembck,remfwhm=removeBadFramesFromList(normalisedBackgroundValues,badFrames_background,badFrames_fwhm)
totalCalibrationFactors_rem,removedFactors_bck,removedFactors_fwhm  = removeBadFramesFromList(totalCalibrationFactors, badFrames_background,badFrames_fwhm)

values_cal = applyCalibrationFactorsToBackgroundValues(values, totalCalibrationFactors_rem)
magnitudesPerArcSecSq = countsToSurfaceBrightnessUnits(values_cal, arcsecPerPx)
if len(rembck)!=0:
    rembck_cal=applyCalibrationFactorsToBackgroundValues(rembck,removedFactors_bck)
    magnitudesPerArcSecSq_rembck=countsToSurfaceBrightnessUnits(rembck_cal,arcsecPerPx)
else:
    magnitudesPerArcSecSq_rembck=[]
if len(remfwhm)!=0:
    remfwhm_cal=applyCalibrationFactorsToBackgroundValues(remfwhm,removedFactors_fwhm)
    magnitudesPerArcSecSq_remfwhm=countsToSurfaceBrightnessUnits(remfwhm_cal,arcsecPerPx)
else:
    magnitudesPerArcSecSq_remfwhm=[]

saveHistogram(np.array(magnitudesPerArcSecSq),np.array(magnitudesPerArcSecSq_rembck),np.array(magnitudesPerArcSecSq_remfwhm), "Distribution of NORMALISED background magnitudes", destinationFolder + "/magnitudeHist.png")
saveScatterFactors(totalCalibrationFactors,removedFactors_bck,removedFactors_fwhm,"Evolution of calibration factors",destinationFolder + "/calibfacEvol.png")
x = [float(i[1]) for i in values]
scatterPlotCountsVsMagnitudes(x, magnitudesPerArcSecSq, destinationFolder + "/countsVsMagnitudes.png")