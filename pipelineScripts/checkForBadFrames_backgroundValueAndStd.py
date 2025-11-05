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

def getBadAstrometrisedFrames(file):
    badAstrometrised = []
    with open(file, 'r') as f:
        for i in f:
            badAstrometrised.append(i)
    return(np.array(badAstrometrised))


def getFilenameWithPattern(folderPath, n):
    pattern1 = re.compile(rf"\bf{n}\b.*\.fits\b", re.IGNORECASE) 
    pattern2 = re.compile(rf"\b{n}\b.*\.fits\b", re.IGNORECASE)  

    for filename in os.listdir(folderPath):
        if pattern1.search(filename):
            return filename  

    for filename in os.listdir(folderPath):
        if pattern2.search(filename):
            return filename  

    return None  

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

        if (numberOfFields == 5):
            return(float(splittedLine[1]))
        elif (numberOfFields == 1):
            return(float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 5 (constant estimation of the background), got " + str(numberOfFields))

def readAirMassesFromFile(fileName):
    values = []
    try:
        with open(fileName, 'r') as file:
            for line in file:
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

def saveHistogram(values, badAstrometrisedIndices, imageName, title, labelX):
    valuesToPlot = values[~np.isnan(values)]
    myBins = calculateFreedmanBins(valuesToPlot)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.tight_layout(pad=7)
    configureAxis(ax, labelX, 'Number of frames', logScale=False)
    ax.set_title(title, fontsize=22, pad=17)

    counts, bins, patches = ax.hist(valuesToPlot, color="teal", bins=myBins)
    max_bin_height = counts.max() + 5
    ax.set_ylim(0, max_bin_height)

    # if (len(badAstrometrisedIndices) != 0):
    #     ax.hist(values[badAstrometrisedIndices], bins=myBins, color='blue', label='Rejected due to Astrometry')

    ax.legend(fontsize=18)
    plt.savefig(imageName)
    return()


def obtainNumberFromFrame(currentFile):
    # I know this code is awful but I use the function in two different places and the names are different
    # Since the names can have multiple numbers it's not trivial to identify the frame number so I have 
    # hardcoded the names used... 
    match = re.search(r'entirecamera_(\d+)\.txt', currentFile)

    if match:
        frameNumber = int(match.group(1))
        return(frameNumber)
    else:
        match = re.search(r'_f(\d+)_', currentFile)
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
    fitsFileNamePattern = getFilenameWithPattern(airMassesFolder, frameNumber)
    fitsFilePath = os.path.join(airMassesFolder, fitsFileNamePattern)
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
            backgroundStd   = float(splittedLine[2])
            backgroundSkew = float(splittedLine[3])
            backgroundKurto = float(splittedLine[4])
        elif (numberOfFields == 1):
            return(float('nan'), float('nan'), float('nan'), float('nan'), float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 5 (constant estimation of the background), got " + str(numberOfFields))

    # Then we read the airmass
    try:
        airmass = obtainAirmassFromFile(currentFile, folderWithAirMasses, airMassKeyWord)
    except:
        print("something went wrong in obtaining the airmass, returning nans (file: " + str(currentFile) + ")")
        return(float('nan'), float('nan'), float('nan'), float('nan'), float('nan')) 

    return(backgroundValue, backgroundValue / airmass, backgroundStd,backgroundSkew,backgroundKurto)

def saveParameterEvolution(files, values, parameter, imageName, astrometryRejectedIndices):
    fig, ax = plt.subplots(2,1, figsize=(20,10))
    configureAxis(ax[0], 'UTC', f'{parameter} (ADU)',logScale=False)
    configureAxis(ax[1], 'Airmass', f'{parameter} (ADU)',logScale=False)
    fig.suptitle(f'{parameter} evolution',fontsize=22)
    pattern=r"entirecamera_(\d+)"

    # Plot all the values
    for i in range(len(files)):
        file = files[i]
        match=re.search(pattern,file)
        frame=match.group(1)
        file=folderWithFramesWithAirmasses+'/'+frame+'.fits'
        date=obtainKeyWordFromFits(file,'DATE-OBS')
        air=obtainKeyWordFromFits(file,'AIRMASS')
        date_ok=datetime.fromisoformat(date)
        bck=values[i]

        ax[0].scatter(date_ok,bck,marker='o',s=50,edgecolor='black',color='teal',zorder=5)
        ax[1].scatter(air,bck,marker='o',s=50,edgecolor='black',color='teal',zorder=5)

    # astrometryRejectedValues        = [x for x in values[astrometryRejectedIndices]]  if len(astrometryRejectedIndices) > 0 else []
    # astrometryRejectedFiles         = [x for x in files[astrometryRejectedIndices]]   if len(astrometryRejectedIndices) > 0 else []

    # for j in range(len(astrometryRejectedFiles)):
    #     match=re.search(pattern, astrometryRejectedFiles[j])
    #     frame=match.group(1)
    #     file=folderWithFramesWithAirmasses+'/'+frame+'.fits'
    #     date=obtainKeyWordFromFits(file,'DATE-OBS')
    #     air=obtainKeyWordFromFits(file,'AIRMASS')
    #     date_ok=datetime.fromisoformat(date)
    #     # ax[0].scatter(date_ok, astrometryRejectedValues[j], facecolors='none', edgecolor='blue', lw=1.5, s=350, zorder=10, label='Rejected astrometry' if (j==0) else "")
    #     # ax[1].scatter(air, astrometryRejectedValues[j], facecolors='none', edgecolor='blue', lw=1.5, s=350, zorder=10, label='Rejected astrometry'if (j==0) else "")
    #     if j==0:
    #         ax[0].legend(fontsize=18)

    for label in ax[0].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')

    plt.tight_layout()
    plt.savefig(imageName)

    return()


def saveValuesVSStats(values, STD, Skew, kurto, backgroundRejectedIndices, stdRejectedIndices, badAstrometrisedIndices, imageName):
    allRejectedIndices = np.append(backgroundRejectedIndices, stdRejectedIndices)

    goodValues = [value for i, value in enumerate(values) if i not in allRejectedIndices]
    goodStd   = [value for i, value in enumerate(STD) if i not in allRejectedIndices]
    goodSkew  = [value for i, value in enumerate(Skew) if i not in allRejectedIndices]
    goodKurto = [value for i, value in enumerate(kurto) if i not in allRejectedIndices]

    valuesRejectedByBack = values[backgroundRejectedIndices]
    valuesRejectedByStd  = values[stdRejectedIndices]

    fig, ax = plt.subplots(1,3,figsize=(30,10))

    configureAxis(ax[0],r'$\sqrt{Background}$','STD',logScale=False)
    configureAxis(ax[1],'Background','Skewness',logScale=False)
    configureAxis(ax[2],'Background','Kurtosis',logScale=False)
    ax[1].set_yscale('log')
    ax[2].set_yscale('log')

    ax[0].scatter(np.sqrt(goodValues), goodStd, s=120,edgecolor='k',color='teal')
    ax[1].scatter(goodValues, goodSkew, s=120,edgecolor='k',color='mediumorchid')
    ax[2].scatter(goodValues, goodKurto,s=120,edgecolor='k',color='red')

    negSkew=[value<0 for value in goodSkew]; negKurto=[value<0 for value in goodKurto]
    if True in negSkew:
        negSkew_vals=[np.abs(value) for value, is_negative in zip(goodSkew,negSkew) if is_negative]
        negSkew_bckvals=[value for value,is_negative in zip(goodValues,negSkew) if is_negative]
        ax[1].scatter(negSkew_bckvals,negSkew_vals,marker='s',s=120,edgecolor='k',color='violet',label='Negative values')
    if True in negKurto:
        negKurto_vals=[np.abs(value) for value, is_negative in zip(goodKurto,negKurto) if is_negative]
        negKurto_bckvals=[value for value,is_negative in zip(goodValues,negKurto) if is_negative]
        ax[2].scatter(negKurto_bckvals,negKurto_vals,marker='s',s=120,edgecolor='k',color='lightsalmon',label='Negative values')
        ax[2].legend(fontsize=20)
    
    if (len(stdRejectedIndices) != 0):
        ax[0].scatter(np.sqrt(values)[stdRejectedIndices], STD[stdRejectedIndices], s=120, marker="D", edgecolor='k',color='gold')
        ax[1].scatter(values[stdRejectedIndices], np.abs(Skew[stdRejectedIndices]), s=120, marker="D", edgecolor='k',color='gold', label="Rejected by background std")
        ax[2].scatter(values[stdRejectedIndices], np.abs(kurto[stdRejectedIndices]), s=120, marker="D", edgecolor='k',color='gold')

    # if (len(badAstrometrisedIndices) != 0):
    #     ax[0].scatter(np.sqrt(values)[badAstrometrisedIndices], STD[badAstrometrisedIndices], s=350, lw=1.5, edgecolor='blue', facecolors='none')
    #     ax[1].scatter(values[badAstrometrisedIndices], np.abs(Skew[badAstrometrisedIndices]), s=350, lw=1.5, edgecolor='blue', facecolors='none', label="Rejected by astrometry")
    #     ax[2].scatter(values[badAstrometrisedIndices], np.abs(kurto[badAstrometrisedIndices]), s=350, lw=1.5, edgecolor='blue', facecolors='none')

    ax[1].legend(fontsize=20)
    plt.tight_layout()
    plt.savefig(imageName)

    return()

def getIndicesOfFiles(files, filesNames):
    indices = []

    for i in filesNames:
        for j in range(len(files)):
            if (i == files[j]):
                indices.append(j)
    return(np.array(indices))

def writeMetricToFile(outputPath, fileNames, values):
    if (len(fileNames) != len(values)):
        raise ValueError(f"Length mismatch: {len(fileNames)} filenames and {len(values)} values")
    
    with open(outputPath, 'w') as file:
        for name, val in zip(fileNames, values):
            file.write(f"{name} {val}\n")


HDU_TO_FIND_AIRMASS = 1

folderWithSkyEstimations      = sys.argv[1]
folderWithFramesWithAirmasses = sys.argv[2] # The airmasses are in the header of the fits files that are in this folder
airMassKeyWord                = sys.argv[3]
outputFolder                  = sys.argv[4]

setMatplotlibConf()

# 0.- Get frames identified as candidates for rejection in previous step (at this point, only astrometrisation)
basAstrometrisedFile = outputFolder + "/identifiedBadFrames_astrometry.txt"
badAstrometrisedFrames = getBadAstrometrisedFrames(basAstrometrisedFile)
badAstrometrisedFrames = [folderWithSkyEstimations + "/entirecamera_" + x[:-1] + ".txt" for x in badAstrometrisedFrames]

# 1.- Obtain the normalised background values and std values ------------------
files = []
originalBackgroundValues   = []
normalisedBackgroundValues = []
backgroundStds             = []
backgroundSkews            = []
backgroundKurtos           = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if fnmatch.fnmatch(currentFile, '*done*.txt'):
        continue
    
    files.append(currentFile)
    originalbackground, normalisedBackground, currentStd, currentSkew, currentKurto = obtainNormalisedBackground(currentFile, folderWithFramesWithAirmasses, airMassKeyWord)
    originalBackgroundValues.append(originalbackground)
    normalisedBackgroundValues.append(normalisedBackground)
    backgroundStds.append(currentStd)
    backgroundSkews.append(currentSkew)
    backgroundKurtos.append(currentKurto)

files = np.array(files)
originalBackgroundValues   = np.array(originalBackgroundValues)
normalisedBackgroundValues = np.array(normalisedBackgroundValues)
backgroundStds             = np.array(backgroundStds)
backgroundSkews            = np.array(backgroundSkews)
backgroundKurtos           = np.array(backgroundKurtos)

writeMetricToFile(outputFolder + "/stdValues.dat", files, backgroundStds)
writeMetricToFile(outputFolder + "/skewnessValues.dat", files, backgroundSkews)
writeMetricToFile(outputFolder + "/kurtosisValues.dat", files, backgroundKurtos)

badAstrometrisedIndices   = getIndicesOfFiles(files, badAstrometrisedFrames)

# 2.- Plots 
saveParameterEvolution(files, originalBackgroundValues, "Background", outputFolder+"/backgroundEvolution_original.png", badAstrometrisedIndices)
saveParameterEvolution(files, normalisedBackgroundValues, "Background", outputFolder+"/backgroundEvolution_normalised.png", badAstrometrisedIndices)
saveParameterEvolution(files, backgroundStds, "STD", outputFolder+"/stdEvolution.png", badAstrometrisedIndices)
saveParameterEvolution(files, backgroundSkews, "Skewness", outputFolder+"/skewnessEvolution.png", badAstrometrisedIndices)
saveParameterEvolution(files, backgroundKurtos, "Kurtosis", outputFolder+"/kurtosisEvolution.png", badAstrometrisedIndices)

saveValuesVSStats(normalisedBackgroundValues, backgroundStds, backgroundSkews, backgroundKurtos, [], [], badAstrometrisedIndices, outputFolder + "/backgroundStats.png")

saveHistogram(normalisedBackgroundValues, badAstrometrisedIndices,   outputFolder + "/backgroundHist.png",  "Background values normalised by the airmass", "Background counts (ADU)")
saveHistogram(backgroundStds, badAstrometrisedIndices,   outputFolder + "/backgroundStdHist.png",  "Background std values", "Background STD (ADU)")