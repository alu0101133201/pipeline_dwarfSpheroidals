import re
import sys
import glob
import math
import fnmatch
import astropy

import numpy as np
import matplotlib.pyplot as plt

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

def saveHistogram(values, median, std, imageName, numOfStd):

    valuesToPlot = values[~np.isnan(values)]
    myBins = calculateFreedmanBins(valuesToPlot)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    counts, bins, patches = ax.hist(valuesToPlot, bins=myBins)
    max_bin_height = counts.max()
    ax.set_ylim(0, max_bin_height)

    ax.vlines(median, ymin=0, ymax=max_bin_height, color="black", linestyle="--", linewidth=2.0)
    ax.vlines(median + numOfStd * std, ymin=0, ymax=max_bin_height, color="red", linestyle="--", linewidth=2.0)
    ax.vlines(median - numOfStd * std, ymin=0, ymax=max_bin_height, color="red", linestyle="--", linewidth=2.0)

    plt.savefig(imageName)
    return()

def removeDoneFileFromList(files):
    validFiles = []
    for f in files:
        currentFile = f.split('/')[-1]
        if re.match(r"done_.*\.txt", currentFile):
            continue
        validFiles.append(currentFile)
    return(validFiles)

folderWithSkyEstimations  = sys.argv[1]
outputFolder              = sys.argv[2]
outputFile                = sys.argv[3]
numberOfStdForRejecting    = float(sys.argv[4])
totalBackgroundValues = []



# Iterate through the sky estimation
# Look for the frame that corresponds in "framesForCommonReduction" and obtain its airmass
# normalise





exit()
# When I manage to get the normalised background the following code should work fine

backgroundValueMean, backgroundValueStd = computeMedianAndStd(totalBackgroundValues)
saveHistogram(totalBackgroundValues, backgroundValueMean, backgroundValueStd, outputFolder + "/backgroundHist.png", numberOfStdForRejecting)

badFiles = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*.txt"):
    if (fnmatch.fnmatch(currentFile, '*/done*')): continue
    backgroundValue = retrieveBackgroundValues(currentFile)

    if (math.isnan(backgroundValue)):
        continue
    if ((backgroundValue > (backgroundValueMean + (numberOfStdForRejecting * backgroundValueStd))) or \
        (backgroundValue < (backgroundValueMean - (numberOfStdForRejecting * backgroundValueStd)))):
        badFiles.append(currentFile.split('/')[-1])

with open(outputFolder + "/" + outputFile, 'w') as file:
    for fileName in badFiles:
        file.write(fileName + '\n')

