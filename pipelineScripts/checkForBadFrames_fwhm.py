import re
import sys
import glob
import math
import astropy

import numpy as np
import matplotlib.pyplot as plt

from astropy.stats import sigma_clipped_stats


def retrieveFWHMValues(currentFile):
    with open(currentFile, 'r') as f:
        lines = f.readlines()
        if( len(lines) != 1):
            raise Exception("File with the FWHM estimation contains more that 1 line. Expected 1 line got " + str(len(lines)))
        
        splittedLine = lines[0].strip().split()
        numberOfFields = len(splittedLine)

        if (numberOfFields == 4):
            return(float(splittedLine[0]))
        elif (numberOfFields == 0):
            return(float('nan')) # Frame which has been lost in reduction (e.g. failed to astrometrise). Just jump to the next iteration
        else:
            raise Exception("Wrong number of fields in the file of background estimation. Expected 4 (constant estimation of the background), got " + str(numberOfFields))

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
    myBins = calculateFreedmanBins(values)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    counts, bins, patches = ax.hist(values, bins=myBins)
    max_bin_height = counts.max()
    ax.set_ylim(0, max_bin_height)

    ax.vlines(median, ymin=0, ymax=max_bin_height, color="black", linestyle="--", linewidth=2.0)
    ax.vlines(median + numOfStd * std, ymin=0, ymax=max_bin_height, color="red", linestyle="--", linewidth=2.0)
    ax.vlines(median - numOfStd * std, ymin=0, ymax=max_bin_height, color="red", linestyle="--", linewidth=2.0)

    plt.savefig(imageName)
    return()


folderWithFWHM            = sys.argv[1]
extensionToLookFor        = sys.argv[2]
outputFolder              = sys.argv[3]
outputFile                = sys.argv[4]
numberOfStdForRejecting    = float(sys.argv[5])


fwhmValues = np.array([])
for currentFile in glob.glob(folderWithFWHM + "/range1_*" + extensionToLookFor):
    fwhmValue = retrieveFWHMValues(currentFile)
    if (not math.isnan(fwhmValue)):
        fwhmValues = np.concatenate((fwhmValues, [fwhmValue]))


fwhmValueMean, fwhmValueStd = computeMedianAndStd(fwhmValues)
saveHistogram(fwhmValues, fwhmValueMean, fwhmValueStd, outputFolder + "/fwhmHist.png", numberOfStdForRejecting)


badFiles = []
for currentFile in glob.glob(folderWithFWHM + "/range1_*" + extensionToLookFor):
    fwhmValue = retrieveFWHMValues(currentFile)

    if (math.isnan(fwhmValue)):
        continue
    if ((fwhmValue > (fwhmValueMean + (numberOfStdForRejecting * fwhmValueStd))) or \
        (fwhmValue < (fwhmValueMean - (numberOfStdForRejecting * fwhmValueStd)))):
        badFiles.append(currentFile.split('/')[-1])

pattern = r"entirecamera_\d+"
with open(outputFolder + "/" + outputFile, 'w') as file:
    for fileName in badFiles:
        match = re.search(pattern, fileName)
        result = match.group()
        file.write(result + '\n')