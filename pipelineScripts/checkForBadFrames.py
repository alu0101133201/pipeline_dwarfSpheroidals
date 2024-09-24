import sys
import glob
import fnmatch

import numpy as np
import matplotlib.pyplot as plt

from astropy.stats import sigma_clipped_stats


def retrieveBackgroundValues(currentFile):
    backgroundValues = []
    backgroundStd    = []
    with open(currentFile, 'r') as f:
        for line in f:
            splittedLine = line.split()
            numberOfFields = len(splittedLine)

            if (numberOfFields == 3):
                backgroundValues.append(float(splittedLine[1]))
                backgroundStd.append(float(splittedLine[2]))
            else:
                raise Exception("Wrong number of fields in the file of background estimation. Expected 3 (constant estimation of the background)")
    return(backgroundValues, backgroundStd)

def computeMedianAndStd(values):
    mean, median, std = sigma_clipped_stats(values)
    return(median, std)

def saveHistogram(values, median, std, imageName, numOfStd):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_ylim(0, 25)
    ax.hist(values)

    ax.vlines(median, ymin=0, ymax=500, color="black", linestyle="--", linewidth=2.0)
    ax.vlines(median + numOfStd * std, ymin=0, ymax=500, color="red", linestyle="--", linewidth=2.0)
    ax.vlines(median - numOfStd * std, ymin=0, ymax=500, color="red", linestyle="--", linewidth=2.0)

    plt.savefig(imageName)
    return()

folderWithSkyEstimations  = sys.argv[1]
extensionToLookFor        = sys.argv[2]
outputFolder              = sys.argv[3]
outputFile                = sys.argv[4]
numberOfStdForRejecting    = float(sys.argv[5])
totalBackgroundValues = []

for currentFile in glob.glob(folderWithSkyEstimations + "/*" + extensionToLookFor):
    if (fnmatch.fnmatch(currentFile, '*/done*')): continue
    backgroundValue, _ = retrieveBackgroundValues(currentFile)

    totalBackgroundValues = np.concatenate((totalBackgroundValues, backgroundValue))


backgroundValueMean, backgroundValueStd = computeMedianAndStd(totalBackgroundValues)

saveHistogram(totalBackgroundValues, backgroundValueMean, backgroundValueStd, outputFolder + "/backgroundHist.png", numberOfStdForRejecting)

badFiles = []
for currentFile in glob.glob(folderWithSkyEstimations + "/*" + extensionToLookFor):
    if (fnmatch.fnmatch(currentFile, '*/done*')): continue
    backgroundValue, backgroundStd = retrieveBackgroundValues(currentFile)

    if ((backgroundValue[0] > (backgroundValueMean + (numberOfStdForRejecting * backgroundValueStd))) or \
        (backgroundValue[0] < (backgroundValueMean - (numberOfStdForRejecting * backgroundValueStd)))):
        badFiles.append(currentFile.split('/')[-1])

with open(outputFolder + "/" + outputFile, 'w') as file:
    for fileName in badFiles:
        file.write(fileName + '\n')

