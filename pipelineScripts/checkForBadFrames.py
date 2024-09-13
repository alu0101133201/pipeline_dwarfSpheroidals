# This script is intented to detect possible bad frames that have been provided to the pipeline

# Since we have a huge number of frames, it could be useful to have some automatic way of detecting possible 
# bad frames. One approach is to check the estimation on the sky and see what frames deviate from the mean values

# When the background estimation is done by a constant, we have available the background value and the std of the sky. We use these two parameters for the analysis
# When the background estimation is done by a polynomial, we have available the coefficients of it. We use these parameters for the analysis

import sys
import glob
import fnmatch

import numpy as np
import matplotlib.pyplot as plt

from astropy.stats import sigma_clipped_stats

# Functions #
def retrieveDataFromFiles(currentFile, dataDict):
    with open(currentFile, 'r') as f:

        for line in f:
            splittedLine = line.split()
            numberOfFields = len(splittedLine)

            # There exists three options:
            #   1- The background is estimated by a constant. In this case the format of the data in the .txt file is:
            #        fileName backgroundValue stdValue
            #      This is the format because the constant background estimation is also used for the weighting of the frames
            #   2.- The background is estimated by a polynomial, and the .txt file contains a set of lines, each of which
            #       has the format key-value
            #   3.- The frame has failed to be astrometrised and there is only the name of the file (1 filed)

            if (numberOfFields == 3):
                if not "backgroundValue" in dataDict: dataDict["backgroundValue"] = []
                if not "backgroundStd" in dataDict: dataDict["backgroundStd"] = []

                dataDict["backgroundValue"].append((currentFile.split('/')[-1], float(splittedLine[1])))
                dataDict["backgroundStd"].append((currentFile.split('/')[-1], float(splittedLine[2])))

            elif (numberOfFields == 2):
                key   = splittedLine[0]
                value = float(splittedLine[1])

                if key in dataDict:
                    dataDict[key].append((currentFile.split('/')[-1], value))
                else:
                    dataDict[splittedLine[0]] = [(currentFile.split('/')[-1], value)]
            elif (numberOfFields == 1):
                return(dataDict) # I do not modify the dictionary, just jump to the next iteration
            else:
                raise Exception("Wrong number of fields in the file of background estimation. Expected 3 (constant estimation), 2 (plane/higher-order estimation) or 1 (frame failed to astrometrise)")
    return(dataDict)

def computeMedianAndStd(dataDict):
    medianValues = {}
    stdValues    = {}

    for currentKey in dataDict:
        currentValues = [value for filename, value in dataDict[currentKey]] 
        mean, median, std = sigma_clipped_stats(currentValues)
        medianValues[currentKey] = median
        stdValues[currentKey]    = std
    return(medianValues, stdValues)

def obtainPotentialBadFiles(dataDict, medianValues, stdValues):
    potentialBadFiles = []
    for currentKey in dataDict:
        median = medianValues[currentKey]
        std = stdValues[currentKey]

        currentValues = dataDict[currentKey]
        potentialBadFiles.append(currentKey + " | median: " + str(median) + " and std: " + str(std))

        for currentImageValue in currentValues:
            if (np.abs(currentImageValue[1] - median) >= (NUMBER_OF_STD_FOR_DETECTING_BAD_FRAMES*std)):
                potentialBadFiles.append(currentImageValue[0] + " due to parameter: " + currentKey)
    return(potentialBadFiles)

def writeArrayToFile(folderToWrite, potentialBadFiles):
    with open(folderToWrite + "/potentialBadFiles.txt", 'w') as file:
        for fileName in potentialBadFiles:
            file.write(fileName + '\n')

def createAndWriteHistograms(dataDict, outputFolder):
    for currentKey in dataDict.keys():
        dataToShow = [value for filename,value in dataDict[currentKey]]
        fig, ax = plt.subplots(1, 1, figsize=(10,10))
        ax.set_xlabel("Value of the parameter " + currentKey)
        ax.set_xlim(np.nanmin(dataToShow), np.nanmax(dataToShow))
        ax.hist(dataToShow)
        plt.savefig(outputFolder + "/" + currentKey + ".png")
# --- #

NUMBER_OF_STD_FOR_DETECTING_BAD_FRAMES = 1.5

folderWithSkyEstimations  = sys.argv[1]
extensionToLookFor        = sys.argv[2]
outputFolder              = sys.argv[3]

dataDictionary = {}

for currentFile in glob.glob(folderWithSkyEstimations + "/*" + extensionToLookFor):
    # We ignore the "done" file of the folder (it is simply a flag)
    if (fnmatch.fnmatch(currentFile, '*/done*')): continue
    dataDictionary = retrieveDataFromFiles(currentFile, dataDictionary)


medianValues, stdValues = computeMedianAndStd(dataDictionary)
potentialBadFiles = obtainPotentialBadFiles(dataDictionary, medianValues, stdValues)

writeArrayToFile(outputFolder, potentialBadFiles)
createAndWriteHistograms(dataDictionary, outputFolder)

