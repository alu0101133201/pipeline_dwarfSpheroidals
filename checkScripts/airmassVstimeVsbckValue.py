import os
import glob
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits
from datetime import datetime
from matplotlib.ticker import ScalarFormatter

from basicPlotConf import *


def readDataAndAirmass(dataFolder):
    dateKeyWord = "DATE-OBS"
    airmassKeyWord = "AIRMASS"

    date_format = "%Y-%m-%dT%H:%M:%S.%f"
    dateAndAirMassPairs = []
    dateArray = []
    massArray = []
    dateAndAirMassPairsTmp = []

    fitsFiles = glob.glob(dataFolder + "/*.fits")

    for j in fitsFiles:
        tmpHeader = fits.open(j)[hduNum].header
        # The conversion to int and the [:26] in the next line are due to excesss of decimals in the seconds
        datetime_obj = datetime.strptime(tmpHeader[dateKeyWord][:26], date_format)

        unix_time = int(datetime_obj.timestamp())

        tmpPair = [unix_time, tmpHeader[airmassKeyWord]]
        dateAndAirMassPairsTmp.append(tmpPair)

    dateAndAirMassPairsTmp.sort(key=lambda x: x[0])
    # dateArray = [(datetime.utcfromtimestamp(i[0])).strftime("%m/%d-%H:%M") for i in dateAndAirMassPairsTmp] # To show the time labels with date format
    dateArray = [pair[0] for pair in dateAndAirMassPairsTmp] # To show the time labels in numerical format
    massArray = [pair[1] for pair in dateAndAirMassPairsTmp]
    return(dateArray, massArray)

def readBackgroundValues(directory):
    values = []

    numOfFiles = len(glob.glob(os.path.join(directory, 'entirecamera_*.txt')))
    for i in range(numOfFiles):
        with open(directory + "/entirecamera_" + str(i + 1) + ".txt", 'r') as f:
            for line in f:
                splittedLine = line.split()
                values.append(float(splittedLine[1]))
    return(np.array(values))

def comparisonPlots(time, airmass, backgroundValue):

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(25, 15))
    plt.tight_layout(pad=10)
    fig.subplots_adjust(left=0.1, right=0.9)

    ax1.set_position([0.1, 0.1, 0.25, 0.8]) 
    ax2.set_position([0.45, 0.1, 0.25, 0.8])
    ax3.set_position([0.7, 0.1, 0.25, 0.8])

    configureAxis(ax1, 'Time', 'Airmass', logScale=False)
    configureAxis(ax2, 'Time', 'Background', logScale=False)
    configureAxis(ax3, 'Airmass', 'Background', logScale=False)
    ax1.yaxis.set_inverted(True)
    plt.xticks(rotation=45, ha='left')

    ax1.xaxis.tick_top()
    ax1.tick_params(axis='x', which='major', labelsize=30, pad=17, rotation=45)

    ax2.xaxis.tick_top()
    ax2.tick_params(axis='x', which='major', labelsize=30, pad=17, rotation=45)
    ax2.tick_params(axis='y', which='minor', labelsize=20, pad=17)

    ax3.xaxis.tick_top()
    ax3.tick_params(axis='x', which='major', labelsize=30, pad=17, rotation=45)
    ax3.tick_params(axis='y', which='minor', labelsize=20, pad=17)

    ax2.set_yscale('log')
    ax3.set_yscale('log')

    ax1.set_xlabel("Time");    ax1.set_ylabel("Airmass")
    ax2.set_xlabel("Time");    ax2.set_ylabel("Background value")
    ax3.set_xlabel("Airmass"); ax3.set_ylabel("")

    ax1.scatter(time, airmass, color="teal", s=150, edgecolor='black', linewidth=2)
    ax2.scatter(time, backgroundValue, color = "teal", s=150, edgecolor='black', linewidth=2)
    ax3.scatter(airmass, backgroundValue, color = "teal", s=150, edgecolor='black', linewidth=2)

    ax2.yaxis.set_major_formatter(ScalarFormatter())
    ax2.yaxis.get_major_formatter().set_scientific(False)
    ax2.yaxis.set_minor_formatter(ScalarFormatter())
    ax2.yaxis.get_minor_formatter().set_scientific(False)

    ax3.tick_params(axis='y', which='both', labelleft=False)

    plt.savefig("airmassVstimeVsBckValue.png")
    return



parser = argparse.ArgumentParser(description='Script that shows the comparisons: (time-airmass; time-backgroundValue; Airmass-backgroundValue)\
 This is done only for one folder (in contrast to the "airmassVsTime" script), so do NOT provide the folder DATA or DATA-OR, but the directly night folder (i.e. DATA/nightX or DATA-or/nightX).\
 The folder with the background values must contain only the values for that night, so THE EASIEST WAY is to only run one night in the pipeline \
(thus having the DATA/night1 folder and the background folder ready to be used)')
parser.add_argument('dataPath', type=str, help='The path to directory that contains the fits files')
parser.add_argument('hdu', type=str, help='The hdu of the fits files that contains the DATE-OBS and AIRMASS')
parser.add_argument('backgroundValuesPath', type=str, help='The path to directory that contains the files with the background values')

args = parser.parse_args()

dataDirectory = args.dataPath
hduNum = int(args.hdu)
backgroundDirectory = args.backgroundValuesPath


if not os.path.exists(dataDirectory):
    raise ValueError(f"The path '{dataDirectory}' is not valid.")
if not os.path.exists(backgroundDirectory):
    raise ValueError(f"The path '{backgroundDirectory}' is not valid.")

timeArray, airMassArray = readDataAndAirmass(dataDirectory)
backgroundValues = readBackgroundValues(backgroundDirectory)

setMatplotlibConf()
comparisonPlots(timeArray, airMassArray, backgroundValues)