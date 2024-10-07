import os
import sys
import glob
import datetime
import argparse
import numpy as np
import matplotlib.pyplot as plt

from astropy    import wcs
from astropy.io import fits
from datetime   import datetime
from astropy    import units as u

from basicPlotConf import *

def dateAirMassPlot(xData, yData):
    for i in range(len(xData)):
        fig, ax = plt.subplots(1, 1, figsize=(15, 15))
        plt.tight_layout(pad=6)
        plt.xticks(rotation=45, ha='left')

        configureAxis(ax, 'Time', 'Airmass', logScale=False)
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_inverted(True)
        ax.xaxis.tick_top()
        ax.tick_params(axis='x', which='major', labelsize=15, pad=17, rotation=45)

        ax.scatter(xData[i], yData[i], s=150, edgecolors='black', color="blue")
        plt.savefig("airMassVsTime_" + str(i) + ".png")
    return()


parser = argparse.ArgumentParser(description='Script that shows the Airmass vs Time.\
 It takes the data path (i.e. the directory DATA or DATA-or) and the hdu to look for the values and \
 reads the fits files looking for DATE-OBS and AIRMASS.\
 The comparison is done for all the "night*" folders that are contained in the data directory')
 
parser.add_argument('path', type=str, help='The path to directory that contains the night folders with the fits files (generally DATA or DATA-or)')
parser.add_argument('hdu', type=str, help='The hdu of the fits files that contains the DATE-OBS and AIRMASS keywords')

args = parser.parse_args()

dataDirectory = args.path
hduNum = int(args.hdu)

if not os.path.exists(dataDirectory):
    raise ValueError(f"The path '{dataDirectory}' is not valid.")

nightFolders = glob.glob(dataDirectory + "/night*")
numOfNights = len(nightFolders)

dateKeyWord = "DATE-OBS"
airmassKeyWord = "AIRMASS"

date_format = "%Y-%m-%dT%H:%M:%S.%f"
dateAndAirMassPairs = []
dateArray = []
massArray = []

for i in nightFolders:
    fitsFiles = glob.glob(i + "/*.fits")
    dateAndAirMassPairsTmp = []
    dateArrayTmp = []
    massArrayTmp = []

    for j in fitsFiles:
        tmpHeader = fits.open(j)[hduNum].header
        # The conversion to int and the [:26] in the next line are due to excesss of decimals in the seconds
        datetime_obj = datetime.strptime(tmpHeader[dateKeyWord][:26], date_format)
        unix_time = int(datetime_obj.timestamp())   
           
        tmpPair = [unix_time, tmpHeader[airmassKeyWord]]
        dateAndAirMassPairsTmp.append(tmpPair)

    dateAndAirMassPairsTmp.sort(key=lambda x: x[0])
    # dateArrayTmp = [(datetime.utcfromtimestamp(i)).strftime("%m/%d-%H:%M") for i in dateArrayTmp] # To show the time labels with date format
    dateArrayTmp = [pair[0] for pair in dateAndAirMassPairsTmp] # To show the time labels in numerical format
    massArrayTmp = [pair[1] for pair in dateAndAirMassPairsTmp]

    dateAndAirMassPairs.append(dateAndAirMassPairsTmp) 
    dateArray.append(dateArrayTmp)
    massArray.append(massArrayTmp) 


setMatplotlibConf()
dateAirMassPlot(dateArray, massArray)
