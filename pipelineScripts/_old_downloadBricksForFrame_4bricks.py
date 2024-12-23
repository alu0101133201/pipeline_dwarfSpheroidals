# This file will receive a frame, and based on the normalisation ring used for the reduction process
# will download 4 decals bricks which are located at the sides of the circle (up, down, left and right)


import sys
import threading 

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

from pipeline.pipelineScripts.GetAndDownloadBricks import *

def getImageData(image, hduNumber):
    hdulist = fits.open(image)
    header = hdulist[hduNumber].header
    imageData_original = hdulist[hduNumber].data
    wcs = WCS(header)
    hdulist.close()
    return(imageData_original, np.shape(imageData_original), wcs)

def getRingRadiusFromFile(path):
    with open(path, 'r') as file:
        firstLine = file.readline().strip()
        values = firstLine.split()
        if (len(values) >= 5):
            return int(values[4])
        else:
            raise Exception ("Error in 'getRingRadiusFromFile'. The ring file is not as expected.")
    return()


if (len(sys.argv) < 5):
    raise Exception("A frame path, ring path, array of filters and download destination has to be provided")

framePath = sys.argv[1]
ringPath = sys.argv[2]
filters = sys.argv[3].split(',')
downloadDestination = sys.argv[4]

radius = getRingRadiusFromFile(ringPath)
data, shape, wcs = getImageData(framePath, 1)
dataCentrePx = (int(shape[1]/2), int(shape[0]/2))

dataLocation1Px = (dataCentrePx[0], dataCentrePx[1] - radius)
dataLocation2Px = (dataCentrePx[0], dataCentrePx[1] + radius)
dataLocation3Px = (dataCentrePx[0] - radius, dataCentrePx[1])
dataLocation4Px = (dataCentrePx[0] + radius, dataCentrePx[1])

dataLocation1Wcs = wcs.pixel_to_world(dataLocation1Px[0], dataLocation1Px[1])
dataLocation2Wcs = wcs.pixel_to_world(dataLocation2Px[0], dataLocation2Px[1])
dataLocation3Wcs = wcs.pixel_to_world(dataLocation3Px[0], dataLocation3Px[1])
dataLocation4Wcs = wcs.pixel_to_world(dataLocation4Px[0], dataLocation4Px[1])

raList = [dataLocation1Wcs.ra.deg, dataLocation2Wcs.ra.deg, dataLocation3Wcs.ra.deg, dataLocation4Wcs.ra.deg]
decList = [dataLocation1Wcs.dec.deg, dataLocation2Wcs.dec.deg, dataLocation3Wcs.dec.deg, dataLocation4Wcs.dec.deg]

brickNames = getBrickNamesFromCoords(raList, decList)

threadList = []
for i in brickNames:
    threadList.append(threading.Thread(target=downloadBrick, args=(i, filters, downloadDestination, False)))

for i in threadList:
    i.start()

for i in threadList:
    i.join()

# This print is to be able to recover the bricknames from the bash pipeline
print(brickNames)
