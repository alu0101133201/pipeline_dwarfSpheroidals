import sys 
import glob

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

def readFits(fitsFile, hdu):
    with fits.open(fitsFile) as hdul:
        header = hdul[hdu].header
        data = hdul[hdu].data
        wcs = WCS(header)
    return(data, header, wcs)

def BrickInRange(raRange, decRange, currentRa, currentDec):
    frameRaMin = raRange[0]
    frameRaMax = raRange[1]
    frameDecMin = decRange[0]
    frameDecMax = decRange[1]

    decalsBrickHalfWidth = DECALS_BRICK_WIDTH_DEG / 2
    if ( (((currentRa - decalsBrickHalfWidth) < frameRaMax) and ((currentRa + decalsBrickHalfWidth) > frameRaMin)) and \
        (((currentDec - decalsBrickHalfWidth) < frameDecMax) and ((currentDec + decalsBrickHalfWidth) > frameDecMin)) ):
        return(True)
    return False

def checkBricksInRange(raRange, decRange, brickIdentificationFile):
    bricksInRange = []

    with open(brickIdentificationFile, 'r') as f:
        next(f)
        for line in f:
            splittedLine = line.split()
            currentBrick = splittedLine[0]
            currentRa    = float(splittedLine[1])
            currentDec   = float(splittedLine[2])

            if ( BrickInRange(raRange, decRange, currentRa, currentDec) ):
                bricksInRange.append(currentBrick)
                
    return(np.array(bricksInRange))

            
def getAssociatedBricksToFrame(frame, brickIdentificationFile, framesDataHdu):
    frameData, _, frameWcs = readFits(frame, framesDataHdu)
    ny, nx = frameData.shape

    corner1_coords = frameWcs.pixel_to_world(0, 0)
    corner2_coords = frameWcs.pixel_to_world(nx, ny)

    ra  = [corner1_coords.ra.deg, corner2_coords.ra.deg]
    dec = [corner1_coords.dec.deg, corner2_coords.dec.deg]

    frameRARange  = (np.min(ra), np.max(ra))
    frameDecRange = (np.min(dec), np.max(dec))

    bricksAssociated = checkBricksInRange(frameRARange, frameDecRange, brickIdentificationFile)
    return(bricksAssociated)

def writeBricksToFile(frame, bricks, outputFile):
    with open(outputFile, 'a') as f:
        f.write(frame + " ")
        for i in bricks:
            f.write(i + " ")
        f.write("\n")

DECALS_BRICK_WIDTH_DEG = 15.5 / 60 

framesDir               = sys.argv[1]
framesDataHdu           = int(sys.argv[2])
brickIdentificationFile = sys.argv[3]
outputFile              = sys.argv[4]
survey                  = sys.argv[5]
if survey=="DECaLS":
    DECALS_BRICK_WIDTH_DEG = 15.5 / 60
elif survey=="PANSTARRS":
    #In Panstarrs we are downloading 3600 pix = 900 asec = 0.25 deg
    DECALS_BRICK_WIDTH_DEG = 0.25

for currentFrame in glob.glob(framesDir + "/*.fits"):
    associatedBricks = getAssociatedBricksToFrame(currentFrame, brickIdentificationFile, framesDataHdu)
    writeBricksToFile(currentFrame, associatedBricks, outputFile)