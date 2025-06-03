# This file will receive a frame, and based on the normalisation ring used for the reduction process
# will download 4 decals bricks which are located at the sides of the circle (up, down, left and right)


import sys
import threading 

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from matplotlibConf import *
from GetAndDownloadBricks import *
import pandas as pd
import os

def getCoordsAndMagFromBrightStars(catalogue, threshold):
    ra = []
    dec = []
    mag_g = []
    with fits.open(catalogue) as hdul:
        data = hdul[1].data
        for row in data:
            if (row[2] < threshold):
                ra.append(row[0])
                dec.append(row[1])
                mag_g.append(row[2])
    return( np.array(ra), np.array(dec) )

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

def getWCSCoordinatesOfFrameCorners(file):
    with fits.open(file) as hdul:
        header = hdul[1].header
        wcs = WCS(header)
        dataShape = hdul[1].data.shape
    corner1 = (0, 0)
    corner2 = (dataShape[0] - 1, dataShape[1] -1)
    cornersCoords = wcs.pixel_to_world_values([corner1, corner2])
    return(cornersCoords)

def brickContainTheGalaxy(ra, dec, galaxyRA, galaxyDec, galaxySMA, galaxyRatio, galaxyPA):
    theta = np.radians(galaxyPA)

    # This is done for flipping the right ascention, the ellipse is defined looking at DS9 (basically the position angle is from west-anticlockwise)
    ra = 360 - ra
    galaxyRA = 360 - galaxyRA

    ra_InEllipseCentre  = ra - galaxyRA
    dec_InEllipseCentre = dec - galaxyDec

    cos_theta = np.cos(-theta)
    sin_theta = np.sin(-theta)
    ra_rot = cos_theta * ra_InEllipseCentre - sin_theta * dec_InEllipseCentre
    dec_rot = sin_theta * ra_InEllipseCentre + cos_theta * dec_InEllipseCentre

    smb = galaxySMA * galaxyRatio
    ellipse_eq = (ra_rot / galaxySMA)**2 + (dec_rot / smb)**2
    return ellipse_eq <= 1

def brickContainBrightStar(bricksNames, bricksRA, bricksDec, brightStarsRa, brighStarsDec, brickWidth):
    mask_bricksWithBrightStars = np.full(len(bricksNames), False)
    bricksWithStarCounter = 0

    for j in range(len(brightStarsRa)):
        for i in range(len(bricksRA)):
            if ( np.abs(bricksRA[i] - brightStarsRa[j]) < (brickWidth/2) and np.abs(bricksDec[i] - brighStarsDec[j]) < (brickWidth/2) ):
                mask_bricksWithBrightStars[i] = True
                bricksWithStarCounter += 1

    if (bricksWithStarCounter > len(brightStarsRa)):
        raise Exception("Error in rejecting the bricks with bright stars. A star has been found in more than one brick, something went wront")
    return(mask_bricksWithBrightStars)

def plotEllipseAndBricks(x0, y0, sma, axis_ratio, pa, x, y, pointsmask_bricksInsideGalaxy, pointsmask_bricksWithBrightStars, maskCombined, mosaicDir):
    t = np.linspace(0, 2 * np.pi, 500)
    ellipse_x = sma * np.cos(t)
    ellipse_y = sma * axis_ratio * np.sin(t)
    
    # Rotate the ellipse by the position angle
    theta = np.radians(180 - pa) # This is for inverting the position angle, because the plot has the x-axis flipped so
                                 # we see the ellipse as we see the galaxy i nDS9
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    ellipse_x_rot = cos_theta * ellipse_x - sin_theta * ellipse_y
    ellipse_y_rot = sin_theta * ellipse_x + cos_theta * ellipse_y
    
    # Translate the ellipse to its center
    ellipse_x_rot += x0
    ellipse_y_rot += y0
    
    fig, ax = plt.subplots(1, 1, figsize=(12 ,12))
    configureAxis(ax, 'RA (deg)', 'Dec (deg)', logScale=False)
    plt.tight_layout(pad=5)
    plt.gca().invert_xaxis()  # Ensure the x-axis reflects the RA convention

    ax.plot(ellipse_x_rot, ellipse_y_rot, label="Galaxy region", color="teal", lw=2)
    ax.scatter(x[pointsmask_bricksInsideGalaxy], y[pointsmask_bricksInsideGalaxy], color="red", marker="^", s=120, label="Rejected - Bricks in galaxy")
    ax.scatter(x[pointsmask_bricksWithBrightStars], y[pointsmask_bricksWithBrightStars], color="orange", marker="v", s=120, label="Rejected - Bricks with Bright stars")
    ax.scatter(x[~maskCombined], y[~maskCombined], color="green", marker="s", s=120, label="Accepted bricks")

    plt.legend(shadow=True, fontsize=17, loc='upper left', bbox_to_anchor=(0.2, 1.1))
    plt.axis("equal")
    plt.savefig(mosaicDir + "/downloadedBricks.png")

def writeBricksAndItsCoordinates(file, brickNames, ra, dec, survey):
    with open(file, 'w') as f:
        f.write("BrickName\tRA_centre\tDec_centre\n")
        for i in range(len(brickNames)):
            if survey=='PANSTARRS':
                #Just to avoid the .fits in the brick identification file, will be useful in the future
                brickName=brickNames[i][:-5]
            elif survey=='DECaLS':
                brickName=brickNames[i]
            f.write(brickName + "\t" + "{:.6f}".format(ra[i]) + "\t" +  "{:.6f}".format(dec[i]) + "\n")


if (len(sys.argv) < 8):
    raise Exception("A frame path, ring path, array of filters and download destination has to be provided")

filters                         = sys.argv[1].split(',')
downloadDestination             = sys.argv[2]
galaxyRA                        = float(sys.argv[3])
galaxyDec                       = float(sys.argv[4])
fieldSize                       = float(sys.argv[5])
mosaicDir                       = sys.argv[6]
bricksIdentificationFile        = sys.argv[7]
gaiaCatalogue                   = sys.argv[8]
survey                          = sys.argv[9]
sizeOfBrick                     = int(sys.argv[10])


setMatplotlibConf()
decalsBrickWidthDeg = 15.5 / 60 
panstarrsBrickWidthDeg=0.25*float(sizeOfBrick) / 3600


# These parameters were used when we were rejecting decals bricks that were on top of the galaxy or with bright stars
# Now we are not using them, so we hide them from the user and just let them here (just in case)
galaxySMA       = 0
galaxyAxisRatio = 0
galaxyPA        = 0
starThresholdForRejectingBricks = -999


# Compute corners of the field
galaxyMinimumRA  = galaxyRA  - (fieldSize / 2)
galaxyMaximumRA  = galaxyRA  + (fieldSize / 2)
galaxyMinimumDec = galaxyDec - (fieldSize / 2)
galaxyMaximumDec = galaxyDec + (fieldSize / 2)
cornersWCSCoords = [(galaxyMinimumRA, galaxyMinimumDec), (galaxyMaximumRA, galaxyMaximumDec)]

if survey=='DECaLS':
    bricksNames, bricksRA, bricksDec = getBrickNamesAndCoordinatesFromRegionDefinedByTwoPoints(cornersWCSCoords[0], cornersWCSCoords[1])
    threadList = []

    for i in bricksNames:
        threadList.append(threading.Thread(target=downloadBrickDecals, args=(i, filters, downloadDestination, False)))
    for i in threadList:
        i.start()
    for i in threadList:
        i.join()

    # This is to keep the format that we are using, so the code is compatible with panstarrs too
    bricksNames = [ i + f".{filters[0]}" for i in bricksNames]



elif survey=='PANSTARRS':
    bricks_fullNames,bricksRA,bricksDec,bricksNames = getPanstarrsBricksFromRegionDefinedByTwoPoints(cornersWCSCoords[0],cornersWCSCoords[1],filters,sizeOfBrick)
    threadList = []
    for i in range(len(bricks_fullNames)):
        b_fname=bricks_fullNames[i]
        b_ra=bricksRA[i]
        b_dec=bricksDec[i]
        b_name=bricksNames[i]
        threadList.append(threading.Thread(target=downloadBrickPanstarrs,args=(b_fname,b_name,b_ra,b_dec,downloadDestination,sizeOfBrick,False)))
    for i in threadList:
        i.start()
    for i in threadList:
        i.join()
else:
    raise Exception (f"Survey {survey} not supported for Photometric calibration")
    
# All this comments concern the fact that we can mask bricks where galaxy or bright stars are allocated
"""
mask_bricksInsideGalaxy = brickContainTheGalaxy(np.array(bricksRA), np.array(bricksDec), galaxyRA, galaxyDec, galaxySMA, galaxyAxisRatio, galaxyPA)

ra, dec = getCoordsAndMagFromBrightStars(gaiaCatalogue, starThresholdForRejectingBricks)
mask_bricksWithBrightStar = brickContainBrightStar(bricksNames, np.array(bricksRA), np.array(bricksDec), ra, dec, decalsBrickWidthDeg)

maskCombined = mask_bricksInsideGalaxy | mask_bricksWithBrightStar

bricksToDownload = bricksNames[~maskCombined]
plotEllipseAndBricks(galaxyRA, galaxyDec, galaxySMA, galaxyAxisRatio, galaxyPA, bricksRA, bricksDec, mask_bricksInsideGalaxy, mask_bricksWithBrightStar, maskCombined, mosaicDir)
"""

writeBricksAndItsCoordinates(bricksIdentificationFile, bricksNames, bricksRA, bricksDec, survey)

"""
For Panstarrs, sometimes the ping does not download the full frame, but the .fits file is created.
Ie, we have a file brick.fits with no data (a size of ~200Kby). Because of that, we will make the following:
we will check the size of each file created and stored in writeBricksAndItsCoordinates. If the file is ok, we
don't do anything, but if the file is not ok, we will re-make the download, and check again. If it still does not work,
we directly remove the brick from the IdentificationFile and the folder.
"""
if survey=='PANSTARRS':
    bricksIdFile=pd.read_csv(bricksIdentificationFile,sep='\t')
    brickSize_mb=float(sizeOfBrick)*float(sizeOfBrick)*50/(3600*3600)
    #1000x1000 pix panstarrs brick has a size of ~4Mb
    badRows=[]
    for row in range(len(bricksIdFile)):
        fname=downloadDestination+'/'+bricksIdFile.loc[row]['BrickName']+'.fits'
        fsize=os.stat(fname).st_size/1e6
        #check test: 0.2 times size of brick in Mb
        if fsize>0.2*brickSize_mb:
            continue
        else:
            #We are assuming a reasonable fit has a size >5Mb
            ra_brick=bricksIdFile.loc[row]['RA_centre']
            dec_brick=bricksIdFile.loc[row]['Dec_centre']
            brick_fullName,brickRA,brickDec,brickName=getPanstarrsBricksFromCentralPoint(ra_brick,dec_brick,filters,sizeOfBrick)
            brick_fullName=brick_fullName[0]; brickRA=brickRA[0]; brickDec=brickDec[0]
            brickName=brickName[0]
            downloadBrickPanstarrs(brick_fullName,brickName,brickRA,brickDec,downloadDestination,sizeOfBrick,overwrite=True)
            #We check again the download
            fname=downloadDestination+'/'+brickName
            fsize=os.stat(fname).st_size/1e6
            #Decission: if download is ok, we change the dataframe and store, if not, we remove the row
            if fsize>0.2*brickSize_mb:
                bricksIdFile.loc[row,'BrickName']=brickName[:-5]
                bricksIdFile.loc[row,'RA_centre']=brickRA
                bricksIdFile.loc[row,'Dec_centre']=brickDec
            else:
                #We remove the fits file and store the bad row
                os.system(f'rm {fname}')
                badRows.append(row)
    ##Now we remove the bad rows and store the new bricksIdentificationFile
    bricksIdFile=bricksIdFile.drop(index=badRows)
    bricksIdFile.to_csv(bricksIdentificationFile,sep='\t',index=False)

