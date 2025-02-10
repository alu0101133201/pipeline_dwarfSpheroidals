# This script performs photometry using photutils on an image in the coordinates provided by a catalogue.

# The output needs to be ID, X, Y, RA, DEC, MAG, SUM

import sys

import numpy as np 

from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.aperture import SkyCircularAperture, SkyCircularAnnulus, aperture_photometry, ApertureStats, CircularAperture, CircularAnnulus

def getImageData(image, hduNumber):
    hdulist = fits.open(image)
    header = hdulist[hduNumber].header
    imageData_original = hdulist[hduNumber].data
    wcs = WCS(header)
    hdulist.close()
    return(imageData_original, np.shape(imageData_original), wcs)

def getCoordinatesFromCatalogue(catalogue, raColNumber, decColNumber,hduNumber,mode):
    ra = []
    dec = []
    if catalogue.endswith(".txt"):
        with open(catalogue, 'r') as file:
            for line in file:
                splittedLine = line.split()
                if (splittedLine[0] == "#"):
                    continue
                ra.append(float(splittedLine[raColNumber]))
                dec.append(float(splittedLine[decColNumber]))

        return (np.array(ra), np.array(dec))
    elif (catalogue.endswith(".fits"))or(catalogue.endswith(".cat")):
        table=fits.open(catalogue)[hduNumber].data
        if mode=='Image':
            x=table['X']; y=table['Y']
            return(x,y)
        elif mode=='WCS':
            ra=table['RA']; dec=table['DEC']
            return(ra,dec)

def writeDataToCatalogue(outputFile, ids, x, y, ra, dec, mag, sums):
    with open(outputFile, 'w') as file:
        file.write("# Column 1: IDs\n")
        file.write("# Column 2: X\n")
        file.write("# Column 3: Y\n")
        file.write("# Column 4: RA\n")
        file.write("# Column 5: DEC\n")
        file.write("# Column 6: MAG\n")
        file.write("# Column 7: sums\n")

        
        for i in range(len(ids)):
            file.write(f"{int(ids[i])} {x[i]:.10f} {y[i]:.10f} {ra[i]:.10f} {dec[i]:.10f} {mag[i]:.4f} {sums[i]:.10f}\n")



catalogue           = sys.argv[1]
image               = sys.argv[2]
aperture_radius_px  = float(sys.argv[3])
output              = sys.argv[4]
zp                  = float(sys.argv[5])
dataHdu             = int(sys.argv[6])
columnWithXCoordPx  = int(sys.argv[7])
columnWithYCoordPx  = int(sys.argv[8])
columnWithXCoordWCS = int(sys.argv[9])
columnWithYCoordWCS = int(sys.argv[10])

imageData, shape, wcs = getImageData(image, dataHdu)

x, y = getCoordinatesFromCatalogue(catalogue, columnWithXCoordPx, columnWithYCoordPx,dataHdu,'Image')
ra, dec = getCoordinatesFromCatalogue(catalogue, columnWithXCoordWCS, columnWithYCoordWCS,dataHdu,'WCS')

sigmaClip = SigmaClip(sigma=3.0, maxiters=3)
innerAnnulus = 3*aperture_radius_px
outerAnnulus = 4*aperture_radius_px

ids = []
mag = []
sums = []

for i in range(len(x)):
    currentAperture = CircularAperture((x[i] - 1, y[i] - 1), r=aperture_radius_px)
    currentAnnulus  =  CircularAnnulus((x[i] - 1, y[i] - 1), r_in=innerAnnulus, r_out=outerAnnulus)

    phot_table_local = aperture_photometry(imageData, currentAperture, mask = np.isnan(imageData), method="exact")
    aperstats = ApertureStats(imageData, currentAnnulus, sigma_clip=sigmaClip)
    bkg_mean = aperstats.median

    aperture_area = currentAperture.area_overlap(imageData)
    total_bkg = bkg_mean * aperture_area
    currentFlux = float(phot_table_local['aperture_sum'][0]) - float(total_bkg)
       
    ids.append(i)

    with np.errstate(invalid='ignore'): # Some values (of bad detections or whatever) are negative and give an error
                                        # since I don't want to exclude them directly because i need that entry in the catalogue I just suppress the warning
        mag.append(-2.5 * np.log10(currentFlux) + zp)
    sums.append(currentFlux)

writeDataToCatalogue(output, ids, x, y, ra, dec, mag, sums)

