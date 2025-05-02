# This script performs photometry using photutils on an image in the coordinates provided by a catalogue.

# The output needs to be ID, X, Y, RA, DEC, MAG, SUM

import sys
import time

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

def getCoordinatesFromCatalogue(catalogue, raColNumber, decColNumber):
    ra = []
    dec = []

    with open(catalogue, 'r') as file:
        for line in file:
            splittedLine = line.split()
            if (splittedLine[0] == "#"):
                continue
            ra.append(float(splittedLine[raColNumber]))
            dec.append(float(splittedLine[decColNumber]))

    return (np.array(ra), np.array(dec))

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

x, y = getCoordinatesFromCatalogue(catalogue, columnWithXCoordPx, columnWithYCoordPx)
z = np.array([(x[i], y[i]) for i in range(len(x))])

ra, dec = getCoordinatesFromCatalogue(catalogue, columnWithXCoordWCS, columnWithYCoordWCS)
zCoord = np.array([(ra[i], dec[i]) for i in range(len(x))])


sigmaClip = SigmaClip(sigma=3.0, maxiters=3)
innerAnnulus = 3*aperture_radius_px
outerAnnulus = 4*aperture_radius_px

ids = np.arange(len(x))
mags = np.full(len(z), np.nan)
sums = np.full(len(z), np.nan)


# The following section aims to reject sources with nan pixels (not accurate flux, we don't want them)
# Since we use really huge apertures, I just impose not to be nans in the half central region of it
aperturesForNanCheck = CircularAperture((z- 1), r=aperture_radius_px/2)
aperturesMask = aperturesForNanCheck.to_mask(method='center')

# validIndices stores the apertures that do NOT contains nans in a radius of aperture/2
validIndices = []
for i, currApertureMask in enumerate(aperturesMask):
    cutout = currApertureMask.cutout(imageData)
    if cutout is None:
        continue

    mask = (currApertureMask.data == 1)
    aperture_pixels = np.where(mask == 1, cutout, np.inf)
    if not np.isnan(aperture_pixels).any():
        validIndices.append(i)

currentApertures = CircularAperture((z[validIndices] - 1), r=aperture_radius_px)
currentAnnuli  =  CircularAnnulus((z[validIndices] - 1), r_in=innerAnnulus, r_out=outerAnnulus)


phot_table_local = aperture_photometry(imageData, currentApertures, mask = np.isnan(imageData), method="exact")
aperstats = ApertureStats(imageData, currentAnnuli, sigma_clip=sigmaClip)
bkg_mean = aperstats.median

aperture_area = currentApertures.area_overlap(imageData)
total_bkg = bkg_mean * aperture_area

fluxes_valid = np.array(phot_table_local['aperture_sum']) - total_bkg
with np.errstate(divide='ignore', invalid='ignore'):
    mag_valid = -2.5 * np.log10(fluxes_valid) + 22.5

for i, idx in enumerate(validIndices):
    sums[idx] = fluxes_valid[i]
    mags[idx]    = mag_valid[i]


writeDataToCatalogue(output, ids, x, y, ra, dec, mags, sums)