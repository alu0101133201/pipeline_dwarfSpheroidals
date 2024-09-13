import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import shift

dataFile = sys.argv[1]
planeFile = sys.argv[2]
hduNum = int(sys.argv[3])
fullGridCoords_x = int(sys.argv[4])
fullGridCoords_y = int(sys.argv[5])
outputPlaneFile = sys.argv[6]

with fits.open(dataFile) as hdul:
    dataImage = hdul[hduNum].data
    dataHeader = hdul[hduNum].header
    dataWcs = WCS(dataHeader)

with fits.open(planeFile) as hdul:
    planeImage = hdul[hduNum].data
    planeHeader = hdul[hduNum].header
    planeWcs = WCS(planeHeader)


# The plane that is fitted to the image is fitted in the small grid, so it begins at 0,0
# But in the full grid the image is somewhere else, so we need to compute the offset and move the image

referencePixel = planeImage[0][0]
worldReferenceCoords = planeWcs.pixel_to_world(0, 0)
imageReferenceCoords = dataWcs.world_to_pixel(worldReferenceCoords)
offset = (int(imageReferenceCoords[1]), int(imageReferenceCoords[0]))
print(offset)

# Create the new image with the plane moved
orig_x, orig_y = planeImage.shape
padded_matrix = np.zeros((fullGridCoords_x, fullGridCoords_y))
start_x = offset[0]
start_y = offset[1]


# I tried to do the subtraction px by px trying to avoid saving a fits but it's much slower
# save a tmp fits and the use astarithmetic. This plane will be deleted after subtracting it so no problem
if (start_x + orig_x <= fullGridCoords_x) and (start_y + orig_y <= fullGridCoords_y):
    padded_matrix[start_x:start_x + orig_x, start_y:start_y + orig_y] = planeImage
else:
    raise ValueError("The original matrix does not fit in the padded matrix with the given shift.")


hdulist_out = fits.HDUList()
hdu_out = fits.ImageHDU(padded_matrix, header=dataHeader)

hdu_empty = fits.PrimaryHDU()
hdu_empty.header['EXTNAME'] = "NODATA"
hdulist_out.append(hdu_empty)
hdulist_out.append(hdu_out)
hdulist_out.writeto(outputPlaneFile, overwrite=True)
