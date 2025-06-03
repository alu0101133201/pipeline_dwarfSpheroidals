####This script is not called inside pipeline, but as a complement if you
#want to manually add an astrometry. The script requires the following:
# 1: RA and DEC of a refference star
# 2: XY of the refference star in the image to astrometrize
# 3: An image with a good WCS whose matrix transformation is useful
#We will only change refference point
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
import astropy.units as u

Ra_star=float(sys.argv[1])
Dec_star=float(sys.argv[2])
X_star=float(sys.argv[3])
Y_star=float(sys.argv[4])
image=sys.argv[5]
ref_image=sys.argv[6]

hdul=fits.open(image,mode='update')
header=hdul[1].header
ref_wcs=WCS(fits.open(ref_image)[1].header)
CD = ref_wcs.wcs.cd if ref_wcs.wcs.has_cd() else ref_wcs.wcs.pc * ref_wcs.wcs.cdelt

wcs = WCS(naxis=2)
wcs.wcs.cd = CD
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

wcs.wcs.crpix = [X_star,Y_star]
coord_star=SkyCoord(Ra_star*u.deg,Dec_star*u.deg)
wcs.wcs.crval = [coord_star.ra.deg, coord_star.dec.deg]

wcs_header = wcs.to_header()
header.update(wcs_header)

# Save and close
hdul.flush()
hdul.close()

print(f"{image} has been updated with {ref_image} wcs.")