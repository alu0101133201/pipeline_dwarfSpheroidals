import sys

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

image = sys.argv[1]
ra  = float(sys.argv[2])
dec = float(sys.argv[3])
radius_arcsec = float(sys.argv[4])
valueToPut = sys.argv[5]


# Load image
hdu = fits.open(image)
data = hdu[1].data
wcs = WCS(hdu[1].header)

# Convert the center coordinates from (RA, Dec) to pixel coordinates
x, y = wcs.world_to_pixel_values(ra, dec)

# Convert radius from arcseconds to pixels
pix_scale = np.abs(wcs.pixel_scale_matrix[1,1]) * 3600  # Pixel scale in arcseconds/pixel
radius_pixels = radius_arcsec / pix_scale


# Create a circular mask
yy, xx = np.ogrid[:data.shape[0], :data.shape[1]]
mask = (xx - x)**2 + (yy - y)**2 <= (radius_pixels)**2

# Apply the mask (set pixels to NaN)
if (valueToPut == "nan"):
    data[mask] = np.nan  # You can also use 0 or another value
else:
    data[mask] = float(valueToPut)

# Save the masked image
# fits.writeto(image, data, hdu[0].header, overwrite=True)
hdu.writeto(image, overwrite=True)