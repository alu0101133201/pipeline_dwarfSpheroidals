import sys

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from astropy.io.fits import PrimaryHDU, ImageHDU, HDUList

image = sys.argv[1]
valueToPut = sys.argv[2]
ra  = float(sys.argv[3])
dec = float(sys.argv[4])
radius_arcsec = float(sys.argv[5])

# Optional arguments. By default (if not provided) a circle with no rotation (axisRatio = 0 and pa=0)
if len(sys.argv) > 6:
    axis_ratio = float(sys.argv[6]) 
else:
    axis_ratio = 1.0  

if len(sys.argv) > 7:
    pa_deg = float(sys.argv[7])
else:
    pa_deg = 0.0  

pa_rad = np.deg2rad(pa_deg) 

# Load image
hdu = fits.open(image)
data = hdu[1].data
wcs = WCS(hdu[1].header)

# Convert the center coordinates from (RA, Dec) to pixel coordinates
x, y = wcs.world_to_pixel_values(ra, dec)

# Convert radius from arcseconds to pixels
pix_scale = np.abs(wcs.pixel_scale_matrix[1,1]) * 3600  # Pixel scale in arcseconds/pixel
radius_pixels = radius_arcsec / pix_scale

# Create a circular or elliptical mask
yy, xx = np.ogrid[:data.shape[0], :data.shape[1]]
dx = xx - x
dy = yy - y

# Rotate coordinates by PA (clockwise)
dx_rot = dx * np.cos(pa_rad) + dy * np.sin(pa_rad)
dy_rot = -dx * np.sin(pa_rad) + dy * np.cos(pa_rad)

# Elliptical mask formula in rotated coordinates
a = radius_pixels
b = a * axis_ratio
mask = (dx_rot / a)**2 + (dy_rot / b)**2 <= 1.0

# Apply the mask (set pixels to NaN)
if (valueToPut == "nan"):
    data[mask] = np.nan  # You can also use 0 or another value
else:
    data[mask] = float(valueToPut)

# Save the masked image
header = hdu[1].header.copy()
if 'BLANK' in header:
    del header['BLANK']

primary_hdu = PrimaryHDU()  
image_hdu = ImageHDU(data=data, header=header)
new_hdul = HDUList([primary_hdu, image_hdu])
new_hdul.writeto(image, overwrite=True)
