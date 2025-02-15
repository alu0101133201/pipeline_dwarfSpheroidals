import sys
from astropy.io import fits
from astropy.wcs import WCS
import os

image=sys.argv[1]
ccd=int(sys.argv[2])
x_camera=float(sys.argv[3])
y_camera=float(sys.argv[4])


image_data=fits.open(image)
image_wcs=WCS(image_data[ccd].header)

x_detector,y_detector=image_wcs.world_to_pixel_values(x_camera,y_camera)

print(f"x_det={x_detector}")
print(f"y_det={y_detector}")