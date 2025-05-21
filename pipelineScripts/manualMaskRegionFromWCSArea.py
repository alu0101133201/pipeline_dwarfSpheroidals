import sys

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

image = sys.argv[1]
ra  = float(sys.argv[2])
dec = float(sys.argv[3])
radius_arcsec = float(sys.argv[4])
valueToPut = sys.argv[5]

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

def is_in_ellipse(px, py):
    dx = px - x 
    dy = py - y
    dx_rot = dx * np.cos(pa_rad) + y * np.sin(pa_rad)
    dy_rot = -dx * np.sin(pa_rad) + dy * np.cos(pa_rad)
    return (dx_rot / a)**2 + (dy_rot / b)**2 <=1.0

# Load image
hdu_list = fits.open(image)
for i in range(len(hdu_list)):
    if isinstance(hdu_list[i], fits.hdu.image.ImageHDU) and hdu_list[i].data is not None:
        try:
            wcs = WCS(hdu_list[1].header)
# Convert the center coordinates from (RA, Dec) to pixel coordinates
            x, y = wcs.world_to_pixel_values(ra, dec)

# Convert radius from arcseconds to pixels
            pix_scale = np.abs(wcs.pixel_scale_matrix[1,1]) * 3600  # Pixel scale in arcseconds/pixel
            radius_pixels = radius_arcsec / pix_scale

            a = radius_pixels
            b = a*axis_ratio

            shape = hdu_list[i].data.shape
            xmin, xmax = 0, shape[1]
            ymin, ymax = 0, shape[0]

            has_overlap = False

            if xmin <= x <= xmax and ymin <= y <= ymax:
                has_overlap = True

            if not has_overlap:
                num_points=20
                for px in np.linspace(xmin,xmax,num_points):
                    if is_in_ellipse(px,ymin) or is_in_ellipse(px,ymax):
                        has_overlap=True
                        break

                if not has_overlap:
                    for py in np.linspace(ymin,ymax,num_points):
                        if is_in_ellipse(xmin,py) or is_in_ellipse(xmax,py):
                            has_overlap=True
                            break

            if not has_overlap:
                corners=[(xmin,ymin),(xmin,ymax),(xmax,ymin),(xmax,ymax)]
                for corner in corners:
                    if is_in_ellipse(*corner):
                        has_overlap=True
                        break

            if not has_overlap:
                if (xmin <= x <= xmax and ymin <= y <= ymax ) and (2*a <= xmax-xmin and 2*b <= ymax-ymin):
                    has_overlap = True

            if has_overlap:
                data = hdu_list[i].data

# Create a circular mask
                yy, xx = np.ogrid[:data.shape[0], :data.shape[1]]
                dx = xx - x
                dy = yy - y

                dx_rot = dx * np.cos(pa_rad) + dy * np.sin(pa_rad)
                dy_rot = -dx * np.sin(pa_rad) + dy * np.cos(pa_rad)
                mask = (dx_rot / a)**2 + (dy_rot / b)**2 <= 1.0

# Apply the mask (set pixels to NaN)
                if (valueToPut == "nan"):
                    data[mask] = np.nan  # You can also use 0 or another value
                else:
                    data[mask] = float(valueToPut)
                if 'BLANK' in hdu_list[i].header:
                    del hdu_list[i].header['BLANK']
                
                print(f"Applied elliptical mask to HDU {i} - center at pixel ({x:.1f}, {y:.1f}), axes {a:.1f}x{b:.1f} pixels, PA {pa_deg} deg")
            else:
                print(f"Skipping HDU {i} - elliptical region does not overlap this layer")
                
        except Exception as e:
            print(f"Skipping HDU {i} - processing error: {str(e)}")

# Save the masked image
hdu_list.writeto(image, overwrite=True)
hdu_list.close()

