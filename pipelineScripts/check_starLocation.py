import sys,os
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from astropy.coordinates import SkyCoord

image=sys.argv[1]
out_file=sys.argv[2]
num_ccd=int(sys.argv[3])
star_ra=float(sys.argv[4])
star_dec=float(sys.argv[5])
circ_r=float(sys.argv[6])
star_coords=SkyCoord(ra=star_ra,dec=star_dec,unit='deg')
base=os.path.basename(image)

###We will store the number of pixels that falls inside the circle
with fits.open(image) as hdul:
    npixels=[]
    for h in np.arange(1,num_ccd+1):
        wcs = WCS(hdul[h].header)
        
        x_star,y_star=wcs.world_to_pixel(star_coords)
        y,x=np.indices(hdul[h].data.shape)
        y=y[::-1]
        distance=np.sqrt((x-x_star)**2+(y-y_star)**2)
        mask = distance <= circ_r
        count=np.sum(mask)
        npixels.append(count)
max_pix=max(npixels)
if max_pix!=0:
    h_max=npixels.index(max_pix)+1
    with open(out_file,'a') as f:
        f.write(f"{image} -h{h_max} \n")
            
        