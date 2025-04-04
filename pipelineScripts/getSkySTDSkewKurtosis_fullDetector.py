import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.stats import skew, kurtosis
import sys, os

input_file=sys.argv[1]
ring_file=sys.argv[2]
output_folder=sys.argv[3]
num_ccd=int(sys.argv[4])

baseName=os.path.basename(input_file)
baseName_txt=os.path.splitext(baseName)[0]+'.txt'
outputFile=os.path.join(output_folder,baseName_txt)

all_pixels=[]
for h in np.arange(1,num_ccd+1):
    input_data=fits.open(input_file)[h].data
    ring_data=fits.open(ring_file)[h].data
    mask = (~np.isnan(input_data)) & (ring_data!=0)
    all_pixels.extend(input_data[mask])

all_pixels=np.array(all_pixels)
mean,median,std=sigma_clipped_stats(all_pixels,sigma=1,maxiters=3)
skewness=skew(all_pixels)
kurt=kurtosis(all_pixels)

output_line=f"{baseName} {mean:.5e} {std:.5e} {skewness:.5e} {kurt:.5e}\n"

with open(outputFile,"a") as f:
    for h in np.arange(1,num_ccd+1):
        f.write(output_line)