import sys
from scipy.stats import skew,kurtosis
from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np

image=sys.argv[1]
value_to_measure=sys.argv[2]

data=np.array(fits.open(image)[1].data)

# If a third argument is provided it is expected to be the ring (estimating the
# background with a ring). Otherwise It must be noisechisel and we ignore the ring parameter
if (len(sys.argv) > 3):
    ring=sys.argv[3]

    mask=fits.open(ring)[1].data
    nonan_ring=(np.logical_not(np.isnan(data))) & (mask!=0)

    maskedData = data[nonan_ring]
else:
    maskedData = data

## We remove the negatives. Here the background has not been substracted so they are all bad pixels
maskedData = maskedData.flatten()
clipped_data = maskedData[maskedData > 0]
##

if value_to_measure=="SKEWNESS":
    print(skew(clipped_data, axis=None, nan_policy='omit'))
elif value_to_measure=="KURTOSIS":
    print(kurtosis(clipped_data, axis=None, nan_policy='omit'))
    
