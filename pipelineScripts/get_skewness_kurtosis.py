import sys
from scipy.stats import skew,kurtosis
from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np

image=sys.argv[1]
value_to_measure=sys.argv[2]
ring=sys.argv[3]
h=int(sys.argv[4])

data=fits.open(image)[h].data
mask=fits.open(ring)[h].data
nonan_ring=(np.logical_not(np.isnan(data))) & (mask!=0)

## We remove the negatives. Here the background has not been substracted so they are all bad pixels
clipped_data = data[nonan_ring]
clipped_data = [x for x in clipped_data if x > 0]
##

if value_to_measure=="SKEWNESS":
    print(skew(data[nonan_ring],axis=None,nan_policy='omit'))
elif value_to_measure=="KURTOSIS":
    print(kurtosis(data[nonan_ring],axis=None,nan_policy='omit'))
    
