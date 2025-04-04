import sys
from scipy.stats import skew,kurtosis
from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np

image=sys.argv[1]
value_to_measure=sys.argv[2]
ring=sys.argv[3]

data=fits.open(image)[1].data
mask=fits.open(ring)[1].data
nonan_ring=(np.logical_not(np.isnan(data))) & (mask!=0)

##
clipped_data = data[nonan_ring]
clipped_data = [x for x in clipped_data if x > 0]
clipped_data = sigma_clip(clipped_data, sigma=30, maxiters=5)  
##

if value_to_measure=="SKEWNESS":
    print(skew(clipped_data,axis=None,nan_policy='omit'))
elif value_to_measure=="KURTOSIS":
    print(kurtosis(clipped_data,axis=None,nan_policy='omit'))
    
