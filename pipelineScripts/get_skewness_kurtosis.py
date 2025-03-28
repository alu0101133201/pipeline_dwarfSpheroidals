import sys
from scipy.stats import skew,kurtosis
from astropy.io import fits
import numpy as np
image=sys.argv[1]
value_to_measure=sys.argv[2]
ring=sys.argv[3]

data=fits.open(image)[1].data
mask=fits.open(ring)[1].data
nonan_ring=(np.logical_not(np.isnan(data))) & (mask!=0)

if value_to_measure=="SKEWNESS":
    print(skew(data[nonan_ring],axis=None,nan_policy='omit'))
elif value_to_measure=="KURTOSIS":
    print(kurtosis(data[nonan_ring],axis=None,nan_policy='omit'))
    
