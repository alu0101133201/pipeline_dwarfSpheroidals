import sys
from scipy.stats import skew,kurtosis
from astropy.io import fits
import numpy as np
image=sys.argv[1]
value_to_measure=sys.argv[2]

data=fits.open(image)[1].data
data_nonan=data[np.logical_not(np.isnan(data))]

if value_to_measure=="SKEWNESS":
    print(skew(data_nonan))
elif value_to_measure=="KURTOSIS":
    print(kurtosis(data_nonan))
    