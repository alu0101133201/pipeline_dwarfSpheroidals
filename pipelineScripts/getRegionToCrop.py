import sys

import numpy as np

from astropy.io import fits

# This script receives an image (and the hdu in which the data is located) and returns four coordinates
# This coordinates correspond to the smaller rectangle that encloses the non-nan data.

# In the pipeline this is needed because we place the data in a huge grid (the grid of the final coadd)
# But, since we do not want to work with this huge grid in all the intermediate steps (would take a long time)
# We crop the region of interest and we use this crop in the intermediate steps to work faster
# This scripts is the one who provides the region when we perform this crop


def find_non_nan_region(matrix):
    isnan_mask = np.isnan(matrix)
    non_nan_indices = np.argwhere(~isnan_mask)

    if non_nan_indices.size > 0:
        row_min, col_min = non_nan_indices.min(axis=0)
        row_max, col_max = non_nan_indices.max(axis=0)
        return (int(row_min), int(row_max + 1), int(col_min), int(col_max + 1))
    else:
        return None  # All values are NaN



imageFile = sys.argv[1]
hduNum = int(sys.argv[2])


with fits.open(imageFile) as hdul:
    data = np.array(hdul[hduNum].data)

region = find_non_nan_region(data)
if region:
    row_min, row_max, col_min, col_max = region
    print(row_min, row_max, col_min, col_max)

else:
    raise Exception("The entire matrix is NaN.")