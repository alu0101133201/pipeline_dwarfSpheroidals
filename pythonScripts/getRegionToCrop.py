import sys

import numpy as np

from astropy.io import fits


def find_non_nan_region(matrix):
    # Create a mask of where NaN values are present
    isnan_mask = np.isnan(matrix)

    # Find the indices of non-NaN values
    non_nan_indices = np.argwhere(~isnan_mask)

    # Get the bounding box of non-NaN values
    if non_nan_indices.size > 0:
        row_min, col_min = non_nan_indices.min(axis=0)
        row_max, col_max = non_nan_indices.max(axis=0)
        return (row_min, row_max + 1, col_min, col_max + 1)
    else:
        return None  # All values are NaN

# Example usage
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