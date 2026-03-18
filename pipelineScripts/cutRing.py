"""
This script processes a FITS image by removing from the normalisation ring some pixels that have border-related defects.
It opens an existing FITS file, creates a writable copy of the image data, optionally enforces binary pixel values,
and sets a specified region of pixels (corresponding to detector borders) to zero. The modified image is then written
back to a FITS file while preserving the original header information.

The region to be masked is provided via command-line
arguments defining the start and end indices for rows and columns. A value of -1 indicates
no cut on that bound (i.e., use the full extent of the array in that direction).

"""


import sys
import numpy as np
from astropy.io import fits


# -----------------------------
# Parse command-line arguments
# -----------------------------
if len(sys.argv) != 5:
    raise ValueError(
        "Usage: python3 file.py row_start row_end col_start col_end\n"
        "Use -1 to indicate no cut on a given bound."
    )

row_start, row_end, col_start, col_end = map(int, sys.argv[1:])

def to_slice(start, end):
    """Convert command-line bounds to a Python slice."""
    s = None if start == -1 else start
    e = None if end == -1 else end
    return slice(s, e)

row_slice = to_slice(row_start, row_end)
col_slice = to_slice(col_start, col_end)

# -----------------------------
# Open FITS file
# -----------------------------
with fits.open("./build/ring/ring.fits") as hdul:
    data = hdul[1].data
    header = hdul[1].header

# Ensure we are working on a writable copy
data = np.array(data, copy=True)

# Optional: enforce binary values (white=0, black=1)
data = (data > 0).astype(data.dtype)

# -----------------------------
# Apply masking
# -----------------------------
if row_slice.start is not None or row_slice.stop is not None:
    data[row_slice, :] = 0

if col_slice.start is not None or col_slice.stop is not None:
    data[:, col_slice] = 0

# -----------------------------
# Write output FITS file
# -----------------------------
new_hdu = fits.ImageHDU(data=data, header=header)

hdul_out = fits.HDUList([
    fits.PrimaryHDU(),
    new_hdu
])

hdul_out.writeto("./build/ring/ring.fits", overwrite=True)