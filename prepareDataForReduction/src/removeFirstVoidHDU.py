import sys
import glob

from astropy.io import fits


folderWithImages = sys.argv[1]

for dataFile in glob.glob(folderWithImages + "/*.fits"):
    with fits.open(dataFile) as hdulData:
        if len(hdulData) == 1:
            print("skipping because it already has 1HDU and seems to have been processed...")
            continue

        # Create a new Primary HDU with data from the second HDU and header from header file
        hdu = fits.PrimaryHDU(data=hdulData[1].data, header=hdulData[1].header)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dataFile, overwrite=True)
