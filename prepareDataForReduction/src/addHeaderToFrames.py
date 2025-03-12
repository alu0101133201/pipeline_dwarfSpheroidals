import os
import re
import sys
import glob

from astropy.io import fits

dataDir=sys.argv[1]
headerDir=sys.argv[2]
dataHdu=int(sys.argv[3])
headerHdu=int(sys.argv[4])
rawDir=sys.argv[5]

# Here is assumed that the fits files have the exact same name in both folders
# i.e. the file nnn.fits in the dataDir is called nnn.fits in headerDir 

for dataFile in glob.glob(dataDir + "/*.fits"):
    print("\nProcessing file: ", dataFile)
    fileName = dataFile.split('/')[-1]

    if (bool(re.search(r"Bias_Dark", fileName))):
        print("Skipping dark...")
        continue

    headerFile = headerDir + "/" + fileName
    try:
        with fits.open(dataFile) as hdulData, fits.open(headerFile) as hdulHeader:
            # Create a new Primary HDU with data from the second HDU and header from header file
            hdu = fits.PrimaryHDU(data=hdulData[dataHdu].data, header=hdulHeader[headerHdu].header)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto(dataFile, overwrite=True)
    except FileNotFoundError:
        print("The reduced file for: " + str(fileName) + " does not exist.\n")
        print("Copying original HDU from raw. Airmass will be updated in pipeline.")
        headerFile=rawDir+"/"+fileName[:-5]+"/"+fileName
        with fits.open(dataFile) as hdulData, fits.open(headerFile) as hdulHeader:
            hdu = fits.PrimaryHDU(data=hdulData[dataHdu].data, header=hdulHeader[headerHdu].header)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto(dataFile, overwrite=True)
        #os.system("mv " + dataFile + " " + dataDir +"/noRedImage_" + fileName)

