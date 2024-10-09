import glob
import argparse

from astropy.io import fits

parser = argparse.ArgumentParser(description='Script that checks the exposures times of the frames.\
 It takes the data path (i.e. the directory DATA or DATA-or) and the hdu to look for the values and \
 reads the fits files looking for EXPOSURE.\
 The comparison is done for all the "night*" folders that are contained in the data directory')
 
parser.add_argument('path', type=str, help='The path to directory that contains the night folders with the fits files (generally DATA or DATA-or)')
parser.add_argument('hdu', type=str, help='The hdu of the fits files that contains the EXPOSURE keyword')

args = parser.parse_args()
dataDirectory = args.path
hduNum = int(args.hdu)


nightFolders = glob.glob(dataDirectory + "/night*")
numOfNights = len(nightFolders)
exposureKeyWord = "EXPTIME"

exposureTimes = set()
framesChecked = 0

for i in nightFolders:
    fitsFiles = glob.glob(i + "/*.fits")
    for j in fitsFiles:
        tmpHeader = fits.open(j)[hduNum].header
        tmpExposureTime = tmpHeader[exposureKeyWord]

        exposureTimes.add(tmpExposureTime)
        framesChecked += 1

print("The exposures times found are: ", exposureTimes)
print("The number of frames checked is: ", framesChecked)
