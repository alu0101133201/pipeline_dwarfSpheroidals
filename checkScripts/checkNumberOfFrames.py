import glob
import argparse

from astropy.io import fits

parser = argparse.ArgumentParser(description='Script that checks the number of frames.\
 It takes the data path (i.e. the directory DATA or DATA-or) \
 The check is done for all the "night*" folders that are contained in the data directory')

parser.add_argument('path', type=str, help='The path to directory that contains the night folders with the fits files (generally DATA or DATA-or)')

args = parser.parse_args()
dataDirectory = args.path


nightFolders = glob.glob(dataDirectory + "/night*")
numOfNights = len(nightFolders)

framesChecked = 0

for i in nightFolders:
    fitsFiles = glob.glob(i + "/*.fits")
    for j in fitsFiles:
        framesChecked += 1

print("The number of nights checked is: ", len(nightFolders))
print("The number of frames checked is: ", framesChecked)
