import sys

import numpy as np

from astropy.io import fits

def getDataFromFitsCatalogue(cat):
    hdul = fits.open(cat)
    data = hdul[1].data
    hdul.close()
    return(data)

def getDataFromPlainCatalogue(cat, raCol, decCol):
    ra = []
    dec = []

    with open(cat, 'r') as file:
        lines = file.readlines()
        for i in lines:
            splittedLine = i.split()

            if (splittedLine[0] == "#"):
                continue

            ra.append(float(splittedLine[raCol]))
            dec.append(float(splittedLine[decCol]))

    dtype = [('ra', float), ('dec', float)]
    data = np.array(list(zip(ra, dec)), dtype=dtype)

    return data


def createDS9RegionFile(data, outputFile):
    firstLine = '# Region file format: DS9 version 4.1\nglobal color=red dashlist=8 3 width=4 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'
    with open(outputFile, 'w') as f:
        f.write(firstLine)
        for i in data:
            f.write("circle(" + str(i['ra']) + "," + str(i['dec']) + ",10\")\n")
    return()

catalogue        = sys.argv[1]
outputRegionFile = sys.argv[2]
catalogueFormat  = sys.argv[3]
raCol            = int(sys.argv[4])   
decCol           = int(sys.argv[5])

if (catalogueFormat == "fits"):
    data = getDataFromFitsCatalogue(catalogue)
elif (catalogueFormat == "plain"):
    data = getDataFromPlainCatalogue(catalogue, raCol, decCol)
else:
    raise Exception(f"Catalofue format {catalogueFormat} not supported")
    exit()

createDS9RegionFile(data, outputRegionFile)

