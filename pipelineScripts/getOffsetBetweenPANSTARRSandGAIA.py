import sys
import glob

import numpy as np

def readCoordsAndMagFromCatalogue(file, colRA, colDec, colMag):
    data = np.loadtxt(file, comments='#', usecols=(colRA, colDec, colMag))

    ra = np.array(data[:, 0].astype(float).tolist())
    dec = np.array(data[:, 1].astype(float).tolist())
    magnitude = np.array(data[:, 2].astype(float).tolist())

    return(ra, dec, magnitude)

panstarrsCataloguesDir=sys.argv[1]
gaiaCatalogueFile=sys.argv[2]
brightLimit=sys.argv[3]
faintLimit=sys.argv[4]



gaiaRA, gaiaDEC, gaiaMag = readCoordsAndMagFromCatalogue(gaiaCatalogueFile, 1, 2, 3)
