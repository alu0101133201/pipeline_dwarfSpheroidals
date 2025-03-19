import sys

import numpy as np

def readCoordsAndMagFromCatalogue(file, colRA, colDec, colMag):
    data = np.loadtxt(file, comments='#', usecols=(colRA, colDec, colMag))

    ra = np.array(data[:, 0].astype(float).tolist())
    dec = np.array(data[:, 1].astype(float).tolist())
    magnitude = np.array(data[:, 2].astype(float).tolist())

    return(ra, dec, magnitude)

def matchSources(ra1, ra2, dec1, dec2, mag1, mag2):
    matchedSources = []
    for i in range(len(ra1)):
        for j in range(len(ra2)):
            if ((np.abs(ra1[i] - ra2[j]) < TOLERANCE_DEC) and (np.abs(dec1[i] - dec2[j]) < TOLERANCE_DEC)):
                matchedSources.append((ra1[i], dec1[i], mag1[i], mag2[j]))

    return(matchedSources)

def getColoursOfMatchedSources(sources):
    colours = []
    for i in sources:
        colours.append((i[0], i[1], i[2] - i[3]))
    return(np.array(colours))


def addColoursToCatalogue(matchedCoordinatesAndcolours, catalogueToAdd, catalogueOutput):

    headers = []
    updatedLines = []

    with open(catalogue3) as infile:
        for i, line in enumerate(infile):
            if line.startswith("#"):  # Preserve headers
                headers.append(line)
            else:
                splittedLine = line.split()
                ra  = float(splittedLine[3])
                dec = float(splittedLine[4])
                colourToAdd = ""

                for currentSourceWithColour in matchedCoordinatesAndcolours:
                    if ((np.abs(currentSourceWithColour[0] - ra) < TOLERANCE_DEC) and (np.abs(currentSourceWithColour[1] - dec) < TOLERANCE_DEC)):
                        colourToAdd = currentSourceWithColour[2]
                        break
                if (colourToAdd == ""):
                    colourToAdd = np.nan
                updatedLines.append(line.strip() + " " + str(colourToAdd) + "\n")

    headers.append("# Column 8: g-r       [None  ,f64,   ] g-r colour from survey (after correcting offset for matching GAIA)\n")

    # I do this in two-steps in order to be able to add the new header
    with open(catalogueOutput, "w") as file:
        file.writelines(headers)
        file.writelines(updatedLines)
    
    return()




# This script receives two catalogues, matches them and compute the colours.
# Then these colours are added to a third catalogue. The colours are mag_catalogue1 - mag_catalogue2 and are added to catalogue3 and written in catalogueOutput


TOLERANCE_DEC = 1/3600

catalogue1 = sys.argv[1]
catalogue2 = sys.argv[2]
catalogue3 = sys.argv[3]
catalogueOutput = sys.argv[4]

ra1, dec1, mag1 = readCoordsAndMagFromCatalogue(catalogue1, 3, 4, 5)
ra2, dec2, mag2 = readCoordsAndMagFromCatalogue(catalogue2, 3, 4, 5)


matchedCoordinatesAndMagnitudes = matchSources(ra1, ra2, dec1, dec2, mag1, mag2)
matchedCoordinatesAndcolours = getColoursOfMatchedSources(matchedCoordinatesAndMagnitudes)

addColoursToCatalogue(matchedCoordinatesAndcolours, catalogue3, catalogueOutput)
