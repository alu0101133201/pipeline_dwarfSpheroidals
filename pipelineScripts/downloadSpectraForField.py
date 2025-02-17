import sys

import concurrent.futures
import threading 

from astroquery.sdss import SDSS

def createRegionFileFromTable(file, table):
    firstLine = '# Region file format: DS9 version 4.1\nglobal color=red dashlist=8 3 width=4 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'
    with open(file, 'w') as f:
        f.write(firstLine)
        for i in table:
            f.write("circle(" + str(i['ra']) + "," + str(i['dec']) + ",10\")\n")
    return()

def downloadSpectrum(result):
    plate      = result["plate"]
    mjd        = result["mjd"]
    fiberID    = result["fiberID"]
    instrument = result["instrument"]
    objectClass = result["class"]

    try:
        sp = SDSS.get_spectra(plate=plate, mjd=mjd, fiberID=fiberID)
        if (sp != None):
            sp[0].writeto(SPECTRA_DIR + f'/spectrum_data_{instrument}_{objectClass}_{plate}_{mjd}_{fiberID}.fits', overwrite=True)
    except:
        pass
    return()

def downloadSpectra(result, max_threads=50):
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_threads) as executor:
        executor.map(downloadSpectrum, result)  # Maps each row to a worker thread

mosaicDir   = sys.argv[1]
spectraDir  = sys.argv[2]
ra          = float(sys.argv[3])
dec         = float(sys.argv[4])
sizeOfField = float(sys.argv[5])

SPECTRA_DIR = spectraDir

halfSizeOfFild = sizeOfField / 2

raMin  = ra - halfSizeOfFild
raMax  = ra + halfSizeOfFild
decMin = dec - halfSizeOfFild
decMax = dec + halfSizeOfFild

sql_query = f"""
    SELECT
        s.bestObjID,
        s.class,  
        s.survey,
        s.instrument,
        s.snMedian_r,
        s.ra, s.dec, 
        s.plate, s.mjd, s.fiberID,
        p.psfMag_g, p.psfMag_r, p.psfMag_i
        
    FROM PhotoObj AS p
    JOIN SpecObj AS s ON s.bestobjid = p.objid
    WHERE 
        class = "STAR" AND
        s.ra BETWEEN {raMin} AND {raMax}
        AND s.dec BETWEEN {decMin} AND {decMax}
"""

result = SDSS.query_sql(sql_query)
createRegionFileFromTable(mosaicDir + "/spectraInField.reg", result)

downloadSpectra(result)
