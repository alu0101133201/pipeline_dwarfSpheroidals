import os
import sys
import glob

import concurrent.futures
import threading 

from astroquery.gaia import Gaia
from astroquery.sdss import SDSS

def createRegionFileFromTable(file, table):
    firstLine = '# Region file format: DS9 version 4.1\nglobal color=red dashlist=8 3 width=4 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'
    with open(file, 'w') as f:
        f.write(firstLine)
        for i in table:
            f.write("circle(" + str(i['ra']) + "," + str(i['dec']) + ",10\")\n")
    return()

def downloadSingleSDSSSpectrum(result):
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

def downloadAllSDSSSpectra(mosaicDir, raMin, raMax, decMin, decMax):
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

    with concurrent.futures.ThreadPoolExecutor(max_workers=50) as executor:
        executor.map(downloadSingleSDSSSpectrum, result)  # Maps each row to a worker thread
    return()

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def downloadAllGaiaSpectra(mosaicDir, raMin, raMax, decMin, decMax):
    query = f"SELECT * FROM gaiadr3.gaia_source \
        WHERE (ra BETWEEN {raMin} AND {raMax}) AND (dec BETWEEN {decMin} AND {decMax}) AND \
        has_xp_sampled = 'True'"

    job     = Gaia.launch_job_async(query)
    results = job.get_results()

    createRegionFileFromTable(mosaicDir + "/spectraInField.reg", results)

    chunk_size   = 10
    ids          = results['source_id']
    ids_chunks   = list(chunks(ids, chunk_size))
    retrieval_type = 'XP_SAMPLED'
    data_structure = 'INDIVIDUAL' 
    data_release   = 'Gaia DR3'   
    datalink_all   = []

    ii = 0
    for chunk in ids_chunks:
        ii = ii + 1
        datalink  = Gaia.load_data(ids=chunk, data_release = data_release, retrieval_type=retrieval_type, format = 'fits', data_structure = data_structure, dump_to_file=True)
        datalink_all.append(datalink) 

    datalink_out = datalink_all[0]
    for inp_dict in datalink_all[1:]:
        datalink_out.update(inp_dict)

    os.system(f"mv ./datalink_output*.zip {SPECTRA_DIR}")
    for file in glob.glob(f"{SPECTRA_DIR}/*.zip"):
        os.system(f"unzip {file} -d {SPECTRA_DIR}")
        os.system(f"rm {file}")
    return()

mosaicDir   = sys.argv[1]
spectraDir  = sys.argv[2]
ra          = float(sys.argv[3])
dec         = float(sys.argv[4])
sizeOfField = float(sys.argv[5])
survey      = sys.argv[6]


SPECTRA_DIR = spectraDir

halfSizeOfFild = sizeOfField / 2

raMin  = ra - halfSizeOfFild
raMax  = ra + halfSizeOfFild
decMin = dec - halfSizeOfFild
decMax = dec + halfSizeOfFild

if (survey == "SDSS"):
    downloadAllSDSSSpectra(mosaicDir, raMin, raMax, decMin, decMax)
elif (survey == "GAIA"):
    downloadAllGaiaSpectra(mosaicDir, raMin, raMax, decMin, decMax)
else:
    raise Exception("Survey for downloading spectra not recognised")