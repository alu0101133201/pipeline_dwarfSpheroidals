# Note:
# Important to have the datalab client in order to be able to perform the queries (https://github.com/astro-datalab/datalab)
# (If you want to install it via pip, simply: pip install --ignore-installed --no-cache-dir astro-datalab )

# The queries are done in ADQL (Astronomical Data Query Language), very similar to SQL

# The table from which the data is being retrieved is "ls_dr10.bricks".
# The rest of the tables can be checked in "https://datalab.noirlab.edu/query.php"

import os
import math
import requests

import numpy as np

from io import StringIO
from astropy.io import fits
from astropy import units as u
from astropy.table import Table
from astroquery.sdss import SDSS
from dl import queryClient as qc
from astropy import coordinates as coords
from concurrent.futures import ThreadPoolExecutor

### Get brick names utilities ###



# DECALS --------------------------------------------------------------------------------

# This function retrieves the brick name to which the coordinate provided belongs
# Arguments:
#   Right ascension [degrees] [float or string]
#   Declination [degrees] [float or string]
# Returns:
#   Block name [string]
###
def getBrickNameFromSingleCoord(ra, dec):
    if ((not (isinstance(ra, (float, str))) or (not (isinstance(dec, (float, str)))))):
        raise Exception("Error in 'getBrickNameFromSingleCoord'. Arguments do not fulfill the requiered type (float or str)")
    ra = str(ra)
    dec = str(dec)

    query = 'SELECT brickname FROM ls_dr10.bricks WHERE \
                ra1 < ' + ra + ' AND ra2 > ' + ra + \
                ' AND dec1 < ' + dec + ' AND dec2 > ' + dec
    result = qc.query(query)
    return(result.split('\n')[1])

# This function retrieves the brick name/s to which the coordinate/s provided belong
# Arguments:
#   Right ascension [degrees] [list/array or float/string]
#   Declination [degrees] [list/array or float/string]
# Returns:
#   Block name/s [list or string]
###
def getBrickNamesFromCoords(ra, dec):
    isRaList = isinstance(ra, (list, np.ndarray))
    isDecList = isinstance(dec, (list, np.ndarray))
    if ((not isRaList) and (not isDecList)):
        return (getBrickNameFromSingleCoord(ra, dec))

    if ((not isRaList) or (not isDecList) or (len(ra) != len(dec))):
        raise Exception ("Error in 'getBrickNamesFromCoords'. Arguments types or lenghts do not match")

    ra = [str(value) for value in ra]
    dec = [str(value) for value in dec]
    result = []

    for currentCoords in zip(ra, dec):
        result.append(getBrickNameFromSingleCoord(currentCoords[0], currentCoords[1]))
    return(result)

# This function retrieves the brick names of the neighbors of the brick provided
# Neighbors includes diagonals, so each brick has 8 neighbors
# Arguments:
#   centreBrickName: Name of the brick [string]
#   includecentreBrickName: Determines if the brick of the centre is included in the output. False by default
# Returns:
#   Brick names of the neighbors [list]
###
def getNeighborsFromBrick(centreBrickName, includecentreBrickName=False):
    if (not isinstance(centreBrickName, str)):
        raise Exception("Error in getNeighborsFromBrick. The brick name must be a string")

    query = "SELECT ra,dec,ra1,ra2,dec1,dec2 FROM ls_dr10.bricks WHERE brickname = '" + centreBrickName + "'"
    result = qc.query(query)

    raAndDec = result.split('\n')[1]
    ra   = float(raAndDec.split(',')[0])
    dec  = float(raAndDec.split(',')[1])
    ra1  = float(raAndDec.split(',')[2])
    ra2  = float(raAndDec.split(',')[3])
    dec1 = float(raAndDec.split(',')[4])
    dec2 = float(raAndDec.split(',')[5])

    raStep = np.abs(ra2 - ra1)
    decStep = np.abs(dec2 - dec1)

    raLowerLimit = ra - raStep
    raUpperLimit = ra + raStep
    decLowerLimit = dec - decStep
    decUpperLimit = dec + decStep

    neighborsBrickNames = []
    currentRa = raLowerLimit
    while (currentRa <= raUpperLimit):
        currentDec = decLowerLimit
        while (currentDec <= decUpperLimit):
            neighborsBrickNames.append(getBrickNameFromSingleCoord(currentRa, currentDec))
            currentDec += decStep
        currentRa += raStep

    if (not includecentreBrickName):
        neighborsBrickNames.remove(centreBrickName)
    return(neighborsBrickNames)

# This function retrieves the brick names of the rectangular region defined by two (opposed) vertices
# Arguments:
#   firstPoint (ra [degrees], dec [degrees]):  First point which defines the region to retrieve bricks
#   secondPoint (ra [degrees], dec [degrees]): Second point which defines the region to retrieve bricks
# Returns:
#   result: Name of the bricks which define the region and its central coordinates
###
def getBrickNamesAndCoordinatesFromRegionDefinedByTwoPoints(firstPoint, secondPoint):
    isRaList = isinstance(firstPoint, (list, np.ndarray, tuple))
    isDecList = isinstance(secondPoint, (list, np.ndarray, tuple))
    if ((not isRaList) and (not isDecList)):
        raise Exception ("Error in 'getBrickNamesFromRegion'. Arguments have to be provided in array or tuple")
    
    firstPointRa = firstPoint[0]
    firstPointDec  = firstPoint[1]
    secondPointRa = secondPoint[0]
    secondPointDec = secondPoint[1]
    raMin, raMax = (firstPointRa, secondPointRa) if (firstPointRa < secondPointRa) else (secondPointRa, firstPointRa)
    decMin, decMax = (firstPointDec, secondPointDec) if (firstPointDec < secondPointDec) else (secondPointDec, firstPointDec)
    
    tmpQuery = 'SELECT brickname, ra, dec FROM ls_dr10.bricks WHERE \
            (ra2 > ' + str(raMin) + ' AND ra1 < ' + str(raMax) + \
            ') AND (dec2 > ' + str(decMin) + ' AND dec1 < ' + str(decMax) + ')'

    result = qc.query(tmpQuery)

    bricksNames = []
    bricksRA    = []
    bricksDec   = []
    for i in result.split()[1:]:
        values = i.split(',')
        if (len(values) > 0):
            bricksNames.append(values[0])
            bricksRA.append(float(values[1]))
            bricksDec.append(float(values[2]))
    
        
    

    return(np.array(bricksNames), np.array(bricksRA), np.array(bricksDec))

# This function retrieves the brick names of the region defined by a central point and the size of the region in ra and dec
# Arguments:
#   Centre (ra [degrees], dec [degrees]): Centre of the region
#   raSize: Size of the region in ra [degrees]
#   decSize: Size of the region in dec [degrees]
# Returns:
#   result: Name of the bricks which define the region
###
def getBrickNamesFromRegionDefinedByCentreAndSides(centre, raSize, decSize, survey):
    if (not isinstance(centre, (list, np.ndarray, tuple))):
        raise Exception ("Error in 'getBrickNamesFromRegionDefinedByCentreAndSides'. The centre must be a tuple or array")
    if ((not isinstance(raSize, (int, float))) or (not isinstance(decSize, (int, float)))):
        raise Exception ("Error in 'getBrickNamesFromRegionDefinedByCentreAndSides'. The ra and dec size must be float or int")
    
    centreRa = centre[0]
    centreDec = centre[1]

    raMin  = centreRa - (raSize / 2)
    raMax  = centreRa + (raSize / 2)
    decMin = centreDec - (decSize / 2)
    decMax = centreDec + (decSize / 2)

    result = getBrickNamesAndCoordinatesFromRegionDefinedByTwoPoints((raMin, decMin), (raMax, decMax))

    return(result)


### Download bricks utilities ###

# This function downloads the images from the brick and filters provided and stores them in the destination folder provided
# Arguments:
#   brickName: Name of the brick [string]
#   filters: Filter/s requested [string / list]
#   destinationFolder: Location where the images will be stored
#   overWrite: Defines if the file will be downloaded in the case that it already exists. True by defaults
# Returns:
#   Nothing
###
def downloadBrickDecals(brickName, filters, destinationFolder, overWrite=True):
    if (isinstance(filters, str)):
        filters = [filters]

    block = getBlockFromBrick(brickName)
    for i in filters:
        
        if (overWrite or (not os.path.exists(f"{destinationFolder}/decal_image_{brickName}_{i}.fits"))):
	    #Decals dr10 south is not complete, we check dr10 north, dr9 south and dr9 north
            base_urls = [
                "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd",
                "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/north/coadd",
                "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/south/coadd",
                "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/north/coadd"]
            success = False
            for base_url in base_urls:
                url = f"{base_url}/{block}/{brickName}/legacysurvey-{brickName}-image-{i}.fits.fz"
                result=os.system(f"wget -O {destinationFolder}/decal_image_{brickName}.{i}.fits '{url}'")
                if result==0:
                    success = True
                    break
            if not success:
                raise Exception(f"Error in 'downloadBrick': Brick {brickName} not present in Decals DR9 or DR10.")	       
    return()

    
# This function obtains the name of the block from the name of the brick
# For what I have seen the block is formed with the first three numbers of the brick name - I have not seen this in documentation, just empirical
# Arguments:
#   brickName: Name of the brick [string]
# Returns:
#   block name [string]
###

def getBlockFromBrick(brickName):
    return(brickName[:3])


# PANSTARRS --------------------------------------------------------------------------------

##Functions for get Panstarrs images
def getPanstarrsQuery(tra, tdec, size=3600, filters="grizy", format="fits", imagetypes="stack"):
    ps1filename="https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    # if imagetypes is a list, convert to a comma-separated string
    if not isinstance(imagetypes,str):
        imagetypes = ",".join(imagetypes)
    # put the positions in an in-memory file object
    cbuf = StringIO()
    cbuf.write('\n'.join(["{} {}".format(ra, dec) for (ra, dec) in zip(tra,tdec)]))
    cbuf.seek(0)
    # use requests.post to pass in positions as a file
    r = requests.post(ps1filename, data=dict(filters=filters, type=imagetypes),
        files=dict(file=cbuf))
    r.raise_for_status()
    tab = Table.read(r.text, format="ascii")
 
    
    return tab


def getPanstarrsBricksFromRegionDefinedByTwoPoints(firstPoint,secondPoint,filters,size):
    #We want panstarrs bricks of 1000 pix = 250 arcsec = 0.07deg. Field of fiew is up to now square
    isRaList = isinstance(firstPoint, (list, np.ndarray, tuple))
    isDecList = isinstance(secondPoint, (list, np.ndarray, tuple))
    if ((not isRaList) and (not isDecList)):
        raise Exception ("Error in 'getBrickNamesFromRegion'. Arguments have to be provided in array or tuple")
    
    firstPointRa = firstPoint[0]
    firstPointDec  = firstPoint[1]
    secondPointRa = secondPoint[0]
    secondPointDec = secondPoint[1]
    raMin, raMax = (firstPointRa, secondPointRa) if (firstPointRa < secondPointRa) else (secondPointRa, firstPointRa)
    decMin, decMax = (firstPointDec, secondPointDec) if (firstPointDec < secondPointDec) else (secondPointDec, firstPointDec)
    dRA=raMax-raMin #Field of fiew is up to now square
    size_deg=float(size)*0.25/3600
    nBricks=dRA // size_deg #Number of bricks per row
    if dRA % size_deg >0:
        nBricks+=1
    overlap_factor=(size_deg*nBricks-dRA)/(nBricks-1) #In a world where dRA//0.25!=0 we need to overlap in order to end in the raMax
    tra=[]; tdec=[]
    for brick in range(int(nBricks)):
        #initial position will be raMin or decMin+(0.25-overlap_factor)*brick
        iniRa=raMin+(size_deg-overlap_factor)*brick
        for brick in range(int(nBricks)):
            iniDec=decMin+(size_deg-overlap_factor)*brick
        #We want to store the central position: ini+0.25/2=ini+0.125
            tra.append(iniRa+size_deg/2); tdec.append(iniDec+size_deg/2)
    tab_panstarrs=getPanstarrsQuery(tra,tdec,size,filters="".join(filters))
    bricks_fullNames=[fname for fname in tab_panstarrs['filename']]
    bricksRA=[ra for ra in tab_panstarrs['ra']]
    bricksDec=[dec for dec in tab_panstarrs['dec']]
    bricksNames=[]
    for i in range(len(bricks_fullNames)):
        filter=tab_panstarrs['filter'][i]
        bricksNames.append("t{:08.4f}{:+07.4f}.{}.fits".format(bricksRA[i],bricksDec[i],filter))
    
    return(np.array(bricks_fullNames),np.array(bricksRA),np.array(bricksDec),np.array(bricksNames))

def getPanstarrsBricksFromCentralPoint(raCen,decCen,filters,size):
    if isinstance(raCen,np.floating):
        raCen=[raCen]; decCen=[decCen]
    tab_panstarrs=getPanstarrsQuery(raCen,decCen,size,filters="".join(filters))
    bricks_fullNames=[fname for fname in tab_panstarrs['filename']]
    bricksRA=[ra for ra in tab_panstarrs['ra']]
    bricksDec=[dec for dec in tab_panstarrs['dec']]
    bricksNames=[]
    for i in range(len(bricks_fullNames)):
        filter=tab_panstarrs['filter'][i]
        bricksNames.append("t{:08.4f}{:+07.4f}.{}.fits".format(bricksRA[i],bricksDec[i],filter))
    
    return(np.array(bricks_fullNames),np.array(bricksRA),np.array(bricksDec),np.array(bricksNames))


def downloadBrickPanstarrs(brick_fullName,brickName,brickRA,brickDEC,destinationFolder,size,overwrite=True):
    
    fitscut="https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    
    urlbase="{}?size={}&format={}".format(fitscut,size,"fits")
    url="{}&ra={}&dec={}&red={}".format(urlbase,brickRA,brickDEC,brick_fullName)
    #r=requests.get(url)
    brick_dir=destinationFolder+'/'+brickName
    #print(url+' - '+ brickName+' - '+str(len(r.content)))
    try:
        os.system(f'wget -O {brick_dir} "{url}"')
    except:
        raise Exception(f"Unable to download brick {brickName}")
    return()



# SDSS --------------------------------------------------------------------------------

def querySDSSbricks(center_ra, center_dec, radius_deg):  
    sql = f"""
    SELECT field, ra, dec
    FROM Field
    WHERE
        ACOS(
            SIN(RADIANS(dec)) * SIN(RADIANS({center_dec})) +
            COS(RADIANS(dec)) * COS(RADIANS({center_dec})) *
            COS(RADIANS(ra - {center_ra}))
        ) < RADIANS({radius_deg})
    """
    result = SDSS.query_sql(sql)
    return result


def download_field_image(field_id, ra, dec, downloadDestination, band='r', data_release='DR17'):
    position = coords.SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
    SDSS.SDSS_MIRROR = f"http://skyserver.sdss.org/{data_release.lower()}"
    
    try:
        images = SDSS.get_images(coordinates=position, radius=0.05 * u.deg, band=band)
    except Exception as e:
        print(f"SDSS query failed for field {field_id} at RA={ra}, Dec={dec}: {e}")
        return
    
    if not images:
        print(f"No image found for field {field_id} at RA={ra}, Dec={dec}")
        return
    
    filename = f"sdss_field{field_id}_{data_release.lower()}.{band}.fits"
    images[0].writeto(downloadDestination + "/" + filename, overwrite=True)
    print(f"Saved field {field_id} image to {filename}")


def download_fields_mosaic(center_ra, center_dec, radius_deg, downloadDestination, bands, data_release=17):
    bricks = querySDSSbricks(center_ra, center_dec, radius_deg)

    def download_task(row, band):
        field_id = row['field']
        ra = row['ra']
        dec = row['dec']
        download_field_image(field_id, ra, dec, downloadDestination, band=band, data_release=f"DR{data_release}")
        return f"sdss_field{field_id}_dr{data_release}.{band}", ra, dec

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(download_task, row, band) for band in bands for row in bricks]
        results = [f.result() for f in futures]

    names, ras, decs = zip(*results) if results else ([], [], [])
    return np.array(names), np.array(ras), np.array(decs)

# This function is needed because even if we are querying to SDSS database for the different fields in the area,
# some of these fields (with different identifiers) download exactly the same brick. So we download everything and
# remove the duplicates (I've no time for thinking/doing a cleaner solution)
def remove_duplicate_bricks(download_folder, tolerance_arcsec=1.0):
    tolerance_deg = tolerance_arcsec / 3600.0

    coord_map = {}
    fits_files = [f for f in os.listdir(download_folder) if f.endswith('.fits')]

    for filename in fits_files:
        filepath = os.path.join(download_folder, filename)
        try:
            with fits.open(filepath) as hdul:
                header = hdul[0].header
                ra = header.get('RA')  # RA center
                dec = header.get('DEC') # Dec center
                if ra is None or dec is None:
                    print(f"Warning: No RA/DEC in {filename}, skipping.")
                    continue
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            continue

        # Remove duplicates
        duplicate_found = False
        for existing_ra, existing_dec in coord_map.keys():
            delta_ra = abs(existing_ra - ra) * 3600.0  # in arcsec (approx)
            delta_dec = abs(existing_dec - dec) * 3600.0
            if delta_ra < tolerance_arcsec and delta_dec < tolerance_arcsec:
                print(f"Duplicate detected: {filename} matches {coord_map[(existing_ra, existing_dec)]}")
                os.remove(filepath)
                duplicate_found = True
                break

        if not duplicate_found:
            coord_map[(ra, dec)] = filename