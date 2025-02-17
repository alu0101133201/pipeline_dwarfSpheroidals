import os
import re
import sys
import glob

import numpy as np

from astropy.io import fits
from scipy.interpolate import interp1d


import matplotlib.pyplot as plt

def getRaAndDecFromSpectrum(spectrumFile):
    with fits.open(spectrumFile) as hdul:
        header = hdul[0].header
        ra = header.get("PLUG_RA")
        dec = header.get("PLUG_DEC")
    return(ra, dec)

def readFilterTransmittance(fileWithFilterTransmittance):
    wavelength = []
    transmittance = []

    with open(fileWithFilterTransmittance, 'r') as f:
        for line in f:
            splittedLine = line.split()
            wavelength.append(float(splittedLine[0]))
            transmittance.append(float(splittedLine[1]))
    
    return(np.array(wavelength), np.array(transmittance))

def getSpectrumData(specFile):
    sp = fits.open(specFile)
    data = sp[1].data

    wavelength = 10 ** data['loglam']  # Convert log(wavelength) to wavelength in Ångstroms
    flux = data['flux']
    return(wavelength, flux)

def obtainFluxFromSpectrumAndFilter(wavelengths, spectrumF,  filterT, plot=False):
    for i in range(len(filterT)):
        spectrumF[i] = spectrumF[i] * (wavelengths[i]**2) / 3e18 # Since the flux is in wavelenths, we need this factor to move it obtain the flux in frequency
                                                                # the 3e18 is the lightspeed in angstroms
        spectrumF[i] = spectrumF[i] / 1e-6   # This factor is because the sloan spectra is divided by 10**(-17) and for going to Yanskins we need to multiply by 10**(-23)

    convolvedSpec = []
    for i in range(len(spectrumF)):
        convolvedSpec.append(spectrumF[i] * filterT[i])

    total_flux = np.trapz(convolvedSpec, wavelengths)
    dividedFlux = total_flux / np.trapz(filterT, wavelengths)
    return(-2.5*np.log10(dividedFlux / 3631))

def extractInstrumentAndObjectType(path):
    pattern = r"spectrum_data_(?P<instrument>[^_]+)_(?P<objectClass>[^_]+)_[^_]+_[^_]+_[^_]+\.fits"
    filename=os.path.basename(path)
    match = re.match(pattern, filename)

    if match:
        return match.group("instrument"), match.group("objectClass")
    else:
        return None, None  # Return None if the pattern doesn't match

def getMagnitudesFromSpectra(spectraFolder, filterFile, transmittanceUnits):
    instruments = []
    objectType  = []
    magnitudes  = []
    ra          = []
    dec         = []

    wavelengths, transmittance_tmp = readFilterTransmittance(filterFile)
    if (transmittanceUnits == "nm"): wavelengths = wavelengths * 10
    interp_func = interp1d(wavelengths, transmittance_tmp, bounds_error=False, fill_value=0)
    transmittance = interp_func(WAVELENGTHS_TO_SAMPLE)

    for i in glob.glob(spectraFolder + "/*.fits"):
        instrument_tmp, objectType_tmp = extractInstrumentAndObjectType(i)

        if (instrument_tmp == "SDSS"):
            instrument_tmp = 0
        elif (instrument_tmp == "BOSS"):
            instrument_tmp = 1
        else:
            instrument_tmp = -1

        if (objectType_tmp == "STAR"):
            objectType_tmp = 0
        elif (objectType_tmp == "QSO"):
            objectType_tmp = 1
        elif (objectType_tmp == "GALAXY"):
            objectType_tmp = 2
        else:
            objectType_tmp = -1
            
        instruments.append(instrument_tmp)
        objectType.append(objectType_tmp)

        ra_tmp, dec_tmp = getRaAndDecFromSpectrum(i)
        ra.append(ra_tmp)
        dec.append(dec_tmp)

        specWavelenth, specFlux_tmp = getSpectrumData(i)
        interp_func = interp1d(specWavelenth, specFlux_tmp, bounds_error=False, fill_value=0)
        specFlux = interp_func(WAVELENGTHS_TO_SAMPLE)

        currentMag = obtainFluxFromSpectrumAndFilter(WAVELENGTHS_TO_SAMPLE, specFlux, transmittance)
        magnitudes.append(currentMag)

    return(np.array(ra), np.array(dec), np.array(magnitudes), np.array(instruments), np.array(objectType))

def createCatalogueFromTable(filename, ra, dec, mag, counts, instruments, objectTypes):
    with open(filename, 'w') as f:

        header1 = f"# Column 1: IDs\n"
        header2 = f"# Column 2: RA\n"
        header3 = f"# Column 3: DEC\n"
        header4 = f"# Column 4: MAGNITUDE\n"
        header5 = f"# Column 5: SUM\n"
        header6 = f"# Column 6: INSTRUMENT\n"
        header7 = f"# Column 7: OBJECT\n"

        f.write(header1 + header2 + header3 + header4 + header5 + header6 + header7)

        for i in range(len(ra)):
            f.write(str(i) + "\t" + str(ra[i]) + "\t" + str(dec[i]) + "\t" + str(mag[i]) + "\t" + str(counts[i]) + "\t" + str(instruments[i]) + "\t" + str(objectTypes[i]) + "\n")
    return()

def getCountsFromMagnitude(mag, zp):
    tmp = (mag - zp) / -2.5
    return(10**tmp)


spectraDir         = sys.argv[1]
transmittanceFile  = sys.argv[2]
transmittanceUnits = sys.argv[3]
catalogueName      = sys.argv[4]

WAVELENGTHS_TO_SAMPLE = np.linspace(3000, 11000, 10000) # Needed in order to have the same wavelengths in filter and spectra
DESIRED_ZP = 22.5

ra, dec, magnitudes, instruments, objectTypes = getMagnitudesFromSpectra(spectraDir, transmittanceFile, transmittanceUnits)


equivalentcountsInZP = getCountsFromMagnitude(magnitudes, DESIRED_ZP)

createCatalogueFromTable(catalogueName, ra, dec, magnitudes, equivalentcountsInZP, instruments, objectTypes)

