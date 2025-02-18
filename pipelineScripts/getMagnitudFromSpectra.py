import os
import re
import sys
import glob

import numpy as np

from astropy.io import fits
from scipy.interpolate import interp1d


import matplotlib.pyplot as plt

def getRaAndDecFromSpectrum(spectrumFile, surveyForSpectra):
    with fits.open(spectrumFile) as hdul:
        if (surveyForSpectra == "SDSS"):
            header = hdul[0].header
            ra = hdul[0].header.get("PLUG_RA")
            dec = hdul[0].header.get("PLUG_DEC")
        elif (surveyForSpectra == "GAIA"):
            header = hdul[1].header
            ra = float(header["POS"].split(',')[0][1:])
            dec = float(header["POS"].split(',')[1][:-1])
        else:
            raise Exception("Survey for spectra not recognised")
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

    if (surveyForSpectra == "SDSS"):
        wavelength = 10 ** data['loglam']  # Convert log(wavelength) to wavelength in Ångstroms
        flux = data['flux'] * 1e-17        # The sdss is given with a factor / 1e-17
    elif (surveyForSpectra == "GAIA"):
        wavelength = 10 * data["wavelength"] # Convert wavelength from nm to Anstroms
        flux = np.array(data["flux"]) * 100  # Convert from W m^(-2) nm^(-1) to erg cm(^-2) s(^-1) A(^-2)
    else:
        raise Exception("Survey for spectra not recognised")
    return(wavelength, flux)

def obtainFluxFromSpectrumAndFilter(wavelengths, spectrumF,  filterT, plot=False):
    for i in range(len(filterT)):
        spectrumF[i] = spectrumF[i] * (wavelengths[i]**2) / 3e18 # Since the flux is in wavelenths, we need this factor to move it obtain the flux in frequency
                                                                # the 3e18 is the lightspeed in angstroms
        spectrumF[i] = spectrumF[i] / 1e-23

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

def getMagnitudesFromSpectra(spectraFolder, filterFile, transmittanceUnits, surveyForSpectra):
    magnitudes  = []
    ra          = []
    dec         = []

    wavelengths, transmittance_tmp = readFilterTransmittance(filterFile)
    if (transmittanceUnits == "nm"): wavelengths = wavelengths * 10
    interp_func = interp1d(wavelengths, transmittance_tmp, bounds_error=False, fill_value=0)
    transmittance = interp_func(WAVELENGTHS_TO_SAMPLE)

    for i in glob.glob(spectraFolder + "/*.fits"):
        instrument_tmp, objectType_tmp = extractInstrumentAndObjectType(i)

        ra_tmp, dec_tmp = getRaAndDecFromSpectrum(i, surveyForSpectra)
        ra.append(ra_tmp)
        dec.append(dec_tmp)

        specWavelenth, specFlux_tmp = getSpectrumData(i)
        interp_func = interp1d(specWavelenth, specFlux_tmp, bounds_error=False, fill_value=0)
        specFlux = interp_func(WAVELENGTHS_TO_SAMPLE)

        currentMag = obtainFluxFromSpectrumAndFilter(WAVELENGTHS_TO_SAMPLE, specFlux, transmittance)
        magnitudes.append(currentMag)

    return(np.array(ra), np.array(dec), np.array(magnitudes))

def createCatalogueFromTable(filename, ra, dec, mag, counts):
    with open(filename, 'w') as f:

        header1 = f"# Column 1: IDs\n"
        header2 = f"# Column 2: RA\n"
        header3 = f"# Column 3: DEC\n"
        header4 = f"# Column 4: MAGNITUDE\n"
        header5 = f"# Column 5: SUM\n"

        f.write(header1 + header2 + header3 + header4 + header5)

        for i in range(len(ra)):
            f.write(str(i) + "\t" + str(ra[i]) + "\t" + str(dec[i]) + "\t" + str(mag[i]) + "\t" + str(counts[i])  + "\n")
    return()

def getCountsFromMagnitude(mag, zp):
    tmp = (mag - zp) / -2.5
    return(10**tmp)


spectraDir         = sys.argv[1]
transmittanceFile  = sys.argv[2]
transmittanceUnits = sys.argv[3]
catalogueName      = sys.argv[4]
surveyForSpectra   = sys.argv[5]


WAVELENGTHS_TO_SAMPLE = np.linspace(3000, 11000, 10000) # Needed in order to have the same wavelengths in filter and spectra
DESIRED_ZP = 22.5

ra, dec, magnitudes = getMagnitudesFromSpectra(spectraDir, transmittanceFile, transmittanceUnits, surveyForSpectra)


equivalentcountsInZP = getCountsFromMagnitude(magnitudes, DESIRED_ZP)

createCatalogueFromTable(catalogueName, ra, dec, magnitudes, equivalentcountsInZP)

