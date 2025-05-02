# This script receives two filters and a field. It downloads the gaia spectra for the field and compares
# the magnitudes obtained for the same stars using both filters, then it computes the offset between magnitudes 
# i.e. the offset due to the differences between filters

# This is later used for applying this offset in the pipeline when calibrating with PANSTARRS as a proxy to the GAIA spectra

import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from astropy.io import fits
from astropy.stats import sigma_clip
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style

def setMatplotlibConf():
    rc_fonts = {
        "font.family": "serif",
        "font.size": 14,
        "font.weight" : "medium",
        # "text.usetex": True,  # laggs a little when generatin plots in my fedora36
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 8.0,
        "xtick.major.width": 2.8,
        "xtick.minor.size": 4.0,
        "xtick.minor.width": 2.5,
        "ytick.major.size": 8.0,
        "ytick.major.width": 1.8,
        "ytick.minor.size": 4.0,
        "ytick.minor.width": 1.8,
        "legend.handlelength": 3.0,
        "axes.linewidth" : 3.5,
        "xtick.major.pad" : 6,
        "ytick.major.pad" : 6,
        "legend.fancybox" : True,
        "mathtext.fontset" : "dejavuserif"
    }
    mpl.rcParams.update(rc_fonts)
    return(rc_fonts)

def configureAxis(ax, xlabel, ylabel, logScale=True):
    ax.xaxis.set_minor_locator(MultipleLocator(1000000))
    ax.yaxis.set_minor_locator(MultipleLocator(1000000))
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis='x', which='major', labelsize=20, pad=17)
    ax.tick_params(axis='y', which='major', labelsize=20, pad=17)
    ax.set_xlabel(xlabel, fontsize=25, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=25, labelpad=10)
    if(logScale): ax.set_yscale('log')

def readFilterTransmittance(fileWithFilterTransmittance, waveUnits, transmittanceUnits):
    wavelengths = []
    transmittance = []

    with open(fileWithFilterTransmittance, 'r') as f:
        for line in f:
            splittedLine = line.split()
            wavelengths.append(float(splittedLine[0]))
            transmittance.append(float(splittedLine[1]))

    wavelengths = np.array(wavelengths)
    transmittance = np.array(transmittance)

    if (waveUnits == "nm"):
        wavelengths = wavelengths * 10
    elif (waveUnits == "A"):
        pass
    else:
        raise Exception(f"Units for wavelength ({waveUnits}) not accepted")

    if (transmittanceUnits == "percentage"):
        transmittance = transmittance / 100
    elif (transmittanceUnits == "normalised"):
        pass
    else:
        raise Exception(f"Units for transmittance ({waveUnits}) not accepted")

    return(wavelengths, transmittance)

def getRaAndDecFromSpectrum(spectrumFile):
    with fits.open(spectrumFile) as hdul:
        header = hdul[1].header
        ra = float(header["POS"].split(',')[0][1:])
        dec = float(header["POS"].split(',')[1][:-1])
    return(ra, dec)

def getSpectrumData(specFile):
    sp = fits.open(specFile)
    data = sp[1].data

    wavelength = 10 * data["wavelength"] # Convert wavelength from nm to Anstroms
    flux = np.array(data["flux"]) * 100  # Convert from W m^(-2) nm^(-1) to erg cm(^-2) s(^-1) A(^-2)
    return(wavelength, flux)

def obtainFluxFromSpectrumAndFilter(wavelengths, spectrumF,  filterT, plot=False):
    for i in range(len(filterT)):
        spectrumF[i] = spectrumF[i] * (wavelengths[i]**2) / 3e18 # Since the flux is in wavelenths, we need this factor to move it obtain the flux in frequency
                                                                # the 3e18 is the lightspeed in angstroms
        spectrumF[i] = spectrumF[i] / 1e-23   # This factor is for going to Yanskins we need to multiply by 10**(-23)

    convolvedSpec = []
    for i in range(len(spectrumF)):
        convolvedSpec.append(spectrumF[i] * filterT[i])

    total_flux = np.trapz(convolvedSpec, wavelengths)
    dividedFlux = total_flux / np.trapz(filterT, wavelengths)
    return(-2.5*np.log10(dividedFlux / 3631))


def getMagnitudesFromSpectra(spectraFolder, wavelengths, transmittance_tmp):
    magnitudes = []
    ra         = []
    dec        = []

    interp_func = interp1d(wavelengths, transmittance_tmp, bounds_error=False, fill_value=0)
    transmittance = interp_func(WAVELENGTHS_TO_SAMPLE)

    for i in glob.glob(spectraFolder + "/*.fits"):
            ra_tmp, dec_tmp = getRaAndDecFromSpectrum(i)
            ra.append(ra_tmp)
            dec.append(dec_tmp)

            specWavelenth, specFlux_tmp = getSpectrumData(i)
            interp_func = interp1d(specWavelenth, specFlux_tmp, bounds_error=False, fill_value=0)
            specFlux = interp_func(WAVELENGTHS_TO_SAMPLE)

            currentMag = obtainFluxFromSpectrumAndFilter(WAVELENGTHS_TO_SAMPLE, specFlux, transmittance)
            magnitudes.append(currentMag)
    return(np.array(ra), np.array(dec), np.array(magnitudes))

def comparisonPlot(mag1, mag2, waveUnits1, transmittanceUnits1, filter1, waveUnits2, transmittanceUnits2, filter2, offset, regionToComputeOffset, field):
    filterName1 = filter1.split("/")[-1].split(".")[0]
    filterName2 = filter2.split("/")[-1].split(".")[0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 12))
    plt.tight_layout(pad=6.0)

    configureAxis(ax1, 'Wavelength (A)', 'Transmittance', logScale=False)
    ax1.set_title("Filters", pad=15, fontsize=20)
    ax1.plot(waveUnits1, transmittanceUnits1, color="red", lw=2.5, label="1 " + filterName1)
    ax1.plot(waveUnits2, transmittanceUnits2, color="blue", lw=2.5, label="2 " + filterName2)

    configureAxis(ax2, f'mag {filterName1}', f'mag {filterName1} - mag {filterName2}', logScale=False)
    ax2.set_title(f"Magnitude comparison - field {field}", pad=15, fontsize=20)
    ax2.set_xlim(10, 16)
    ax2.set_ylim(-0.5, 0.5)
    ax2.scatter(mag1, mag1 - mag2, s=80, color="teal", edgecolors="black", linewidths=1.5)
    ax2.plot([7, 16], [0, 0], lw=2, ls="--", color="grey")

    ax2.text(10.5, 0.2, f"Offset ({filterName1} - {filterName2}): " + "{:.2}".format(offset), fontsize=18, color="blue")
    ax2.vlines(x=regionToComputeOffset[0], ymin=-0.5, ymax=0.5, color="black", lw=1.5, ls="--", label="Region to compute offset")
    ax2.vlines(x=regionToComputeOffset[1], ymin=-0.5, ymax=0.5, color="black", lw=1.5, ls="--")
    ax2.hlines(y=offset, xmin=7, xmax=16, color="red", lw=2, ls="-.", label="offset")

    ax1.legend(fontsize=15, loc="upper right")
    ax2.legend(fontsize=15, loc="upper right")

    plt.savefig("./figure.png")

def magInRange(mag, brighLimit, faintLimit):
    if ( (mag > brighLimit) and (mag < faintLimit)):
        return(True)
    return False

def computeOffset(mag1, mag2, limits):
    diffs = []
    for i in range(len(mag1)):
        if (magInRange(mag1[i], limits[0], limits[1]) and magInRange(mag2[i], limits[0], limits[1])):
            diffs.append(mag1[i] - mag2[i])

    clippedDiffs = sigma_clip(diffs, sigma=3, masked=False)
    return(np.nanmean(clippedDiffs))


setMatplotlibConf()

field = "NGC5308"
filterName1 = "./filters/panstarrs_i.dat"; waveUnits1 = "A";  transmittanceUnits1 = "normalised"
filterName2 = "./filters/WHT_TWFC_i.dat";       waveUnits2 = "A"; transmittanceUnits2 = "normalised"

regionToComputeOffset = [12, 15]

# Here the magnitudes using both filters are computed
WAVELENGTHS_TO_SAMPLE = np.linspace(3000, 11000, 10000) # Needed in order to have the same wavelengths in filter and spectra
spectraFolder = f"./gaiaSpectra_{field}"

wavelengths1, transmittance1 = readFilterTransmittance(filterName1, waveUnits1, transmittanceUnits1)
wavelengths2, transmittance2 = readFilterTransmittance(filterName2, waveUnits2, transmittanceUnits2)

ra1, dec1, magnitudes1 = getMagnitudesFromSpectra(spectraFolder, wavelengths1, transmittance1)  
ra2, dec2, magnitudes2 = getMagnitudesFromSpectra(spectraFolder, wavelengths2, transmittance2)

offset = computeOffset(magnitudes1, magnitudes2, regionToComputeOffset)
comparisonPlot(magnitudes1, magnitudes2, wavelengths1, transmittance1, filterName1, wavelengths2, transmittance2, filterName2, offset, regionToComputeOffset, field)
