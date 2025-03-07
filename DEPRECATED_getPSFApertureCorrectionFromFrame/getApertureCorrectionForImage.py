import os
import sys
import glob
import math

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy.integrate as integrate

from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import SigmaClip, sigma_clip
from scipy.interpolate import interp1d
from matplotlib.colors import LogNorm
from photutils.profiles import RadialProfile
from photutils.centroids import centroid_2dg
from photutils.aperture import aperture_photometry, ApertureStats, CircularAperture, CircularAnnulus
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

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

def createFitsFunc(data, fitsName, wcs = None):
    if (wcs):
        header = wcs.to_header()
        hdu = fits.PrimaryHDU(data, header=header)
    else:
        hdu = fits.PrimaryHDU(data)
    hdul = fits.HDUList([hdu])
    hdul.writeto(fitsName, overwrite=True)
    
def readCoordsAndMagFromCatalogue(file, colRA, colDec, colMag):
    data = np.loadtxt(file, comments='#', usecols=(colRA, colDec, colMag))

    ra = np.array(data[:, 0].astype(float).tolist())
    dec = np.array(data[:, 1].astype(float).tolist())
    magnitude = np.array(data[:, 2].astype(float).tolist())

    return(ra, dec, magnitude)

def getPixelCoordinates(wcs, ra, dec):
    x = []
    y = []

    for i in range(len(ra)):
        currentRa  = ra[i]
        currentDec = dec[i]
        pixel_coords = wcs.all_world2pix(currentRa, currentDec, 1)  

        y.append(int(pixel_coords[0]))
        x.append(int(pixel_coords[1]))
    return(x, y)

def getDataAndWCSFromImage(image, hdu):
    hdul = fits.open(image)
    wcs = WCS(hdul[hdu].header)
    data = hdul[hdu].data
    hdul.close()
    return(data, wcs)

def calculateRadiusOfHalfFWHM(radius, profile):
    maxValue = profile[0] 
    halfOfMaxValue = maxValue / 2


    lasDiffValue = np.inf
    lastDiffRad = -1
    for i in range(len(radius)):
        currentDiff = np.abs(halfOfMaxValue - profile[i])
        if (currentDiff <= lasDiffValue):
            lastDiffRad = radius[i]
            lasDiffValue = currentDiff
        if (currentDiff > lasDiffValue):
            break
    return(lastDiffRad)

def plotMedianProf(profiles, medianProfRad, weightedMeanProfile, fwhm, corr, numOfStarsUsed, galaxy):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.tight_layout(pad=5.5)
    configureAxis(ax, "Radius (px)", "Flux (ADU)", logScale=False)
    ax.set_title(f"Combining {numOfStarsUsed} stars between magnitude {magnitudeRangeToUse[0]} and {magnitudeRangeToUse[1]}")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-6, 1.5)
    ax.set_xlim(0.1, 40)

    for i in range(len(profiles)):
        ax.plot(profiles[i].radius, profiles[i].profile / np.nanmax(profiles[i].profile), color='grey')
    ax.plot(medianProfRad, weightedMeanProfile / np.nanmax(weightedMeanProfile), color='blue', linewidth=3, label = "Weighted mean Prof")

    currentAnnulus = patches.Rectangle((INNER_ANNULUS, 1.5), (OUTER_ANNULUS - INNER_ANNULUS), -2, linewidth=0.25, facecolor="teal", alpha=0.5, label="Background annulus estimation")
    ax.add_patch(currentAnnulus)
    
    correctionMag = "{:.3}".format(-2.5*np.log10(corr))
    ax.text(0.7, 1e-4, f"Correction {correctionMag} mag", fontsize=20, color="blue")
    ax.vlines(x=(NUM_OF_FWHM * fwhm) / 2, ymin=1e-6, ymax=1.5, color="red", ls="dotted", label="Aperture used")
    ax.vlines(x=LIMIT_OF_THE_STAR_PROFILES, ymin=1e-6, ymax=1.5, color="black", ls="dotted", label=f"Limit of the prof ({LIMIT_OF_THE_STAR_PROFILES} px)")

    ax.legend(fontsize=18, shadow=True)
    plt.savefig(f"./{galaxy}/medianProf.png")

def getEECorrectionFromAperture(medianProfileRad, medianProfile, integrationLimit, aperture):
    interp_rp = interp1d(medianProfileRad, medianProfile, kind='linear', fill_value="extrapolate", bounds_error=False)
    profileSize = len(medianProfileRad)
    x, y = np.meshgrid(np.linspace(-medianProfileRad[-1], medianProfileRad[-1], 1 + (2 * profileSize)), np.linspace(-medianProfileRad[-1], medianProfileRad[-1], 1 + (2 * profileSize)))
    radiusGrid = np.sqrt(x**2 + y**2)
    profile2D = interp_rp(radiusGrid)
    # createFitsFunc(profile2D, "./test.fits")


    # We have to take into account the resolution of the profile generated
    # tmpCenter = (((REGION_SIZE - 1) * RADIAL_PROF_RESOLUTION) - 1) / 2 # profileSize
    tmpCenter = profileSize

    currentAper = CircularAperture((tmpCenter, tmpCenter), r = aperture*RADIAL_PROF_RESOLUTION)
    currentFullAperture = CircularAperture((tmpCenter, tmpCenter), r = integrationLimit*RADIAL_PROF_RESOLUTION)

    apertureFlux = aperture_photometry(profile2D, currentAper)
    totalFlux = aperture_photometry(profile2D, currentFullAperture)
    return(apertureFlux['aperture_sum'][0] / totalFlux['aperture_sum'][0])



def filterSourcesByMagRange(raOfMatches, decOfMatches, magMatched, magnitudeRangeToUse):
    filteredRa  = []
    filteredDec = []
    filteredMag = []
    for i in range(len(raOfMatches)):
        if ((magMatched[i] >= magnitudeRangeToUse[0]) and (magMatched[i] <= magnitudeRangeToUse[1])):
            filteredRa.append(raOfMatches[i])
            filteredDec.append(decOfMatches[i])
            filteredMag.append(magMatched[i])

    return(filteredRa, filteredDec, filteredMag)

def sortSourcesByMagnitude(ra, dec, mag):
    sorted_indices = np.argsort(mag)

    raSorted  = ra[sorted_indices]
    decSorted = dec[sorted_indices]
    magSorted = mag[sorted_indices]

    return(raSorted, decSorted, magSorted)

def trimImage(imageData, x, y, regionSize):
    xStart = int(x - int(regionSize/2)); xEnd = int(x + int(regionSize/2))
    yStart = int(y - int(regionSize/2)); yEnd = int(y + int(regionSize/2))

    if ( xStart < 0 ): xStart = 0
    if ( yStart < 0 ): yStart = 0

    trimmedData = imageData[xStart:xEnd, yStart:yEnd]
    return(trimmedData, xStart, yStart)

def radialProf(data, centroid, numberOfRad, index, plot=False):
    resolution = RADIAL_PROF_RESOLUTION
    edge_radii = np.arange(0, numberOfRad, (1/resolution))

    rp = RadialProfile(data, centroid, edge_radii)
    
    if (plot):
        fig, ax = plt.subplots(figsize=(15,15))
        ax.set_yscale('log')
        ax.plot(rp.radius, rp.profile)
        plt.savefig(f"./prof_{index}.png")
    return(rp)

def getRadialProfileFromCoordinates(imageOriginal, x, y, innerAnnulus, outerAnnulus):
    radialProfs = []
    totalFluxes = []

    sigmaClip = SigmaClip(sigma=3.0, maxiters=1)

    for i in range(len(x)):
        imageData = np.copy(imageOriginal)
        imageDataMask = np.isnan(imageOriginal)

        starData, xStart, yStart = trimImage(imageData, x[i], y[i], CENTROID_REGION_SIZE)

        try:
            centroid = centroid_2dg(starData)
        except:
            continue
        xCentroidGlobal = centroid[1] + xStart
        yCentroidGlobal = centroid[0] + yStart

        starDataExtended, xStartExtended, yStartExtended = trimImage(imageData, x[i], y[i], REGION_SIZE)
        xCentroidExtended = xCentroidGlobal - xStartExtended
        yCentroidExtended = yCentroidGlobal - yStartExtended

        # fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        # ax.scatter(yCentroidExtended, xCentroidExtended, s=150, color="red", marker="+", label="Computed centroid")
        # ax.imshow(starDataExtended,  norm=LogNorm(vmin=0.01, vmax=1000), cmap="viridis")
        # ax.legend(fontsize=22, shadow=True)
        # plt.savefig(f"./{i}.png")

        annulus = CircularAnnulus((yCentroidGlobal, xCentroidGlobal), r_in = innerAnnulus, r_out=  outerAnnulus)
        aperstats = ApertureStats(imageData, annulus, mask = imageDataMask, sigma_clip=sigmaClip)
        bkg_mean = aperstats.mean

        for index in range(len(starDataExtended)):
            for j in range(len(starDataExtended[0])):
                starDataExtended[index][j] = starDataExtended[index][j] - bkg_mean

        radialProfs.append(radialProf(starDataExtended, (yCentroidExtended, xCentroidExtended), int(REGION_SIZE/2), i, plot=False))

        currentAper = CircularAperture((yCentroidExtended, xCentroidExtended), r = LIMIT_OF_THE_STAR_PROFILES)
        aperstatsNormalizationAper = ApertureStats(starDataExtended, currentAper, sigma_clip=sigmaClip, mask = np.isnan(starDataExtended))
        totalFlux = aperstatsNormalizationAper.sum
        totalFluxes.append(totalFlux)
    return(radialProfs, totalFluxes)

def getMeanProfileSNRWeighted(profiles, totalFluxes):
    meanProfRadius = profiles[0].radius
    weightedMeanProf = []
    weights = []
    maximumFlux = np.nanmax(totalFluxes)

    for i in totalFluxes:
        weights.append(i/maximumFlux)
    weights = np.array(weights)

    for i in range(len(meanProfRadius)):
        currentProfValues = []
        for j in range(len(profiles)):
            tmp = profiles[j].profile[i]
            currentProfValues.append(tmp)

        currentProfValues = np.array(currentProfValues)
        nansMask = ~np.isnan(currentProfValues)
        currentProfValues = currentProfValues[nansMask]
        weights_tmp = weights[nansMask]

        currentProfValues = sigma_clip(currentProfValues, sigma = 3, maxiters=3)
        tmpMask = currentProfValues.mask
        average = np.average(currentProfValues.data[~tmpMask], weights=np.array(weights_tmp)[~tmpMask])
        weightedMeanProf.append(average)
    return(meanProfRadius, weightedMeanProf)

def getFWHMFromMedianProf(radii, medianProfile, verbose=False):
    halfFWHM = calculateRadiusOfHalfFWHM(radii, medianProfile)
    if (verbose): print("The FWHM is: ", 2 * halfFWHM)
    if (verbose): print("The 2FWHM is: ", 4*halfFWHM)

    if (verbose):
        fig, ax = plt.subplots(1, 1, figsize=(15 ,15))
        configureAxis(ax, 'radius', 'flux', logScale=False)
        ax.set_xlim(0.1, 5)
        ax.plot(radii, medianProfile, color="black")
        ax.vlines(x=halfFWHM, ymin=0, ymax=np.max(medianProfile), color='red', lw=2, label="Half-FWHM")
        # ax.vlines(x=2*halfFWHM, ymin=0, ymax=np.max(medianProfile), color='red', lw=2, ls="--", label="FWHM")
        plt.legend(fontsize=18)
        plt.savefig("fwhmCalculus.png")
    return(2*halfFWHM)

def getFWHMAndEeCorrection(medianProfileRad, medianProfile, integrationLimit, verbose=False):
    fwhm = getFWHMFromMedianProf(medianProfileRad, medianProfile, verbose)
    print("\n\nTHE FWHM MEASURED IS: " + str(fwhm))

    # apertureRadius = fwhm # 2FWHM of diameter
    apertureRadius = ((NUM_OF_FWHM) * fwhm) / 2 # 

    
    print("\n\n")
    eeAperture = getEECorrectionFromAperture(medianProfileRad, medianProfile, integrationLimit, apertureRadius)
    return(fwhm, eeAperture, apertureRadius)

def line(x, a, b):
    return(a*x + b)

def exponential(x, a, b):
    return((10**b) * (x**a))

def plotProfileAndFit(x, y, indices, slopeFitted, constantCoeff):

    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    configureAxis(ax, '', '', logScale=False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-6, 1.5)
    ax.set_xlim(0.1, 40)
    ax.plot(x, y / np.nanmax(y), color="blue", linewidth=2, label="Mean profile")

    ax.vlines(x=x[indices[0]], ymin=1e-8, ymax=5, color="black", linestyle="--", label="Region to fit")
    ax.vlines(x=x[indices[1]], ymin=1e-8, ymax=5, color="black", linestyle="--")
    ax.axvspan(x[indices[0]], x[indices[1]], alpha=0.5, color='orange')

    xValues = np.linspace(x[indices[0]], 100, 50)
    ax.plot(xValues, exponential(xValues, slopeFitted, constantCoeff) / np.nanmax(y), color='red', linewidth=2, label='Fitted exponential')

    ax.text(0.5, 1e-4, s="Slope: " + "{:.2f}".format(slopeFitted), fontsize=20)
    ax.legend(fontsize=20)
    plt.savefig("./testExponential.png")
    return()

def summaryPlot(profiles, medianProfRad, weightedMeanProfile, FWHMInitial, corr, galaxy, finalProfRad, finalProf, apertureCorrectionFinal, FWHMFinal):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 12))
    plt.tight_layout(pad=5.5)

    configureAxis(ax1, "Radius (px)", "Flux (ADU)", logScale=False)
    ax1.set_title(f"Field of {galaxy}")
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim(1e-6, 1.5)
    ax1.set_xlim(0.1, 40)

    for i in range(len(profiles)):
        ax1.plot(profiles[i].radius, profiles[i].profile / np.nanmax(profiles[i].profile), color='grey')
    ax1.plot(medianProfRad, weightedMeanProfile / np.nanmax(weightedMeanProfile), color='blue', linewidth=3, label = "Weighted mean Prof")

    currentAnnulus = patches.Rectangle((INNER_ANNULUS, 1.5), (OUTER_ANNULUS - INNER_ANNULUS), -2, linewidth=0.25, facecolor="teal", alpha=0.5, label="Background annulus estimation")
    ax1.add_patch(currentAnnulus)
    
    correctionMag = "{:.3}".format(-2.5*np.log10(corr))
    ax1.text(0.4, 1e-4, f"Correction {correctionMag} mag", fontsize=20, color="blue")
    ax1.vlines(x=(NUM_OF_FWHM * FWHMInitial) / 2, ymin=1e-6, ymax=1.5, color="red", ls="dotted", label="Aperture used")
    ax1.vlines(x=LIMIT_OF_THE_STAR_PROFILES, ymin=1e-6, ymax=1.5, color="black", ls="dotted", label=f"Limit of the prof ({LIMIT_OF_THE_STAR_PROFILES} px)")


    configureAxis(ax2, "Radius (px)", "Flux (ADU)", logScale=False)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(1e-6, 1.5)
    ax2.set_xlim(0.1, 40)

    ax2.plot(finalProfRad, finalProf / np.nanmax(finalProf), color="blue", lw=1.5)
    ax2.vlines(x=(NUM_OF_FWHM * FWHMFinal) / 2, ymin=1e-6, ymax=1.5, color="red", ls="dotted")

    correction = "{:.3}".format(apertureCorrectionFinal)
    ax2.text(0.4, 5e-4, f"Correction {correction}", fontsize=20, color="blue")

    correctionMag = "{:.3}".format(-2.5*np.log10(apertureCorrectionFinal))
    ax2.text(0.4, 1e-4, f"Correction {correctionMag} mag", fontsize=20, color="blue")

    ax1.legend(fontsize=18, shadow=True)
    plt.savefig(f"./{galaxy}/summaryPlot.png")

NUM_OF_FWHM = 4
INNER_ANNULUS = 35
OUTER_ANNULUS = 40

CENTROID_REGION_SIZE = 20
REGION_SIZE = 50
# LIMIT_OF_THE_STAR_PROFILES = 8
LIMIT_OF_THE_STAR_PROFILES = 5

RADIAL_PROF_RESOLUTION = 10

setMatplotlibConf()


GALAXY = "UGC5364"

# magnitudeRangeToUse = [15.5, 16]
magnitudeRangeToUse = [14.5, 15]

print("GALAXY USED: ", GALAXY)

# The calibrated catalogues are needed for selecting stars in a magnitude range 
if (GALAXY == "IC1613"):
    # image   = "/home/sguerra/IC1613/build_panstarrs_till18_withCoadd/coadds/IC1613_coadd_g.fits"
    # image = "/home/sguerra/IC1613/build_onlyStars_till18_withCoadd/coadds/IC1613_coadd_g.fits"

    image = "/home/sguerra/IC1613/build_onlyStars_till18_withCoadd/photCorrSmallGrid-dir_it1/entirecamera_1.fits"

    imageDataHdu = 1
    catalogueDir = "/home/sguerra/IC1613/build_panstarrs_till18_withCoadd/calibratedCatalogues"

elif (GALAXY == "UGC5364"):
    image = "/home/sguerra/UGC5364/build_panstarrs/coadds/UGC5364_coadd_g.fits"
    # image = "/home/sguerra/UGC5364/build_panstarrs/photCorrSmallGrid-dir_it1/entirecamera_1.fits"
    
    imageDataHdu = 1
    catalogueDir = "/home/sguerra/UGC5364/build_panstarrs/calibratedCatalogues"

elif (GALAXY == "NGC598"):
    image = "/home/sguerra/NGC598/build_panstarrs_18/coadds/NGC598_coadd_g.fits"
    imageDataHdu = 1
    catalogueDir = "/home/sguerra/NGC598/build_panstarrs_18/calibratedCatalogues"
else:
    raise Exception("Galaxy not supported")

data, wcs = getDataAndWCSFromImage(image, imageDataHdu)



mergedCatalogue = catalogueDir + "/mergedCatalogued.cat"
os.system(f"rm {mergedCatalogue}")
file_list = glob.glob(catalogueDir + "/*.cat")
df_list = [pd.read_csv(f, delim_whitespace=True, header=None) for f in file_list]
df = pd.concat(df_list, ignore_index=True)
df = df.drop_duplicates()
df.to_csv(mergedCatalogue, sep=" ", index=False, header=False)

ra, dec, mag = readCoordsAndMagFromCatalogue(mergedCatalogue, 0, 1, 2)
ra, dec, mag = filterSourcesByMagRange(ra, dec, mag,magnitudeRangeToUse)
numOfStarsUsed = len(ra)

print(numOfStarsUsed)

xCoords, yCoords = getPixelCoordinates(wcs, ra, dec)
radialProfs, fluxes = getRadialProfileFromCoordinates(data, xCoords, yCoords, INNER_ANNULUS, OUTER_ANNULUS)

meanProfRad, meanProf = getMeanProfileSNRWeighted(radialProfs, fluxes)


FWHMInitial, apertureCorrectionInitial, apertureUsed = getFWHMAndEeCorrection(meanProfRad, meanProf, LIMIT_OF_THE_STAR_PROFILES, verbose=False)
print(f"Using an aperture of {apertureUsed}")
print(FWHMInitial, apertureCorrectionInitial)
print("correction in magnitudes: ", -2.5*np.log10(apertureCorrectionInitial))

# plotMedianProf(radialProfs, meanProfRad, meanProf, FWHMInitial, apertureCorrectionInitial, numOfStarsUsed, GALAXY)

# Correction with exponential
# fitRegion = (4, 8)
fitRegion = (2, 6)


fitIndices = []
for count, i in enumerate(meanProfRad):
    if (i > fitRegion[0]):
        fitIndices.append(count)
        break
for count, i in enumerate(meanProfRad):
    if (i > fitRegion[1]):
        fitIndices.append(count)
        break


logX = np.log10(meanProfRad)
logY = np.log10(meanProf)

popt, pcov = curve_fit(line, logX[fitIndices[0]:fitIndices[1]], logY[fitIndices[0]:fitIndices[1]])
slopeFitted = popt[0]
constantCoeff = popt[1]

print("slope: ", slopeFitted)
print("coeff: ", constantCoeff)
plotProfileAndFit(10**logX, 10**logY, fitIndices, slopeFitted, constantCoeff)

# Combine profile
LIMIT_OF_THE_FULL_PROF = 100

newProfRad = np.linspace(0.05, LIMIT_OF_THE_STAR_PROFILES, LIMIT_OF_THE_STAR_PROFILES*RADIAL_PROF_RESOLUTION)
interp_Prof = interp1d(meanProfRad, meanProf, kind='linear', fill_value=0, bounds_error=False)
newProf = interp_Prof(newProfRad)

extensionRad = np.linspace(LIMIT_OF_THE_STAR_PROFILES, LIMIT_OF_THE_FULL_PROF, (LIMIT_OF_THE_FULL_PROF-LIMIT_OF_THE_STAR_PROFILES)*RADIAL_PROF_RESOLUTION)
extensionProf = exponential(extensionRad, slopeFitted, constantCoeff)

finalProfRad = np.concatenate((newProfRad, extensionRad))
finalProf    = np.concatenate((newProf, extensionProf))

if (finalProf[-1] < 0):
    raise Exception("The exponential extension goes to negative values. Check this.")

# FWHMFinal, apertureCorrectionFinal, apertureUsed = getFWHMAndEeCorrection(finalProfRad, finalProf, LIMIT_OF_THE_STAR_PROFILES, verbose=False)
FWHMFinal, apertureCorrectionFinal, apertureUsed = getFWHMAndEeCorrection(finalProfRad, finalProf, LIMIT_OF_THE_FULL_PROF, verbose=False)
print(f"Using an aperture of {apertureUsed}")
print(FWHMFinal, apertureCorrectionFinal)
print("correction in magnitudes: ", -2.5*np.log10(apertureCorrectionFinal))

summaryPlot(radialProfs, meanProfRad, meanProf, FWHMInitial, apertureCorrectionInitial, GALAXY, \
            finalProfRad, finalProf, apertureCorrectionFinal, FWHMFinal)