import sys
import glob

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style
from astropy.stats import sigma_clip

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
                matchedSources.append((mag1[i], mag2[j]))

    return(matchedSources)

def magnitudeComparisonPlot(magnitudeMatched, outputName, brightLimit, faintLimit, offset=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 10))

    if (offset != None ):
        fig.suptitle(f"Offset computed between " + "{:.1f}".format(brightLimit) + " and " + "{:.1f}".format(faintLimit) + ": " + "{:.4f}".format(offset), fontsize=22)
    else:
        fig.suptitle(f"Comparison before correcting", fontsize=22)

    configureAxis(ax1, 'GAIA', 'PANSTARRS', logScale=False)
    ax1.set_xlim(12, 16)
    ax1.set_ylim(12, 16)

    ax1.plot([12, 16], [12, 16], lw=2, ls="dotted", color="black")
    for i in magnitudeMatched:
        ax1.scatter(i[0], i[1], s=50, color="teal", edgecolor="black", lw=1.5)

    ax1.vlines(x=brightLimit, ymin=12, ymax=16, color="gray", lw=1.5, ls="-.")
    ax1.vlines(x=faintLimit, ymin=12, ymax=16,  color="gray", lw=1.5, ls="-.")

    configureAxis(ax2,  'GAIA', 'GAIA - PANSTARRS', logScale=False)
    ax2.hlines(y=0, xmin=12, xmax=16, lw=2, ls="dotted", color="black")
    ax2.set_xlim(12, 16)
    ax2.set_ylim(-1, 1)

    for i in magnitudeMatched:
        ax2.scatter(i[0], i[0] - i[1], s=50, color="teal", edgecolor="black", lw=1.5)

    ax2.vlines(brightLimit, ymin=-1, ymax=1, color="gray", lw=1.5, ls="-.")
    ax2.vlines(faintLimit, ymin=-1,  ymax=1, color="gray", lw=1.5, ls="-.")

    plt.tight_layout(pad=2.0)
    plt.savefig(outputName)

def getOffsetToCorrect(matchedSources, brightLimit, faintLimit):
    diffs = np.array([i[0] - i[1] for i in matchedSources])
    clippedDiffs = sigma_clip(diffs, sigma=3, masked=False)
    return(np.nanmean(clippedDiffs))

def correctPanstarrsMag(magnitudeMatched, offset):
    correctedMags = []
    for i in matchedSources:
        correctedMags.append((i[0], i[1] + offset))
    return(correctedMags)

setMatplotlibConf()

panstarrsCataloguesFile = sys.argv[1]
gaiaCatalogueFile       = sys.argv[2]
brightLimit             = float(sys.argv[3])
faintLimit              = float(sys.argv[4])
mosaicDir               = sys.argv[5]

TOLERANCE_DEC = 1.5/3600


gaiaRA, gaiaDEC, gaiaMag = readCoordsAndMagFromCatalogue(gaiaCatalogueFile, 1, 2, 3)
panstarrsRA, panstarrsDEC, panstarrsMag = readCoordsAndMagFromCatalogue(panstarrsCataloguesFile, 3, 4, 5)

matchedSources = matchSources(gaiaRA, panstarrsRA, gaiaDEC, panstarrsDEC, gaiaMag, panstarrsMag)

outputName = mosaicDir + "/beforeOffsetCorrection.png"
magnitudeComparisonPlot(matchedSources, outputName, brightLimit, faintLimit)

offset = getOffsetToCorrect(matchedSources, brightLimit, faintLimit)
correctedMagnitudes = correctPanstarrsMag(matchedSources, offset)

outputName = mosaicDir + "/afterOffsetCorrection.png"
magnitudeComparisonPlot(correctedMagnitudes, outputName, brightLimit, faintLimit, offset)

factorToApplyToCounts = 10**(-offset / 2.5)
print(offset, factorToApplyToCounts)