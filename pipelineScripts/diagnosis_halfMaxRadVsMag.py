# This script receives a table, extracts the half-max-radius and magnitudes and
# saves the plot in the diagnosis folder

import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style

import numpy as np

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
    ax.tick_params(axis='x', which='major', labelsize=25, pad=17)
    ax.tick_params(axis='y', which='major', labelsize=25, pad=17)
    ax.set_xlabel(xlabel, fontsize=30, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=30, labelpad=10)
    if(logScale): ax.set_yscale('log')

def readHalfMaxRadAndMag(file, halfMaxRadCol, magnitudeColumn):
    data = np.loadtxt(file, comments='#', usecols=(halfMaxRadCol, magnitudeColumn))
    
    half_max_radius = data[:, 0]
    magnitude = data[:, 1]
    
    return half_max_radius, magnitude

wholeTableFile   = sys.argv[1]
matchedtableFile = sys.argv[2]
meanRad          = float(sys.argv[3])
minRad           = float(sys.argv[4])
maxRad           = float(sys.argv[5])
imageOutputName  = sys.argv[6]
plotXLowerLimit  = float(sys.argv[7])
plotXHigherLimit = float(sys.argv[8])
plotYLowerLimit  = float(sys.argv[9])
plotYHigherLimit = float(sys.argv[10])
apertureUnits    = sys.argv[11]
minRange         = float(sys.argv[12])
maxRange         = float(sys.argv[13])
brightCalLimit   = float(sys.argv[14])
faintCalLimit    = float(sys.argv[15])

setMatplotlibConf()
halfMaxRadAll, magnitudeAll = readHalfMaxRadAndMag(wholeTableFile, 5, 4)
halfMaxRadMatched, magnitudeMatched = readHalfMaxRadAndMag(matchedtableFile, 4, 5)

x, y =  readHalfMaxRadAndMag(matchedtableFile, 0, 1)

fig, ax = plt.subplots(1, 1, figsize=(12, 12))
plt.tight_layout(pad=7.0)
configureAxis(ax, f'{apertureUnits} (px)', 'Magnitude (mag)', logScale=False)
ax.set_xscale('log')
ax.set_xlim(plotXLowerLimit, plotXHigherLimit)
ax.set_ylim(plotYLowerLimit, plotYHigherLimit)
ax.scatter(halfMaxRadAll, magnitudeAll, s=30, color="crimson", alpha=0.85, linewidths=1.5, edgecolor="black", label="All sources")
ax.scatter(halfMaxRadMatched, magnitudeMatched, s=40, color="blue", linewidths=1.5, edgecolor="black", label="Matched Gaia")

if (minRange and maxRange):
    ax.vlines(minRange, lw=2.5, color="black", ymin=plotYLowerLimit, ymax=plotYHigherLimit)
    ax.vlines(maxRange, lw=2.5, color="black", ymin=plotYLowerLimit, ymax=plotYHigherLimit)

if (brightCalLimit and faintCalLimit):
    ax.hlines(brightCalLimit, lw=2.5, color="grey", xmin=plotXLowerLimit, xmax=plotXHigherLimit)
    ax.hlines(faintCalLimit, lw=2.5, color="grey", xmin=plotXLowerLimit, xmax=plotXHigherLimit)

if (meanRad > 0): 
    ax.vlines(x=meanRad, ymin = 10, ymax = 26, color="black", lw=2.5, ls="--")
if (minRad > 0): 
    ax.vlines(x=minRad, ymin = 10, ymax = 26, color="black", lw=1.5, ls="--", label=r"Point-like region")
if (maxRad > 0): 
    ax.vlines(x=maxRad, ymin = 10, ymax = 26, color="black", lw=1.5, ls="--")
ax.legend(fontsize=18, shadow=True)
plt.savefig(imageOutputName)