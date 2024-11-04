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


print("Table file with all sources: ", wholeTableFile)
print("Table file with matched sources with gaia (stars): ", matchedtableFile)
print("imageoutputName: ", imageOutputName)


setMatplotlibConf()
halfMaxRadAll, magnitudeAll = readHalfMaxRadAndMag(wholeTableFile, 3, 2)
halfMaxRadMatched, magnitudeMatched = readHalfMaxRadAndMag(matchedtableFile, 2, 3)



fig, ax = plt.subplots(1, 1, figsize=(12, 12))
plt.tight_layout(pad=7.0)
ax.set_title("DECaLS objects", fontsize=20, pad=17)
configureAxis(ax, 'HALF_MAX_RADIUS (px)', 'Magnitude (mag)', logScale=False)
ax.set_xscale('log')
ax.set_ylim(11, 28)
ax.set_xlim(0.5, 40)
ax.scatter(halfMaxRadAll, magnitudeAll, s=30, color="crimson", alpha=0.85, linewidths=1.5, edgecolor="black", label="All sources")
ax.scatter(halfMaxRadMatched, magnitudeMatched, s=60, color="blue", linewidths=1.5, edgecolor="black", label="Matched Gaia")

if (meanRad > 0): 
    ax.vlines(x=meanRad, ymin = 10, ymax = 26, color="black", lw=2.5, ls="--")
if (minRad > 0): 
    ax.vlines(x=minRad, ymin = 10, ymax = 26, color="black", lw=1.5, ls="--", label=r"Point-like region $(1\sigma)$")
if (maxRad > 0): 
    ax.vlines(x=maxRad, ymin = 10, ymax = 26, color="black", lw=1.5, ls="--")
ax.legend(fontsize=18, shadow=True)
plt.savefig(imageOutputName)