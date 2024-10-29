import sys
import glob

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
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
    ax.tick_params(axis='x', which='major', labelsize=25, pad=17)
    ax.tick_params(axis='y', which='major', labelsize=25, pad=17)
    ax.set_xlabel(xlabel, fontsize=30, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=30, labelpad=10)
    if(logScale): ax.set_yscale('log')

def read_columns_from_file(file_path):
    data = np.loadtxt(file_path, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    mag1  = np.array(data[:, 0].astype(float).tolist())
    mag2 = np.array(data[:, 1].astype(float).tolist())
    return mag1, mag2

def getMagnitudeDiffScatterInMagnitudeRange(mag, magDiff, faintLimit, brightLimit):
    diffMagInRange = []

    for i in range(len(mag)):
        if ( (mag[i] > brightLimit) and (mag[i] < faintLimit) ):
            diffMagInRange.append(magDiff[i])

    return(np.nanstd(np.array(diffMagInRange)))


directoryWithTheCatalogues = sys.argv[1]
imageName = sys.argv[2]
calibrationBrightLimit = float(sys.argv[3])
calibrationFaintLimit  = float(sys.argv[4])

magDiff = np.array([])
mag1Total = np.array([])
frameNumber = np.array([])

currentFrame = 1
for index, file in enumerate(glob.glob(directoryWithTheCatalogues + "/*.cat")):
    mag1, mag2 = read_columns_from_file(file)
    mag1Total = np.append(mag1Total, mag1)
    magDiff = np.append(magDiff, np.array(mag1 - mag2))
    frameNumber = np.append(frameNumber, np.repeat(index, len(mag1)))


totalScatter = np.nanstd(magDiff)
scatterInRange = getMagnitudeDiffScatterInMagnitudeRange(mag1Total, magDiff, calibrationFaintLimit, calibrationBrightLimit)

setMatplotlibConf()

fig, ax = plt.subplots(1, 1, figsize=(15, 15))
configureAxis(ax, "DECaLS mag (mag)", "DECaLS - reduced_Data (mag)", logScale=False)
ax.set_ylim(-1, 2)
ax.set_xlim(15, 23)

ax.vlines(x=calibrationBrightLimit, ymin = -1, ymax = 2, lw=2, color="grey", linestyle="--", label="Calibration region")
ax.vlines(x=calibrationFaintLimit,  ymin = -1, ymax = 2, lw=2, color="grey", linestyle="--")

ax.hlines(y=0, xmin=15, xmax=23, color="black", lw=1.5, ls="--")
scatter = ax.scatter(mag1Total, magDiff, s=25, c=frameNumber, cmap='viridis', edgecolor="black")
cbar = plt.colorbar(scatter)
cbar.set_label('Frame Number', rotation=270, labelpad=20, fontsize=22)

ax.text(0.08, 0.925, r"Total $\sigma$: " + "{:.2f}".format(totalScatter) + " mag", transform=ax.transAxes, 
    fontsize=24, verticalalignment='top', horizontalalignment='left')

ax.text(0.08, 0.875, r"Calibration region $\sigma$: " + "{:.2f}".format((scatterInRange)) + " mag", transform=ax.transAxes, 
    fontsize=24, verticalalignment='top', horizontalalignment='left')

ax.legend(fontsize=22)
plt.savefig(imageName + ".jpg")