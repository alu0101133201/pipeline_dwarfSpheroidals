import glob

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style
from astropy.io import fits
import os
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

def read_columns_from_file(file_path,h):
    """
    try:
        data = np.loadtxt(file_path, comments='#')
        if (data.size == 0):
            raise Exception("File: " + file_path + " is empty")

        ra1  = np.array(data[:, 0].astype(float).tolist())
        dec1 = np.array(data[:, 1].astype(float).tolist())
        ra2  = np.array(data[:, 2].astype(float).tolist())
        dec2 = np.array(data[:, 3].astype(float).tolist())
        return ra1, dec1, ra2, dec2

    except Exception as e:
        print(f"Tried to load an empty file: {e}")
        return np.array([]), np.array([]), np.array([]), np.array([]), 
    """
    try:
        #Data is in .cat and we want only the .fits
        
        data = fits.open(file_path)[h].data
        
        if (data.size == 0):
            raise Exception("File: " + file_path + " is empty")
        ra1=data['RA_CALIBRATED']
        dec1=data['DEC_CALIBRATED']
        ra2=data['RA_NONCALIBRATED']
        dec2=data['DEC_NONCALIBRATED']
        return ra1, dec1, ra2, dec2
    except Exception as e:
        print(f"Tried to load an empty file: {e}")
        return np.array([]), np.array([]), np.array([]), np.array([]) 



cataloguesDir = sys.argv[1]
imageName     = sys.argv[2]
pixelScale    = sys.argv[3]
h             =int(sys.argv[4])
raArrays = []
decArrays = []

count=0

for i in glob.glob(cataloguesDir + "/*.cat"):
    ra1, dec1, ra2, dec2 = read_columns_from_file(i,h)

    raArrays.append(ra1-ra2)
    decArrays.append(dec1-dec2)


for i in range(len(raArrays)):
    raArrays[i] = raArrays[i]*3600
    decArrays[i] = decArrays[i]*3600


setMatplotlibConf()

fig, ax = plt.subplots(1, 1, figsize=(15, 15))
ax.set_title("Checking astrometry. The data pixel scale is " + str(pixelScale) + " (arcsec/px)", fontsize=22, pad=17)
plt.tight_layout(pad=8)
configureAxis(ax, r'$\delta$ ra (arcsec)', r'$\delta$ dec (arcsec)', logScale=False)
for i in range(len(raArrays)):
    ax.scatter(raArrays[i], decArrays[i], s=50, linewidths=1.5, edgecolor="black")
ax.set_xlim(-float(pixelScale)*3, float(pixelScale)*3)
ax.set_ylim(-float(pixelScale)*3, float(pixelScale)*3)

ax.hlines(y=-float(pixelScale), xmin=-float(pixelScale)*3, xmax=float(pixelScale)*3, color="black", linestyle="--", label="Pixelscale")
ax.hlines(y=+float(pixelScale), xmin=-float(pixelScale)*3, xmax=float(pixelScale)*3, color="black", linestyle="--")

ax.vlines(x=-float(pixelScale), ymin=-float(pixelScale)*3, ymax=float(pixelScale)*3, color="black", linestyle="--")
ax.vlines(x=float(pixelScale), ymin=-float(pixelScale)*3, ymax=float(pixelScale)*3, color="black", linestyle="--")

plt.legend(fontsize=20)
plt.savefig(imageName)


# for i in range(len(raArrays)):
#     fig, ax = plt.subplots(1, 1, figsize=(15, 15))
#     ax.set_title("Checking astrometry. The data pixel scale is " + str(pixelScale) + " (arcsec/px)", fontsize=22, pad=17)
#     plt.tight_layout(pad=8)
#     configureAxis(ax, r'$\delta$ ra (arcsec)', r'$\delta$ dec (arcsec)', logScale=False)
#     ax.scatter(raArrays[i], decArrays[i], color="teal", s=50, linewidths=1.5, edgecolor="black")
#     ax.set_xlim(-1.5, 1.5)
#     ax.set_ylim(-1.5, 1.5)
#     ax.hlines(y=0, xmin=-1.5, xmax=1.5, color="black", linestyle="--")
#     ax.vlines(x=0, ymin=-1.5, ymax=1.5, color="black", linestyle="--")
#     plt.savefig(str(i) + "_astrometry.png")