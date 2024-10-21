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



directoryWithTheCatalogues = sys.argv[1]
imageName = sys.argv[2]

magDiff = np.array([])
mag1Total = np.array([])

for i in glob.glob(directoryWithTheCatalogues + "/*.cat"):
    print("Processing ", i)
    mag1, mag2 = read_columns_from_file(i)
    mag1Total = np.append(mag1Total, mag1)
    magDiff = np.append(magDiff, np.array(mag1 - mag2))



setMatplotlibConf()

fig, ax = plt.subplots(1, 1, figsize=(15, 15))
configureAxis(ax, "DECaLS mag (mag)", "DECaLS - TST (mag)", logScale=False)
ax.set_ylim(-1, 4)
ax.scatter(mag1Total, magDiff, color="blue", s=5)
plt.savefig(imageName + ".jpg")
