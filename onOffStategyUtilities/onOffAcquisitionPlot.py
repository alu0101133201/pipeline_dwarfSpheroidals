import sys
import datetime

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.dates as mdates
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

setMatplotlibConf()

objectData = sys.argv[1].split()
flatFieldData = sys.argv[2].split() # This is the data from the field used for building the flat field
flats = sys.argv[3].split()  # This is the flats themselves
imageName = sys.argv[4]

objectDates = np.array(objectData[1::2]).astype(float)
objectDatesY = np.full(len(objectDates), 0.5)

flatFieldDates = np.array(flatFieldData[1::2]).astype(float)
flatFieldDatesY = np.full(len(flatFieldDates), 1)

flatDates   = np.array(flats[1::2]).astype(float)
flatDatesY = np.full(len(flatDates), 0.75)


fig, ax = plt.subplots(1, 1, figsize=(20, 10))
plt.tight_layout(pad=12.0)
configureAxis(ax, 'Time', '', logScale=False)
ax.set_ylim(-0.5, 1.5)

ax.scatter([datetime.datetime.fromtimestamp(ts) for ts in flatFieldDates], flatFieldDatesY, s=150, color="green", label="flat field ima")
ax.scatter([datetime.datetime.fromtimestamp(ts) for ts in objectDates], objectDatesY, s=150, color="red", label="Object field Ima")

ax.scatter([datetime.datetime.fromtimestamp(ts) for ts in flatDates], flatDatesY, s=150, color="teal", marker="D", label="Flats")

ax.scatter(datetime.datetime.fromtimestamp(flatDates[0]), 0.75, s=150, color="orange", marker="s", zorder=15, label="Border flats")
ax.scatter(datetime.datetime.fromtimestamp(flatDates[-1]), 0.75, s=150, color="orange", marker="s", zorder=15)
ax.vlines([datetime.datetime.fromtimestamp(ts) for ts in flatDates], ymin=-0.5, ymax=1.5, ls="--", color="teal", lw=1.5)
ax.legend(fontsize=20)

for label in ax.get_xticklabels():
    label.set_rotation(45)
    label.set_horizontalalignment('right')

plt.savefig(imageName)