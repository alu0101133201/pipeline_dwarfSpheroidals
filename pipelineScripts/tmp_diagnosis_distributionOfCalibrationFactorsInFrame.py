import os
import re
import sys
import glob
import math
import fnmatch
import astropy

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from datetime import datetime
import time

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

tableWithValues = sys.argv[1]
index = sys.argv[2]

factors = np.loadtxt(tableWithValues, comments='#', usecols=2)
# instruments = np.loadtxt(tableWithValues, comments='#', usecols=3)
# objects = np.loadtxt(tableWithValues, comments='#', usecols=4)




fig, ax = plt.subplots(1, 1, figsize=(12, 12))
configureAxis(ax, 'Factors for each star', '', logScale=False)
counts1, bins1, patches1 = ax.hist(factors, color="teal")
max_bin_height = np.max(counts1)
ax.vlines(x=np.median(factors[0]), ymin=0, ymax=max_bin_height, color="red", lw=2, ls="--", label="median")
ax.vlines(x=np.mean(factors[0]), ymin=0, ymax=max_bin_height, color="red", lw=2, ls="-", label="mean")
ax.legend(fontsize=18)

plt.savefig(f"./tmpDistributions/distributionOfFactorsInFrame_{index}.png")


exit()


sortedFactorsIns = [[], []]
for i in range(len(factors)):
    if (instruments[i] == 0):
        sortedFactorsIns[0].append(factors[i])
    else:
        sortedFactorsIns[1].append(factors[i])


sortedFactorsObj = [[], [], []]
for i in range(len(factors)):
    if (objects[i] == 0):
        sortedFactorsObj[0].append(factors[i])
    if (objects[i] == 1):
        sortedFactorsObj[1].append(factors[i])
    else:
        sortedFactorsObj[2].append(factors[i])


setMatplotlibConf()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 12))

configureAxis(ax1, 'Factors for each star', '', logScale=False)
configureAxis(ax2, 'Factors for each star', '', logScale=False)

ax1.set_xlim(np.nanmin(factors)/1.5, np.nanmax(factors)*1.2)
ax2.set_xlim(np.nanmin(factors)/1.5, np.nanmax(factors)*1.2)

counts1, bins1, patches1 = ax1.hist(sortedFactorsIns[0], alpha=0.8, color="red", label="Spectrograph SDSS")
counts2, bins2, patches2 = ax1.hist(sortedFactorsIns[1], alpha=0.8, color="teal", label="Spectrograph BOSS")
max_bin_height1 = np.nanmax(counts1)
max_bin_height2 = np.nanmax(counts2)
max_bin_height = np.nanmax([max_bin_height1, max_bin_height2])

ax1.vlines(x=np.median(sortedFactorsIns[0]), ymin=0, ymax=max_bin_height, color="red", lw=2, ls="--")
ax1.vlines(x=np.mean(sortedFactorsIns[0]),   ymin=0, ymax=max_bin_height, color="red", lw=2, ls="-")
ax1.vlines(x=np.median(sortedFactorsIns[1]), ymin=0, ymax=max_bin_height, color="teal", lw=2, ls="--")
ax1.vlines(x=np.mean(sortedFactorsIns[1]),   ymin=0, ymax=max_bin_height, color="teal", lw=2, ls="-")

counts2, bins2, patches2 = ax2.hist(sortedFactorsObj[2], alpha=0.8, color="gold", label="Galaxy")
counts3, bins3, patches3 = ax2.hist(sortedFactorsObj[0], alpha=0.8, color="teal", label="Stars")
counts1, bins1, patches1 = ax2.hist(sortedFactorsObj[1], alpha=0.8, color="red", label="QSO")
max_bin_height1 = np.nanmax(counts1)
max_bin_height2 = np.nanmax(counts2)
max_bin_height3 = np.nanmax(counts3)

max_bin_height = np.nanmax([max_bin_height1, max_bin_height2, max_bin_height3])

ax2.vlines(x=np.nanmedian(sortedFactorsObj[1]), ymin=0, ymax=max_bin_height, lw=2, ls="--", color="red")
ax2.vlines(x=np.nanmean(sortedFactorsObj[1]),   ymin=0, ymax=max_bin_height, lw=2, ls="-", color="red")
ax2.vlines(x=np.nanmedian(sortedFactorsObj[2]), ymin=0, ymax=max_bin_height, lw=2, ls="--", color="gold")
ax2.vlines(x=np.nanmean(sortedFactorsObj[2]),   ymin=0, ymax=max_bin_height, lw=2, ls="-", color="gold")
ax2.vlines(x=np.nanmedian(sortedFactorsObj[0]), ymin=0, ymax=max_bin_height, lw=2, ls="--", color="teal")
ax2.vlines(x=np.nanmean(sortedFactorsObj[0]),   ymin=0, ymax=max_bin_height, lw=2, ls="-", color="teal")

ax1.legend(fontsize=18)
ax2.legend(fontsize=18)

plt.savefig(f"./tmpDistributions/distributionOfFactorsInFrame_{index}.png")

