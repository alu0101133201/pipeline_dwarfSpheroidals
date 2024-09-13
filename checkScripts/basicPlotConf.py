import sys
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
        "xtick.major.size": 16.0,
        "xtick.major.width": 4,
        "xtick.minor.size": 8.0,
        "xtick.minor.width": 4,
        "ytick.major.size": 16,
        "ytick.major.width": 4,
        "ytick.minor.size": 8.0,
        "ytick.minor.width": 4,
        "legend.handlelength": 3.0,
        "axes.linewidth" : 3.5,
        "xtick.major.pad" : 6,
        "ytick.major.pad" : 6,
        "legend.fancybox" : True,
        "mathtext.fontset" : "dejavuserif"
    }
    mpl.rcParams.update(rc_fonts)


def configureAxis(ax, xlabel, ylabel, logScale=True):
    ax.xaxis.set_minor_locator(MultipleLocator(1000000))
    ax.yaxis.set_minor_locator(MultipleLocator(1000000))
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis='x', which='major', labelsize=30, pad=17)
    ax.tick_params(axis='y', which='major', labelsize=30, pad=17)
    ax.set_xlabel(xlabel, fontsize=35, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=35, labelpad=10)
    if(logScale): ax.set_yscale('log')
