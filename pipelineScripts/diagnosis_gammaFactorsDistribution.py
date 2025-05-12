import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator

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

def readGammaFactors(file):
    gammaFactors = []
    with open(file, 'r') as f:
        for line in f:
            splittedLine = line.split(":")
            gammaFactors.append(float(splittedLine[1]))
    return(np.array(gammaFactors))

def plotGammaValues(values, outputDir):
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_title("Gamma factors from air mass map correction", fontsize=25)
    configureAxis(ax, r'$\gamma$ distribution', '', logScale=False)
    ax.hist(values, color='teal')
    ax.text(0.90, 0.90, r'Correction map = ($\frac{X}{X_0})^{\gamma}$', transform=ax.transAxes,
        fontsize=25, ha='right', va='top', color="blue")
    plt.savefig(outputDir + "/gammaFactorsHist.png")

setMatplotlibConf()

diagnosisAndBadFilesDir = sys.argv[1]

gammaValuesFile = diagnosisAndBadFilesDir + "/gammaFactorsFromAirMassMapCorr.txt"
gammaValues = readGammaFactors(gammaValuesFile)

plotGammaValues(gammaValues, diagnosisAndBadFilesDir)