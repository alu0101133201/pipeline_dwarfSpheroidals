import sys 
import astropy

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip
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

def calculateFreedmanBins(data, initialValue = None):
    if (initialValue == None):
        bins = [min(data)]
    else:
        bins = [initialValue]

    binWidht = astropy.stats.freedman_bin_width(data)
    while(bins[-1] <= max(data)):
        bins.append(bins[-1] + binWidht)

    return(bins)

def savePlot(col1Name, table, diagnosisFolder, threshold):
    myBins = calculateFreedmanBins(df[col1_name])

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_title("Relevant parameter to check astrometrisation; Rejecting frames with < " + str(threshold), pad=17)
    configureAxis(ax, col1_name, '', logScale=False)
    ax.set_xlim(0, 65)
    n, _, _ = ax.hist(df[col1_name], bins=myBins, color="teal")
    ax.hist([x for x in df[col1_name] if x < threshold], bins=myBins, color="darkred")

    ax.vlines(x=threshold, ymin=0, ymax=(2*n.max()), color="black", ls="--", lw=2.5)
    ax.set_ylim(0, n.max() + 0.5)
    plt.savefig(diagnosisFolder + "/scamp_contrastParameters_hist.png")
    return()

def identifyBadFrames(df, col1_name, threshold1):
    potentiallyBadAstrometrised = []
    numberOfRows = len(df[col1_name])

    for i in range(numberOfRows):
        if (df[col1_name][i] < threshold1): # 
            potentiallyBadAstrometrised.append(i + 1) # Index of the frames starts in 1 not in 0
    return(potentiallyBadAstrometrised)

def writeBadAstrometrisedFramesToFile(diagnosisFolder, warningFile, badAstrometryFrames):
    with open(diagnosisFolder + "/" + warningFile, 'w') as file:
        for fileName in badAstrometryFrames:
            file.write("entirecamera_" + str(fileName) + '\n')

setMatplotlibConf()

diagnosisFolder = sys.argv[1]
xmlFile         = sys.argv[2]
warningFile     = sys.argv[3]

tree = ET.parse(xmlFile)
root = tree.getroot()
ns = {'vo': 'http://www.ivoa.net/xml/VOTable/v1.1'}
rows = root.findall('.//vo:DATA/vo:TABLEDATA/vo:TR', namespaces=ns)
fields = root.findall('.//vo:FIELD', namespaces=ns)
headers = [field.get('name') for field in fields]

# Column that I found to be key for astrometry. 
# The threshold comes by observing the frames and the documentation. In the docs the say that a threshold of 2 or lower
# is not trustable; I have defined 2.5 to be conservative
col1_name = 'XY_Contrast'
threshold = 5
data = []

for row in rows:
    values = [td.text for td in row.findall('./vo:TD', namespaces=ns)]
    data.append(values)
header_indices = {name: idx for idx, name in enumerate(headers)}

if col1_name in header_indices :
    col1_idx = header_indices[col1_name]

    data = []
    for row in rows:
        values = [td.text for td in row.findall('./vo:TD', namespaces=ns)]

        # There are some rows with different data (not frame information) so I skip them
        # I just store the values with 43 columns which are the proper data rows
        if (len(values) == 43):
            data.append((float(values[col1_idx])))
    
    df = pd.DataFrame(data, columns=[col1_name])
    potentiallyBadAstrometrised = identifyBadFrames(df, col1_name, threshold)
    writeBadAstrometrisedFramesToFile(diagnosisFolder, warningFile, potentiallyBadAstrometrised)
    savePlot(col1_name, df, diagnosisFolder, threshold)
else:
    print(f"The specified column ('{col1_name}') does not exist in the XML.")


