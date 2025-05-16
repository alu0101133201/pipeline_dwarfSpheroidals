import sys

import numpy as np
import matplotlib as mpl
import astropy.units as u
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from filelock import FileLock
from astropy.time import Time
from scipy.optimize import minimize
from matplotlib.ticker import MultipleLocator
from astropy.visualization import astropy_mpl_style
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

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
        "axes.linewidth" : 2,
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

def savePlot(airmass_map):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.tight_layout(pad=8.0)
    configureAxis(ax, 'x (px)', 'y (px)', logScale=False)
    ax.set_title("Airmass map")
    im = ax.imshow(airmass_map)
    ax.invert_yaxis()  
    fig.colorbar(im, ax=ax)  
    plt.savefig(airMassDir + f"/{imageNumber}.png")

def readFitsImage(image, hdu):
    hdul = fits.open(image)
    header = hdul[hdu].header
    data = hdul[hdu].data
    wcs = WCS(header)
    return(data, wcs)

def sec(z_deg):
    z = z_deg
    value = 1 / np.cos(z * (np.pi/180))
    return value

def objective_function(gamma, masked_image, airmass_map, X0):
    corrected_image = masked_image / ((airmass_map / X0) ** gamma)
    return np.nanstd(corrected_image)

def writeGammaToFile(diagnosis_and_badFilesDir, imageNumber, gamma):
    lock = FileLock(diagnosis_and_badFilesDir + "/gammaFactorsFromAirMassMapCorr.txt.lock")
    with lock:
        with open(diagnosis_and_badFilesDir + "/gammaFactorsFromAirMassMapCorr.txt", 'a') as f:
            f.write(f"{imageNumber}: {gamma}\n")

def createMassMapCorrectionFits(airMap, wcs, imageNumber):
    header = wcs.to_header()

    hdulist_out = fits.HDUList()
    hdu_out = fits.ImageHDU(airMap, header=header)
    hdu_empty = fits.PrimaryHDU()
    hdu_empty.header['EXTNAME'] = "NODATA"
    hdulist_out.append(hdu_empty)
    hdulist_out.append(hdu_out)
    hdulist_out.writeto(f'{airMassDir}/airMassMapCorr_smallGrid_entirecamera_{imageNumber}.fits', overwrite=True)


setMatplotlibConf()

maskedImage   = sys.argv[1]
airMassDir    = sys.argv[2]
telescopeLat  = float(sys.argv[3])
telescopeLong = float(sys.argv[4])
telescopeElevation = float(sys.argv[5])
dateObs = sys.argv[6]
diagnosis_and_badFilesDir = sys.argv[7]

imageNumber = maskedImage.split("/")[-1].split("_")[1].split("_")[0]
data, wcs = readFitsImage(maskedImage, 1)

location = EarthLocation(lat=telescopeLat*u.deg, lon=telescopeLong*u.deg, height=telescopeElevation*u.m)
time = Time(dateObs)   

ny, nx = data.shape
y, x = np.mgrid[0:ny, 0:nx]
sky_coords = wcs.pixel_to_world(x, y)

altaz_frame = AltAz(obstime=time, location=location)
altaz_coords = sky_coords.transform_to(altaz_frame)

zenith_angle = 90 * u.deg - altaz_coords.alt
z_map = zenith_angle.to(u.deg).value
airmass_map = sec(z_map)

savePlot(airmass_map)

center_y, center_x = np.array(data.shape) // 2
X0 = airmass_map[center_y, center_x]

gamma_init = 0.5
result = minimize(objective_function, gamma_init, args=(data, airmass_map, X0), bounds=[(0.0, 1.5)])
optimal_gamma = result.x[0]

writeGammaToFile(diagnosis_and_badFilesDir, imageNumber, optimal_gamma)

finalMapForCorrection = (airmass_map / X0) ** optimal_gamma
createMassMapCorrectionFits(finalMapForCorrection, wcs, imageNumber)

