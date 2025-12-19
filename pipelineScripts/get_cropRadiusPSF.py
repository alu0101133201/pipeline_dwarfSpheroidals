import sys
import numpy as np
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
starMag=float(sys.argv[1])
psfProf=sys.argv[2]
limitMag=float(sys.argv[3])
pixelScale=float(sys.argv[4])

with fits.open(psfProf) as hdul:
    I_psf=hdul[1].data['MEAN']
    R_psf=hdul[1].data['RADIUS']
    I_tot=2*np.pi*cumulative_trapezoid(I_psf*R_psf,R_psf,initial=0)
    alfa=(10**(-0.4*starMag))/I_tot[-1]
    mag_psf=-2.5*np.log10(I_psf*alfa)+5*np.log10(pixelScale) 
    indexes=np.where(mag_psf<=limitMag)
    if indexes[0][-1]==hdul[1].data['RADIUS'][-1]:
        print(0)
    else:
        id=indexes[0][-1]
        print(int(hdul[1].data['RADIUS'][id]))
