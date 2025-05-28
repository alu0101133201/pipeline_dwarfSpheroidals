###This script will try to measure a theoretical Airmass and get the DATA-obs based on Pickering 2002 interpolation formulae
from astropy.coordinates import EarthLocation,SkyCoord,AltAz
from astropy.time import Time
from astropy import units as u
import sys
from astropy.io import fits

###Variables: file, dateheaderkeyword, ra, dec
img = sys.argv[1]
datK = sys.argv[2]
ra=float(sys.argv[3])
dec=float(sys.argv[4])
lat=float(sys.argv[5])
long=float(sys.argv[6])
elev=float(sys.argv[7])

#Observing site: Izaña
observing_location = EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=elev*u.m)  

#Header
hed = fits.open(img)[1].header
if datK.startswith("DATE"):
    date_obs = Time(hed[datK],format='isot',scaleç='utc')
elif datK.startswith("MJD"):
    date_obs = Time(hed[datK],format='mjd',scale='utc')
else:
    raise Exception("Non supported Time format.") 
aa = AltAz(location=observing_location, obstime=date_obs)

coord = SkyCoord(ra=ra*u.deg, dec = dec*u.deg)
azalt = coord.transform_to(aa)

airmass = azalt.secz.value
print(airmass)