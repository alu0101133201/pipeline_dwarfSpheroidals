###This script will try to measure a theoretical Airmass and get the DATA-obs based on Pickering 2002 interpolation formulae
from astropy.coordinates import EarthLocation,SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import AltAz
from astropy.io import fits

###Variables: file, dateheaderkeyword, ra, dec
img = sys.argv[1]
datK = sys.argv[2]
ra=float(sys.argv[3])
dec=float(sys.argv[4])

#Observing site: Izaña
observing_location = EarthLocation(lat='28d18m04s', lon='16d30m38s', height=2390*u.m)  

#Header
hed = fits.open(img)[1].header
date_obs = Time(hed[datK][:10]+' '+hed[datK][12:])
 
aa = AltAz(location=observing_location, obstime=date_obs)

coord = SkyCoord(ra=ra*u.deg, dec = dec*u.deg)
coord.transform_to(aa)
