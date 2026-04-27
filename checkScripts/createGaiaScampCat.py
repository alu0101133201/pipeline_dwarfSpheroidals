import numpy as np
from astropy.io import fits
from astropy.table import Table
import sys
input = sys.argv[1]
out = sys.argv[2]
gaia_data=Table.read(input, format='csv')
scamp_table=Table()
scamp_table['RA']=gaia_data['ra']
scamp_table['DEC']=gaia_data['dec']
scamp_table['ERRA_WORLD']=gaia_data['ra_error'] / 3600000.0
scamp_table['ERRB_WORLD']=gaia_data['dec_error'] / 3600000.0
scamp_table['ERRTHETA_WORLD']=np.zeros(len(gaia_data))
scamp_table['PMRA']=np.nan_to_num(gaia_data['pmra'], nan=0.0)
scamp_table['PMDEC']=np.nan_to_num(gaia_data['pmdec'], nan=0.0)
scamp_table['PMRA_ERR']=np.nan_to_num(gaia_data['pmra_error'], nan=0.0)
scamp_table['PMDEC_ERR']=np.nan_to_num(gaia_data['pmdec_error'], nan=0.0)
scamp_table['MAG'] = gaia_data['phot_g_mean_mag']
flux_str = gaia_data['phot_g_mean_flux'].astype(str)
flux_err_str = gaia_data['phot_g_mean_flux_error'].astype(str)
for bad_value in ['NULL', 'null', 'None', '']:
    flux_str = np.where(flux_str == bad_value, 'NaN', flux_str)
    flux_err_str = np.where(flux_err_str == bad_value, 'NaN', flux_err_str)
flux = flux_str.astype(float)
flux_err = flux_err_str.astype(float)
scamp_table['MAG_ERR']=2.5*np.abs(flux_err/flux/np.log(10))
primary_hdu = fits.PrimaryHDU()
dummy_header = fits.Header()
dummy_header['NAXIS']=0
header_str=dummy_header.tostring(endcard=True)
col = fits.Column(name='Field Header Card',format='2880A',array=[header_str])
imhead_hdu = fits.BinTableHDU.from_columns([col], name='LDAC_IMHEAD')

objects_hdu = fits.table_to_hdu(scamp_table)
objects_hdu.name = 'LDAC_OBJECTS'
hdul = fits.HDUList([primary_hdu, imhead_hdu, objects_hdu])
hdul.writeto(out, overwrite=True)