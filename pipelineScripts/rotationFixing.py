from astropy.io import fits
import sys

def fix_rotation_to_0_inplace(fits_file):
    with fits.open(fits_file,mode='update') as hdul:
        hdr=hdul[1].header
        hdr['PC1_1'] = 1.0
        hdr['PC1_2'] = 0.0
        hdr['PC2_1'] = 0.0
        hdr['PC2_2'] = 1.0
        if 'CDELT1' in hdr:
            hdr['CDELT1'] = abs(hdr['CDELT1'])
        if 'CDELT2' in hdr:
            hdr['CDELT2'] = abs(hdr['CDELT2'])

        # Save changes in place
        hdul.flush()

image=sys.argv[1]
fix_rotation_to_0_inplace(image)