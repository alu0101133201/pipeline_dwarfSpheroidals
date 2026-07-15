import os
import sys
import glob
import pandas as pd
import concurrent.futures
import threading 
from datetime import datetime
from astropy import units as u
from astropy_healpix import HEALPix
from astroquery.gaia import Gaia
import numpy as np
from astropy.io import fits

def getSpectraCSV(ra,dec,sizeOfField,spectraDir,mosaicDir):
    DR3=True
    target_table='Spectroscopy/xp_sampled_mean_spectrum'
    hpx_level = 6
    lon=ra*u.deg
    lat=dec*u.deg
    radius=(sizeOfField/2)*u.deg

    output_file=mosaicDir+"/gaia_spectra.txt"
    output_dir=spectraDir
    if not os.path.isdir(f'{output_dir}'):
        os.system(f'mkdir {output_dir}')
    gaia_dr_flag='DR3'
    url_prefix      = f'http://cdn.gea.esac.esa.int/Gaia/g{gaia_dr_flag.lower()}/{target_table}/'
    md5sum_file_url = url_prefix + '_MD5SUM.txt'
    md5sum_file     = pd.read_csv(md5sum_file_url, header=None, delim_whitespace=True, names=['md5Sum', 'file'])
    md5sum_file.drop(md5sum_file.tail(1).index,inplace=True)
    # Extract HEALPix level-8 from file name ======================================
    healpix_8_min  = [int(file[file.find('_')+1:file.rfind('-')])     for file in md5sum_file['file']]
    healpix_8_max  = [int(file[file.rfind('-')+1:file.rfind('.csv')]) for file in md5sum_file['file']]
    reference_file = pd.DataFrame({'file':md5sum_file['file'], 'healpix8_min':healpix_8_min, 'healpix8_max':healpix_8_max}).reset_index(drop=True)

    # Compute HEALPix levels 6,7, and 9 ===========================================
    reference_file['healpix7_min'] = [inp >> 2 for inp in reference_file['healpix8_min']]
    reference_file['healpix7_max'] = [inp >> 2 for inp in reference_file['healpix8_max']]

    reference_file['healpix6_min'] = [inp >> 2 for inp in reference_file['healpix7_min']]
    reference_file['healpix6_max'] = [inp >> 2 for inp in reference_file['healpix7_max']]

    reference_file['healpix9_min'] = [inp << 2       for inp in reference_file['healpix8_min']]
    reference_file['healpix9_max'] = [(inp << 2) + 3 for inp in reference_file['healpix8_max']]

    # Generate reference file =====================================================
    ncols          = ['file', 'healpix6_min', 'healpix6_max', 'healpix7_min', 'healpix7_max', 'healpix8_min', 'healpix8_max', 'healpix9_min', 'healpix9_max']
    reference_file = reference_file[ncols]
    hp             = HEALPix(nside=2**hpx_level, order='nested')
    hp_cone_search = hp.cone_search_lonlat(lon, lat, radius = radius) 

    f = open(output_file,"w")
    subset=[]
    for index in reference_file.index:
        row = reference_file.iloc[index]
        hp_min, hp_max = row[f'healpix{hpx_level}_min'], row[f'healpix{hpx_level}_max']
        if np.any(np.logical_and(hp_min <= hp_cone_search, hp_cone_search <= hp_max)):
            bulk_file = url_prefix + row['file'] + '\n'
            f.write(bulk_file)
            subset.append(bulk_file)
    f.close()
    os.system(f'wget -i {output_file} -P {output_dir}/ -q  --show-progress --progress=bar:force 2>&1')
    os.system(f'gunzip {output_dir}/*.gz')
    return()
def parse_gaia_array(s):
    # 1. s[1:-1] elimina el '[' del inicio y el ']' del final
    # 2. sep=',' le dice a numpy que corte por comas
    return np.fromstring(s[1:-1], sep=',')

def transformCSVintoFITS(spectraDir, ra, dec, sizeOfField):
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

    # Do the positional query ONCE for the whole field, not per file/chunk
    radius_deg = sizeOfField / 2
    query = f"""
    SELECT source_id, ra, dec
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {ra}, {dec}, {radius_deg})
    )
    """
    job = Gaia.launch_job_async(query)  # async avoids polling limits on big results
    field_sources = job.get_results().to_pandas()
    valid_ids_set = set(field_sources['source_id'])
    coords_lookup = field_sources.set_index('source_id')[['ra', 'dec']]

    sampling_grid = np.linspace(336.0, 1020.0, 343)

    for file in glob.glob(spectraDir + '/*.csv'):
        df = pd.read_csv(file, skiprows=63, sep=',')
        df_filtered = df[df['source_id'].isin(valid_ids_set)].copy()
        df_filtered['flux_arr'] = df_filtered['flux'].apply(parse_gaia_array)
        df_filtered['error_arr'] = df_filtered['flux_error'].apply(parse_gaia_array)

        for _, row in df_filtered.iterrows():
            source_id = row['source_id']
            flux = row['flux_arr']
            flux_error = row['error_arr']
            ra_i, dec_i = coords_lookup.loc[source_id, ['ra', 'dec']]
            solution_id = row['solution_id']

            col_wave = fits.Column(name='wavelength', format='D', unit='nm', array=sampling_grid)
            col_flux = fits.Column(name='flux', format='E', unit='W.m**-2.nm**-1', array=flux)
            col_error = fits.Column(name='flux_error', format='E', unit='W.m**-2.nm**-1', array=flux_error)
            hdu_table = fits.BinTableHDU.from_columns([col_wave, col_flux, col_error])
            head = hdu_table.header
            head['EXTNAME'] = 'BINTABLE'
            head['SOURCEID'] = str(source_id)
            head['SOLUTION'] = str(solution_id)
            head['POS'] = f'({ra_i}, {dec_i})'
            head['REFEPOCH'] = '2016.0'
            head['EPOCHEXT'] = '2.83'
            head['WAVEERRO'] = '0.0'
            head['SPECTRAL'] = '0.0'
            head['WAVEEXTE'] = '684.0'
            head['WAVESTAR'] = '336.0'
            head['WAVEEND'] = '1020.0'
            head['APERTURE'] = '5.8932666E-4'
            head['DATAMODE'] = 'Spectrum 1.01'
            head['PUBLISHE'] = 'ESA/Gaia/DPAC'
            head['TITLE'] = 'Spectrum'
            head['DATE-HDU'] = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')

            prihdu = fits.PrimaryHDU()
            prihdu.header.add_comment("Dummy header; see following table extension")

            hdul = fits.HDUList([prihdu, hdu_table])
            filename = f"XP_SAMPLED-Gaia DR3 {source_id}.fits"
            outfile = os.path.join(spectraDir, filename)
            hdul.writeto(outfile, overwrite=True)

        os.system(f'rm {file}')
        
fieldName   = "UGC00180"
ra          = 4.7840
dec         = 15.7492

sizeOfField = 1


retrieval_type = 'XP_SAMPLED'
data_structure = 'INDIVIDUAL' 
data_release   = 'Gaia DR3'   
datalink_all   = []


if not os.path.exists(f'./gaiaSpectra_{fieldName}'):
	os.mkdir(f'./gaiaSpectra_{fieldName}')

spectraDir  = f'./gaiaSpectra_{fieldName}'
mosaicDir   = f'./'
getSpectraCSV(ra,dec,sizeOfField,spectraDir,mosaicDir)
transformCSVintoFITS(spectraDir,ra,dec,sizeOfField)
