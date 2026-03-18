# This script receives 

# off_images_folder = sys.argv[1] # The DATA-or or DATA of the off field images
# off_noise_sky_folder = sys.argv[2] # The noise-sky folder of the off field images
# on_images_folder = sys.argv[3] # The DATA-or or DATA of the on field images
# on_noise_sky_folder = sys.argv[4] # The noise-sky folder TO BE CREATED of the off images
# filt = sys.argv[5] # The filter (for the done.txt)
# exposureTimeRatioOff = float(sys.argv[6]) # exposure times to scale the sky level
# exposureTimeRatioOn  = float(sys.argv[7]) # exposure times to scale the sky level

# And creates the on_noise_sky_folder. For this it takes the DATE-OBS and sky level of the off images,
# interpolates with the DATE-OBS of the on images and creates the noise-sky folder needed


#!/usr/bin/env python3

import os
import sys
from astropy.io import fits
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def collect_off_sky_data(off_images_folder, off_noise_sky_folder):
    """Collect DATE-OBS (as MJD) and sky value pairs from off images."""
    sky_data_pairs = []

    for f_name in os.listdir(off_images_folder):
        if not f_name.endswith(".fits"):
            continue
        
        # extract number between _f and _ccd
        try:
            number = f_name.split("_f")[1].split("_ccd")[0]
        except IndexError:
            print(f"Warning: filename {f_name} does not match expected format")
            continue

        fits_path = os.path.join(off_images_folder, f_name)
        txt_path = os.path.join(off_noise_sky_folder, f"entirecamera_{number}.txt")
        if not os.path.isfile(txt_path):
            print(f"Warning: {txt_path} does not exist")
            continue

        # read first number (sky) from txt
        with open(txt_path, "r") as f:
            line = f.readline().strip()
            parts = line.split()
            if len(parts) < 2:
                print(f"Warning: {txt_path} has unexpected format")
                continue
            sky = float(parts[1])

        # read DATE-OBS from HDU 1
        with fits.open(fits_path) as hdul:
            header = hdul[1].header
            date_obs = header.get("DATE-OBS", None)
            if date_obs is None:
                print(f"Warning: {fits_path} missing DATE-OBS")
                continue

        # convert DATE-OBS to MJD
        time_mjd = Time(date_obs, format='isot', scale='utc').mjd

        sky_data_pairs.append((time_mjd, sky))

    return sky_data_pairs

def interpolate_sky(sky_data_pairs):
  """
  Create an interpolation function for sky vs time.
  Outside the time range, returns the first or last sky value.
  """
  # Extract times and skies
  times = np.array([t for t, s in sky_data_pairs])
  skies = np.array([s for t, s in sky_data_pairs])

  # Sort by time
  sorted_indices = np.argsort(times)
  times_sorted = times[sorted_indices]
  skies_sorted = skies[sorted_indices]

  # Create interpolation function with boundary clipping
  interp_func = interp1d(
      times_sorted,
      skies_sorted,
      kind='linear',
      bounds_error=False,
      fill_value=(skies_sorted[0], skies_sorted[-1])
  )
  return interp_func, times_sorted, skies_sorted


def plot_sky_with_on_images(times, skies, interp_func, on_image_sky):
    """
    Plot original off-sky data, interpolation, and on-image interpolated sky.
    """
    times_fine = np.linspace(times.min(), times.max(), 500)
    skies_interp = interp_func(times_fine)

    plt.figure(figsize=(10, 5))
    plt.scatter(times, skies, color='blue', label='Off-sky data')
    plt.plot(times_fine, skies_interp, color='red', linestyle='--', label='Interpolation')

    # Plot on-images interpolated sky values
    if on_image_sky:
        on_times = np.array([t_mjd for (_, _, t_mjd, _) in on_image_sky])
        on_skies = np.array([sky for (_, _, _, sky) in on_image_sky])
        plt.scatter(on_times, on_skies, color='green', marker='x', s=80, label='On-images interpolated')

    plt.xlabel('MJD')
    plt.ylabel('Sky value')
    plt.title('Sky values vs Time')
    plt.legend()
    plt.grid(True)
    plt.savefig("./interpolateSkyCheck.png")


def assign_sky_to_on_images(on_images_folder, interp_func):
    """
    For each on-image, extract number and DATE-OBS,
    convert to MJD, and evaluate interpolated sky value.
    Returns:
        list of tuples: [(number, date_obs, time_mjd, sky_interp), ...]
    """
    results = []

    for f_name in os.listdir(on_images_folder):
        if not f_name.endswith(".fits"):
            continue
        
        try:
            number = f_name.split("_f")[1].split("_ccd")[0]
        except IndexError:
            print(f"Warning: filename {f_name} does not match expected format")
            continue

        fits_path = os.path.join(on_images_folder, f_name)
        with fits.open(fits_path) as hdul:
            header = hdul[1].header
            date_obs = header.get("DATE-OBS", None)
            if date_obs is None:
                print(f"Warning: {fits_path} missing DATE-OBS")
                continue
            time_mjd = Time(date_obs, format='isot', scale='utc').mjd

        # get interpolated sky
        sky_interp = float(interp_func(time_mjd))
        results.append((number, date_obs, time_mjd, sky_interp))

        print(f"File number: {number}, DATE-OBS: {date_obs}, "
              f"MJD: {time_mjd:.5f}, Interpolated sky: {sky_interp:.3f}")

    return results


def save_interpolated_sky(on_image_sky, dest_folder, filt, exposureTimeRatio):
    """
    Save interpolated sky values into destination folder.
    Each file will be named like 'entirecamera_<number>.txt'
    and contain just the sky value.
    
    Args:
        on_image_sky: list of tuples [(number, date_obs, time_mjd, sky_interp), ...]
        dest_folder: path to folder where files will be saved
    """
    os.makedirs(dest_folder, exist_ok=True)

    for number, date_obs, time_mjd, sky_interp in on_image_sky:
        scaled_sky = sky_interp * exposureTimeRatio

        filename = f"entirecamera_{number}.txt"
        file_path = os.path.join(dest_folder, filename)
        with open(file_path, "w") as f:
            f.write(f"entirecamera_{number}.fits {scaled_sky:.6e} 1 \n")  # scientific notation
        print(f"Saved interpolated sky for number {number} -> {file_path}")

    done_file_path = os.path.join(dest_folder, f"done_{filt}.txt")
    with open(done_file_path, "w") as f:
        f.write("")  # empty file
    print(f"Created completion marker: {done_file_path}")
# -----------------------
# Example usage:


off_images_folder = sys.argv[1]
off_noise_sky_folder = sys.argv[2]
on_images_folder = sys.argv[3]
on_noise_sky_folder = sys.argv[4]
filt = sys.argv[5]
exposureTimeRatioOff = float(sys.argv[6])
exposureTimeRatioOn  = float(sys.argv[7])

# 1. Collect off-sky data
sky_data_pairs = collect_off_sky_data(off_images_folder, off_noise_sky_folder)

# 2. Interpolation
interp_func, times, skies = interpolate_sky(sky_data_pairs)

# 3. Assign interpolated sky values to on-images
on_image_sky = assign_sky_to_on_images(on_images_folder, interp_func)

plot_sky_with_on_images(times, skies, interp_func, on_image_sky)

save_interpolated_sky(on_image_sky, on_noise_sky_folder, filt, float(exposureTimeRatioOn / exposureTimeRatioOff))
