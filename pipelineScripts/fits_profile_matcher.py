#!/usr/bin/env python3
"""
Script to find the scaling factor between mean values of two FITS radial profile tables.
Matches profiles within a specified radius range after converting to common units (arcsec).
"""

import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import argparse
import sys

def read_fits_table(filename):
    """
    Read FITS file and extract radius, mean, and surface brightness data.
    
    Parameters:
    -----------
    filename : str
        Path to the FITS file
        
    Returns:
    --------
    dict : Dictionary containing radius, mean, surface_brightness arrays
    """
    try:
        with fits.open(filename) as hdul:
            # Find the table extension
            table_hdu = None
            for hdu in hdul:
                if hasattr(hdu, 'data') and hdu.data is not None:
                    if hasattr(hdu.data, 'dtype') and len(hdu.data.dtype.names or []) > 0:
                        table_hdu = hdu
                        break
            
            if table_hdu is None:
                raise ValueError(f"No table data found in {filename}")
            
            data = table_hdu.data
            
            # Extract columns (handle different possible column names)
            radius_cols = ['RADIUS', 'radius', 'R', 'r']
            mean_cols = ['MEAN', 'mean', 'FLUX', 'flux']
            sb_cols = ['SURFACE_BRIGHTNESS', 'surface_brightness', 'SB', 'sb']
            
            radius = None
            mean = None
            surface_brightness = None
            
            for col in radius_cols:
                if col in data.dtype.names:
                    radius = data[col]
                    break
                    
            for col in mean_cols:
                if col in data.dtype.names:
                    mean = data[col]
                    break
                    
            for col in sb_cols:
                if col in data.dtype.names:
                    surface_brightness = data[col]
                    break
            
            if radius is None or mean is None:
                print(f"Available columns in {filename}: {data.dtype.names}")
                raise ValueError("Could not find required RADIUS and MEAN columns")
            
            return {
                'radius': radius,
                'mean': mean,
                'surface_brightness': surface_brightness,
                'header': table_hdu.header
            }
            
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        sys.exit(1)

def get_pixel_scale(header, default_scale=1.0):
    """
    Extract pixel scale from FITS header.
    
    Parameters:
    -----------
    header : astropy.io.fits.Header
        FITS header
    default_scale : float
        Default pixel scale if not found in header
        
    Returns:
    --------
    float : Pixel scale in arcsec/pixel
    """
    # Common keywords for pixel scale
    scale_keywords = ['CDELT1', 'CD1_1', 'PIXSCALE', 'SCALE', 'PLTSCALE']
    
    for keyword in scale_keywords:
        if keyword in header:
            scale = abs(float(header[keyword]))
            # Convert from degrees to arcseconds if necessary
            if scale < 0.01:  # Likely in degrees
                scale *= 3600
            return scale
    
    print(f"Warning: Pixel scale not found in header. Using default: {default_scale} arcsec/pixel")
    return default_scale

def calculate_surface_brightness(mean, pixel_scale):
    """
    Calculate surface brightness from mean values.
    SB = -2.5*log10(mean) + 22.5 + 5*log10(pixel_scale)
    
    Parameters:
    -----------
    mean : array
        Mean flux values
    pixel_scale : float
        Pixel scale in arcsec/pixel
        
    Returns:
    --------
    array : Surface brightness in mag/arcsec²
    """
    # Handle zero or negative values
    mean_clean = np.where(mean <= 0, np.nan, mean)
    sb = -2.5 * np.log10(mean_clean) + 22.5 + 5 * np.log10(pixel_scale)
    return sb
    """
    Convert radius from pixels to arcseconds.
    
    Parameters:
    -----------
    radius_pix : array
        Radius in pixels
    pixel_scale : float
        Pixel scale in arcsec/pixel
        
    Returns:
    --------
    array : Radius in arcseconds
    """
    return radius_pix * pixel_scale

def interpolate_profile(radius, mean, target_radius):
    """
    Interpolate profile to target radius points.
    
    Parameters:
    -----------
    radius : array
        Original radius points
    mean : array
        Original mean values
    target_radius : array
        Target radius points for interpolation
        
    Returns:
    --------
    array : Interpolated mean values
    """
    # Remove any NaN or infinite values
    mask = np.isfinite(radius) & np.isfinite(mean)
    radius_clean = radius[mask]
    mean_clean = mean[mask]
    
    if len(radius_clean) < 2:
        raise ValueError("Not enough valid data points for interpolation")
    
    # Sort by radius
    sort_idx = np.argsort(radius_clean)
    radius_clean = radius_clean[sort_idx]
    mean_clean = mean_clean[sort_idx]
    
    # Create interpolation function
    interp_func = interp1d(radius_clean, mean_clean, kind='linear', 
                          bounds_error=False, fill_value=np.nan)
    
    return interp_func(target_radius)

def find_scaling_factor_from_sb(sb1, sb2):
    """
    Find the surface brightness offset between two profiles.
    
    Parameters:
    -----------
    sb1 : array
        Reference surface brightness values (mag/arcsec²)
    sb2 : array
        Target surface brightness values (mag/arcsec²)
        
    Returns:
    --------
    tuple : (sb_offset, mean_scale_factor)
        sb_offset: Surface brightness offset (sb2 - sb1) in mag/arcsec²
        mean_scale_factor: Equivalent scaling factor for mean values
    """
    # Remove NaN values
    mask = np.isfinite(sb1) & np.isfinite(sb2)
    sb1_clean = sb1[mask]
    sb2_clean = sb2[mask]
    
    if len(sb1_clean) < 2:
        raise ValueError("Not enough overlapping data points for fitting")
    
    # Calculate the mean offset in surface brightness
    sb_offset = np.mean(sb2_clean - sb1_clean)
    
    # Convert surface brightness offset to equivalent mean scaling factor
    # If SB2 = SB1 + offset, then:
    # -2.5*log10(mean2) = -2.5*log10(mean1) + offset
    # log10(mean2) = log10(mean1) - offset/2.5
    # mean2 = mean1 * 10^(-offset/2.5)
    # So scaling factor = 10^(-offset/2.5)
    mean_scale_factor = 10**(-sb_offset / 2.5)
    
    return sb_offset, mean_scale_factor

def find_scaling_factor(mean1, mean2, scale1, scale2):
    """
    Find the scaling factor to match surface brightness profiles (fallback method).
    This is used when surface brightness columns are not available.
    """
    # Remove NaN values and ensure positive values
    mask = np.isfinite(mean1) & np.isfinite(mean2) & (mean1 > 0) & (mean2 > 0)
    mean1_clean = mean1[mask]
    mean2_clean = mean2[mask]
    
    if len(mean1_clean) < 2:
        raise ValueError("Not enough overlapping data points for fitting")
    
    # Calculate surface brightness for both profiles
    sb1 = calculate_surface_brightness(mean1_clean, scale1)
    sb2 = calculate_surface_brightness(mean2_clean, scale2)
    
    # Use the surface brightness method
    sb_offset, mean_scale_factor = find_scaling_factor_from_sb(sb1, sb2)
    
    return mean_scale_factor

def convert_to_arcsec(radius,scale):
    return radius*scale

def plot_comparison(radius1, mean1, sb1, radius2, mean2, sb2, scale_factor, r_min, r_max, save_plot=False):
    """
    Plot the original and scaled surface brightness profiles for comparison.
    """
    plt.figure(figsize=(12, 8))
    
    # Calculate scaled surface brightness: new_sb = old_sb - 2.5*log10(scale_factor)
    # If scale_factor > 1, we're making it brighter, so SB should decrease (become more negative)
    sb_correction = -2.5 * np.log10(scale_factor)
    sb1_scaled = sb1 + sb_correction
    
    # Plot original profiles
    plt.subplot(2, 1, 1)
    if sb1 is not None and sb2 is not None:
        plt.plot(radius1, sb1, 'b-', label='Profile 1 (reference)', linewidth=2)
        plt.plot(radius2, sb2, 'r-', label='Profile 2 (target)', linewidth=2)
        plt.ylabel('Surface Brightness (mag/arcsec²)')
        plt.gca().invert_yaxis()  # Invert y-axis for magnitudes (brighter = lower values)
    else:
        # Fallback to mean if surface brightness not available
        plt.plot(radius1, mean1, 'b-', label='Profile 1 (reference)', linewidth=2)
        plt.plot(radius2, mean2, 'r-', label='Profile 2 (target)', linewidth=2)
        plt.ylabel('Mean')
        plt.yscale('log')
    
    plt.axvline(r_min, color='gray', linestyle='--', alpha=0.7, label=f'R range: {r_min}-{r_max}"')
    plt.axvline(r_max, color='gray', linestyle='--', alpha=0.7)
    plt.xlabel('Radius (arcsec)')
    plt.title('Original Profiles')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot scaled comparison
    plt.subplot(2, 1, 2)
    if sb1 is not None and sb2 is not None:
        plt.plot(radius1, sb1_scaled, 'b-', label=f'Profile 1 scaled (SF={scale_factor:.4f})', linewidth=2)
        plt.plot(radius2, sb2, 'r-', label='Profile 2 (target)', linewidth=2)
        plt.ylabel('Surface Brightness (mag/arcsec²)')
        plt.gca().invert_yaxis()  # Invert y-axis for magnitudes
    else:
        # Fallback to mean if surface brightness not available
        plt.plot(radius1, mean1 * scale_factor, 'b-', label=f'Profile 1 × {scale_factor:.4f}', linewidth=2)
        plt.plot(radius2, mean2, 'r-', label='Profile 2 (target)', linewidth=2)
        plt.ylabel('Mean')
        plt.yscale('log')
    
    plt.axvline(r_min, color='gray', linestyle='--', alpha=0.7, label=f'R range: {r_min}-{r_max}"')
    plt.axvline(r_max, color='gray', linestyle='--', alpha=0.7)
    plt.xlabel('Radius (arcsec)')
    plt.title('Scaled Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig('profile_comparison.png', dpi=300, bbox_inches='tight')
        print("Plot saved as 'profile_comparison.png'")
    
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Find scaling factor between FITS radial profiles')
    parser.add_argument('file1', help='First FITS file (reference)')
    parser.add_argument('file2', help='Second FITS file (target)')
    parser.add_argument('r_min', type=float, help='Minimum radius for matching (arcsec)')
    parser.add_argument('r_max', type=float, help='Maximum radius for matching (arcsec)')
    parser.add_argument('--scale1', type=float, default=None, 
                       help='Pixel scale for file1 (arcsec/pixel). If not provided, will try to read from header')
    parser.add_argument('--scale2', type=float, default=None,
                       help='Pixel scale for file2 (arcsec/pixel). If not provided, will try to read from header')
    parser.add_argument('--plot', action='store_true', help='Show comparison plot')
    parser.add_argument('--save-plot', action='store_true', help='Save comparison plot')
    
    args = parser.parse_args()
    
    print(f"Reading FITS files...")
    print(f"File 1: {args.file1}")
    print(f"File 2: {args.file2}")
    print(f"Radius range: {args.r_min} - {args.r_max} arcsec")
    print("-" * 50)
    
    # Read FITS files
    data1 = read_fits_table(args.file1)
    data2 = read_fits_table(args.file2)
    
    # Get pixel scales
    if args.scale1 is not None:
        scale1 = args.scale1
    else:
        scale1 = get_pixel_scale(data1['header'])
    
    if args.scale2 is not None:
        scale2 = args.scale2
    else:
        scale2 = get_pixel_scale(data2['header'])
    
    print(f"Pixel scale 1: {scale1} arcsec/pixel")
    print(f"Pixel scale 2: {scale2} arcsec/pixel")
    
    # Convert to arcseconds
    radius1_arcsec = convert_to_arcsec(data1['radius'], scale1)
    radius2_arcsec = convert_to_arcsec(data2['radius'], scale2)
    
    # Create common radius grid for comparison
    r_min_common = max(np.min(radius1_arcsec), np.min(radius2_arcsec), args.r_min)
    r_max_common = min(np.max(radius1_arcsec), np.max(radius2_arcsec), args.r_max)
    
    if r_min_common >= r_max_common:
        print("Error: No overlapping radius range between the two profiles!")
        print(f"Profile 1 range: {np.min(radius1_arcsec):.2f} - {np.max(radius1_arcsec):.2f} arcsec")
        print(f"Profile 2 range: {np.min(radius2_arcsec):.2f} - {np.max(radius2_arcsec):.2f} arcsec")
        print(f"Requested range: {args.r_min} - {args.r_max} arcsec")
        sys.exit(1)
    
    print(f"Effective matching range: {r_min_common:.2f} - {r_max_common:.2f} arcsec")
    
    # Find actual data points within the range for both profiles
    has_sb1 = data1['surface_brightness'] is not None
    has_sb2 = data2['surface_brightness'] is not None
    
    print(f"Surface brightness column in file 1: {'Yes' if has_sb1 else 'No'}")
    print(f"Surface brightness column in file 2: {'Yes' if has_sb2 else 'No'}")
    
    # Find actual data points within the range for both profiles
    mask1 = (radius1_arcsec >= r_min_common) & (radius1_arcsec <= r_max_common)
    mask2 = (radius2_arcsec >= r_min_common) & (radius2_arcsec <= r_max_common)
    
    radius1_range = radius1_arcsec[mask1]
    mean1_range = data1['mean'][mask1]
    radius2_range = radius2_arcsec[mask2]
    mean2_range = data2['mean'][mask2]
    
    if has_sb1:
        sb1_range = data1['surface_brightness'][mask1]
    if has_sb2:
        sb2_range = data2['surface_brightness'][mask2]
    
    print(f"Profile 1 points in range: {len(radius1_range)}")
    print(f"Profile 2 points in range: {len(radius2_range)}")
    
    # Use the profile with fewer points as the reference for common radii
    if len(radius1_range) <= len(radius2_range):
        common_radius = radius1_range
        mean1_interp = mean1_range
        mean2_interp = interpolate_profile(radius2_arcsec, data2['mean'], common_radius)
        
        if has_sb1:
            sb1_interp = sb1_range
        else:
            sb1_interp = calculate_surface_brightness(mean1_interp, scale1)
            
        if has_sb2:
            sb2_interp = interpolate_profile(radius2_arcsec, data2['surface_brightness'], common_radius)
        else:
            sb2_interp = calculate_surface_brightness(mean2_interp, scale2)
    else:
        common_radius = radius2_range
        mean1_interp = interpolate_profile(radius1_arcsec, data1['mean'], common_radius)
        mean2_interp = mean2_range
        
        if has_sb1:
            sb1_interp = interpolate_profile(radius1_arcsec, data1['surface_brightness'], common_radius)
        else:
            sb1_interp = calculate_surface_brightness(mean1_interp, scale1)
            
        if has_sb2:
            sb2_interp = sb2_range  
        else:
            sb2_interp = calculate_surface_brightness(mean2_interp, scale2)
    
    try:
        # Use surface brightness directly if available, otherwise calculate from mean
        if has_sb1 and has_sb2:
            print("Using surface brightness columns directly...")
            sb_offset, scale_factor = find_scaling_factor_from_sb(sb1_interp, sb2_interp)
        else:
            print("Calculating surface brightness from mean values...")
            scale_factor = find_scaling_factor(mean1_interp, mean2_interp, scale1, scale2)
            sb_offset = sb2_interp.mean() - sb1_interp.mean()  # Approximate offset
        
        # Calculate goodness of fit in surface brightness space
        mask = np.isfinite(sb1_interp) & np.isfinite(sb2_interp)
        if np.sum(mask) > 0:
            residuals = sb2_interp[mask] - sb1_interp[mask]
            rms_residual = np.sqrt(np.mean(residuals**2))
            relative_rms = rms_residual / np.std(sb2_interp[mask]) * 100 if np.std(sb2_interp[mask]) > 0 else np.nan
        else:
            rms_residual = np.nan
            relative_rms = np.nan
        
        print("-" * 50)
        print("RESULTS:")
        
        if has_sb1 and has_sb2:
            print(f"Surface brightness offset (Profile2 - Profile1): {sb_offset:.4f} mag/arcsec²")
            if sb_offset > 0:
                print(f"  Profile 2 is {sb_offset:.4f} mag fainter than Profile 1")
            else:
                print(f"  Profile 2 is {-sb_offset:.4f} mag brighter than Profile 1")
            print(f"Equivalent mean scaling factor: {scale_factor:.6f}")
        else:
            print(f"Mean scaling factor: {scale_factor:.6f}")
            print(f"Surface brightness offset: {sb_offset:.4f} mag/arcsec²")
        
        print(f"To match profile 1 to profile 2, multiply profile 1 mean by: {scale_factor:.6f}")
        
        # Show the surface brightness comparison
        if np.sum(mask) > 0:
            sb1_mean = np.nanmean(sb1_interp[mask])
            sb2_mean = np.nanmean(sb2_interp[mask])
            
            print(f"\nSurface brightness comparison:")
            print(f"  Profile 1 mean SB: {sb1_mean:.3f} ± {np.nanstd(sb1_interp[mask]):.3f} mag/arcsec²")
            print(f"  Profile 2 mean SB: {sb2_mean:.3f} ± {np.nanstd(sb2_interp[mask]):.3f} mag/arcsec²")
            print(f"  Difference (P2-P1): {sb2_mean - sb1_mean:.3f} mag/arcsec²")
            
            print(f"\nFit quality:")
            print(f"  RMS residual in SB: {rms_residual:.4f} mag/arcsec²")
            print(f"  Relative RMS: {relative_rms:.2f}%")
            print(f"  Number of matched points: {np.sum(mask)}")
        
        # Show plot if requested
        if args.plot or args.save_plot:
            # Use existing or calculated surface brightness for plotting
            if has_sb1:
                sb1_full = data1['surface_brightness']
            else:
                sb1_full = calculate_surface_brightness(data1['mean'], scale1)
                
            if has_sb2:
                sb2_full = data2['surface_brightness']  
            else:
                sb2_full = calculate_surface_brightness(data2['mean'], scale2)
                
            plot_comparison(radius1_arcsec, data1['mean'], sb1_full,
                          radius2_arcsec, data2['mean'], sb2_full,
                          scale_factor, args.r_min, args.r_max, args.save_plot)
        
    except Exception as e:
        print(f"Error finding scaling factor: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()