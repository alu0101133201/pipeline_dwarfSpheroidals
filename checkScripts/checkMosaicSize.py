"""
With this python script we aim to compute the minimum size required for the final mosaic
Taking into account that our observation strategy is based on a dithering pattern, the 
mosaic must be larger than the individual frames. How much it deppends on the dithering.

Usage:
    python3 checkMosaicSize.py --pixel-scale 0.598 --center ra dec files*.fits

Optional parameters:
    --margin add extra pixels around the bounding box (default=2)
    --output specify output text filename (default: mosaic_size_summary.txt)
"""
import argparse
import math
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord,SkyOffsetFrame
import astropy.units as u
from concurrent.futures import ProcessPoolExecutor, as_completed
import os,sys

def corners_radec_from_fits(fname):
    """Return (filename, SkyCoord) or (filename, None, error_message)."""
    try:
        try:
            hdul = fits.open(fname, memmap=True)
            header, data = None, None
            for hdu in hdul:
                if hasattr(hdu, "data") and hdu.data is not None:
                    header = hdu.header
                    data = hdu.data
                    break
        except Exception:
            # fallback: open without memmap if scaling keywords prevent memmap
            hdul = fits.open(fname, memmap=False)
            header, data = None, None
            for hdu in hdul:
                if hasattr(hdu, "data") and hdu.data is not None:
                    header = hdu.header
                    data = hdu.data
                    break
        
        if header is None or data is None:
            hdul.close()
            raise ValueError("No image HDU with data")

        w = WCS(header)
        if data.ndim < 2:
            hdul.close()
            raise ValueError(f"Unexpected data shape: {data.shape}")

        ny, nx = data.shape[-2], data.shape[-1]
        pix = np.array(
            [[0, 0],
             [nx - 1, 0],
             [0, ny - 1],
             [nx - 1, ny - 1]],
            dtype=float
        )

        try:
            world = w.pixel_to_world(pix[:, 0], pix[:, 1])
        except:
            # fallback: use low-level transformation
            ra_dec = w.wcs_pix2world(pix, 0)
            world = SkyCoord(ra_dec[:, 0] * u.deg, ra_dec[:, 1] * u.deg, frame="icrs")
        if not isinstance(world, SkyCoord):
            raise ValueError(f"Invalid WCS conversion in {fname}")
        hdul.close()
        return (fname, world)

    except Exception as e:
        return (fname, None, str(e))


def main(file_patterns,pixel_scale_arcsec,center_ra_deg,center_dec_deg,margin=2,outfile="mosaic_size_summary.txt",nproc=None):
    files = []
    for p in file_patterns:
        files.extend(sorted(glob.glob(p)))
    if not files:
        raise SystemExit("No input files found")
    
    all_coords = []
    failed = []
    with ProcessPoolExecutor(max_workers=nproc) as executor:
        futures = [executor.submit(corners_radec_from_fits, fn) for fn in files]
        for fut in as_completed(futures):
            result = fut.result()
            if len(result) == 2:
                fn, sky = result
                all_coords.extend(sky)
            elif len(result) == 3:
                fn, _, msg = result
                failed.append((fn, f"ERROR: {msg}"))
    if failed:
        print(f"⚠️  {len(failed)} files failed to read corners:")
        for fn, msg in failed:
            print(f"   {fn}: {msg}")

    if not all_coords:
        raise SystemExit("No valid coordinates extracted from FITS files.")

    
    #Define center from user-supplied
    center=SkyCoord(center_ra_deg*u.deg,center_dec_deg*u.deg,frame='icrs')

    #Convert corners to offsets in the tangent/offset frame centered at center
    offset_frame=SkyOffsetFrame(origin=center)
    xs_arcsec = []
    ys_arcsec = []
    for c in all_coords:
        off = c.transform_to(offset_frame)
        xs_arcsec.append(off.lon.to(u.arcsec).value) # east offset
        ys_arcsec.append(off.lat.to(u.arcsec).value) # north offset
    
    xs = np.array(xs_arcsec)
    ys = np.array(ys_arcsec)

    xmin,xmax = xs.min(), xs.max()
    ymin,ymax = ys.min(), ys.max()
    width_arcsec = xmax-xmin
    height_arcsec = ymax-ymin

    #Convert to pixels and add margin
    nx = math.ceil(width_arcsec / pixel_scale_arcsec) + margin*2
    ny = math.ceil(height_arcsec / pixel_scale_arcsec) + margin*2

    lines=[]
    lines.append(f"Using supplied projection center (CRVAL): RA={center.ra.deg:.8f} deg Dec={center.dec.deg:.8f}\"")
    lines.append(f"Span in arcsec: X [{xmin:.1f}, {xmax:.1f}] -> width {width_arcsec:.1f}\"")
    lines.append(f"              : Y [{ymin:.1f}, {ymax:.1f}] -> height {height_arcsec:.1f}\"")
    lines.append(f"Pixel scale: {pixel_scale_arcsec:.6f} arcsec/pix")
    lines.append(f"Required mosaic size (including margin {margin}px): {nx} x {ny} pixels (W x H)")
    max_span_deg = max(width_arcsec, height_arcsec) / 3600.0
    if max_span_deg > 2.0:
        lines.append("")
        lines.append("Warning: mosaic spans >~2 degrees. Tangent-plane approximations may be poor;")
        lines.append("consider verifying with a projection-aware reprojection tool (reproject or montage).")

    output_text="\n".join(lines)
    print(output_text)

    with open(outfile,"w") as f:
        f.write(output_text + "\n")
    print(f"\nResults saved to: {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute minimal mosaic size using tangent-plane offsets (handles rotation).")
    parser.add_argument("files", nargs="+", help="input FITS files or glob patterns")
    parser.add_argument("--center", nargs=2, type=float, metavar=("RA_DEG", "DEC_DEG"), required=True,
                        help="projection center (CRVAL) in degrees")
    parser.add_argument("--pixel-scale", type=float, required=True, help="desired output pixel scale in arcsec/pix")
    parser.add_argument("--margin", type=int, default=2, help="extra pixel margin (default 2)")
    parser.add_argument("--output", default="mosaic_size_summary.txt", help="output text file name (default mosaic_size_summary.txt)")
    parser.add_argument("--nproc", type=int, default=None, help="number of parallel processes (default = CPU count)")
    args = parser.parse_args()

    # Expand globs
    filelist = []
    for pattern in args.files:
        filelist.extend(sorted(glob.glob(pattern)))
    if not filelist:
        raise SystemExit("No FITS files found — check your input patterns.")

    main(filelist, args.pixel_scale, args.center[0], args.center[1], args.margin, args.output)    