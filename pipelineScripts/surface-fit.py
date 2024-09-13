#!/usr/bin/env python3

# Fit 2D surface model to a given image
#
# This is a Python script to obtain a 2D surface fit of a given image.  It
# considers an input image and uses it as the basis for obtaining a surface
# fitting.  By default, the script will save an image of the same
# dimensions than the input image as the fited model. If the user provides
# a new grid dimension, this will be used to evaluate the fited model. The
# models are analytical and comes from the Astropy 2D functions
# (Polynomial2D, Hermite2D, Legendre2D, Chebyshev2D):
# https://docs.astropy.org/en/stable/modeling/fitting.html
#
# Original author:
#     Raul Infante-Sainz <infantesainz@gmail.com>
# Contributing author(s):
# Copyright (C) 2019-2021, Raul Infante-Sainz.
#
# This Python script is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This Python script is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details. See <http://www.gnu.org/licenses/>.




# Usage
# -----
#
# python3 polynomial-fit.py -i image.fits -h N -d D \
#                           -g x y -o model.fits

# python3 polynomial-fit.py --input=image.fits --hdu=N --degree=D \
#                           --grid x y --output=model.fits





# Import modules
import sys
import argparse
import warnings
import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting





# Available Astropy 2D models
surfacemodels = ['Polynomial2D', 'Hermite2D', 'Legendre2D', 'Chebyshev2D']






# Read input arguments from the command line
# ------------------------------------------
#

parser = argparse.ArgumentParser(description='Arguments from the command line')
# By default, argparse uses -h, --help as arguments to provide help. To be
# able to supress this behaviour and use -h for the extension number (as
# all Gnuastro programs do), it is necessary to disable this option. Here,
# the options to obtain help are: -H, --help
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-H', '--help', action='help', default=argparse.SUPPRESS,
                    help='Show this help message and exit.')

parser.add_argument('-i', '--image', type=str,
                    help='Input image name.')
parser.add_argument('-f', '--fileToWrite', type=str, 
                    help='File to write the coefficients of the plane')
parser.add_argument('-o', '--output', type=str, default="default",
                    help='Output image name.')
parser.add_argument('-h', '--hdu', type=int, default=1,
                    help='Input HDU number.')
parser.add_argument('-d', '--degree', type=int, default=1,
                    help='Degree of the polynomial model.')
parser.add_argument('-g', '--grid', type=str,  default="0,0",
                    help='x,y dimensions for the output.')
parser.add_argument('-m', '--model', type=str,  default=surfacemodels[0],
                    help='Available models: ' + str(surfacemodels).replace("'","") + \
                        '. Default: ' + surfacemodels[0])

# Define variables from command line arguments
args = parser.parse_args()
nhdu = args.hdu
image = args.image
model = args.model
degree = args.degree
output = args.output
grid = args.grid.split(',')
fileToWriteCoeff = args.fileToWrite
# Obtain xdim and ydim from --grid parameter
xdim = int(grid[0])
ydim = int(grid[1])





# Set the model to be fitted
# --------------------------
#
# Only some models already availables from Astropy are considered here.
# The user has to specify one of them. Default one will be 'Polynomial2D'
if model not in surfacemodels:
    raise ValueError('WARNING! -m, --model has to be one of:', surfacemodels)

if model == surfacemodels[0]:
    p_init = models.Polynomial2D(degree=degree)
if model == surfacemodels[1]:
    p_init = models.Hermite2D(x_degree=degree, y_degree=degree)
if model == surfacemodels[2]:
    p_init = models.Legendre2D(x_degree=degree, y_degree=degree)
if model == surfacemodels[3]:
    p_init = models.Chebyshev2D(x_degree=degree, y_degree=degree)
fit_p = fitting.LevMarLSQFitter()





# Fit the polynomial function
# ---------------------------
#
# Open the input image
with fits.open(image) as hdul:

    # Define HDU list and read input (data and header)
    hdulist_out = fits.HDUList()
    data_orig = hdul[nhdu].data
    header_orig = hdul[nhdu].header

    # Original x,y grid
    y_max = data_orig.shape[0]
    x_max = data_orig.shape[1]
    y_orig, x_orig = np.mgrid[:y_max, :x_max]

    # Mask of nan values (True --> not nan)
    nanmask = ~(np.isnan(data_orig))

    # Do the fitting ignoring the warnings
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')

        print("Fitting ", model, " model to: ", image)
        polyfit = fit_p(p_init,
                        x_orig[nanmask],
                        y_orig[nanmask],
                        data_orig[nanmask])


    # Obtain the dimension of the output array image, if they are zero,
    # then use the original dimensions (same shape as the input image).
    if ydim == 0: ydim = data_orig.shape[0]
    if xdim == 0: xdim = data_orig.shape[1]

    # Generate the output (new) grid. Use linspace with the number of
    # elements equal to the desired dimension. In practice, this generate a
    # vector of (over)sampled values.
    x_sampled = np.linspace(0, x_max, xdim)
    y_sampled = np.linspace(0, y_max, ydim)

    # Mix the two dimension (x,y) to generate the grid over which the
    # function will be evaluated.
    x_out, y_out = np.meshgrid(x_sampled, y_sampled)

    # Generate the output HDU with the polynomial model by evaluating the
    # function obtained over the wanted pixel grid values.
    hdu_out = fits.ImageHDU(polyfit(x_out, y_out).astype(np.float32))

    # Metadata: original header and general keywords.
    hdu_out.header = header_orig
    hdu_out.header['SFITDEG'] = degree
    hdu_out.header['EXTNAME'] = "SURFIT"
    hdu_out.header['COMMENT'] = model + " fit of degree " + str(degree)

    # Inject coefficients of the fit.
    for pname, pvalue in zip(polyfit.param_names, polyfit.parameters):
        hdu_out.header[pname] = pvalue

    # Write coefficients in a file
    with open(fileToWriteCoeff, 'w') as file:
        for pname, pvalue in zip(polyfit.param_names, polyfit.parameters):
            file.write(str(pname) + " " + str(pvalue) + " ")
            file.write("\n")
            
    # Construct the output HDU list:
    #  0 NODATA
    #  1 SURFIT
    hdu_empty = fits.PrimaryHDU()
    hdu_empty.header['EXTNAME'] = "NODATA"
    hdulist_out.append(hdu_empty)
    hdulist_out.append(hdu_out)

    # Save the HDU list as FITS image.  If the user don't specify an
    # output, then use the same as the input but adding the model before
    # the '.fits'.
    if output == 'default':
        output_name = image.replace('.fits','_'+model+'.fits')
    else:
        output_name = output
    hdulist_out.writeto(output_name, overwrite=True)
    print("Output written to: ", output_name)


