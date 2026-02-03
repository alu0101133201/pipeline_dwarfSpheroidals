# Decode original raw data from HiPERCAM instrument at GTC telescope
#
# This Python script makes the decodification of the original raw data taken
# with HiPERCAM camera. This camera takes data in 5 filters at the same
# time, so it should has 5 HDU extensions. Because of it has 4 channels
# (windows) for each filter, then it would be 5 bands x 4 channels = 20 HDU
# extensions. This is the wanted structure of the images, but this is not
# the structure in which the very raw data images came. The original raw
# data cames in a really strange shape. Each single `.fits' image is a table
# that needs some treatment, and this is what the actual Python script does!
#
# It takes the input image and make multiple outputs from them. Renamming
# according to the internal convention of HiPERCAM pipeline, the filter
# name, the channel notation, and so on.
#
# IMPORTANT NOTE: in order to be able to run this script, it is necessary to
# have the HIPERCAM pipeline installed. Information about how to install and
# run it can be found at:
# http://deneb.astro.warwick.ac.uk/phsaap/hipercam/docs/html/index.html
#
# Original author:
#     Raul Infante-Sainz <infantesainz@gmail.com>
# Contributing author(s):
# Copyright (C) 2019-2021, Raul Infante-Sainz.
#
# This Python script is free software: you can redistribute it and/or modify
# it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This Python script is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details. See <http://www.gnu.org/licenses/>.




# Usage:
# -----

# python3 raw-to-fits.py input_name path_out
#
#
# input_name = 'XXXXXXXXXX-YYYYYYYY-HIPERCAM-Imaging.fits'
# path_out = '/path/where/you/want/to/save/the/images'


# Importing modules
import os
import sys
import hipercam
from hipercam.hcam import Rdata


# Defining arguments passed from command line
input_name = sys.argv[1]
path_out   = sys.argv[2]


# Select only the basename, so, delete path
basename_extension = os.path.basename(input_name)
basename = os.path.splitext(basename_extension)[0]


# Main loop, it will see how many images are stored in the input raw data,
# and then iterate over them to construct the wanted images. Operation are
# done in place using `mccd' object.
for n, mccd in enumerate(Rdata(input_name)):
    ## Comment the following lines to NOT compute the median
    ## subtract median from each window of each CCD
    #for ccd in mccd.values():
    #    for wind in ccd.values():
    #        wind -= wind.median()


    # Rename the basename according with the counter value
    counter_renamed = '{:s}-{:03d}.hcm'.format(basename,n+1)

    # Add output path to the output name
    output_file = path_out + '/' + counter_renamed

    # If output file exists, continue. If not, save it.
    if os.path.isfile(output_file):
        print('FILE ALREADY EXISTS: ', output_file)
        continue

    else:
        # Decode the input file
        mccd.write(output_file)

        # Print on screen the saved imate
        print("SAVED: ", output_file)
