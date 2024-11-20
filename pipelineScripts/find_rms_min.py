# Code to find the minimum std value of the sky

# Original author:
# Giulia Golini <giulia.golini@gmail.com>
# Contributing author(s)
# Copyright (C) 2020, Giulia Golini.
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

# System imports

import os
import pdb
import sys
import warnings
import glob

# 3rd parties


import numpy as np
import matplotlib.pyplot as plt




filter         = sys.argv[1]   
start          = int(sys.argv[2])
end            = int(sys.argv[3])
h              = str(sys.argv[4])   
noiseFilesPath = str(sys.argv[5])
rootDir        = str(sys.argv[6])
iteration      = str(sys.argv[7])
outputFile     = str(sys.argv[8])

rms = []

noiseFiles = glob.glob(noiseFilesPath + "/entirecamera_*.txt")

for file in noiseFiles:
    try:
        with open(file, 'r') as currentFile:
            data = currentFile.read().split(' ')
            rms.append(float(data[2]))

    except Exception as e:
        print(f'Error processing file {file}: {e}')



rms_min = np.nanmin(rms)

file = open(rootDir + '/build/' + outputFile, "w")
file.write(str(rms_min))

