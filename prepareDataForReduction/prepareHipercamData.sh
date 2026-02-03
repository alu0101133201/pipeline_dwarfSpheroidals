#!/bin/bash
# This script decompress HiPERCAM data and prepares it for reduction with the pipeline. It needs: folder with raw HiPERCAM data, output folder and, if needed, a block for warping
# Usage: prepareHipercamData.sh <inputHipercamFolder> <outputFolder> [<warpBlockFile>]
inputHipercamFolder=$1
outputFolder=$2
warpBlock=$3 #Integer number
if [ -z "$warpBlock" ]; then
  warpBlock=0
fi
#Create output folder if it does not exist
if ! [ -d "$outputFolder" ]; then
  mkdir -p "$outputFolder"
fi

#Decoding raw data
for file in $inputHipercamFolder/*.fits; do
    python3 ./src/decoding-rawdata.py $file $outputFolder
done

#Re-distributing into different folders per filter
for file in $outputFolder/*.hcm; do
    a=1
    name=$( basename $file )
    name_ok="${name%.hcm}.hcm"
    for filter in u g r i z; do
        if ! [ -d "$outputFolder/$filter" ]; then
            mkdir -p "$outputFolder/$filter"
        fi
        out= "$outputFolder/$filter/$name_ok"
        astfits $file --copy=0 --primaryimghdu -o$out
        for ccd in $(seq 1 4); do
            h=$(((a - 1) * 4 + ccd))
            if [ "$warpBlock" -gt 0 ]; then
                astwarp $file --hdu=$h --scale=$warpBlock -otemp.fits
                astfits temp.fits --copy=1 -o$out
                rm temp.fits
            else
                astfits $file --copy=$h -o$out
            fi
        done
        ((a++))
    done
done



