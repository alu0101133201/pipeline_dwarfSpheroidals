#!/bin/bash

# This script will take some data and change its px scale

objectName=$1
dataDir=$2
destinationDir=$3
binningFactor=$4
filter=$5


if ! [ -d $destinationDir ]; then mkdir $destinationDir; fi

for i in $dataDir/*$objectName*$filter/*.fits; do
    fileName=$( basename $i)
    astwarp $i -h0 --scale=1/$binningFactor -o $destinationDir/$fileName
done

