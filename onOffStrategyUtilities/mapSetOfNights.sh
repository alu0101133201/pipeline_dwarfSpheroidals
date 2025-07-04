#!/bin/bash

#Example of usage:
# ./mapSetOfNights.sh 10 ~/Malin1/DATA-or/night ~/Malin1/DATA-or_offg/night ~/Malin1/offFieldBuilds/build_g_final/flat-it3-Running_n ./correctedFlats_n


numberOfNights=$1
objectImagesDir=$2
flatImagesDir=$3
flatsDir=$4
newFlatDir=$5

for ((i = 1; i <= numberOfNights; i++)); do
    ./prepareMappingBetweenObjectAndFlat.sh $objectImagesDir$i $flatImagesDir$i $flatsDir$i $newFlatDir$i $i
done