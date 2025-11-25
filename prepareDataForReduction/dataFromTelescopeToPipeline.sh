#!/bin/bash

# This script performs several steps in order to prepare the data taken from the TST to the pipeline
# 1.- It rebins the data using the $rebinFactor argument
# 2.- It adds the headers of the reduced data to the raw data (we need to do this for reducing the data)
# 3.- It sorts all the data based on the filter
# 4.- It sorts all the data based on the date of the data

# If a filter ($6) is given it will only process data from that filter. Otherwise it will process
# the data of any filter


objectName=$1
rawImages=$2
redImages=$3
destination=$4
rebinFactor=$5
filter=$6
numCpus=$7

echo "Preparing data for the pipeline"
echo -e "\n\tObject: $objectName"
echo -e "\n\tRaw images path: $rawImages"
echo -e "\n\tRed images path (for the headers): $redImages"
echo -e "\n\tDestionation path: $destination"
echo -e "\n\tRebinn factor: $rebinFactor"
echo -e "\n\tFilter (if given): $filter"


if [ -z $filter ]; then
    filter=""
else
    filter=${filter}"*"
    # The following command goes through all the folders of $destination/$filter (if exists) in order to
    # take out all the files. This is simply not to compute again a frame already rebinned
    find $destination/$filter -mindepth 1 -type f -exec mv {} $destination \; 
    rm -rf $destination/$filter
fi


# # If $rebinFactor is greater than 0 that factor is applied
# # Otherwise the images are left at its original resolution
if (( $(echo "$rebinFactor > 0" | bc -l) )); then
    ./src/changePxScale.sh $objectName $rawImages $destination $rebinFactor $filter $numCpus
else
    if ! [ -d $destination ]; then mkdir $destination; fi
    for i in $rawImages/*$objectName*$filter/*.fits; do
        cp $i $destination
    done
fi

# Add the headers of the reduced data to the raw data
# The data HDU depends on if it has been warped (then hdu 0 is the gnuastro header) or not
headerHdu=0
if (( $(echo "$rebinFactor > 0" | bc -l) )); then
    dataHdu=1
else
    dataHdu=0
fi
python3 ./src/addHeaderToFrames.py $destination $redImages $dataHdu $headerHdu $rawImages $rebinFactor

# Ordering by filter
./src/orderFitsByFilters.sh $destination

# Moving data from hdu1 to hdu0 for the darks (this is needed because darks do not have reduced headers)
# But the pipeline expects the data in the hdu0, so we have to move the hdu (due to having warped gnuastro introduces an hdu0)
python3 ./src/removeFirstVoidHDU.py $destination/Dark

# Order by nights
for currentFilter in $destination/*; do
    ./src/orderFitsByNights.sh $currentFilter DATE-OBS
done
