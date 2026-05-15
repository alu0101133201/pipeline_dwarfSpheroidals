#!/bin/bash

mosaicDir=$1
catDir=$2
indexesDir=$3
objectName=$4
ra=$5
dec=$6
sizeOfField=$7
survey=$8

mkdir -p $mosaicDir
mkdir -p $catDir
mkdir -p $indexesDir

spectraDir=$mosaicDir/spectra
done=$spectraDir/done.txt
# python3 getSpectraFromBulk.py $mosaicDir $spectraDir $ra $dec $sizeOfField $survey
echo "done" > $done
 
source getGaiaCatalogue.sh
catName=$catDir/"$objectName"_gaia.fits
# downloadCatalogue "gaia" $ra $dec $sizeOfField $catDir $catName

python3 createScampCat.py ./fullCatalogueForScamp.fits ./"$objectName"_scampref.cat
rm ./fullCatalogueForScamp.fits


source getIndex.sh
indexes=()
for re in $(seq 4 10 ); do
  downloadIndex $re $catName $indexesDir
done
