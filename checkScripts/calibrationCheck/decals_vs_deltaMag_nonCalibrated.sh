module load gnuastro/0.22

GALAXY="IC1613"

myCatalogue_nonCalibrated="../../../$GALAXY/build/ourData-catalogs-apertures_it1"
decalsCatalogue="../../../$GALAXY/build/decals-aperture-catalog_it1"

numOfCatalogues1=$( ls $myCatalogue_nonCalibrated/*.cat | wc -l)
numOfCatalogues2=$( ls $decalsCatalogue/*.cat | wc -l)

if [ "$numOfCatalogues1" != "$numOfCatalogues2" ]; then
    echo "Error - The number of catalogues must be the same" >&2
    exit 1
fi

tmpDir="./NonCalibratedCatalogues"
if ! [ -d $tmpDir ]; then mkdir $tmpDir; fi

for i in $myCatalogue_nonCalibrated/*.cat; do
    myFrame=$i
    frameNumber=$(echo "$i" | awk -F '[/]' '{print $(NF)}' | awk -F '[.]' '{print $(1)}' | awk -F '[_]' '{print $(NF)}')

    decalsFrame=$decalsCatalogue/*_$frameNumber.*

    astmatch $decalsFrame --hdu=1 $myFrame --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=1/3600 --outcols=aMAGNITUDE,bMAGNITUDE -o./$tmpDir/$frameNumber.cat
done

python3 decals_vs_deltaMag.py $tmpDir "nonCalibrated"
rm $tmpDir/*.cat
rmdir $tmpDir