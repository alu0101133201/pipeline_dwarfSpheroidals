
# This folders are the ones that the pipeline produces in this moment.
# If things change in the future or you want to use this script in another context just modify how to access the data
GALAXY="IC1613"

myCatalogue="../../../$GALAXY/build_g/my-catalog-halfmaxradius_it1"
decalsCatalogue="../../../$GALAXY/build_g/decals-aperture-catalog_it1"

numOfCatalogues1=$( ls $myCatalogue/match*.txt | wc -l)
numOfCatalogues2=$( ls $decalsCatalogue/decals*.cat | wc -l)

if [ "$numOfCatalogues1" != "$numOfCatalogues2" ]; then
    echo "Error - The number of catalogues must be the same" >&2
    exit 1
fi

for i in $myCatalogue/match*.txt; do
    myFrame=$i
    frameNumber=$(echo "$i" | awk -F '[/]' '{print $(NF)}' | awk -F '[.]' '{print $(1)}' | awk -F '[_]' '{print $(NF)}')
    decalsFrame=$decalsCatalogue/*_$frameNumber.*

    astmatch $decalsFrame --hdu=1 $myFrame --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=1/3600 --outcols=aRA,aDEC,bRA,bDEC -o./tmp$frameNumber.cat
done

python3 deltaRAdeltaDecPlot.py

rm ./tmp*
