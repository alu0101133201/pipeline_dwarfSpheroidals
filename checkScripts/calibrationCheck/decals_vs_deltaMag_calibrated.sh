module load gnuastro/0.22

GALAXY="IC1613"

decalsCatalogue="../../../$GALAXY/build_i/decals-aperture-catalog_it1"

myCatalogue_nonCalibrated="../../../$GALAXY/build_i/ourData-catalogs-apertures_it1"
myFrames_calibrated="../../../$GALAXY/build_i/photCorrSmallGrid-dir_it1"
aperturesForMyData_dir="../../../$GALAXY/build_i/my-catalog-halfmaxradius_it1"


numOfCatalogues1=$( ls $myCatalogue_nonCalibrated/*.cat | wc -l)
numOfCatalogues2=$( ls $decalsCatalogue/*.cat | wc -l)

if [ "$numOfCatalogues1" != "$numOfCatalogues2" ]; then
    echo "Error - The number of catalogues must be the same" >&2
    exit 1
fi


tmpDir="./calibratedCatalogues"
if ! [ -d $tmpDir ]; then mkdir $tmpDir; fi

for i in $myCatalogue_nonCalibrated/*.cat; do
    myFrame=$i
    frameNumber=$(echo "$i" | awk -F '[/]' '{print $(NF)}' | awk -F '[.]' '{print $(1)}' | awk -F '[_]' '{print $(NF)}')
    currentDecalsCatalogue=$decalsCatalogue/*_$frameNumber.*

    myCalibratedFrame=$myFrames_calibrated/entirecamera_$frameNumber.fits
    myNonCalibratedCatalogue=$myCatalogue_nonCalibrated/entirecamera_$frameNumber.fits*
    fileWithMyApertureData=$aperturesForMyData_dir/range1_entirecamera_$frameNumber*

    r_myData_pix_=$(awk 'NR==1 {printf $1}' $fileWithMyApertureData)
    r_myData_pix=$(astarithmetic $r_myData_pix_ 2. x -q )

    asttable $myNonCalibratedCatalogue -hSOURCE_ID -cRA,DEC | awk '!/^#/{print NR, $1, $2, 5, '$r_myData_pix', 0, 0, 1, NR, 1}' > $tmpDir/apertures.txt
    astmkprof $tmpDir/apertures.txt --background=$myCalibratedFrame --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=$tmpDir/aperture_myData.fits
        
    astmkcatalog $tmpDir/aperture_myData.fits -h1 --zeropoint=22.5 \
            --valuesfile=$myCalibratedFrame --valueshdu=1 \
            --ids --ra --dec --magnitude --sum \
            --output=$tmpDir/$frameNumber.cat


    astmatch $currentDecalsCatalogue --hdu=1 $tmpDir/$frameNumber.cat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=1/3600 --outcols=aMAGNITUDE,bMAGNITUDE -o$tmpDir/"$frameNumber"_matched.cat
    rm $tmpDir/apertures.txt
    rm $tmpDir/aperture_myData.fits
    rm $tmpDir/$frameNumber.cat
done

python3 decals_vs_deltaMag.py $tmpDir "calibrated"
rm -rf $tmpDir
