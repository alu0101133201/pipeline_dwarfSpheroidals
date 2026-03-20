#!/bin/bash
### This script performs photometry on images using the surface brightness profiles of galaxies. The script requires: an image of hipercam, and a .conf file containing the following data:
# 1) RA and DEC of the target galaxy
# 2) PA and axis ratio of the target galaxy.
# 3) Azimuth (in degrees) if wanted (default is 0)
# 4) Inner and outer radius (in arcsec) to perform the fitting
# 5) Survey of reference and filter
# 6) Directory of working
# The script will output a file with the photometry results, and the photometrized image.
# Usage: photometryWithGalaxies.sh  <config.conf> <image.fits> <output_file>
scriptPath=`dirname "$0"`
scriptPath=`( cd "$scriptPath" && pwd )`
pipelinePath=$(dirname "$scriptPath")
pythonScriptsPath=$pipelinePath/pipelineScripts
scriptPath=$pipelinePath/checkScripts
export pipelinePath
export pythonScriptsPath
export scriptPath

# Load functions of pipeline
source "$pipelinePath/pipeline_functions.sh"

#Loading gnuastro
gnuastroModuleName="gnuastro/0.22"
load_module $gnuastroModuleName

OPTSTRING=":h"
while getopts ${OPTSTRING} opt; do
  case ${opt} in
    h)
      help
      exit 0;;
    \?)
      echo "Invalid option: -${OPTARG}"
      exit 1;;
  esac
done


configFile=$1
imageFile=$2
outputPrefix=$3

if [ -z "$configFile" ] || [ ! -f "$configFile" ]; then
    echo "Error: Config file not provided or does not exist"
    echo "Usage: photometryWithGalaxies.sh <config.conf> <image.fits> <output_file>"
    exit 1
fi

loadVariablesFromFile $configFile

##First: download bricks for frame
mosaicDir=$DIR
surveyImagesDir=$mosaicDir/surveyImages
bricksIdentificationFile=$mosaicDir/bricksForFrame.txt
sizeOfOurFieldDegrees=0.005 # 18 arcsec to avoid edge effects
sizeOfBrick=100 #pixels
## We need a gaia catalogue to use pipeline functions, unless it wont be used
downloadCatalogue gaia $ra_gal $dec_gal $sizeOfOurFieldDegrees $mosaicDir $mosaicDir/gaia_catalgue.fits
downloadSurveyData $mosaicDir $surveyImagesDir $bricksIdentificationFile $filter $ra_gal $dec_gal $sizeOfOurFieldDegrees $surveyForPhotometry $sizeOfBrick

#If survey is DECALS we need to decompress, and if survey is PANSTARRS we need to recalibrate
if [ "$surveyForPhotometry" = "DECaLS" ]; then
    decompressBricks $surveyImagesDir
elif [ "$surveyForPhotometry" = "PANSTARRS" ]; then
    divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDir
elif [ "$surveyForPhotometry" = "SDSS" ]; then
    #No need to do anything for SDSS images
    echo "SDSS images do not require preprocessing"
fi

#Name of brick
case "$surveyForPhotometry" in
    DECaLS)
        brickName=$(ls $surveyImagesDir/decompressed_*.fits)
        pixel_scale=0.262
        ;;
    PANSTARRS)
        brickName=$(ls $surveyImagesDir/cal_*.fits)
        pixel_scale=0.25
        ;;
    SDSS)
        brickName=$(ls $surveyImagesDir/sdss*.fits) 
        pixel_scale=0.396
        ;;
    *)
        echo "Error: Survey $surveyForPhotometry not supported for photometry with galaxies"
        exit 1;;
esac

# Perform radial profile on image and on brick
radProfParams=" --mode=wcs --center=$ra_gal,$dec_gal --position-angle=$pa_gal --axis-ratio=$ar_gal --measure=mean,std --zeroisnotblank"
if [ -n "$azimuth_gal" ]; then
    radProfParams+=" --azimuth=$azimuth_gal"
fi
#Radial profile of hipercam
astscript-radial-profile $imageFile $radProfParams --rmax=50 --undersample=4 -o $DIR/RP_"$filter"_hipercam.fits
#Radial profile of survey brick
if [ "$surveyForPhotometry" = "SDSS" ]; then
    astscript-radial-profile $brickName $radProfParams --rmax=13 --undersample=1 -o $DIR/RP_"$filter"_"$surveyForPhotometry".fits
else
    astscript-radial-profile $brickName $radProfParams --rmax=25 --undersample=2 -o $DIR/RP_"$filter"_"$surveyForPhotometry".fits
fi

#Match radial profiles
python3 $scriptsPath/matchRadialProfiles.py $DIR/RP_"$filter"_hipercam.fits $DIR/RP_"$filter"_"$surveyForPhotometry".fits $innerRadiusArcsec $outerRadiusArcsec $filter $surveyForPhotometry $DIR/alfa_"$surveyForPhotometry"_"$filter".txt --scale1=0.16 --scale2=$pixel_scale --plot_name=${output_file%.fits}_profile_comparison.pdf --save-plot

#Apply photometric calibration to hipercam image
scalingFactor=$(cat $DIR/alfa_"$surveyForPhotometry"_"$filter".txt)
astarithmetic $imageFile $scalingFactor x float32 -o $output_file

echo "Photometrized image saved as $output_file"