#!/bin/bash

# Script for recalibrating the coadd, for whatever reason (e.g. you used an incorrect camera response, not based on personal experience...)
# The recalibrated coadd with updated header keywords will be stored in build_reCalibration/recalibratedCoadd



# Load the functions from the pipeline
pipelinePath=`dirname "$0"`
pipelinePath=`( cd "$pipelinePath" && pwd )`
pythonScriptsPath=$pipelinePath/pipelineScripts
export pipelinePath
export pythonScriptsPath

source "$pipelinePath/pipeline_functions.sh"

# Load modules

echo -e "\n ${GREEN} ---Loading Modules--- ${NOCOLOUR}"

gnuastroModuleName="gnuastro/0.22"
load_module $gnuastroModuleName

astrometryModuleName="astrometry.net/0.94"
load_module $astrometryModuleName

# Load config -----------------------------------------

confFile=$1
loadVariablesFromFile $confFile

filterCorrectionCoeff=$( checkIfNeededFilterCorrectionIsGiven $telescope $filter $surveyForPhotometry $ROOTDIR/"$objectName"/config )
if [[ $filterCorrectionCoeff == 11 ]]; then
  echo "The filter corrections for the filter $filter, telescope $telescope and survey $survey were not found"
  echo "Exiting with error code 11"
  exit 11
fi


ra=$ra_gal
dec=$dec_gal
export ra
export dec

num_cpus=$SLURM_CPUS_ON_NODE
if [ -z $num_cpus ]; then
  num_cpus=$defaultNumOfCPUs
fi
export num_cpus

DIR=$ROOTDIR/$objectName/build_reCalibration_$filter
BDIR=$DIR
mkdir -p $DIR
export DIR
export BDIR

tileSize=35
export tileSize


# Files to use ------------------------------------------------------------------

coaddToRecalibrateName=$objectName"_coadd_"$filter"_it2.fits"
coaddMask=$objectName"_coadd_"$filter"_mask.fits"
folderWithInputToRecalibrate=$ROOTDIR/$objectName/coadd_"$filter"_ToRecalibrate


exposureMap=$folderWithInputToRecalibrate/exposureMap.fits
coaddToRecalibrate=$folderWithInputToRecalibrate/$coaddToRecalibrateName
coaddMask=$folderWithInputToRecalibrate/$coaddMask

if [[ -f "$exposureMap" && -f "$coaddToRecalibrate" && -f "$coaddMask" ]]; then
    echo "Input files found, they are the following:"
    echo $exposureMap
    echo $coaddToRecalibrate
    echo $coaddMask
else
    echo "One or more files are missing"
    [[ ! -f "$exposureMap" ]] && echo "Missing: $exposureMap"
    [[ ! -f "$coaddToRecalibrate" ]] && echo "Missing: $coaddToRecalibrate"
    [[ ! -f "$coaddMask" ]] && echo "Missing: $coaddMask"
    exit 1
fi

# Gaia catalogue ---------------------------------------

if ((  $(echo "$sizeOfOurFieldDegrees > 1.0" | bc -l) )); then
  radiusToDownloadCatalogue=$( echo "$sizeOfOurFieldDegrees + 0.25" | bc -l | awk '{printf "%.1f", $0}' ) #The awk part is to avoiod problems when R<1
else
  radiusToDownloadCatalogue=$( echo "$sizeOfOurFieldDegrees + 0.5" | bc -l | awk '{printf "%.1f", $0}' ) #The awk part is to avoiod problems when R<1
fi

catdir=$BDIR/catalogs
catdone=$catdir/done.txt

surveys_to_download=("gaia")

if ! [ -d $catdir ]; then mkdir $catdir; fi
if [ -f $catdone ]; then
  echo -e "\n\tCatalogue is already downloaded\n"
else
  for survey in "${surveys_to_download[@]}"; do
      catName=$catdir/"$objectName"_"$survey".fits
      catRegionName=$catdir/"$objectName"_"$survey"_regions.reg
      
      downloadCatalogue $survey $ra_gal $dec_gal $radiusToDownloadCatalogue $catdir $catName
  done
  echo "done" > $catdone
fi

# We leave the variables of the catalogue names with the selected survey. This is will be used in next steps
catName=$catdir/"$objectName"_"$surveyToUseInSolveField".fits
catRegionName=$catdir/"$objectName"_"$surveyToUseInSolveField"_regions.reg

# We get rid of bright stars (based on the calibration range) because saturated stars might make our estimation of point-like region wrong
if ! [ "$surveyToUseInSolveField" = "gaia" ]; then
  mv $BDIR/catalogs/"$objectName"_gaia.fits $BDIR/catalogs/"$objectName"_gaia_tmp.fits
  asttable $BDIR/catalogs/"$objectName"_gaia_tmp.fits --range=phot_g_mean_mag,$calibrationBrightLimitIndividualFrames,30 -o $BDIR/catalogs/"$objectName"_gaia.fits
  rm $BDIR/catalogs/"$objectName"_gaia_tmp.fits
fi



# Step 2.- Calibration data ------------------------------------------------

# Variables needed

toleranceForMatching=1.5 #arcsec
sigmaForPLRegion=3 # Parameter for deciding the selection region (half-max-rad region)
export toleranceForMatching
export sigmaForPLRegion

# Parameters for performing the sigma clipping to the different samples in aststatistics
sigmaForStdSigclip=2
iterationsForStdSigClip=3
export sigmaForStdSigclip
export iterationsForStdSigClip

mosaicDir=$DIR/reCalibrationMosaic
selectedCalibrationStarsDir=$mosaicDir/automaticallySelectedStarsForCalibration
rangeUsedCalibrationDir=$mosaicDir/rangesUsedForCalibration
aperturePhotDir=$mosaicDir/aperturePhotometryCatalogues # This is the final product that "prepareCalibrationData" produces and will be used in "computeCalibrationFactors"
mosaicDone=$mosaicDir/done_prep.txt

# We get the survey data needed for perform the calibration
referenceImagesForMosaic="-" # This parameter is used for associating individual images to bricks. Since we are doing the calibration of the coadd we don't need it
catName="-" # This was used when avoiding the download of bricks with bright stars. It is not used anymore so
sizeOfBrick=1000

prepareCalibrationData $surveyForPhotometry $referenceImagesForMosaic $aperturePhotDir $filter $ra $dec $mosaicDir $selectedCalibrationStarsDir $rangeUsedCalibrationDir \
                                            $pixelScale $sizeOfOurFieldDegrees $catName $surveyForSpectra $apertureUnits $folderWithTransmittances "$filterCorrectionCoeff" \
                                            $surveyCalibrationToGaiaBrightLimit $surveyCalibrationToGaiaFaintLimit $mosaicDone $sizeOfBrick



# 3.- Coadd recalibration ---------------------------------------------------------


# Computing calibration factor

iteration=1
alphatruedir=$BDIR/alpha-stars-true_coaddPrephot_it$iteration
matchdir=$BDIR/match-decals-myData_coaddPrephot_it$iteration
ourDataCatalogueDir=$BDIR/ourData-aperture-photometry_coaddPrephot_it$iteration
prepareCalibrationCataloguePerFrame=$BDIR/survey-aperture-photometry_perBrick_coaddPrephot_it$iteration
mycatdir=$BDIR/my-catalog-halfmaxradius_coaddPrephot_it$iteration
calibratingMosaic=true
imagesForCalibration=$BDIR/coaddForCalibration_it$iteration

mkdir -p $BDIR/coaddForCalibration_it$iteration
# Calibrate the coadd only in the high snr section 
expMax=$(aststatistics $exposureMap --maximum -q)
exp_fr=$(astarithmetic $expMax 0.5 x -q)
astarithmetic $coaddToRecalibrate $exposureMap -g1 $exp_fr lt nan where --output=$imagesForCalibration/entirecamera_1.fits


computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $mycatdir $rangeUsedCalibrationDir \
                          $mosaicDir $alphatruedir $calibrationBrightLimitCoaddPrephot $calibrationFaintLimitCoaddPrephot $tileSize $apertureUnits $numberOfApertureUnitsForCalibration $calibratingMosaic 


# Apply calibration factor

recalibratedCoadd=$BDIR/recalibratedCoadd
applyCalibrationFactors $imagesForCalibration $alphatruedir $recalibratedCoadd $iteration False


# Calibration plot
aperturesFolder=$BDIR/my-catalog-halfmaxradius_coaddPrephot_it$iteration
calibrationPlotName=$BDIR/calibrationPlot_coaddPrephot.png
onlyPointLikeCat=$BDIR/my-catalog-halfmaxradius_coaddPrephot_it$iteration

if [ -f $calibrationPlotName ]; then
    echo -e "\nCalibration diagnosis plot for coadd prephot already done\n"
else
    if [[ "$surveyForCalibration" == "SPECTRA" ]]; then
      dirWithReferenceCat=$mosaicDir
    else
      dirWithReferenceCat=$BDIR/survey-aperture-photometry_perBrick_coaddPrephot_it1
    fi
    mosaicPlot=true
    produceCalibrationCheckPlot $BDIR/ourData-aperture-photometry_coaddPrephot_it1 $recalibratedCoadd $aperturesFolder $dirWithReferenceCat \
                                  $pythonScriptsPath $calibrationPlotName $calibrationBrightLimitCoaddPrephot $calibrationFaintLimitCoaddPrephot $numberOfApertureUnitsForCalibration $BDIR $surveyForPhotometry $BDIR $mosaicPlot $BDIR/calibratedCatalogue_prehot_it$iteration $onlyPointLikeCat
fi




# Since we calibrated the coadd in the high snr region, we need to restore the whole image
# I don't do that before so the calibration plot is with the calibrated area

# We clean the input and output of the calibrated high-snr crop of the coadd and we get it again of the compelte field
rm $recalibratedCoadd/*
rm $imagesForCalibration/*
cp $coaddToRecalibrate $imagesForCalibration/entirecamera_1.fits
applyCalibrationFactors $imagesForCalibration $alphatruedir $recalibratedCoadd $iteration False
mv $recalibratedCoadd/entirecamera_1.fits $recalibratedCoadd/$coaddToRecalibrateName



# Compute surface brightness limit

sblimitFile=$recalibratedCoadd/"$objectName"_"$filter"_sblimit.txt
exposuremapName=$exposureMap
maskName=$coaddMask
surfaceBrightnessLimit=$( limitingSurfaceBrightness $recalibratedCoadd/$coaddToRecalibrateName $maskName $exposuremapName $recalibratedCoadd $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )

# Propagate variables


keyWords=("FRAMES_COMBINED" \
          "NUMBER_OF_DIFFERENT_NIGHTS" \
          "INITIAL_DATE_OBS" \
          "MEAN_DATE_OBS" \
          "FINAL_DATE_OBS" \
          "FILTER" \
          "LOWER_VIGNETTING_THRESHOLD" \
          "UPPER_VIGNETTING_THRESHOLD" \
          "SATURATION_THRESHOLD" \
          "CALIBRATED_USING" \
          "CALIBRATION_BRIGHTLIMIT" \
          "CALIBRATION_FAINTLIMIT" \
          "RUNNING_FLAT" \
          "WINDOW_SIZE" \
          "SURFACE_BRIGHTNESS_LIMIT")
          
comments=("" "" "" "" "" "" "" "" "" "" "" "" "" "Running flat built with +-N frames" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")

# values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$lowerVignettingThreshold" "$upperVignettingThreshold" "$saturationThreshold" "$surveyForPhotometry" "$calibrationBrightLimitCoaddPrephot" "$calibrationFaintLimitCoaddPrephot" "$RUNNING_FLAT" "$halfWindowSize" "$surfaceBrightnessLimit")


values=()
for i in "${!keyWords[@]}"; do
  currentKeyWord="${keyWords[$i]}"
  values+=($(astfits $coaddToRecalibrate -h1 --keyvalue=$currentKeyWord -q))
done

addkeywords $recalibratedCoadd/$coaddToRecalibrateName keyWords values comments
astfits $recalibratedCoadd/$coaddToRecalibrateName --update=SURFACE_BRIGHTNESS_LIMIT,$surfaceBrightnessLimit,"[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)"


echo "Depth of the previous coadd: " $( astfits $coaddToRecalibrate -h1 --keyvalue=surface_brightness_limit -q)
echo "Depth of the new coadd: " $( astfits $recalibratedCoadd/$coaddToRecalibrateName -h1 --keyvalue=surface_brightness_limit -q)

