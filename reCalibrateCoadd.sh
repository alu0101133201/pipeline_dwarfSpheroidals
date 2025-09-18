#!/bin/bash

# Script for recalibrating the coadd, for whatever reason (e.g. you used an incorrect camera response, not based on personal experience...)


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

# Load config file
confFile=$1
loadVariablesFromFile $confFile

ra=$ra_gal
dec=$dec_gal
export ra
export dec

filterCorrectionCoeff=$( checkIfNeededFilterCorrectionIsGiven $telescope $filter $surveyForPhotometry $ROOTDIR/"$objectName"/config )


# Step 1.- Calibration data

mosaicDir=$DIR/reCalibrationMosaic
selectedCalibrationStarsDir=$mosaicDir/automaticallySelectedStarsForCalibration
rangeUsedCalibrationDir=$mosaicDir/rangesUsedForCalibration
aperturePhotDir=$mosaicDir/aperturePhotometryCatalogues # This is the final product that "prepareCalibrationData" produces and will be used in "computeCalibrationFactors"
mosaicDone=$mosaicDir/done_prep.txt

echo -e "\n"
echo $surveyForPhotometry
echo $referenceImagesForMosaic
echo $aperturePhotDir
echo $filter
echo $ra
echo $dec
echo $mosaicDir
echo $selectedCalibrationStarsDir
echo $rangeUsedCalibrationDir
echo $pixelScale
echo $sizeOfOurFieldDegrees
echo $catName
echo $surveyForSpectra
echo $apertureUnits
echo $folderWithTransmittances
echo "$filterCorrectionCoeff"
echo $surveyCalibrationToGaiaBrightLimit 
echo $surveyCalibrationToGaiaFaintLimit 
echo $mosaicDone 
echo $sizeOfBrick
# We get the survey data needed for perform the calibration

# referenceImagesForMosaic = "" # This parameter is used for associating individual images to bricks. Since we are doing the calibration of the coadd we don't need it
# catName = "" # This was used when avoiding the download of bricks with bright stars. It is not used anymore so
# sizeOfBrick=1000

# prepareCalibrationData $surveyForPhotometry $referenceImagesForMosaic $aperturePhotDir $filter $ra $dec $mosaicDir $selectedCalibrationStarsDir $rangeUsedCalibrationDir \
#                                             $pixelScale $sizeOfOurFieldDegrees $catName $surveyForSpectra $apertureUnits $folderWithTransmittances "$filterCorrectionCoeff" \
#                                             $surveyCalibrationToGaiaBrightLimit $surveyCalibrationToGaiaFaintLimit $mosaicDone $sizeOfBrick
