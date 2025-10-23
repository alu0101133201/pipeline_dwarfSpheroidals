#!/bin/bash
########## Functions ##########

# Functions of general sue
help() {
    echo -e "Pipeline for reducing data from small telescopes"
    echo -e ""
    echo -e "Syntax: pipelineName [-h] configurationFile.conf"
    echo -e "Options:"
    echo -e "    h    Print this help"
    echo -e "Arguments:"
    echo -e "    Configuration file"
    echo -e "\n"
}
export -f help

load_module() {
    local moduleName="$1"  
    errorNumber=8

    if [[ -z "$moduleName" ]]; then
        echo "Error: No module name provided"
        return $errorNumber
    fi

    echo -e "\nLoading $moduleName"
    module load "$moduleName"

    if [[ $? -eq 0 ]]; then
        echo -e "$moduleName loaded successfully"
    else
        echo -e "Failed to load $moduleName"
        # return 1 # I comment this because in the ICR that I'm running they are already loaded so...
    fi
}
export -f load_module

writeTimeOfStepToFile() {
    local step=$1
    local file=$2
    echo "Step: $step. Start time:  $(date +%D-%T)" >> $file
}
export -f writeTimeOfStepToFile

loadVariablesFromFile() {
  local file=$1
  initialShellVariables=$(compgen -v)

  if [[ -f $confFile ]]; then 
    source $confFile
    echo -e "\nVariables loaded from $confFile file"
  else
    errorNumber=1
    echo -e "\nA configuration file has to be provided in order to run the pipeline"  >&2
    echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
    exit $errorNumber
  fi

  allShellVariablesAfterLoadingConfFile=$(compgen -v)
  
  # This code exports only the variables of the configuration file
  for var in $allShellVariablesAfterLoadingConfFile; do
    if ! grep -q "^$var$" <<< "$initialShellVariables"; then
      export "$var"
    fi
  done
  return 0
}
export -f loadVariablesFromFile

outputConfigurationVariablesInformation() {
    data=(
        "·Telescope:$telescope"
        "·Object name:$objectName"
        "·Right ascension:$ra_gal:[deg]"
        "·Declination:$dec_gal:[deg]"
        "·Latitud of the telescope:$telescopeLat:[deg]"
        "·Longitude of the telescope:$telescopeLong:[deg]"
        "·Elevation of the telescope:$telescopeElevation:[deg]"
        ""
        "·Keyword for the airmass:$airMassKeyWord"
        "·Keyword for date:$dateHeaderKey"
        "·Keyword for RA of the pointings:$pointingRA:[$pointingRAUnits]"
        "·Keyword for DEC of the pointings:$pointingDEC:[$pointingDECUnits]"
        "·Root directory to perform the reduction:$ROOTDIR"
        ""
        "·Calibration range"
        "  Bright limit individual frames:$calibrationBrightLimitIndividualFrames:[mag]"
        "  Faint limit individual frames:$calibrationFaintLimitIndividualFrames:[mag]"
        "  Bright limit coadd prephot:$calibrationBrightLimitCoaddPrephot:[mag]"
        "  Faint limit coadd prephot:$calibrationFaintLimitCoaddPrephot:[mag]"

        "·The aperture photometry for calibrating our data is in units of:$apertureUnits"
        "·The aperture photometry will be done with an aperture of:$numberOfApertureUnitsForCalibration:[$apertureUnits]"
        "·The calibration will be done using data from survey:$surveyForPhotometry"
        "·If the calibration is done with survey, the range for calibrating the survey to GAIA is:"

        "·If the calibration is done with spectra, the survey to use is:$surveyForSpectra"
        "·The transmittances of the filters are in the folder:$folderWithTransmittances"
        "  Bright limit:$surveyCalibrationToGaiaBrightLimit:[mag]"
        "  Faint limit:$surveyCalibrationToGaiaFaintLimit:[mag]"
        "·Saturation threshold:$saturationThreshold:[ADU]"
        "·Gain:$gain:[e-/ADU]"
        "·Approximately size of the field:$sizeOfOurFieldDegrees:[deg]"
        "·Size of the coadd:$coaddSizePx:[px]"
        "·Lower Vignetting threshold (mask every px with flat value below this one):$lowerVignettingThreshold "
        "·Upper Vignetting threshold (mask every px with flat value above this one):$upperVignettingThreshold "

        " "
        "·Width of the normalisation ring:$ringWidth:[px]"
        "·Common normalisation ring:$USE_COMMON_RING"
        "  If so, the file with the ring specification is:$commonRingDefinitionFile"
        ""
        " Otherwise, parameters for using multiple normalisation rings"
        "  KeyWord to decide the ring to use:$keyWordToDecideRing"
        "  Threshold to decide the ring to use:$keyWordThreshold"
        "  File with the first ring specification:$firstRingDefinitionFile"
        "  Value of the keyword for using the first ring:$keyWordValueForFirstRing"

        "  File with the second ring specification:$secondRingDefinitionFile"
        "  Value of the keyword for using the second ring:$keyWordValueForSecondRing"
        " "
        "·Running flat:$RUNNING_FLAT"
        "  If so, half of the window size is:$halfWindowSize:[frames]"
        " "
        "·The background is modelled as a constant:$MODEL_SKY_AS_CONSTANT"
        "  If so, the sky estimation method is:$sky_estimation_method"
        "  Otherwise, the polynomial degree is:$polynomialDegree"
        "  Noisechisel will be run with the following params: $noisechisel_param"
        "  Prior to noisechisel, a block will be applied with value: $blockScale (1=No block)"
        " "
        "·Indices scales for astrometrisation"
        "  Lowest index:$lowestScaleForIndex"
        "  Highest index:$highestScaleForIndex"
        " "
        "·Scale-low parameters for solve-field (astrometry.net):$solve_field_L_Param"
        "·Scale-high parameters for solve-field (astrometry.net):$solve_field_H_Param"
        "·Scale units for parameters for solve-field (astrometry.net):$solve_field_u_Param"
        " "
        "·Filter:$filter"
        "·Pixel scale:$pixelScale:[arcsec/px]"
        "·Detector width:$detectorWidth:[px]"
        "·Detector height:$detectorHeight:[px]"
        "·Is there overscan:$overscan"
        "·Keyword for illuminated section:$trimsecKey"
        "·Number of CCDs of the camers:$num_ccd"
         " "
        "Parameters for rejecting frames"
        "·Maxmimum background brightness:$maximumBackgroundBrightness:[mag/arcsec²]"
        "·Maxmimum seeing:$maximumSeeing:[FWHM]"
        " "
        "Parameters for measuring the surface brightness limit"
        "·Exp map fraction:$fractionExpMap"
        "·Area of the SB limit metric:$areaSBlimit: [arcsec]"
        " "
        "·Produce coadd prephot:$produceCoaddPrephot"
    )

    echo -e "Summary of the configuration variables provided for the reduction\n"
    for entry in "${data[@]}"; do
    IFS=":" read -r text value unit <<< "$entry" 
    printf "\t%-60s $ORANGE %-20s $GREEN %-10s $NOCOLOUR\n" "$text" "$value" "$unit"
    done
}
export -f outputConfigurationVariablesInformation

escapeSpacesFromString() {
    local input="$1"
    escaped_string="${input// /\\ }"
    echo $escaped_string
}
export -f escapeSpacesFromString

checkIfExist_DATEOBS() {
    local DATEOBSValue=$1

    if [ "$DATEOBSValue" = "n/a" ]; then
        errorNumber=3
        echo -e "The file $i do not has the $dateHeaderKey, used for sorting the raw files for the pipeline"  >&2
        echo -e "Exiting with error number: $errorNumber"  >&2
        exit $errorNumber
    fi
}
export -f checkIfExist_DATEOBS

getHighestNumberFromFilesInFolder() {
    local folderToCheck=$1
    # Find all files, remove non-numeric parts, sort numerically, get the highest number
    highest=$(ls "$folderToCheck" | grep -oE '[0-9]+' | sort -n | tail -1)
    if [ -z "$highest" ]; then
            highest=0
    fi
    echo $highest
}
export -f getHighestNumberFromFilesInFolder

checkIfAllVariablesAreSet() {
    errorNumber=2
    flagToExit=""
    variablesToCheck=(objectName \
                ra_gal \
                dec_gal \
                telescope \
                defaultNumOfCPUs \
                ROOTDIR \
                airMassKeyWord \ 
                dateHeaderKey \
                saturationThreshold \
                gain \
                sizeOfOurFieldDegrees \
                coaddSizePx \
                lowerVignettingThreshold \
		upperVignettingThreshold \
                calibrationBrightLimitIndividualFrames \
                calibrationFaintLimitIndividualFrames \
                calibrationBrightLimitCoaddPrephot \
                calibrationFaintLimitCoaddPrephot \
                apertureUnits \
                numberOfApertureUnitsForCalibration \
                surveyForPhotometry \
                folderWithTransmittances \
		surveyCalibrationToGaiaBrightLimit \
                surveyCalibrationToGaiaFaintLimit \
                surveyForSpectra \
                ringWidth \
                USE_COMMON_RING \
                commonRingDefinitionFile \
                keyWordToDecideRing
                keyWordThreshold
                firstRingDefinitionFile
                keyWordValueForFirstRing
                secondRingDefinitionFile
                keyWordValueForSecondRing
                RUNNING_FLAT \
                halfWindowSize \
                MODEL_SKY_AS_CONSTANT \
                sky_estimation_method \
                polynomialDegree \
                noisechisel_param \
                blockScale \
                filter \
                pixelScale \
                detectorWidth \
                detectorHeight \ 
                num_ccd \
                overscan \
                trimsecKey \
                lowestScaleForIndex \
                highestScaleForIndex \ 
                solve_field_L_Param \
                solve_field_H_Param \
                solve_field_u_Param \ 
                maximumBackgroundBrightness \
                maximumSeeing \
                fractionExpMap\
                areaSBlimit \
		produceCoaddPrephot)

    echo -e "\n"
    for currentVar in ${variablesToCheck[@]}; do
        [[ -z ${!currentVar} ]] && echo "${currentVar} variable not defined" && flagToExit=true
    done

    # I exit here and not when I find the variable missing because I want to show all the messages of "___ variable not defined", so the user knows all the variables that are needed
    [[ $flagToExit ]] && echo -e "Exiting with error number: $errorNumber" && exit $errorNumber
}
export -f checkIfAllVariablesAreSet

checkIfStringVariablesHaveValidValues() {
    errorCode=10

    if [[ ("$apertureUnits" != "FWHM") && ("$apertureUnits" != "Re") ]]; then
        echo "Error. The variable apertureUnits has a value ($apertureUnits) which is not accepted"
        exit $errorCode
    fi

    if [[ ("$surveyForPhotometry" != "PANSTARRS") && ("$surveyForPhotometry" != "DECaLS") && ("$surveyForPhotometry" != "SPECTRA")]]; then
        echo "Error. The variable surveyForPhotometry has a value ($surveyForPhotometry) which is not accepted"
        exit $errorCode
    fi

    if [[ ("$surveyForPhotometry" == "SPECTRA") && ("$surveyForSpectra" != "SDSS" && "$surveyForSpectra" != "GAIA") ]]; then
        echo "Error. The variable surveyForSpectra has a value ($surveyForSpectra) which is not accepted"
        exit $errorCode
    fi

}
export -f checkIfStringVariablesHaveValidValues

checkTransmittanceFilterAndItsUnits() {
    local telescope=$1
    local survey=$2
    local filterFolder=$3
    local filterToUse=$4

    filterFileNeeded=$( checkIfAllTheTransmittancesNeededAreGiven $telescope $surveyForPhotometry $folderWithTransmittances $filter )
    checkUnitsAndConvertToCommonUnitsIfNeeded $filterFileNeeded
}
export -f checkTransmittanceFilterAndItsUnits

checkIfNeededFilterCorrectionIsGiven() {
    local telescope=$1
    local filter=$2
    local survey=$3
    local configDir=$4

    if [[ ("$survey" == "SPECTRA") ]]; then
        return # If we calibrate with spectra then we don't have to correct between filters
    else
        fileWithFilterCorrections=$configDir/filterCorrections.dat

        # Get the coefficients of the correction
        coefficients=$( awk -v tel="$telescope" -v surv="$survey" -v fil="$filter" '
            BEGIN { found=0 }
            /^telescope:/ { found=(tolower($2) == tolower(tel)) ? 1 : 0 }
            /^reference survey:/ { if (found && tolower($3) != tolower(surv)) found=0 }
            found && tolower($1) == tolower(fil)":" { print $2, $3, $4 }
        ' $fileWithFilterCorrections )

        if [[ -z "$coefficients" ]]; then
            errorCode=11
            echo $errorCode
        else
            echo $coefficients
        fi
    fi
}
export -f checkIfNeededFilterCorrectionIsGiven

checkUnitsAndConvertToCommonUnitsIfNeeded() {
    local filterFile=$1
    errorCode=12

    # Since we may change the file (depending on the units) I create a backup in order not to loose the original file
    cp $filterFile $filterFile.original

    read units transmittanceFormat < <(head -n 1 $filterFile)
    transmittanceFormat="${transmittanceFormat//$'\r'/}" # Remove the final endline

    if [[ ("$units" != "A") && ("$units" != "nm" ) ]]; then
        echo "Units ($units) for transmittance wavelengths not accepted. It is expected either nm or A"
        exit $errorCode
    fi

    if [[ ("$transmittanceFormat" != "normalised") && ("$transmittanceFormat" != "percentage") ]]; then
        echo "format ($transmittanceFormat) for transmittance not accepted. It is expected either normalised or percentage"
        exit $errorCode
    fi

    awk 'NR>1 {print $1, $2}' $filterFile.original |  
    {
        if [[ "$units" == "nm" ]]; then
            awk '{print $1 * 10, $2}'  # Convert nm to Å
        else
            cat  # Keep as is
        fi
    } | {
        if [[ "$transmittanceFormat" == "percentage" ]]; then
            awk '{print $1, $2 / 100}'  # Convert percentage to normalised
        else
            cat  # Keep as is
        fi
    } > "$filterFile"

    sed -i '1s/^/A normalised\n/' "$filterFile"
}
export -f checkUnitsAndConvertToCommonUnitsIfNeeded

checkIfAllTheTransmittancesNeededAreGiven() {
    local telescope=$1
    local survey=$2
    local filterFolder=$3
    local filterToUse=$4

    errroCode=11

    # If we calibrate with spectra (survey="SPECTRA") then only the transmittance of the filter to use is needed
    if [[ "$survey" == "SPECTRA" ]]; then
        filterPath=$filterFolder/"$telescope"_"$filterToUse".dat
        if ! [ -e $filterPath ]; then
            echo "Error. Calibrating with SPECTRA option but not found the transmittance of the filter to calibrate (expected file: $filterPath)"
            exit $errorCode  
        fi
    else   
        # If calibrating with imaging survey then the filter from the survey is needed (to calibrate it to GAIA). 
        # No filter from the telescope of the data to reduce is needed because the imaging survey will be used
        filterPath=$filterFolder/"$survey"_"$filterToUse".dat
        if ! [ -e $filterPath ]; then
            echo "Error. Calibrating with IMAGING SURVEY ($survey) option but not found the transmittance of the filter to calibrate this survey to our reference (GAIA). expected file: $filterPath"
            exit $errorCode  
        fi
    fi
    echo $filterPath
}
export -f checkIfAllTheTransmittancesNeededAreGiven

subtractBiasFromFrame(){
    local base=$1
    local dark=$2
    local satThres=$3
    local inputDir=$4
    local outDir=$5

    i=$inputDir/$base
    out=$outDir/$base

    for h in $(seq 0 $num_ccd); do
        if [ $h -eq 0 ]; then
            astfits $1 --copy=$h --primaryimghdu -o $out
            if [ "$USE_COMMON_RING" = false ]; then
                propagateKeyword $i $keyWordToDecideRing $out 0
            fi
        else
            astarithmetic $i -h$h set-i $dark -h$h set-m \
                  i i $satThres gt i isblank or 2 dilate nan where m - float32 \
                  -o $outDir/temp_$base
            astfits $outDir/temp_$base --copy=1 -o $out
            rm $outDir/temp_$base
            propagateKeyword $i $gain $out $H
        fi
    done
}
export -f subtractBiasFromFrame

# Functions used in Flat
maskImages() {
    local inputDirectory=$1
    local masksDirectory=$2
    local outputDirectory=$3
    local useCommonRing=$4
    local keyWordToDecideRing=$5
    imagesToMask=()
    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
        imagesToMask+=("$base")
    done
    printf "%s\n" "${imagesToMask[@]}" | parallel -j "$num_cpus" maskIndividualImage {} $inputDirectory $masksDirectory $outputDirectory $useCommonRing $keyWordToDecideRing
}
export -f maskImages

maskIndividualImage() {
    local base=$1
    local inputDirectory=$2
    local masksDirectory=$3 
    local outputDirectory=$4
    local useCommonRing=$5
    local keyWordToDecideRing=$6
    i=$inputDirectory/$base
    out=$outputDirectory/$base
    astfits $i --copy=0 --primaryimghdu -o $out
    for h in $(seq 1 $num_ccd); do
        astarithmetic $i -h$h $masksDirectory/$base -h$h 1 eq nan where float32 -o $outputDirectory/temp_"$base" -q
        astfits $outputDirectory/temp_"$base" --copy=1 -o $out
        rm $outputDirectory/temp_"$base"
    done
        #propagateKeyword $i $airMassKeyWord $out 
        # If we are not doing a normalisation with a common ring we propagate the keyword that will be used to decide
        # which ring is to be used. This way we can check this value in a comfortable way in the normalisation section
    if [ "$useCommonRing" = false ]; then
        propagateKeyword $i $keyWordToDecideRing $out 0
    fi
}
export -f maskIndividualImage


getInitialMidAndFinalFrameTimes() {
  local directoryWithNights=$1
  local dateKey=$2
  declare -a date_obs_array
  
  while IFS= read -r -d '' file; do
    currentDateObs=$( gethead $file $dateKey)

    if [[ -n "$currentDateObs" ]]; then
        if [[ $dateKey == "MJD-OBS" ]]; then
          unixTime=$(astarithmetic $currentDateObs 40587 - 86400 x -q)
          unixTime=$(printf "%.0f" "$unixTimeInSeconds")
          
	    else
	        ## MACOS does not support -d in date, so it is better to use coreutils:gdata
	        if [[ $OSTYPE == 'darwin'* ]]; then
	        	unixTime=$(gdate -d "$currentDateObs" +"%s")
	        else
	        	unixTime=$(date -d "$currentDateObs" +"%s")
	        fi
        fi
        date_obs_array+=("$unixTime")
    fi
  done < <(find "$directoryWithNights" -type f -name "*.fits" -print0)

  sortedDateObsArray=($(for date in "${date_obs_array[@]}"; do echo "$date"; done | sort))
  arrayLength=$( echo "${#sortedDateObsArray[@]}" )

  initialTime=${sortedDateObsArray[0]}
  meanTime=${sortedDateObsArray[(( (( $arrayLength - 1)) / 2))]}
  finalTime=${sortedDateObsArray[(( $arrayLength - 1))]}
  echo "$initialTime $meanTime $finalTime"
}
export -f getInitialMidAndFinalFrameTimes


writeKeywordToFits() {
    local fitsFile=$1
    local header=$2
    local keyWord=$3
    local value=$4
    local comment=$5

    astfits --write=$keyWord,$value,"$comment" $fitsFile -h$header
}
export -f writeKeywordToFits

propagateKeyword() {
    local image=$1
    local keyWordToPropagate=$2
    local out=$3
    local h=$4
    variableToDecideRingToNormalise=$(gethead $image $keyWordToPropagate -x $h)
    eval "astfits --write=$keyWordToPropagate,$variableToDecideRingToNormalise $out -h$h" 
}
export -f propagateKeyword

addkeywords() {
    local fits_file=$1
    shift
    local -n keys_array=$1
    local -n values_array=$2
    local -n comments_array=$3

    if [[ -z "$fits_file" || ${#keys_array[@]} -eq 0 || ${#values_array[@]} -eq 0 ]]; then
        errorNumber=7
        echo -e "Error in 'addkeywords', some argument is empty"
        echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
        exit $errorNumber 
    fi

    if [[ ${#keys_array[@]} -ne ${#values_array[@]} ]]; then
        echo -e "Error in 'addkeywords', the length of keys and values does not match"
        echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
        exit $errorNumber   
    fi

    for i in "${!keys_array[@]}"; do
        local key="${keys_array[$i]}"
        local value="${values_array[$i]}"
        local comment="${comments_array[$i]}"

        writeKeywordToFits $fits_file 1 "$key" "$value" "$comment"
    done
}
export -f addkeywords


getMedianValueInsideRing() {
    local i=$1
    local commonRing=$2
    local doubleRing_first=$3
    local doubleRing_second=$4
    local useCommonRing=$5
    local keyWordToDecideRing=$6
    local keyWordThreshold=$7
    local keyWordValueForFirstRing=$8
    local keyWordValueForSecondRing=$9
    local h=${10}

    if [ "$useCommonRing" = true ]; then
            # Case when we have one common normalisation ring
            me=$(astarithmetic $i -h$h $commonRing -h$h 0 eq nan where medianvalue --quiet)
    else
        # Case when we do NOT have one common normalisation ring
        # All the following logic is to decide which normalisation ring apply
        variableToDecideRingToNormalise=$(gethead $i $keyWordToDecideRing)
        firstRingLowerBound=$(echo "$keyWordValueForFirstRing - $keyWordThreshold" | bc)
        firstRingUpperBound=$(echo "$keyWordValueForFirstRing + $keyWordThreshold" | bc)
        secondRingLowerBound=$(echo "$keyWordValueForSecondRing - $keyWordThreshold" | bc)
        secondRingUpperBound=$(echo "$keyWordValueForSecondRing + $keyWordThreshold" | bc)

        if (( $(echo "$variableToDecideRingToNormalise >= $firstRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $firstRingUpperBound" | bc -l) )); then
            me=$(astarithmetic $i -h$h $doubleRing_first -h$h 0 eq nan where medianvalue --quiet)
        elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
            me=$(astarithmetic $i -h$h $doubleRing_second -h$h 0 eq nan where medianvalue --quiet)
        else
            errorNumber=4
            echo -e "\nMultiple normalisation ring have been tried to be used. The keyword selection value of one has not matched with the ranges provided" >&2
            echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
            exit $errorNumber 
        fi
    fi
    echo $me # This is for "returning" the value
}
export -f getMedianValueInsideRing

getStdValueInsideRing() {
    local i=$1
    local commonRing=$2
    local doubleRing_first=$3
    local doubleRing_second=$4
    local useCommonRing=$5
    local keyWordToDecideRing=$6
    local keyWordThreshold=$7
    local keyWordValueForFirstRing=$8
    local keyWordValueForSecondRing=$9
    local h=${10}

    if [ "$useCommonRing" = true ]; then
            # Case when we have one common normalisation ring
            std=$(astarithmetic $i -h$h $commonRing -h$h 0 eq nan where stdvalue --quiet)
    else
        # Case when we do NOT have one common normalisation ring
        # All the following logic is to decide which normalisation ring apply
        variableToDecideRingToNormalise=$(gethead $i $keyWordToDecideRing)
        firstRingLowerBound=$(echo "$keyWordValueForFirstRing - $keyWordThreshold" | bc)
        firstRingUpperBound=$(echo "$keyWordValueForFirstRing + $keyWordThreshold" | bc)
        secondRingLowerBound=$(echo "$keyWordValueForSecondRing - $keyWordThreshold" | bc)
        secondRingUpperBound=$(echo "$keyWordValueForSecondRing + $keyWordThreshold" | bc)

        if (( $(echo "$variableToDecideRingToNormalise >= $firstRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $firstRingUpperBound" | bc -l) )); then
            std=$(astarithmetic $i -h$h $doubleRing_first -h$h 0 eq nan where stdvalue --quiet)
        elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
            std=$(astarithmetic $i -h$h $doubleRing_second -h$h 0 eq nan where stdvalue --quiet)
        else
            errorNumber=5
            echo -e "\nMultiple normalisation ring have been tried to be used. The keyword selection value of one has not matched with the ranges provided" >&2
            echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
            exit $errorNumber 
        fi
    fi

    echo $std # This is for "returning" the value
}
export -f getStdValueInsideRing

getSkewKurtoValueInsideRing(){
    local i=$1
    local commonRing=$2
    local doubleRing_first=$3
    local doubleRing_second=$4
    local useCommonRing=$5
    local keyWordToDecideRing=$6
    local keyWordThreshold=$7
    local keyWordValueForFirstRing=$8
    local keyWordValueForSecondRing=$9
    local h=${10}

    if [ "$useCommonRing" = true ]; then
            # Case when we have one common normalisation ring
            #astarithmetic $i -h1 $commonRing -h1 0 eq nan where -q -o ring_masked.fits
            skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS $commonRing $h)
            kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS $commonRing $h)
            #rm ring_masked.fits
    else
        # Case when we do NOT have one common normalisation ring
        # All the following logic is to decide which normalisation ring apply
        variableToDecideRingToNormalise=$(gethead $i $keyWordToDecideRing)
        firstRingLowerBound=$(echo "$keyWordValueForFirstRing - $keyWordThreshold" | bc)
        firstRingUpperBound=$(echo "$keyWordValueForFirstRing + $keyWordThreshold" | bc)
        secondRingLowerBound=$(echo "$keyWordValueForSecondRing - $keyWordThreshold" | bc)
        secondRingUpperBound=$(echo "$keyWordValueForSecondRing + $keyWordThreshold" | bc)

        if (( $(echo "$variableToDecideRingToNormalise >= $firstRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $firstRingUpperBound" | bc -l) )); then
            #astarithmetic $i -h1 $doubleRing_first -h1 0 eq nan where -q -o ring_masked.fits
            skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS $doubleRing_first $h)
            kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS $doubleRing_first $h)
            #rm ring_masked.fits
        elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
            #astarithmetic $i -h1 $doubleRing_second -h1 0 eq nan where -q -o ring_masked.fits
            skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS $doubleRing_second $h)
            kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS $doubleRing_second $h)
            #rm ring_masked.fits
        else
            errorNumber=5
            echo -e "\nMultiple normalisation ring have been tried to be used. The keyword selection value of one has not matched with the ranges provided" >&2
            echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
            exit $errorNumber 
        fi
    fi
    echo "$skew $kurto"
}
export -f getSkewKurtoValueInsideRing

normaliseImagesWithRing() {
    local imageDir=$1
    local outputDir=$2
    local useCommonRing=$3
    local commonRing=$4
    local doubleRing_first=$5
    local doubleRing_second=$6
    local keyWordToDecideRing=$7
    local keyWordThreshold=$8
    local keyWordValueForFirstRing=$9
    local keyWordValueForSecondRing=${10}
    imagesToNormalise=()
    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
        imagesToNormalise+=("$base")
    done
    printf "%s\n" "${imagesToNormalise[@]}" | parallel -j "$num_cpus" normaliseIndividualImageWithRing {} $imageDir $outputDir $useCommonRing $commonRing $doubleRing_first $doubleRing_second $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing
}
export -f normaliseImagesWithRing
normaliseIndividualImageWithRing() {
    local base=$1
    local imageDir=$2
    local outputDir=$3
    local useCommonRing=$4
    local commonRing=$5
    local doubleRing_first=$6
    local doubleRing_second=$7
    local keyWordToDecideRing=$8
    local keyWordThreshold=$9
    local keyWordValueForFirstRing=${10}
    local keyWordValueForSecondRing=${11}
    i=$imageDir/$base
    out=$outputDir/$base
    astfits $i --copy=0 --primaryimghdu -o $out
    if [ "$USE_COMMON_RING" = false ]; then
        propagateKeyword $i $keyWordToDecideRing $out 0
    fi
    for h in $(seq 1 $num_ccd); do

        me=$(getMedianValueInsideRing $i $commonRing $doubleRing_first $doubleRing_second $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $h)
        astarithmetic $i -h$h $me / -o $outputDir/temp_$base
        astfits $outputDir/temp_$base --copy=1 -o $out
        rm $outputDir/temp_$base
        propagateKeyword $i $gain $out $h
    done
}
export -f normaliseIndividualImageWithRing


calculateFlat() {
    local flatName="$1"
    local flatIteration="$2"
    shift 2
    local filesToUse="$@"
    numberOfFiles=$#
    
    # ****** Decision note *******
    # The rejection parameters for the construction of the flat has been chosen to be 2 sigmas
    # The running flat implies that we only have fewer frames for the flat (in our case 11 for example)
    # So we have to be a little bit aggresive in order to be able to remove the outliers
    sigmaValue=2
    iterations=10
    gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
    first_file=$(echo "$filesToUse" | awk '{print $1}')
    astfits $first_file --copy=0 --primaryimghdu -o $flatName
    for h in $(seq 1 $num_ccd); do
        if [ "$flatIteration" -ne "3" ]; then
            if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
                astarithmetic $filesToUse $numberOfFiles $sigmaValue $iterations sigclip-median -g$h --writeall -o "${flatName%.fits}_temp.fits"
            else
                astarithmetic $filesToUse $numberOfFiles $sigmaValue $iterations sigclip-median -g$h -o "${flatName%.fits}_temp.fits"
            fi
        else
            badFilesWarningFile=$BDIR/diagnosis_and_badFiles/CCD"$h"/identifiedBadFrames_preFlat_onlyStd.txt
            declare -A badFramesMap
            unset badFramesMap
            declare -A badFramesMap
            if [ -s $badFilesWarningFile ]; then    
                while IFS= read -r line; do
                    baseName=$(basename "$line" .txt)
                    badFramesMap["$baseName"]=1
                done < $badFilesWarningFile
            fi
            filesToUseOk=()
            for file in $filesToUse; do
                baseNameToCheck=$(basename "$file" .fits)
                if [[ -z "${badFramesMap[$baseNameToCheck]}" ]]; then
                    filesToUseOk+=("$file")
                fi
            done
            if [ ${#filesToUseOk[@]} -eq 0 ]; then 
                filesToUseOk=($filesToUse)
            fi
            numberOfFilesOk=${#filesToUseOk[@]}
            if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
                astarithmetic "${filesToUseOk[@]}" $numberOfFilesOk $sigmaValue $iterations sigclip-median -g$h --writeall -o "${flatName%.fits}_temp.fits"
            else
                astarithmetic "${filesToUseOk[@]}" $numberOfFilesOk $sigmaValue $iterations sigclip-median -g$h -o "${flatName%.fits}_temp.fits"
            fi
        fi
        astfits "${flatName%.fits}_temp.fits" --copy=1 -o $flatName
        rm "${flatName%.fits}_temp.fits"
    done
}
export -f calculateFlat

calculateRunningFlat() {
    local normalisedDir=$1
    local outputDir=$2
    local doneFile=$3
    local iteration=$4
    windowSize=$(( (halfWindowSize * 2) + 1 ))
    fileArray=()
    fileArray=( $(ls -v $normalisedDir/*Decals-"$filter"_n*_f*.fits) )
    fileArrayLength=( $(ls -v $normalisedDir/*Decals-"$filter"_n*_f*.fits | wc -l) )

    lefFlatFiles=("${fileArray[@]:0:$windowSize}")
    
    echo "Computing left flat - iteration $iteration"
    calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_left.fits" "$iteration" "${lefFlatFiles[@]}"
    
    rightFlatFiles=("${fileArray[@]:(fileArrayLength-$windowSize):fileArrayLength}")
    echo "Computing right flat - iteration $iteration"
    calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_right.fits" "$iteration" "${rightFlatFiles[@]}"

    echo "Computing non-common flats - iteration $iteration"
    for a in $(seq 1 $n_exp); do
        if [ "$a" -gt "$((halfWindowSize + 1))" ] && [ "$((a))" -lt "$(($n_exp - $halfWindowSize))" ]; then
            leftLimit=$(( a - $halfWindowSize - 1))
            calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_f"$a".fits" "$iteration" "${fileArray[@]:$leftLimit:$windowSize}"
        fi
    done
    echo done > $doneFile
}
export -f calculateRunningFlat

divideImagesByRunningFlats(){
    local imageDir=$1
    local outputDir=$2
    local flatDir=$3
    local flatDone=$4
    imagesToDivide=()
    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
        imagesToDivide+=("$base")
    done
    printf "%s\n" "${imagesToDivide[@]}" | parallel -j "$num_cpus" divideIndividualImageByRunningFlats {} $imageDir $outputDir $flatDir $n_exp $currentNight $iteration 
    echo done > $flatDone
}
export -f divideImagesByRunningFlats
divideIndividualImageByRunningFlats(){
    local base=$1
    local imageDir=$2
    local outputDir=$3
    local flatDir=$4
    local n_exp=$5
    local currentNight=$6
    local iteration=$7

    i=$imageDir/$base
    out=$outputDir/$base
    a="${base#*_f}"
    a="${a%.fits}"
    astfits $i --copy=0 --primaryimghdu -o $out
    if [ "$a" -le "$((halfWindowSize + 1))" ]; then
        flatToUse=$flatDir/flat-it*_"$filter"_n"$currentNight"_left.fits
    elif [ "$a" -ge "$((n_exp - halfWindowSize))" ]; then
        flatToUse=$flatDir/flat-it"$iteration"_"$filter"_n"$currentNight"_right.fits
    else
        flatToUse=$flatDir/flat-it"$iteration"_"$filter"_n"$currentNight"_f"$a".fits
    fi
    for h in $(seq 1 $num_ccd); do
        astarithmetic $i -h$h $flatToUse -h$h / -o $outputDir/temp_$base
        astfits $outputDir/temp_$base --copy=1 -o$out
        rm $outputDir/temp_$base
        propagateKeyword $i $gain $out $h
    done
}
export -f divideIndividualImageByRunningFlats
    

divideImagesByWholeNightFlat(){
    local imageDir=$1
    local outputDir=$2
    local flatToUse=$3
    local flatDone=$4
    imagesToDivide=()
    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
        imagesToDivide+=("$base")
    done
    printf "%s\n" "${imagesToDivide[@]}" | parallel -j "$num_cpus" divideIndividualImageByWholeNightFlat {} $imageDir $outputDir $flatToUse
    echo done > $flatDone
}
export -f divideImagesByWholeNightFlat
divideIndividualImageByWholeNightFlat(){
    local base=$1
    local imageDir=$2
    local outputDir=$3
    local flatToUse=$4
    
    i=$imageDir/$base
    out=$outputDir/$base
    astfits $i --copy=0 --primaryimghdu -o $out
    for h in $(seq 1 $num_ccd); do
        astarithmetic $i -h$h $flatToUse -h$h / -o $outputDir/temp_$base
        astfits $outputDir/temp_$base --copy=1 -o$out
        rm $outputDir/temp_$base
        propagateKeyword $i $gain $out $h
    done
}
export -f divideIndividualImageByWholeNightFlat

correctRunningFlatWithWholeNightFlat() {
    local base=$1
    local beforeDir=$2
    local wholeNightFlat=$3 
    local outputDir=$4

    i=$beforeDir/$base
    out=$outputDir/$base
    astfits $i --copy=0 --primaryimghdu -o $out
    for h in $(seq 1 $num_ccd); do
        tmpRatio=$outputDir/tmpRatio_$base
        astarithmetic $wholeNightFlat -h$h $i -h$h / -o $tmpRatio
        tmpCorrected=$outputDir/tmpCorrected_$base
        astarithmetic $i -h$h set-m $tmpRatio -h1 set-f m f 0.85 lt nan where -o $tmpCorrected
        astfits $tmpCorrected --copy=1 -o $out
        rm $tmpRatio $tmpCorrected
    done
}    
export -f correctRunningFlatWithWholeNightFlat

maskVignettingOnImages() {
    local base=$1
    local imaDir=$2 
    local outDir=$3
    local flatDir=$4
    local wholeFlatDir=$5
    local runningFlat=$6
    local n_exp=$7
    local currentNight=$8
    local lowerVignettingThreshold=$9
    local upperVignettingThreshold=${10}
    a="${base#*_f}"
    a="${a%.fits}"
    if $runningFlat; then
      if [ "$a" -le "$((halfWindowSize + 1))" ]; then
        currentFlatImage=$flatDir/flat-it3_"$filter"_n"$currentNight"_left.fits
      elif [ "$a" -ge "$((n_exp - halfWindowSize))" ]; then
        currentFlatImage=$flatDir/flat-it3_"$filter"_n"$currentNight"_right.fits
      else
        currentFlatImage=$flatDir/flat-it3_"$filter"_n"$currentNight"_f"$a".fits
      fi
    else
      currentFlatImage=$wholeFlatDir/flat-it3_wholeNight_n$currentNight.fits
    fi 

    i=$imaDir/$base
    out=$outDir/$base
    tempStep=$outDir/temp_$base
    astfits $i --copy=0 --primaryimghdu -o $out
    for h in $(seq 1 $num_ccd); do
        astarithmetic $i -h$h set-m $currentFlatImage -h$h set-f m f $lowerVignettingThreshold lt  nan where set-n n f $upperVignettingThreshold gt nan where -o $tempStep
        astfits $tempStep --copy=1 -o $out
        rm -f $tempStep
    done
    propagateKeyword $i $gain $out 0
}
export -f maskVignettingOnImages

runNoiseChiselOnFrame() {
    local baseName=$1
    local inputFileDir=$2
    local outputDir=$3
    local blockScale=$4
    local noiseChiselParams=$5
    
    ###If a block scale is given, we will block, highlighting LSB regions, detect, and un-block the mask
    
    imageToUse=$inputFileDir/$baseName
    output=$outputDir/$baseName
    
    for h in $(seq 1 $num_ccd); do
        if [ "$blockScale" -eq 1 ]; then
            astnoisechisel $imageToUse -h$h $noiseChiselParams --numthreads=$num_threads -o $outputDir/temp_"$baseName"
        else
            wFile=$outputDir/imW_$baseName
            wMask=$outputDir/mkW_$baseName
            wMask2=$outputDir/mkW2_$baseName
            astwarp $imageToUse -h$h --scale=1/$blockScale --numthreads=$num_threads -o $wFile
            astnoisechisel $wFile -h1 $noiseChiselParams --numthreads=$num_threads -o $wMask
            astwarp $wMask -h1 --gridfile=$imageToUse --gridhdu=$h --numthreads=$num_threads -o $wMask2
            astarithmetic $wMask2 -h1 set-i i i 0 gt 1 where -q float32 -o $outputDir/temp_"$baseName"
            rm $wFile $wMask $wMask2
        fi
        astfits $outputDir/temp_"$baseName" --copy=1 -o $output
        rm $outputDir/temp_"$baseName"
    done
}
export -f runNoiseChiselOnFrame

# Functions for Warping the frames
getCentralCoordinate(){
    local image=$1
    local hdu=$2

    NAXIS1=$(fitsheader $image -e $hdu | grep "NAXIS1" | awk '{print $3'})
    NAXIS2=$(fitsheader $image -e $hdu | grep "NAXIS2" | awk '{print $3'})
                   

    # Calculate the center pixel coordinates
    center_x=$((NAXIS1 / 2))
    center_y=$((NAXIS2 / 2))

    # Use xy2sky to get the celestial coordinates of the center pixel
    imageCentre=$( xy2sky $image,$hdu $center_x $center_y )
    echo $imageCentre
}
export -f getCentralCoordinate

warpImage() {
    local imageToSwarp=$1
    local entiredir=$2
    local ra=$3
    local dec=$4
    local coaddSizePx=$5

    # ****** Decision note *******
    # We need to regrid the frames into the final coadd grid. But if we do this right now we will be processing
    # frames too big (which are mostly Nans) and the noisechisel routine takes a looot of time
    # The approach taken is to move the frame to that grid, and then crop it to the dimension of the data itself
    # We need to store both. I have tried to store the small one and then warp it again to the big grid, it's more time consuming
    # and the nan wholes grow so we end up with less light in the final coadd.

    # Parameters for identifing our frame in the full grid
    currentIndex=$(basename $imageToSwarp .fits)

    
    

    ##Multiple layers treatement: we want to preserve the multi-layer structure of the .fits file in the output, something swarp apparently doesn't like. 
    #The idea here will be to create a new folder called astro-ima-single with the .fits and the .head broken into each ccd, in order to run swarp
    
    #if ! [ -d $tmpDir ]; then mkdir $tmpDir; fi
    #for h in $(seq 1 $num_ccd); do
    #    astfits $imageToSwarp --copy=$h -o $tmpDir/"$currentIndex"_ccd"$h".fits
    #done
   
    #h=1
    #header=$astroimadir/"$currentIndex".head
    #awk -v h="$h" -v a="$currentIndex" -v dir="$tmpDir" '
    #/^HISTORY/ {
    #    if (h > 1) close(output);
    #    output = dir "/" a "_ccd" h ".head";
    #    h++;
    #}
    #{ print > output }
    #' "$header.head"
    

    # Resample into the final grid
    # Be careful with how do you have to call this package, because in the SIE sofware is "SWarp" and in the TST-ICR is "swarp"
    #for h in $(seq 1 $num_ccd); do

    #frameFullGrid=$entiredir/entirecamera_"$currentIndex"_fullGrid.fits
    frameSmallGrid=$entiredir/entirecamera_$currentIndex.fits
    # Resample into the final grid
    # Be careful with how do you have to call this package, because in the SIE sofware is "SWarp" and in the TST-ICR is "swarp"
    detect_swarp() {
        for cmd in swarp SWarp; do
            if command -v "$cmd" >/dev/null 2>&1; then
                echo "$cmd"
                return
            fi
        done
        echo "Error: SWarp not found" >&2
        exit 1
    }

    SWARP_CMD=$(detect_swarp)
    $SWARP_CMD -c $swarpcfg $imageToSwarp -CENTER $ra,$dec -IMAGE_SIZE $coaddSizePx,$coaddSizePx -IMAGEOUT_NAME $entiredir/"$currentIndex"_swarp1.fits -WEIGHTOUT_NAME $entiredir/"$currentIndex"_swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $pixelScale -PIXELSCALE_TYPE MANUAL -DELETE_TMPFILES N
    
    # Mask bad pixels

    ##Temporary files are of the tipe $currentIndex.000$h.resamp(.weight).fits
    for h in $(seq 1 $num_ccd); do
        tmpFile1=$entiredir/"$currentIndex"_ccd"$h"_temp.fits
        tmpFile2=$entiredir/"$currentIndex"_ccd"$h"_temp2.fits
        tmpFile3=$entiredir/"$currentIndex"_ccd"$h"_temp3.fits
        
        astarithmetic "$currentIndex".000"$h".resamp.weight.fits -h0 set-i i i 0 lt nan where -o $tmpFile1
        astarithmetic "$currentIndex".000"$h".resamp.fits -h0 $tmpFile1 -h1 0 eq nan where -o $tmpFile2
        #astfits $tmpFile2 --copy=1 -o $frameSmallGrid #Would be easier but will make the maskPointings of 2nd iteration coadd not work
        astcrop $tmpFile2 --mode=wcs --center=$ra,$dec --widthinpix --width=$coaddSizePx,$coaddSizePx --zeroisnotblank -o $tmpFile3
        #astfits $tmpFile3 --copy=1 -o $frameFullGrid
        #Now we use tmpFIle3 to generate the frameSmallGrid
        regionOfDataInFullGrid=$(python3 $pythonScriptsPath/getRegionToCrop.py $tmpFile3 1)
        read row_min row_max col_min col_max <<< "$regionOfDataInFullGrid"
        echo $row_min $row_max $col_min $col_max >> $entiredir/entirecamera_"$currentIndex"_cropRegion.txt
        tmpFile4=$entiredir/entirecamera_sg_"$currentIndex"_ccd"$h".fits 
        astcrop $tmpFile3 --polygon=$col_min,$row_min:$col_max,$row_min:$col_max,$row_max:$col_min,$row_max --mode=img  -o $tmpFile4 --quiet
        astfits $tmpFile4 --copy=1 -o $frameSmallGrid
        rm $tmpFile1 $tmpFile2 $tmpFile3  "$currentIndex".000"$h"*.fits $tmpFile4
        #propagateKeyword $imageToSwarp $gain $frameFullGrid $h 
        propagateKeyword $imageToSwarp $gain $frameSmallGrid $h
    done
    #propagateKeyword $imageToSwarp $airMassKeyWord $frameFullGrid 0 
    propagateKeyword $imageToSwarp $airMassKeyWord $frameSmallGrid 0
    #propagateKeyword $imageToSwarp $dateHeaderKey $frameFullGrid 0 
    propagateKeyword $imageToSwarp $dateHeaderKey $frameSmallGrid 0
    #Swarp temporary files are named as
        # Mask bad pixels
        #tmpFile1=$tmpDir"/$currentIndex"_temp1_ccd"$h".fits
        #astarithmetic $tmpDir/"$currentIndex"_swarp_w1_ccd"$h".fits -h0 set-i i i 0 lt nan where -o$tmpFile1
        #astarithmetic $tmpDir/"$currentIndex"_swarp1_ccd"$h".fits -h0 $tmpFile1 -h1 0 eq nan where -o$tmpDir/entirecamera_fg_"$currentIndex"_ccd"$h".fits

        #regionOfDataInFullGrid=$(python3 $pythonScriptsPath/getRegionToCrop.py $tmpDir/entirecamera_fg_"$currentIndex"_ccd"$h".fits 1)
        #read row_min row_max col_min col_max <<< "$regionOfDataInFullGrid"
        #astcrop $tmpDir/entirecamera_fg_"$currentIndex"_ccd"$h".fits --polygon=$col_min,$row_min:$col_max,$row_min:$col_max,$row_max:$col_min,$row_max --mode=img  -o $tmpDir/entirecamera_sg_"$currentIndex"_ccd"$h".fits --quiet
        #astfits $tmpDir/entirecamera_fg_"$currentIndex"_ccd"$h".fits --copy=1 -o$frameFullGrid
        #astfits $tmpDir/entirecamera_sg_"$currentIndex"_ccd"$h".fits --copy=1 -o$entiredir/entirecamera_"$currentIndex".fits
        #propagateKeyword $imageToSwarp $gain $frameFullGrid $h 
        #propagateKeyword $imageToSwarp $gain $entiredir/entirecamera_"$currentIndex".fits $h
    #done
    
    rm -rf $entiredir/"$currentIndex"_swarp*.fits
}
export -f warpImage

removeBadFramesFromReduction() {
    local sourceToRemoveFiles=$1
    local destinationDir=$2
    local badFilesWarningDir=$3
    local badFilesWarningFile=$4

    filePath=$badFilesWarningDir/$badFilesWarningFile


    while IFS= read -r file_name; do
        file_name=$(basename "$file_name")
        fileName=$prefixOfFilesToRemove"${file_name%.*}".fits
        if [ -f $sourceToRemoveFiles/$fileName ]; then
            mv $sourceToRemoveFiles/$fileName $destinationDir/$fileName
        fi

        txtFileName=$prefixOfFilesToRemove"${file_name%.*}".txt
        if [ -f $sourceToRemoveFiles/$txtFileName ]; then
            mv $sourceToRemoveFiles/$txtFileName $destinationDir/$txtFileName
        fi

        maskFileName=$prefixOfFilesToRemove"${file_name%.*}"_masked.fits
        if [ -f $sourceToRemoveFiles/$maskFileName ]; then
            mv $sourceToRemoveFiles/$maskFileName $destinationDir/$maskFileName
        fi

    done < "$filePath"
}
export -f removeBadFramesFromReduction

# Functions for compute and subtract sky from frames
computeSkyForFrame(){
    local base=$1
    local entiredir=$2
    local noiseskydir=$3
    local constantSky=$4
    local constantSkyMethod=$5
    local polyDegree=$6
    local inputImagesAreMasked=$7
    local ringDir=$8
    local useCommonRing=$9
    local keyWordToDecideRing=${10}
    local keyWordThreshold=${11}
    local keyWordValueForFirstRing=${12}
    local keyWordValueForSecondRing=${13}
    local ringWidth=${14}
    local swarped=${15}
    local blockScale=${16}
    local noisechisel_param=${17}
    local maskParams=${18}
    i=$entiredir/$1

    # ****** Decision note *******
    # Here we have implemented two possibilities. Either the background is estimated by a constant or by a polynomial.
    # If it is a constant we store the fileName, the background value and the std. This is implemented this way because we
    # need the background to subtract but also later the std for weighing the frames
    # If it is a polynomial we only use it to subtract the background (the weighing is always with a constat) so we only store
    # the coefficients of the polynomial.
    # 
    # Storing this values is also relevant for checking for potential bad frames
    out=$(echo $base | sed 's/.fits/.txt/')

    if [ "$constantSky" = true ]; then  # Case when we subtract a constant
        # Here we have two possibilities
        # Estimate the background within the normalisation ring or using noisechisel

        # The problem is that I can't use the same ring/s as in the normalisation because here we have warped and cropped the images... So I create a new normalisation ring from the centre of the images
        # I cannot even create a common ring for all, because they are cropped based on the number of non-nan (depending on the vignetting and how the NAN are distributed), so i create a ring per image
        # For that reason the subtraction of the background using the ring is always using a ring centered in the frame
        # More logic should be implemented to use the normalisation ring(s) and recover them after the warping and cropping
        if [ "$constantSkyMethod" = "ring" ]; then

            # Mask the image if they are not already masked
            if ! [ "$inputImagesAreMasked" = true ]; then
                tmpMask=$(echo $base | sed 's/.fits/_mask.fits/')
                tmpMaskedImage=$(echo $base | sed 's/.fits/_masked.fits/')
                tmpMaskedImage_single=$(echo $base | sed 's/.fits/_masked_ccd.fits/')
                runNoiseChiselOnFrame $1 $entiredir $noiseskydir $blockScale "$noisechisel_param"
                
                mv $noiseskydir/$1 $noiseskydir/$tmpMask
                for h in $(seq 1 $num_ccd); do
                    astarithmetic $i -h$h $noiseskydir/$tmpMask -h$h 1 eq nan where float32 -o $noiseskydir/$tmpMaskedImage_single -q
                    astfits $noiseskydir/$tmpMaskedImage_single --copy=1 -o$noiseskydir/$tmpMaskedImage
                    rm -f $noiseskydir/$tmpMaskedImage_single
                    
                done    
                rm -f $noiseskydir/$tmpMask
                imageToUse=$noiseskydir/$tmpMaskedImage
                
                ##Aply manual mask defined by user
                valueToPut=nan
                read -r -a maskArray <<< "$maskParams"
                for ((i=0; i<${#maskArray[@]}; i+=5)); do
                    ra="${maskArray[i]}"
                    dec="${maskArray[i+1]}"
                    r="${maskArray[i+2]}"
                    axisRatio="${maskArray[i+3]}"
                    pa="${maskArray[i+4]}"
                    python3 $pythonScriptsPath/manualMaskRegionFromWCSArea.py $imageToUse $valueToPut $ra $dec $r $axisRatio $pa
                done 
            else
                imageToUse=$i
            fi

            # We generate the ring (we cannot use the normalisation ring because we have warped and cropped) and compute the background value within it
            #tmpRingDefinition=$(echo $base | sed 's/.fits/_ring.txt/')
            tmpRingFits=$(echo $base | sed 's/.fits/_ring.fits/')
            tmpRingFits_single=$(echo $base | sed 's/.fits/_ring_single.fits/')
            ##We get the reference from the first non rotated ccd
            
            #half_naxis1=$(echo "$naxis1 / 2" | bc)
            #half_naxis2=$(echo "$naxis2 / 2" | bc)

            #ringRadius=$( awk '{print $5}' $ringDir/ring.txt )
            #echo "1 $half_naxis1 $half_naxis2 6 $ringRadius 1 1 1 1 1" > $ringDir/$tmpRingDefinition
            for h in $(seq 1 $num_ccd); do
                if [ "$swarped" = "YES" ]; then
                
                    x_ring=$( awk ' {print $2}' $ringDir/ring_ccd"$h".txt )
                    y_ring=$( awk ' {print $3}' $ringDir/ring_ccd"$h".txt )
                    tmpRingDefinition=$(echo $base | sed 's/.fits/_ring_ccd.txt/')
                    ringRadius=$( awk '{print $5}' $ringDir/ring_ccd"$h".txt )
                    ##We first check if a rotation is needed by comparing naxis of the ring (which is already created) and naxis of the frame
                    naxis1=$(fitsheader $imageToUse -e $h | grep "NAXIS1" | awk '{print $3'})
                    naxis2=$(fitsheader $imageToUse -e $h | grep "NAXIS2" | awk '{print $3'})
                    naxis1_r=$(fitsheader $ringDir/ring.fits -e $h | grep "NAXIS1" | awk '{print $3'})
                    naxis2_r=$(fitsheader $ringDir/ring.fits -e $h | grep "NAXIS2" | awk '{print $3'})

                    #If the axis on ring and on image keeps the comparison, we don't need to do anything
                    if [[ $naxis1 -gt $naxis2 && $naxis1_r -gt $naxis2_r ]] || [[ $naxis1 -lt $naxis2 && $naxis1_r -lt $naxis2_r ]]; then
                        echo "1 $x_ring $y_ring 6 $ringRadius 1 1 1 1 1" > $ringDir/$tmpRingDefinition
                    else
                        #If image new is rotated, we look for the astrometrized, not rotated in order to get the correct position
                        image_astro=${base#entirecamera_}
                        ringCentre=$( xy2sky $BDIR/astro-ima/$image_astro,$h $x_ring $y_ring )
                        ringRa=$(echo "$ringCentre" | awk '{print $1}')
                        ringDec=$(echo "$ringCentre" | awk '{print $2}')
                        newringCentre=$( sky2xy $imageToUse,$h $ringRa $ringDec )
                        x_new=$(echo "$newringCentre" | awk '{print $5}')
                        y_new=$(echo "$newringCentre" | awk '{print $6}')
                        echo "1 $x_new $y_new 6 $ringRadius 1 1 1 1 1" > $ringDir/$tmpRingDefinition
                   
                        
                    fi
             
                    astmkprof --background=$imageToUse --backhdu=$h --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas --quiet -o $ringDir/$tmpRingFits_single $ringDir/$tmpRingDefinition
                    rm -f $ringDir/$tmpRingDefinition
                else
                    astmkprof --background=$imageToUse  --backhdu=$h --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas --quiet -o $ringDir/$tmpRingFits_single $ringDir/ring_ccd"$h".txt
                fi
                astfits $ringDir/$tmpRingFits_single --copy=1 -o $ringDir/$tmpRingFits
                rm -f $ringDir/$tmpRingFits_single
               
            
            done
            #if [ "$swarped" = "YES" ]; then rm $ringDir/$tmpRingDefinition; fi
            #Since getMedianValueInsideRing is modified to treat with a multiple layer ring we are forced to split the loops :(
            for h in $(seq 1 $num_ccd); do
                me=$(getMedianValueInsideRing $imageToUse  $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $h)
                std=$(getStdValueInsideRing $imageToUse $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $h)
                read skew kurto < <(getSkewKurtoValueInsideRing $imageToUse $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $h)
                #This is introduced in order to make a jump line unless file is not written
                echo "$base $me $std $skew $kurto" >> $noiseskydir/$out
            done
            #rm $ringDir/$tmpRingDefinition
            
            rm $ringDir/$tmpRingFits

        elif [ "$constantSkyMethod" = "noisechisel" ]; then
            for h in $(seq 1 $num_ccd); do
                sky=$(echo $base | sed 's/.fits/_sky.fits/')

                # The sky substraction is done by using the --checksky option in noisechisel
                astnoisechisel $i -h$h --tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2 --checksky $noisechisel_param --numthreads=$num_cpus -o $noiseskydir/$base
                mean=$(aststatistics $noiseskydir/$sky -hSKY --sigclip-mean)
                std=$(aststatistics $noiseskydir/$sky -hSTD --sigclip-mean)
                skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $noiseskydir/$sky SKEWNESS NO $h)
                kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $noiseskydir/$sky KURTOSIS NO $h)
                echo "$base $mean $std $skew $kurto" >> $noiseskydir/$out
                rm -f $noiseskydir/$sky
            done
        elif [ "$constantSkyMethod" = "wholeImage" ]; then
            for h in $(seq 1 $num_ccd); do
                mean=$(aststatistics $i -h$h --sigclip-mean -q)
                std=$(aststatistics $i -h$h --sigclip-std -q)
                skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS NO $h)
                kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS NO $h)
                echo "$base $mean $std $skew $kurto" >> $noiseskydir/$out
            done
        elif [ "$constantSkyMethod" = "ringAndDetector" ]; then
            tmpRingFits=$(echo $base | sed 's/.fits/_ring.fits/')
            tmpRingFits_single=$(echo $base | sed 's/.fits/_ring_single.fits/')
            imageToUse=$i
            if [ "$swarped" = "YES" ]; then
                for h in $(seq 1 $num_ccd); do
                    x_ring=$( awk ' {print $2}' $ringDir/ring_ccd"$h".txt )
                    y_ring=$( awk ' {print $3}' $ringDir/ring_ccd"$h".txt )
                    tmpRingDefinition=$(echo $base | sed 's/.fits/_ring_ccd.txt/')
                    ringRadius=$( awk '{print $5}' $ringDir/ring_ccd"$h".txt )
                    ##We first check if a rotation is needed by comparing naxis of the ring (which is already created) and naxis of the frame
                    naxis1=$(fitsheader $imageToUse -e $h | grep "NAXIS1" | awk '{print $3'})
                    naxis2=$(fitsheader $imageToUse -e $h | grep "NAXIS2" | awk '{print $3'})
                    naxis1_r=$(fitsheader $ringDir/ring.fits -e $h | grep "NAXIS1" | awk '{print $3'})
                    naxis2_r=$(fitsheader $ringDir/ring.fits -e $h | grep "NAXIS2" | awk '{print $3'})

                    #If the axis on ring and on image keeps the comparison, we don't need to do anything
                    if [[ $naxis1 -gt $naxis2 && $naxis1_r -gt $naxis2_r ]] || [[ $naxis1 -lt $naxis2 && $naxis1_r -lt $naxis2_r ]]; then
                        echo "1 $x_ring $y_ring 6 $ringRadius 1 1 1 1 1" > $ringDir/$tmpRingDefinition
                    else
                        #If image new is rotated, we look for the astrometrized, not rotated in order to get the correct position
                        image_astro=${base#entirecamera_}
                        ringCentre=$( xy2sky $BDIR/astro-ima/$image_astro,$h $x_ring $y_ring )
                        ringRa=$(echo "$ringCentre" | awk '{print $1}')
                        ringDec=$(echo "$ringCentre" | awk '{print $2}')
                        newringCentre=$( sky2xy $imageToUse,$h $ringRa $ringDec )
                        x_new=$(echo "$newringCentre" | awk '{print $5}')
                        y_new=$(echo "$newringCentre" | awk '{print $6}')
                        echo "1 $x_new $y_new 6 $ringRadius 1 1 1 1 1" > $ringDir/$tmpRingDefinition
                   
                        
                    fi
             
                    astmkprof --background=$imageToUse --backhdu=$h --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas --quiet -o $ringDir/$tmpRingFits_single $ringDir/$tmpRingDefinition
                    rm -f $ringDir/$tmpRingDefinition
                    astfits $ringDir/$tmpRingFits_single --copy=1 -o $ringDir/$tmpRingFits
                    rm -f $ringDir/$tmpRingFits_single
                done
            else
                for h in $(seq 1 $num_ccd); do
                    astmkprof --background=$imageToUse  --backhdu=$h --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas --quiet -o $ringDir/$tmpRingFits_single $ringDir/ring_ccd"$h".txt
                    astfits $ringDir/$tmpRingFits_single --copy=1 -o $ringDir/$tmpRingFits
                    rm -f $ringDir/$tmpRingFits_single
                done
            fi
            
            python3 $pythonScriptsPath/getSkySTDSkewKurtosis_fullDetector.py $imageToUse $ringDir/$tmpRingFits $noiseskydir $num_ccd
            rm $ringDir/$tmpRingFits

        else
            errorNumber=6
            echo -e "\nAn invalid value for the sky_estimation_method was provided" >&2
            echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
            exit $errorNumber
        fi

    else
        echo "\n\tMultiDetector pipeline has not been prepared yet to work with polynomial fitting of sky"
        exit 35
        # Case when we model a plane
        noiseOutTmp=$(echo $base | sed 's/.fits/_sky.fits/')
        maskTmp=$(echo $base | sed 's/.fits/_masked.fits/')
        planeOutput=$(echo $base | sed 's/.fits/_poly.fits/')
        planeCoeffFile=$(echo $base | sed 's/.fits/.txt/')

        # This conditional allows us to introduce the images already masked (masked with the mask of the coadd) in the second and next iterations
        if ! [ "$inputImagesAreMasked" = true ]; then
            astnoisechisel $i --tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2 --checksky $noisechisel_param --numthreads=$num_cpus -o $noiseskydir/$base
            astarithmetic $i -h1 $noiseskydir/$noiseOutTmp -hDETECTED 1 eq nan where -q float32 -o $noiseskydir/$maskTmp
            python3 $pythonScriptsPath/surface-fit.py -i $noiseskydir/$maskTmp -o $noiseskydir/$planeOutput -d $polyDegree -f $noiseskydir/$planeCoeffFile
        else
            python3 $pythonScriptsPath/surface-fit.py -i $i -o $noiseskydir/$planeOutput -d $polyDegree -f $noiseskydir/$planeCoeffFile
        fi

        rm -f $noiseskydir/$noiseOutTmp
        rm -f $noiseskydir/$maskTmp
    fi
}
export -f computeSkyForFrame

computeSky() {
    local framesToUseDir=$1
    local noiseskydir=$2
    local noiseskydone=$3
    local constantSky=$4
    local constantSkyMethod=$5
    local polyDegree=$6
    local inputImagesAreMasked=$7
    local ringDir=$8
    local useCommonRing=$9
    local keyWordToDecideRing=${10}
    local keyWordThreshold=${11}
    local keyWordValueForFirstRing=${12}
    local keyWordValueForSecondRing=${13}
    local ringWidth=${14}
    local swarped=${15}
    local blockScale=${16}
    local noisechisel_param=${17}
    local maskParams=${18}
    
    if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
    if [ -f $noiseskydone ]; then
        echo -e "\n\tScience images are 'noisechiseled' for constant sky substraction\n"
    else
        framesToComputeSky=()
        for a in $( ls $framesToUseDir/*.fits ); do
            base=$( basename $a )
            framesToComputeSky+=("$base")
        done
        
        printf "%s\n" "${framesToComputeSky[@]}" | parallel -j "$num_cpus" computeSkyForFrame {} $framesToUseDir $noiseskydir $constantSky $constantSkyMethod $polyDegree $inputImagesAreMasked $ringDir $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth $swarped $blockScale $noisechisel_param $maskParams
        echo done > $noiseskydone
    fi
}
export -f computeSky

subtractSkyForFrame() {
    local a=$1
    local directoryWithSkyValues=$2
    local framesToSubtract=$3
    local directoryToStoreSkySubtracted=$4
    local constantSky=$5

    base="entirecamera_"$a.fits
    input=$framesToSubtract/$base
    output=$directoryToStoreSkySubtracted/$base

    if [ "$constantSky" = true ]; then
        i=$directoryWithSkyValues/"entirecamera_"$a.txt
        astfits $input --copy=0 --primaryimghdu -o $output
        for h in $(seq 1 $num_ccd); do
            temp_file=$directoryToStoreSkySubtracted/temp_$base
            me=$(awk 'NR=='$h'{print $2}' $i)
            astarithmetic $input -h$h $me - -o$temp_file
            astfits $temp_file --copy=1 -o $output
            rm $temp_file
        done
    else
        i=$directoryWithSkyValues/"entirecamera_"$a"_poly.fits"
        for h in $(seq 1 $num_ccd); do
            temp_file=$directoryToStoreSkySubtracted/temp_$base
            NAXIS1_image=$(gethead $input -x $h NAXIS1); NAXIS2_image=$(gethead $input -x $h NAXIS2)
            NAXIS1_plane=$(gethead $i -x $h NAXIS1); NAXIS2_plane=$(gethead $i -x $h NAXIS2)

            if [[ "$NAXIS1_image" == "$NAXIS1_plane" ]] && [[ "$NAXIS2_image" == "$NAXIS2_plane" ]]; then
                astarithmetic $input -h$h $i -h$h - -o$temp_file
                astfits $temp_file --copy=1 -o $output
                rm $temp_file
            else
                python3 $pythonScriptsPath/moveSurfaceFitToFullGrid.py $input $i $h $NAXIS1_image $NAXIS2_image $directoryToStoreSkySubtracted/"planeToSubtract_"$a".fits"
                astarithmetic $input -h$h $directoryToStoreSkySubtracted/"planeToSubtract_"$a".fits" -h1 - -o$temp_file
                astfits $temp_file --copy=1 -o $output
                rm $temp_file
                rm $directoryToStoreSkySubtracted/"planeToSubtract_"$a".fits"
            fi
        done
    fi
}
export -f subtractSkyForFrame

subtractSky() {
    local framesToSubtract=$1
    local directoryToStoreSkySubtracted=$2
    local directoryToStoreSkySubtracteddone=$3
    local directoryWithSkyValues=$4
    local constantSky=$5


    if ! [ -d $directoryToStoreSkySubtracted ]; then mkdir $directoryToStoreSkySubtracted; fi
    if [ -f $directoryToStoreSkySubtracteddone ]; then
        echo -e "\n\tSky substraction is already done for the science images\n"
    else
    framesToSubtractSky=()
    for a in $(seq 1 $totalNumberOfFrames); do
            framesToSubtractSky+=("$a")            
    done
    printf "%s\n" "${framesToSubtractSky[@]}" | parallel -j "$num_cpus" subtractSkyForFrame {} $directoryWithSkyValues $framesToSubtract $directoryToStoreSkySubtracted $constantSky
    echo done > $directoryToStoreSkySubtracteddone
    fi
}
export -f subtractSky

# Functions for decals data
# The function that is to be used (the 'public' function using OOP terminology) is 'prepareSurveyDataForPhotometricCalibration'
getBricksWhichCorrespondToFrame() {
    local frame=$1
    local frameBrickMapFile=$2

    bricks=$( awk -v var=$(basename $frame) '$1==var { match($0, /\[([^]]+)\]/, arr); print arr[1] }' $frameBrickMapFile )
    IFS=", "
    read -r -a array <<< $bricks

    # Remove the single quotes from elements
    for ((i=0; i<${#array[@]}; i++)); do
            array[$i]=${array[$i]//\'/}
    done
    echo "${array[@]}"
}
export -f getBricksWhichCorrespondToFrame

getParametersFromHalfMaxRadius() {
    local image=$1
    local gaiaCatalogue=$2
    local kernel=$3
    local tmpFolder=$4

    # The output of the commands are redirected to /dev/null because otherwise I cannot return the median and std.
    # Quite uncomfortable the return way of bash. Nevertheless, the error output is not modified so if an instruction fails we still get the error message.
    astconvolve $image --kernel=$kernel --domain=spatial --output=$tmpFolder/convolved.fits 1>/dev/null
    astnoisechisel $image -h1 -o $tmpFolder/det.fits --convolved=$tmpFolder/convolved.fits --tilesize=20,20 --detgrowquant=0.95 --erode=4 --numthreads=$num_cpus 1>/dev/null
    astsegment $tmpFolder/det.fits -o $tmpFolder/seg.fits --snquant=0.1 --gthresh=-10 --objbordersn=0    --minriverlength=3 1>/dev/null
    astmkcatalog $tmpFolder/seg.fits --ra --dec --magnitude --half-max-radius --sum --clumpscat -o $tmpFolder/decals.txt --zeropoint=22.5 1>/dev/null
    astmatch $tmpFolder/decals_c.txt --hdu=1    $BDIR/catalogs/"$objectName"_Gaia_eDR3.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=bRA,bDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o $tmpFolder/match_decals_gaia.txt 1>/dev/null

    numOfStars=$( cat $tmpFolder/match_decals_gaia.txt | wc -l )
    median=$( asttable $tmpFolder/match_decals_gaia.txt -h1 -c3 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median )
    std=$( asttable $tmpFolder/match_decals_gaia.txt -h1 -c3 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std )
    rm $tmpFolder/*
    echo $median $std $numOfStars
}
export -f getParametersFromHalfMaxRadius


downloadSurveyData() {
    local mosaicDir=$1
    local surveyImagesDir=$2
    local bricksIdentificationFile=$3
    local filters=$4
    local ra=$5
    local dec=$6
    local fieldSizeDeg=$7
    local gaiaCatalogue=$8
    local survey=$9
    local sizeOfBrick=${10}
    
    echo -e "\n·Downloading ${survey} bricks"
    donwloadMosaicDone=$surveyImagesDir/done_downloads.txt

    if ! [ -d $mosaicDir ]; then mkdir $mosaicDir; fi
    if ! [ -d $surveyImagesDir ]; then mkdir $surveyImagesDir; fi
    if [ -f $donwloadMosaicDone ]; then
        echo -e "\n\tMosaic images already downloaded\n"
    else
        rm $bricksIdentificationFile # Remove the brick indentification file. This is done to avoid problems with files of previous executions        
        echo "Downloading $survey bricks for field centered at ($ra, $dec) and size $fieldSizeDeg deg; filters: " $filters
        python3 $pythonScriptsPath/downloadBricksForFrame.py $filters $surveyImagesDir $ra $dec $fieldSizeDeg $mosaicDir $bricksIdentificationFile $gaiaCatalogue $survey $sizeOfBrick
        echo "done" > $donwloadMosaicDone
    fi
}
export -f downloadSurveyData

downloadSpectra() {
    local mosaicDir=$1
    local spectraDir=$2
    local ra=$3
    local dec=$4
    local sizeOfOurFieldDegrees=$5
    local surveyForSpectra=$6

    spectraDone=$spectraDir/done.txt
    if [ -f $spectraDone ]; then
        echo -e "\nSpectra already downloaded\n"
    else
        python3 $pythonScriptsPath/downloadSpectraForField.py $mosaicDir $spectraDir $ra $dec $sizeOfOurFieldDegrees $surveyForSpectra
        echo "done" > $spectraDone
    fi
}
export -f downloadSpectra

addTwoFiltersAndDivideByTwo() {
    local decalsImagesDir=$1
    local filter1=$2
    local filter2=$3

    addBricksDone=$mosaicDir/decalsImages/done_adding.txt
    if [ -f $addBricksDone ]; then
        echo -e "\nDecals '$filter1' and '$filter2' bricks are already added\n"
    else
        for file in "$decalsImagesDir"/*_$filter1.fits; do
            # The following lines depend on the name of the decals images, which is defined in the python script "./decals_GetAndDownloadBricks.py"
            brickName=$(basename "$file" | cut -d '_' -f3)
            filter1File="decal_image_"$brickName"_$filter1.fits"
            filter2File="decal_image_"$brickName"_$filter2.fits"

            echo -e "Adding the files " $filter1File " and the file " $filter2File
            astarithmetic $decalsImagesDir/$filter1File -h1 $decalsImagesDir/$filter2File -h1 + 2 / -o$decalsImagesDir/"decal_image_"$brickName"_"$filter1"+"$filter2"_div2.fits"
        done
        echo done > $addBricksDone
    fi
}
export -f addTwoFiltersAndDivideByTwo

gainCorrection() {
    local base=$1
    local inputDir=$2
    local ringDir=$3
    local outDir=$4
    local blockScale=$5
    local noisechisel_param=$6
    image=$inputDir/$base
    output=$outDir/$base
    astfits $image --copy=0 --primaryimghdu -o$output
    ringFile=$ringDir/ring.fits
    noiseOut=$outDir/noise_$base
    maskOut=$outDir/mask_$base
    gainOut=$outDir/gain_$base
    runNoiseChiselOnFrame $base $inputDir $outDir $blockScale $noisechisel_param
    mv $outDir/$base $noiseOut
    for h in $(seq 1 $num_ccd); do
        astarithmetic $image -h$h $noiseOut -h$h 0 ne nan where -q -o$outDir/temp_$base
        astarithmetic $outDir/temp_$base -h1 $ringFile -h$h 0 eq nan where -q -o$maskOut
        gain_h=$(aststatistics $maskOut --sigclip-median -q)
        if [ $h -eq 1 ]; then
            gain_ref=$gain_h
            astfits $image --copy=$h -o $output
        else
            astarithmetic $image -h$h $gain_ref x $gain_h / -o$gainOut
            astfits $gainOut --copy=1 -o$output
            rm $gainOut
        fi
        rm $outDir/temp_$base $maskOut
    done
    rm $noiseOut

}
export -f gainCorrection

downloadGaiaCatalogue() {
    local query=$1
    local catdir=$2
    local catName=$3

    astquery $query -o $catdir/"$objectName"_Gaia_eDR3_tmp.fits
    asttable $catdir/"$objectName"_Gaia_eDR3_tmp.fits -c1,2,3 -c'arith $4 abs' -c'arith $5 3 x' -c'arith $6 abs' -c'arith $7 3 x' -c'arith $8 abs' -c'arith $9 3 x' --noblank=4 -o$catdir/tmp.txt

    # I have explored 3 different ways of selecting good stars. 
    # From the most restrictive to the less restrictive:

    # # Here I demand that the gaia object fulfills simultaneously that:
    # # 1.- Parallax > 3 times its error
    # # 2.- Proper motion (ra) > 3 times its error
    # # 3.- Proper motion (dec) > 3 times its error
    # asttable $catdir/tmp.txt -c1,2,3 -c'arith $4 $4 $5 gt 1000 where' -c'arith $6 $6 $7 gt 1000 where' -c'arith $8 $8 $9 gt 1000 where' -o$catdir/test_.txt
    # asttable $catdir/test_.txt -c1,2,3 -c'arith $4 $5 + $6 +' -o$catdir/test1.txt
    # asttable $catdir/test1.txt -c1,2,3 --range=ARITH_2,2999,3001 -o $catName

    # # Here I only demand that the parallax is > 3 times its error
    # asttable $catdir/tmp.txt -c1,2,3 -c'arith $4 $4 $5 gt 1000 where' -o$catdir/test_.txt
    # asttable $catdir/test_.txt -c1,2,3 --range=ARITH_2,999,1001 -o $catName

    # Here I  demand that the parallax OR a proper motion is > 3 times its error
    asttable $catdir/tmp.txt -c1,2,3 -c'arith $4 $4 $5 gt 1000 where' -c'arith $6 $6 $7 gt 1000 where' -c'arith $8 $8 $9 gt 1000 where' -o$catdir/test_.txt
    asttable $catdir/test_.txt -c1,2,3 -c'arith $4 $5 + $6 +' -o$catdir/test1.txt
    asttable $catdir/test1.txt -c1,2,3 --range=ARITH_2,999,3001 -o $catName

    rm $catdir/test1.txt $catdir/tmp.txt $catdir/"$objectName"_Gaia_eDR3_tmp.fits $catdir/test_.txt
}
export -f downloadGaiaCatalogue

downloadPanstarrsCatalogue() {
    local query=$1
    local catdir=$2
    local catName=$3

    astquery $query -o $catdir/"$objectName"_Panstarrs_S1_tmp.fits
    asttable $catdir/"$objectName"_Panstarrs_S1_tmp.fits -c1,2,3  --colmetadata=1,RA,deg --colmetadata=2,DEC,deg --colmetadata=3,phot_g_mean_mag,mag -o$catName

    # We are downloading Panstarrs for hipercam, this will stay as it is. We change the metadata to mimic that of 

    rm $catdir/"$objectName"_Panstarrs_S1_tmp.fits 
}
export -f downloadPanstarrsCatalogue

downloadIndex() {
    local re=$1
    local catName=$2
    local indexdir=$3

    build-astrometry-index -i $catName -e1 \
                            -P $re \
                            -S phot_g_mean_mag \
                            -E -A RA -D  DEC \
                            -o $indexdir/index_$re.fits;
}
export -f downloadIndex

solveField() {
    local i=$1
    local solve_field_L_Param=$2
    local solve_field_H_Param=$3
    local solve_field_u_Param=$4
    local ra_gal=$5
    local dec_gal=$6
    local confFile=$7
    local astroimadir_layer=$8
    local sexcfg_sf=$9
    local sizeOfOurFieldDegrees=${10}
    base=$( basename $i)
    LC_NUMERIC=C  # Format to get rid of scientific notation if needed

    pointingRAValue=$( astfits $i -h0 --keyvalue=$pointingRA --quiet)
    #pointingRAValue=$( printf "%.8f\n" " $pointingRAValue")
    pointingDecValue=$( astfits $i -h0 --keyvalue=$pointingDEC --quiet)
    if [[ "$pointingRAUnits" == "hours" ]]; then
        pointingRAValue=$( printf "%.8f\n" " $pointingRAValue")
        pointRA=$(echo "$pointingRAValue * 15" | bc -l)
    elif [[ "$pointingRAUnits" == "deg" || "$pointingRAUnits" == "degrees" ]]; then
        pointRA="$( printf "%.8f\n" " $pointingRAValue")"
    elif [[ "$pointingRAUnits" == "hms" ]]; then
            ra_dec=$(skycoor -d "$pointingRAValue" "$pointingDecValue")
            pointRA=$(echo "$ra_dec" | awk '{print $1}')
    else
        echo "Error: Unsupported RA units: $pointingRAUnits"
        exit 888
    fi

    #pointingDecValue=$( astfits $i -h0 --keyvalue=$pointingDEC --quiet)
    #pointingDecValue=$( printf "%.8f\n" " $pointingDecValue")
    if [[ "$pointingDECUnits" == "hours" ]]; then
        pointDec=$(echo "$pointingDecValue * 15" | bc -l)
    elif [[ "$pointingDECUnits" == "deg" || "$pointingDECUnits" == "degrees" ]]; then
        pointDec="$( printf "%.8f\n" " $pointingDecValue")"
    elif [[ "$pointingDECUnits" == "dms" ]]; then
            ra_dec=$(skycoor -d "$pointingRAValue" "$pointingDecValue")
            pointDec=$(echo "$ra_dec" | awk '{print $2}')
    else
        echo "Error: Unsupported RA units: $pointingDECUnits"
        exit 888
    fi
    
    
    # The default sextractor parameter file is used.
    # I tried to use the one of the config directory (which is used in other steps), but even using the default one, it fails
    # Maybe a bug? I have not managed to make it work
    ### Multi-layer problem: solve-field does not work with multiple layers. Because of that, we run solve-field into each of the layers and then store them into a single .fits with multiple layers
    layer_temp=$astroimadir_layer/layer_$base
    
    for h in $(seq 1 $num_ccd); do
        image_temp=image"$h"_$base
        astfits $i --copy=$h -o $image_temp
        layer_temp=$astroimadir_layer/layer"$h"_$base
        max_attempts=4
        attempt=1
        while [ $attempt -le $max_attempts ]; do
            #Sometimes the output of solve-field is not properly writen in the computer (.i.e, size of file=0). 
            #Because of that, we iterate solve-field in a maximum of 4 times until file is properly saved
            solve-field $image_temp --no-plots \
            -L $solve_field_L_Param -H $solve_field_H_Param -u $solve_field_u_Param \
            --ra $pointRA --dec=$pointDec --radius $sizeOfOurFieldDegrees \
            --overwrite --extension 1 --config $confFile/astrometry_$objectName.cfg --no-verify \
            --use-source-extractor --source-extractor-path=/usr/bin/source-extractor \
            --source-extractor-config=$sexcfg_sf --x-column X_IMAGE --y-column Y_IMAGE \
            --sort-column MAG_AUTO --sort-ascending  \
            -Unone --temp-axy  -Snone -Mnone -Rnone -Bnone -N$layer_temp ;
            if [ -s "$layer_temp" ]; then
                attempt=$max_attempts
            fi
            
            ((attempt++))
        done
        
        rm $image_temp image"$h"_*.wcs
        #rm $layer_temp
    done
}
export -f solveField

runSextractorOnImage() {
    local a=$1
    local sexcfg=$2
    local sexparam=$3
    local sexconv=$4
    local astroimadir=$5
    local sexdir=$6
    local saturationThreshold=$7
     

    # Here I put the saturation threshold and the gain directly.
    # This is because it's likely that we end up forgetting about tuning the sextractor configuration file but we will be more careful with the configuration file of the reductions
    # These two values (saturation level and gain) are key for astrometrising correctly, they are used by scamp for identifying saturated sources and weighting the sources
    # I was, in fact, having frames bad astrometrised due to this parameters.
    i=$astroimadir/"$a".fits
    
    source-extractor $i   -c $sexcfg -PARAMETERS_NAME $sexparam -FILTER_NAME $sexconv -CATALOG_NAME $sexdir/"$a".cat -SATUR_LEVEL=$saturationThreshold 

}
export -f runSextractorOnImage


# warpDecalsBrick() {
#     a=$1
#     swarpedImagesDir=$2
#     decalsImagesDir=$3
#     scaleFactor=$4
#     swarpcfg=$5
#     ra=$6
#     dec=$7
#     mosaicSize=$8
#     decalsPxScale=$9

#     decalsImage=$decalsImagesDir/$a
#     downSampledImages="$swarpedImagesDir/originalGrid_$(basename $a)"

#     astwarp $decalsImage --scale=$scaleFactor -o $downSampledImages
#     swarp -c $swarpcfg $downSampledImages -CENTER $ra,$dec -IMAGE_SIZE $mosaicSize,$mosaicSize -IMAGEOUT_NAME $swarpedImagesDir/"$a"_swarp1.fits \
#                         -WEIGHTOUT_NAME $swarpedImagesDir/"$a"_swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $decalsPxScale -PIXELSCALE_TYPE MANUAL
#     astarithmetic $swarpedImagesDir/"$a"_swarp_w1.fits -h0 set-i i i 0 lt nan where -o$swarpedImagesDir/"$a"_temp1.fits
#     astarithmetic $swarpedImagesDir/"$a"_swarp1.fits -h0 $swarpedImagesDir/"$a"_temp1.fits -h1 0 eq nan where -o$swarpedImagesDir/commonGrid_"$(basename $a)"
#     rm $swarpedImagesDir/"$a"_swarp_w1.fits $swarpedImagesDir/"$a"_swarp1.fits $swarpedImagesDir/"$a"_temp1.fits
# }
# export -f warpDecalsBrick

# buildDecalsMosaic() {
#     # We only need the mosaic in order to download the gaia catalogue. That's why downgrade the bricks
#     # Values for original decals resolution. As a reminder decals original pxScale 0.2626 arcsec/px
#     mosaicDir=$1
#     decalsImagesDir=$2
#     swarpcfg=$3
#     ra=$4
#     dec=$5
#     filter=$6
#     swarpedImagesDir=$7
#     dataPixelScale=$8 # Pixel scale of our data. In order to do realisticly we should downgrade decals data to the same resolution as our data
#     sizeOfOurFieldDegrees=$9 # Estimation of how big is our field

#     originalDecalsPxScale=0.262 # arcsec/px

#     buildMosaicDone=$swarpedImagesDir/done_t.xt
#     if ! [ -d $swarpedImagesDir ]; then mkdir $swarpedImagesDir; fi
#     if [ -f $buildMosaicDone ]; then
#         echo -e "\n\tMosaic already built\n"
#     else


#         # This depends if you want to calibrate with the original resolution (accurate but slower) or downgrade it to your data resolution
#         scaleFactor=1 # Original resolution
#         # scaleFactor=$(awk "BEGIN {print $originalDecalsPxScale / $dataPixelScale}") # Your data resolution

#         # We dinamically compute the new pixelscale and the grid for the full mosaic of the field
#         newDecalsPxScale=$( echo "($originalDecalsPxScale / $scaleFactor)" | bc -l)
#         mosaicSize=$(echo "($sizeOfOurFieldDegrees * 3600) / ($originalDecalsPxScale / $scaleFactor)" | bc -l)

#         if [ "$filter" = "lum" ]; then
#             bricks=$(ls -v $decalsImagesDir/*_g+r_div2.fits)
#         elif [ "$filter" = "i" ]; then
#             bricks=$(ls -v $decalsImagesDir/*_r+z_div2.fits)
#         else
#             bricks=$(ls -v $decalsImagesDir/*$filter.fits)
#         fi

#         brickList=()
#         for a in $bricks; do
#             brickList+=("$( basename $a )")
#         done
#         printf "%s\n" "${brickList[@]}" | parallel -j "$num_cpus" warpDecalsBrick {} $swarpedImagesDir $decalsImagesDir $scaleFactor $swarpcfg $ra $dec $mosaicSize $newDecalsPxScale

#         sigma=2
#         astarithmetic $(ls -v $swarpedImagesDir/commonGrid*.fits) $(ls -v $swarpedImagesDir/commonGrid*.fits | wc -l) -g1 $sigma 0.2 sigclip-median -o $mosaicDir/mosaic.fits
#         echo done > $buildMosaicDone
#     fi
# }
# export -f buildDecalsMosaic

decompressDecalsFrame() {
    local frameName=$1
    local dirWithbricks=$2

    funpack -O $dirWithbricks/decompressed_$frameName $dirWithbricks/$frameName 
}
export -f decompressDecalsFrame

decompressBricks() {
    local dirWithBricks=$1

    decompressedBricks=$dirWithBricks/decompressed_done.txt
    if [ -f $decompressedBricks ]; then
        echo -e "\n\tBricks already decompressed\n"
    else
        brickList=()
        for i in $( ls $dirWithBricks/*.fits); do
            currentName=$( basename $i )
            brickList+=("$currentName")
        done
        printf "%s\n" "${brickList[@]}" | parallel -j "$num_cpus" decompressDecalsFrame {} $dirWithBricks
        echo "done" > $decompressedBricks
    fi
}
export -f decompressBricks

divideByExpTimeAndMoveZPForPanstarrsFrame() {
    local brickName=$1
    local dirWithBricks=$2

    cal_tm="cal_tm_"$brickName
    ##We need to divide /10 to get into zp=22.5, and to divide between exposure time
    texp=$(astfits $dirWithBricks/$brickName -h0 --keyvalue=EXPTIME -q )
    astarithmetic $dirWithBricks/$brickName -h0 10 / $texp / -o $dirWithBricks/$cal_tm
    cal="cal_"$brickName
    astfits $dirWithBricks/$cal_tm --copy=1 --primaryimghdu -o $dirWithBricks/$cal
    astfits $dirWithBricks/$cal -h0 --write=ZP,22.5
    rm $dirWithBricks/$cal_tm
}
export -f divideByExpTimeAndMoveZPForPanstarrsFrame

divideByExpTimeAndMoveZPForPanstarrs() {
    local dirWithBricks=$1

    calibratedBricks=$dirWithBricks/recalibrated_done.txt
    if [ -f $calibratedBricks ]; then
        echo -e "\n\tPanstarrs bricks already re-calibrated\n"
    else
        brickList=()
        for i in $( ls $dirWithBricks/t*.fits); do
            currentName=$( basename $i )
            brickList+=("$currentName")
        done
        printf "%s\n" "${brickList[@]}" | parallel -j "$num_cpus" divideByExpTimeAndMoveZPForPanstarrsFrame {} $dirWithBricks
        echo "done" > $calibratedBricks
    fi
}
export -f divideByExpTimeAndMoveZPForPanstarrs

selectStarsAndSelectionRangeSurvey() {
    local dirWithBricks=$1
    local cataloguedir=$2
    local methodToUse=$3
    local survey=$4
    local apertureUnits=$5

    starSelectionDone=$cataloguedir/done.txt
    if ! [ -d $cataloguedir ]; then mkdir $cataloguedir; fi
    if [ -f $starSelectionDone ]; then
        echo -e "\n\tStar and range selection for calibration already done\n"
    else
        brickList=()
        if [ "$survey" = "DECaLS" ]; then
            for i in $( ls $dirWithBricks/decompressed_*.fits); do
                currentName=$( basename $i )
                brickList+=("$currentName")
            done
        elif [ "$survey" = "PANSTARRS" ]; then
            for i in $( ls $dirWithBricks/cal_*.fits); do
                currentName=$( basename $i )
                brickList+=("$currentName")
            done
        fi

        headerWithData=0 # After decompressing the data ends up in the hdu 0
        noisechiselTileSize=50
        printf "%s\n" "${brickList[@]}" | parallel -j "$num_cpus" selectStarsAndRangeForCalibrateSingleFrame {} $dirWithBricks $cataloguedir $headerWithData $methodToUse $noisechiselTileSize YES $apertureUnits
        echo "done" > $starSelectionDone
    fi
}
export -f selectStarsAndSelectionRangeSurvey

produceHalfMaxRadiusPlotsForDecals() {
    local folderWithCatalogues=$1
    local outputDir=$2
    local filter=$3

    decalsHalfMaxRadiusPlotsDone=$outputDir/done.txt

    if ! [ -d $outputDir ]; then mkdir $outputDir; fi
    if [ -f $decalsHalfMaxRadiusPlotsDone ]; then
        echo -e "\n\tHalf-max radius vs magnitude plots already done\n"
    else
        for currentFullCatalogue in $( ls $folderWithCatalogues/catalogue_*.cat ); do
            brickName=$(echo "$currentFullCatalogue" | awk -F'_image_' '{print $2}' | awk -F'_' '{print $1}')
            matchedCatalogue=$folderWithCatalogues/match_decompressed_decal_image_"$brickName"_"$filter".fits_my_gaia.txt
            rangeFile=$folderWithCatalogues/range_decompressed_decal_image_"$brickName"_"$filter".fits.txt
            read minr maxr < <(awk '{print $3, $4}' $rangeFile)

            plotXLowerLimit=1
            plotXHigherLimit=10
            plotYLowerLimit=13
            plotYHigherLimit=24
            python3 $pythonScriptsPath/diagnosis_halfMaxRadVsMag.py $currentFullCatalogue $matchedCatalogue -1 $minr $maxr $outputDir/$brickName.png  \
                $plotXLowerLimit $plotXHigherLimit $plotYLowerLimit $plotYHigherLimit
        done
        
        echo "done" > $decalsHalfMaxRadiusPlotsDone
    fi
}
export -f produceHalfMaxRadiusPlotsForDecals

produceHalfMaxRadiusPlotsForPanstarrs() {
    local folderWithCatalogues=$1
    local outputDir=$2
    local filter=$3
    #ps==panstarrs
    psHalfMaxRadiusPlotsDone=$outputDir/done.txt
    if ! [ -d $outputDir ]; then mkdir $outputDir; fi
    if [ -f $psHalfMaxRadiusPlotsDone ]; then
        echo -e "\n\tHalf-max radius vs magnitude plots already done\n"
    else
        for currentFullCatalogue in $( ls $folderWithCatalogues/catalogue_*.cat ); do
            brickName=$(echo "$currentFullCatalogue" | awk -F'catalogue_' '{print $2}' | awk -F'.cat' '{print $1}')
            matchedCatalogue=$folderWithCatalogues/match_"$brickName"_my_gaia.txt
            rangeFile=$folderWithCatalogues/range_"$brickName".txt
            
            read minr maxr < <(awk '{print $3, $4}' $rangeFile)
            plotXLowerLimit=1
            plotXHigherLimit=10
            plotYLowerLimit=13
            plotYHigherLimit=24
            python3 $pythonScriptsPath/diagnosis_halfMaxRadVsMag.py $currentFullCatalogue $matchedCatalogue -1 $minr $maxr $outputDir/$brickName.png  \
                $plotXLowerLimit $plotXHigherLimit $plotYLowerLimit $plotYHigherLimit
        done
        
        echo "done" > $psHalfMaxRadiusPlotsDone
    fi
}
export -f produceHalfMaxRadiusPlotsForPanstarrs

performAperturePhotometryToSingleBrick() {
    local brick=$1
    local brickDir=$2
    local automaticallySelectedDir=$3
    local outputCat=$4
    local filter=$5
    local numberOfApertureForRecuperateGAIA=$6
    local survey=$7
    if [[ "$survey" = "DECaLS" ]]; then
        brickName=decompressed_decal_image_"$brick"_"$filter".fits
    elif [[ "$survey" = "PANSTARRS" ]]; then
        brickName=cal_"$brick".fits
    fi
    brickImage=$brickDir/$brickName
    fileWithAperture=$automaticallySelectedDir/range_$brickName.txt
    automaticCatalogue=$automaticallySelectedDir/selected_"$brickName"_automatic.txt

    r_decals_pix_=$(awk 'NR==1 {printf $1}' $fileWithAperture)
    r_decals_pix=$(astarithmetic $r_decals_pix_ $numberOfApertureForRecuperateGAIA. x -q )

    dataHdu=0

    # raColumnName=RA
    # decColumnName=DEC
    # photometryOnImage_noisechisel $brick $outputCat $automaticCatalogue $brickImage $r_decals_pix $outputCat/$brick.cat \
    #                             22.5 $dataHdu $raColumnName $decColumnName

    # This code is only for checking the same aperture as the sloan SDSS spectrograph fiber
    # r_decals_pix=4
    # echo "Realizando fotometría en brick $brick con apertura: $r_decals_pix"


    columnWithXCoordForDecalsPx=0 # These numbers come from how the catalogue of the matches stars is built. This is not very clear right now, should be improved
    columnWithYCoordForDecalsPx=1
    columnWithXCoordForDecalsWCS=2
    columnWithYCoordForDecalsWCS=3 
    photometryOnImage_photutils $brick $outputCat $automaticCatalogue $brickImage $r_decals_pix $outputCat/$brick.cat 22.5 $dataHdu \
                                $columnWithXCoordForDecalsPx $columnWithYCoordForDecalsPx $columnWithXCoordForDecalsWCS $columnWithYCoordForDecalsWCS
}
export -f performAperturePhotometryToSingleBrick

performAperturePhotometryToBricks() {
    local brickDir=$1
    local automaticallySelectedDir=$2
    local outputCat=$3
    local filter=$4
    local survey=$5
    local nomberOfApertureForRecuperateGaia=$6
   
    outputDone=$outputCat/done.txt
    if ! [ -d $outputCat ]; then mkdir $outputCat; fi
    if [ -f $outputDone ]; then
        echo -e "\n\tCatalogues done with aperture photometry already done\n"
    else
        brickList=()
        #If aperture photometry is done in Panstarrs, brickNames are different
        
        if [[ "$survey" = "DECaLS" ]]; then
            for a in $( ls $brickDir/decompressed*.fits); do
                brickName=$(echo "$a" | awk -F'_image_' '{print $2}' | awk -F'_' '{print $1}')
                brickList+=("$brickName")
            done
        elif [[ "$survey" = "PANSTARRS" ]]; then
            for a in $( ls $brickDir/cal_*.fits); do
                brickName=$(echo "$a" | awk -F'cal_' '{print $2}' | awk -F'.fits' '{print $1}')
                brickList+=("$brickName")
            done
        else
            echo "Error: "$survey" not supported for photometric calibration!"
            exit 1
        fi
        
        
        printf "%s\n" "${brickList[@]}" | parallel -j "$num_cpus" performAperturePhotometryToSingleBrick {}  $brickDir $automaticallySelectedDir $outputCat $filter $numberOfApertureForRecuperateGAIA $survey
        echo "done" > $outputDone
    fi
}
export -f performAperturePhotometryToBricks

prepareCalibrationData() {
    local surveyForCalibration=$1
    local referenceImagesForMosaic=$2
    local aperturePhotDir=$3
    local filter=$4
    local ra=$5
    local dec=$6
    local mosaicDir=$7
    local selectedSurveyStarsDir=$8
    local rangeUsedSurveyDir=$9
    local dataPixelScale=${10}
    local sizeOfOurFieldDegrees=${11}
    local gaiaCatalogue=${12}
    local surveyForSpectra=${13}
    local apertureUnits=${14}
    local folderWithTransmittances=${15}
    local filterCorrectionCoeff=${16}
    local calibrationBrightLimit=${17}
    local calibrationFaintLimit=${18}
    local mosaicDone=${19}
    local sizeOfBrick=${20}


    if ! [ -d $mosaicDir ]; then mkdir $mosaicDir; fi
    if [ -f $mosaicDone ]; then
        echo -e "\nSurvey data already prepared for photometric calibration\n"
    else
        if [[ "$surveyForCalibration" == "SPECTRA" ]]; then
            spectraDir=$mosaicDir/spectra

            writeTimeOfStepToFile "Spectra data processing" $fileForTimeStamps
            transmittanceCurveFile=$folderWithTransmittances/"$telescope"_"$filter".dat
            prepareSpectraDataForPhotometricCalibration $spectraDir $filter $ra $dec $mosaicDir $aperturePhotDir $sizeOfOurFieldDegrees $surveyForSpectra $transmittanceCurveFile

        else
            surveyImagesDir=$mosaicDir/surveyImages
            writeTimeOfStepToFile "Survey data processing" $fileForTimeStamps

            prepareSurveyDataForPhotometricCalibration $referenceImagesForMosaic $surveyImagesDir $filter $ra $dec $mosaicDir $selectedSurveyStarsDir $rangeUsedSurveyDir \
                                                $dataPixelScale $surveyForCalibration $sizeOfOurFieldDegrees $gaiaCatalogue $aperturePhotDir $apertureUnits $folderWithTransmittances "$filterCorrectionCoeff" \
                                                $calibrationBrightLimit $calibrationFaintLimit $sizeOfBrick
        fi
        
        echo done > $mosaicDone
    fi
}
export -f prepareCalibrationData

prepareSpectraDataForPhotometricCalibration() {
    local spectraDir=$1
    local filter=$2
    local ra=$3
    local dec=$4
    local mosaicDir=$5
    local aperturePhotDir=$6
    local sizeOfOurFieldDegrees=$7
    local surveyForSpectra=$8
    local transmittanceCurveFile=$9

    if ! [ -d $spectraDir ]; then mkdir $spectraDir; fi
    downloadSpectra $mosaicDir $spectraDir $ra $dec $sizeOfOurFieldDegrees $surveyForSpectra
    
    aperturePhotDone=$aperturePhotDir/done.txt
    if ! [ -d $aperturePhotDir ]; then mkdir $aperturePhotDir; fi
    if [ -f $aperturePhotDone ]; then
        echo -e "\nThe catalogue with the magnitudes of the spectra already built\n"
    else
        output_tmpCat=$aperturePhotDir/wholeFieldPhotometricCatalogue_tmp.cat
        outputCat=$aperturePhotDir/wholeFieldPhotometricCatalogue.cat
        transmittanceWavelengthUnits=A # At the beginning of the pipeline everything was transformed to A. Nevertheless the python scripts handles transmittances in nm 
        python3 $pythonScriptsPath/getMagnitudFromSpectra.py $spectraDir $transmittanceCurveFile $transmittanceWavelengthUnits $output_tmpCat $surveyForSpectra

        asttable $output_tmpCat -p4 --colmetadata=2,RA,deg,"Right ascension" \
                        --colmetadata=3,DEC,none,"Declination" \
                        --colmetadata=4,MAGNITUDE,none,"Magnitude" \
                        --colmetadata=5,SUM,none,"sum" \
                        --output=$outputCat
        rm $output_tmpCat
        echo "done" > $aperturePhotDone
    fi
}
export -f prepareSpectraDataForPhotometricCalibration

prepareSurveyDataForPhotometricCalibration() {
    local referenceImagesForMosaic=$1
    local surveyImagesDir=$2
    local filter=$3
    local ra=$4
    local dec=$5
    local mosaicDir=$6
    local selectedSurveyStarsDir=$7
    local rangeUsedSurveyDir=$8
    local dataPixelScale=$9
    local survey=${10}
    local sizeOfOurFieldDegrees=${11} 
    local gaiaCatalogue=${12}
    local aperturePhotDir=${13}
    local apertureUnits=${14}
    local folderWithTransmittances=${15}
    local filterCorrectionCoeff=${16}
    local calibrationBrightLimit=${17}
    local calibrationFaintLimit=${18}
    local sizeOfBrick=${19}
   
    sizeOfFieldForCalibratingPANSTARRStoGAIA=1.5
    sizeOfBrick_gaia=3600

    echo -e "\n ${GREEN} ---Preparing ${survey} data--- ${NOCOLOUR}"
    bricksIdentificationFile=$surveyImagesDir/brickIdentification.txt

    surveyImagesDirForGaiaCalibration=$surveyImagesDir"ForGAIACalibration"
    bricksIdentificationFileForGaiaCalibration=$surveyImagesDirForGaiaCalibration/brickIdentification.txt

    # When calibrating, we are taking into account the difference in the filters between our filter and the filter of the survey. For correcting this we apply a
    # colour correction. So we always need to download panstarrs 'g' and panstarrs 'r' band.

    # Additionally, we work with two downloads of survey bricks. One of the field to calibrate the data to reduce and another one
    # (independent of the field to reduce) to calibrate the survey used for calibration (PANSTARRS or DECaLS) to our reference framework of GAIA
    # This is useful because regardless of the size of our field to reduce we can compute the offset in a constant size field
    # Could this be done more efficient and not have duplicity? yes. Would it need quite extra logic? Yes
    # Because sometimes the field for reducing the data is bigger than the field for calibrating the imaging survey, but other times
    # it is not. So it would need to do checks and stuff that we are not doing right now. Either way this is not bottle neck in the pipeline.
    surveyImagesDir_g="$surveyImagesDir"_g
    surveyImagesDirForGaiaCalibration_g="$surveyImagesDirForGaiaCalibration"_g
    surveyImagesDir_r="$surveyImagesDir"_r
    surveyImagesDirForGaiaCalibration_r="$surveyImagesDirForGaiaCalibration"_r

    bricksIdentificationFile_g=$surveyImagesDir_g/brickIdentification.txt
    bricksIdentificationFileForGaiaCalibration_g=$surveyImagesDirForGaiaCalibration_g/brickIdentification.txt
    bricksIdentificationFile_r=$surveyImagesDir_r/brickIdentification.txt
    bricksIdentificationFileForGaiaCalibration_r=$surveyImagesDirForGaiaCalibration_r/brickIdentification.txt

    downloadSurveyData $mosaicDir $surveyImagesDir_g $bricksIdentificationFile_g "g" $ra $dec $sizeOfOurFieldDegrees $gaiaCatalogue $survey $sizeOfBrick
    downloadSurveyData $mosaicDir $surveyImagesDirForGaiaCalibration_g $bricksIdentificationFileForGaiaCalibration_g "g" $ra $dec $sizeOfFieldForCalibratingPANSTARRStoGAIA $gaiaCatalogue $survey $sizeOfBrick_gaia
    downloadSurveyData $mosaicDir $surveyImagesDir_r $bricksIdentificationFile_r "r" $ra $dec $sizeOfOurFieldDegrees $gaiaCatalogue $survey $sizeOfBrick
    downloadSurveyData $mosaicDir $surveyImagesDirForGaiaCalibration_r $bricksIdentificationFileForGaiaCalibration_r "r" $ra $dec $sizeOfFieldForCalibratingPANSTARRStoGAIA $gaiaCatalogue $survey $sizeOfBrick_gaia

    
    # If the images are donwloaded but the done.txt file is no present, the images won't be donwloaded again but
    # the step takes a while even if the images are already downloaded because we have to do a query to the database
    # in order to obtain the brickname and check if it is already downloaded or no
    if [ "$filter" = "lum" ]; then
        # This was used at the beginning as a sort of proxy to lum. This is not used anymore (we calibrate lum with spectra) but we prefer not to
        # remove the code. Anyway this is not well tested because a lot of things have changed and should need some checks.

        if ! [ -d $surveyImagesDir ]; then mkdir $surveyImagesDir; fi
        cp $surveyImagesDir_g/* $surveyImagesDir
        cp $surveyImagesDir_r/* $surveyImagesDir

        if ! [ -d $surveyImagesDirForGaiaCalibration ]; then mkdir $surveyImagesDirForGaiaCalibration; fi
        cp $surveyImagesDirForGaiaCalibration_g/* $surveyImagesDirForGaiaCalibration
        cp $surveyImagesDirForGaiaCalibration_r/* $surveyImagesDirForGaiaCalibration

        # This step creates the images (g+r)/2. This is needed because we are using a luminance filter which is a sort of (g+r)
        # The division by 2 is because in AB system we work with Janskys, which are W Hz^-1 m^-2. So we have to give a flux per wavelenght
        # So, when we add two filters we have to take into account that we are increasing the wavelength rage. In our case, 'g' and 'r' have
        # practically the same wavelenght width, so dividing by 2 is enough
        addTwoFiltersAndDivideByTwo $surveyImagesDirForGaiaCalibration "g" "r"
        addTwoFiltersAndDivideByTwo $surveyImagesDir "g" "r"
    else 
        if [[ ("$filter" == "g") || ("$filter" == "r") ]]; then
            [ -L $surveyImagesDir  ] || ln -s "$surveyImagesDir"_$filter $surveyImagesDir
            [ -L $surveyImagesDirForGaiaCalibration  ] || ln -s "$surveyImagesDirForGaiaCalibration"_$filter $surveyImagesDirForGaiaCalibration
        else
            downloadSurveyData $mosaicDir $surveyImagesDirForGaiaCalibration $bricksIdentificationFileForGaiaCalibration $filter $ra $dec $sizeOfFieldForCalibratingPANSTARRStoGAIA $gaiaCatalogue $survey $sizeOfBrick_gaia
            downloadSurveyData $mosaicDir $surveyImagesDir $bricksIdentificationFile $filter $ra $dec $sizeOfOurFieldDegrees $gaiaCatalogue $survey $sizeOfBrick
        fi
    fi


    # The photometric calibration is frame by frame, so we are not going to use the mosaic for calibration. 
    # This was implemented due to the LBT origin of the pipeline, but for doing it at original resolution takes time and memory so
    # it's not worth. I dont' delete it to let the option of using it
    # resampledDecalsBricks=$mosaicDir/resampled
    # buildDecalsMosaic $mosaicDir $decalsImagesDir $swarpcfg $ra $dec $filter $resampledDecalsBricks $dataPixelScale $sizeOfOurFieldDegrees   
    
    if [ "$survey" = "DECaLS" ]; then 
        decompressBricks $surveyImagesDir_g
        decompressBricks $surveyImagesDirForGaiaCalibration_g
        decompressBricks $surveyImagesDir_r
        decompressBricks $surveyImagesDirForGaiaCalibration_r
        decompressBricks $surveyImagesDir
        decompressBricks $surveyImagesDirForGaiaCalibration                                   # I can run noisechisel, but since it is quickly to decompress I prefer to simplify the logic of the code and decompress always
    elif [ "$survey" = "PANSTARRS" ]; then
        # Panstarrs is calibrated at ZP=25. In order to be the most general possible, we transform panstarrs data so it is at zp=22.5 and we work always with the same zp
        divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDir_g
        divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDirForGaiaCalibration_g
        divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDir_r
        divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDirForGaiaCalibration_r
        divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDir
        divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDirForGaiaCalibration
    fi
    
    methodToUse="sextractor"
    selectStarsAndSelectionRangeSurvey $surveyImagesDir_g "$selectedSurveyStarsDir"_g $methodToUse $survey $apertureUnits
    selectStarsAndSelectionRangeSurvey $surveyImagesDirForGaiaCalibration_g $selectedSurveyStarsDir"ForGAIACalibration_g" $methodToUse $survey $apertureUnits
    selectStarsAndSelectionRangeSurvey $surveyImagesDir_r "$selectedSurveyStarsDir"_r $methodToUse $survey $apertureUnits
    selectStarsAndSelectionRangeSurvey $surveyImagesDirForGaiaCalibration_r $selectedSurveyStarsDir"ForGAIACalibration_r" $methodToUse $survey $apertureUnits
    if [[ ("$filter" == "g") || ("$filter" == "r") ]]; then
        [ -L $selectedSurveyStarsDir ] || ln -s "$selectedSurveyStarsDir"_$filter $selectedSurveyStarsDir
        [ -L $selectedSurveyStarsDir"ForGAIACalibration"  ] || ln -s $selectedSurveyStarsDir"ForGAIACalibration_"$filter $selectedSurveyStarsDir"ForGAIACalibration"
    else
        selectStarsAndSelectionRangeSurvey $surveyImagesDir $selectedSurveyStarsDir $methodToUse $survey $apertureUnits
        selectStarsAndSelectionRangeSurvey $surveyImagesDirForGaiaCalibration $selectedSurveyStarsDir"ForGAIACalibration" $methodToUse $survey $apertureUnits
    fi
    
    #halfMaxRad_Mag_plots=$mosaicDir/halfMaxradius_Magnitude_plots
    #if [ "$survey" = "DECalS" ]; then
    #    produceHalfMaxRadiusPlotsForDecals $selectedSurveyStarsDir $halfMaxRad_Mag_plots $filter
    #elif [ "$survey" = "PANSTARRS" ]; then
    #    produceHalfMaxRadiusPlotsForPanstarrs $selectedSurveyStarsDir $halfMaxRad_Mag_plots $filter
    #fi
    
    if [ $apertureUnits == "FWHM" ]; then
        numberOfApertureForRecuperateGAIA=10
    elif [ $apertureUnits == "Re" ]; then 
        numberOfApertureForRecuperateGAIA=9
    else
        echo "Error. Aperture Units not recognised. We should not get there never"
    fi
    performAperturePhotometryToBricks $surveyImagesDir_g "$selectedSurveyStarsDir"_g "$aperturePhotDir"_g "g" $survey $numberOfApertureForRecuperateGAIA
    performAperturePhotometryToBricks $surveyImagesDirForGaiaCalibration_g $selectedSurveyStarsDir"ForGAIACalibration_g" $aperturePhotDir"ForGAIACalibration_g" "g" $survey $numberOfApertureForRecuperateGAIA
    performAperturePhotometryToBricks $surveyImagesDir_r "$selectedSurveyStarsDir"_r "$aperturePhotDir"_r "r" $survey $numberOfApertureForRecuperateGAIA
    performAperturePhotometryToBricks $surveyImagesDirForGaiaCalibration_r $selectedSurveyStarsDir"ForGAIACalibration_r" $aperturePhotDir"ForGAIACalibration_r" "r" $survey $numberOfApertureForRecuperateGAIA
    if [[ ("$filter" == "g") || ("$filter" == "r") ]]; then
        [ -L $aperturePhotDir ] || ln -s "$aperturePhotDir"_$filter $aperturePhotDir
        [ -L $aperturePhotDir"ForGAIACalibration" ] || ln -s $aperturePhotDir"ForGAIACalibration_"$filter $aperturePhotDir"ForGAIACalibration"
    else
        performAperturePhotometryToBricks $surveyImagesDir $selectedSurveyStarsDir $aperturePhotDir $filter $survey $numberOfApertureForRecuperateGAIA
        performAperturePhotometryToBricks $surveyImagesDirForGaiaCalibration $selectedSurveyStarsDir"ForGAIACalibration" $aperturePhotDir"ForGAIACalibration" $filter $survey $numberOfApertureForRecuperateGAIA
    fi
    
    # As mentioned in other comments and in the README, our reference framework is gaia, so I compute any offset to the 
    # photometry of PANSTARRS and GAIA magnitudes (from spectra) and correct it
    spectraDir=$mosaicDir/gaiaSpectra
    magFromSpectraDir_g=$mosaicDir/magnitudesFromGaiaSpectra_g
    magFromSpectraDir_r=$mosaicDir/magnitudesFromGaiaSpectra_r

    # These two ranges (14.5-15.5 for g and 13.65-15 for r) are tested that work for calibrating panstarrs to gaia in these bands. 
    calibrationToGAIA $spectraDir $folderWithTransmittances "g" $ra $dec $mosaicDir $sizeOfFieldForCalibratingPANSTARRStoGAIA $magFromSpectraDir_g $aperturePhotDir"ForGAIACalibration_g" 14.5 15.5
    calibrationToGAIA $spectraDir $folderWithTransmittances "r" $ra $dec $mosaicDir $sizeOfFieldForCalibratingPANSTARRStoGAIA $magFromSpectraDir_r $aperturePhotDir"ForGAIACalibration_r" 14 15
    
    read offset_g factorToApplyToCounts_g < "$mosaicDir/offsetToCorrectSurveyToGaia_g.txt"
    read offset_r factorToApplyToCounts_r < "$mosaicDir/offsetToCorrectSurveyToGaia_r.txt"

    correctOffsetFromCatalogues $aperturePhotDir"_g" $offset_g $factorToApplyToCounts_g "beforeCorrectingPanstarrsGAIAOffset"
    correctOffsetFromCatalogues $aperturePhotDir"ForGAIACalibration_g" $offset_g $factorToApplyToCounts_g "beforeCorrectingPanstarrsGAIAOffset"
    correctOffsetFromCatalogues $aperturePhotDir"_r" $offset_r $factorToApplyToCounts_r "beforeCorrectingPanstarrsGAIAOffset"
    correctOffsetFromCatalogues $aperturePhotDir"ForGAIACalibration_r" $offset_r $factorToApplyToCounts_r "beforeCorrectingPanstarrsGAIAOffset"

    if [[ ("$filter" == "g") || ("$filter" == "r") ]]; then
        : # Since the correct offset happens in the aperturePhotDir, this soft link has already been done
    else
        magFromSpectraDir=$mosaicDir/magnitudesFromGaiaSpectra

        calibrationToGAIA $spectraDir $folderWithTransmittances $filter $ra $dec $mosaicDir $sizeOfFieldForCalibratingPANSTARRStoGAIA $magFromSpectraDir $aperturePhotDir"ForGAIACalibration" $calibrationBrightLimit $calibrationFaintLimit
        read offset factorToApplyToCounts < "$mosaicDir/offsetToCorrectSurveyToGaia_"$filter".txt"
        correctOffsetFromCatalogues $aperturePhotDir $offset $factorToApplyToCounts "beforeCorrectingPanstarrsGAIAOffset"
        correctOffsetFromCatalogues $aperturePhotDir"ForGAIACalibration" $offset $factorToApplyToCounts "beforeCorrectingPanstarrsGAIAOffset"
    fi

    computeColoursAndAddThemToCatalogues $aperturePhotDir $aperturePhotDir"_g" $aperturePhotDir"_r" $filter
    applyColourcorrectionToAllCatalogues $aperturePhotDir "$filterCorrectionCoeff"

    ##Decision note: I'm gonna create multiple frame_bricks_association_ccd$h.txt for each ccd, i think it will be easier
    for imagesHdu in $(seq 1 $num_ccd); do
    
        brickDecalsAssociationFile=$mosaicDir/frames_bricks_association_ccd"$imagesHdu".txt
        if [[ -f $brickDecalsAssociationFile ]]; then
            rm $brickDecalsAssociationFile
        fi
        python3 $pythonScriptsPath/associateDecalsBricksToFrames.py $referenceImagesForMosaic $imagesHdu $bricksIdentificationFile $brickDecalsAssociationFile $survey
    done
}
export -f prepareSurveyDataForPhotometricCalibration

applyColourcorrectionToAllCatalogues() {
    local catDir=$1
    local filterCorrectionCoeff=$2

    echo -e "\n·Computing and applying magnitude and flux correction based on g-r colour\n"
    if [ -f "$catDir/colourcorrection_done.txt" ]; then
        echo -e "\n\tThe magnitude and flux correction based on the flux already applied\n"
    else
        parallel applyColourcorrectionToSingleCatalogue ::: $catDir/*.cat ::: "$filterCorrectionCoeff"

        # for i in $catDir/*.cat; do
        #     applyColourcorrectionToSingleCatalogue $i "$filterCorrectionCoeff"
        # done
        echo "done" > $catDir/colourcorrection_done.txt
    fi
}
export -f applyColourcorrectionToAllCatalogues

applyColourcorrectionToSingleCatalogue() {
    local currentCatalogue=$1
    local coeffs=$2

    while read -r line; do
        if [[ "$line" =~ ^# ]]; then
            echo "$line"
            continue
        fi

        colour=$(echo "$line" | awk '{if (NF>=8) print $8; else print "nan"}')
        flux=$(echo "$line" | awk '{print $7}')
        magnitude=$(echo "$line" | awk '{print $6}')

        if [[ "$colour" == "nan" || -z "$colour" ]]; then
            echo "$line"
        else
            magnitudeCorrection=$( getColourCorrectionFromPolynomial "$coeffs" $colour )
            factorToApply=$( awk -v magCorr="$magnitudeCorrection" 'BEGIN {print 10^(-magCorr/2.5)}' )
            newMag=$( awk -v currentMag="$magnitude" -v magCorr="$magnitudeCorrection" 'BEGIN {print currentMag - magCorr}')
            newFlux=$( awk -v currentFlux="$flux" -v factor="$factorToApply" 'BEGIN {print currentFlux / factor}')
            new_line=$(echo "$line" | awk -v newMag="$newMag" -v newFlux="$newFlux" '{
                        $6 = newMag;  # Replace magnitude column
                        $7 = newFlux
                        print $0;
                    }')
            echo $new_line
        fi
    done < "$currentCatalogue" > "$currentCatalogue"_tmp
    mv "$currentCatalogue" "$currentCatalogue"_beforeColourcorrection
    mv "$currentCatalogue"_tmp $currentCatalogue
}
export -f applyColourcorrectionToSingleCatalogue

getColourCorrectionFromPolynomial() {
    local coeffs=$1
    local x=$2
    local power=0
    local result=0

    IFS=',' read -r -a coeffs_array <<< "$coeffs"
    x=$(printf "%f" "$x") # This is needed in case $x is given for transforming $x from scientific notation to normal, so bc can handle it

    for ((i=${#coeffs_array[@]}-1; i>=0; i--)); do
        coeff=${coeffs_array[i]}
        term=$(echo "$coeff * ($x^$power)" | bc -l)
        result=$(echo "$result + $term" | bc -l)
        ((power++))
    done
    echo $result
}

export -f getColourCorrectionFromPolynomial

computeColoursAndAddThemToCatalogues() {
    local cataloguesToUseDir=$1
    local cataloguesDir_g=$2
    local cataloguesDir_r=$3
    local filt=$4

    echo -e "\n·Computing and adding g-r colour to the catalogues\n"

    if [ -f "$cataloguesToUseDir/colourAdded_done.txt" ]; then
        echo -e "\n\tThe g-r colour has been already added to the catalogues\n"
    else
        for i in "$cataloguesDir_g"/*.g.cat; do
            fileName=$( basename "$i" .g.cat )

            if [ ! -f "$cataloguesToUseDir/$fileName.$filt.cat" ] || [ ! -f "$cataloguesDir_r/$fileName.r.cat" ]; then
                errorCode=999
                echo "ERROR - Catalogues in the different filters for calibration are not equal (missing a catalogue in some filter). We should never get there. Exiting with errorCode $errorCode"
                exit 999
            fi

            fileNameInitial=$cataloguesToUseDir/$fileName.$filt.cat
            fileNameWithColour=$cataloguesToUseDir/$fileName.$filt.cat_withColour

            fileName_g=$cataloguesDir_g/$fileName.g.cat
            fileName_r=$cataloguesDir_r/$fileName.r.cat
            python3 $pythonScriptsPath/getColoursFromTwoCataloguesAndAddItToThirdCatalogue.py $fileName_g $fileName_r $fileNameInitial $fileNameWithColour

            mv $cataloguesToUseDir/$fileName.$filt.cat $cataloguesToUseDir/$fileName.$filt.cat_beforeAddingColour
            mv $cataloguesToUseDir/$fileName.$filt.cat_withColour $cataloguesToUseDir/$fileName.$filt.cat
        done
        echo "done" > $cataloguesToUseDir/colourAdded_done.txt
    fi
}
export -f computeColoursAndAddThemToCatalogues


correctOffsetFromCatalogues() {
    local dirWithCatalogues=$1
    local offset=$2
    local factorToApplyToCounts=$3

    # The following awk command corrects the magnitude and the counts of all the catalogues
    for i in $( ls $dirWithCatalogues/*.cat ); do
        awk -v offset="$offset" -v factorToApplyToCounts="$factorToApplyToCounts" '
        /^#/ {print; next} 
        {
            if ($6 == "nan") {
                $6 = sprintf("%-9s", "nan")  # This is needed in order to have the nan and maintain the format (the 9 is for the widht of the field) 
                print;
            } else {
                $6 = sprintf("%.6f", $6 + offset); 
                $7 = sprintf("%.6f", $7 * factorToApplyToCounts)
                print;
            }
        }' $i > "$i"_tmp

        mv $i "${i%.txt}"_beforeGAIACalibration
        mv "$i"_tmp $i
    done
}
export -f correctOffsetFromCatalogues

calibrationToGAIA() {
    local spectraDir=$1
    local folderWithTransmittances=$2
    local filter=$3
    local ra=$4
    local dec=$5
    local mosaicDir=$6
    local sizeOfFieldForCalibratingPANSTARRStoGAIA=$7
    local magFromSpectraDir=$8
    local panstarrsCatalogueDir=$9

    local brightLimitToCompareGAIAandPANSTARRS=${10}
    local faintLimitToCompareGAIAandPANSTARRS=${11}


    if [ -f "$panstarrsCatalogueDir/mergedCatalogue.cat" ]; then
            echo -e "\n\tSurvey catalogues already merged for calibration with GAIA\n"
    else    
        listOfCatalogues=()
        for i in "$panstarrsCatalogueDir"/*.cat; do
            tmpName=$( basename $i )
            listOfCatalogues+=("${tmpName%.cat}")
        done 
        combineCatalogues $panstarrsCatalogueDir $panstarrsCatalogueDir "mergedCatalogue.cat" "${listOfCatalogues[@]}"
    fi

    transmittanceCurveFile="$folderWithTransmittances"/PANSTARRS_$filter.dat
    prepareSpectraDataForPhotometricCalibration $spectraDir $filter $ra $dec $mosaicDir $magFromSpectraDir $sizeOfFieldForCalibratingPANSTARRStoGAIA "GAIA" $transmittanceCurveFile
    offsetValues=$( python3 $pythonScriptsPath/getOffsetBetweenPANSTARRSandGAIA.py $panstarrsCatalogueDir/mergedCatalogue.cat $magFromSpectraDir/wholeFieldPhotometricCatalogue.cat $brightLimitToCompareGAIAandPANSTARRS $faintLimitToCompareGAIAandPANSTARRS $mosaicDir $filter )
    
    offset=$(echo $offsetValues | awk '{print $1}')
    factorToApplyToCounts=$(echo $offsetValues | awk '{print $2}')
    
    #rm $panstarrsCatalogueDir/mergedCatalogue.cat
    echo $offset $factorToApplyToCounts > $mosaicDir/offsetToCorrectSurveyToGaia_"$filter".txt
}
export -f calibrationToGAIA
    
# Photometric calibration functions
# The function that is to be used (the 'public' function using OOP terminology)
# Is 'computeCalibrationFactors' and 'applyCalibrationFactors'
selectStarsAndRangeForCalibrateSingleFrame(){
    local a=$1
    local framesForCalibrationDir=$2
    local mycatdir=$3
    local headerToUse=$4
    local methodToUse=$5         # This parameter will only be used if the catalogue is being generated with noisechisel
    local survey=$6
    local apertureUnits=$7
    local noisechisel_param=$8
    i=$framesForCalibrationDir/$a
    ##In the case of using it for Decals or Panstarrs, we need the variable survey
    
    if [[ "$methodToUse" == "sextractor" ]]; then
        outputCatalogue=$( generateCatalogueFromImage_sextractor $i $mycatdir $a $survey )
    elif [[ "$methodToUse" == "noisechisel" ]]; then
        outputCatalogue=$( generateCatalogueFromImage_noisechisel $i $mycatdir $a $headerToUse "'$noisechisel_param"  )
    else
        errorNumber=9
        echo "Error, method for selecting stars and the range in the calibration not recognised"
        echo "Exiting with error number: $erroNumber"
        exit $erroNumber
    fi
    
    if [[ "$survey" == "YES" ]]; then
        astmatch $outputCatalogue --hdu=1 $BDIR/catalogs/"$objectName"_Gaia_DR3.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aMAGNITUDE,a$apertureUnits -o$mycatdir/match_"$a"_my_gaia.txt
        s=$(asttable $mycatdir/match_"$a"_my_gaia.txt -h1 -c6 --noblank=MAGNITUDE   | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
        std=$(asttable $mycatdir/match_"$a"_my_gaia.txt -h1 -c6 --noblank=MAGNITUDE | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
        minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
        maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)
        echo $s $std $minr $maxr > $mycatdir/range_"$a".txt
        asttable $outputCatalogue --range=$apertureUnits,$minr,$maxr -o $mycatdir/selected_"$a"_automatic.txt
    else
        for h in $(seq 1 $num_ccd); do
            tmp_cat=$mycatdir/match_"$a"_ccd"$h"_my_gaia.fits
            astmatch $outputCatalogue --hdu=$h $BDIR/catalogs/"$objectName"_Gaia_DR3.fits --hdu2=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aMAGNITUDE,a$apertureUnits -o$tmp_cat
            astfits $tmp_cat --copy=1 -o $mycatdir/match_"$a"_my_gaia.fits
            rm $tmp_cat
        done
    
    
    # The intermediate step with awk is because I have come across an Inf value which make the std calculus fail
    # Maybe there is some beautiful way of ignoring it in gnuastro. I didn't find int, I just clean de inf fields.
        for h in $(seq 1 $num_ccd); do
            s=$(asttable $mycatdir/match_"$a"_my_gaia.fits -h$h -c6 --noblank=MAGNITUDE   | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=1,$iterationsForStdSigClip --sigclip-median)
            std=$(asttable $mycatdir/match_"$a"_my_gaia.fits -h$h -c6 --noblank=MAGNITUDE | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=1,$iterationsForStdSigClip --sigclip-std)
            
            minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
            maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)
            echo "$s $std $minr $maxr" >> $mycatdir/range_"$a".txt
            asttable $outputCatalogue -h$h --range=$apertureUnits,$minr,$maxr -o $mycatdir/selected_"$a"_ccd"$h"_automatic.fits
            astfits $mycatdir/selected_"$a"_ccd"$h"_automatic.fits --copy=1 -o$mycatdir/selected_"$a"_automatic.fits.cat
            rm $mycatdir/selected_"$a"_ccd"$h"_automatic.fits
            #mv $mycatdir/selected_"$a"_automatic.fits #$mycatdir/selected_"$a"_automatic.fits.cat
        done 
    fi
}
export -f selectStarsAndRangeForCalibrateSingleFrame

selectStarsAndSelectionRangeOurData() {
    local iteration=$1
    local framesForCalibrationDir=$2
    local mycatdir=$3
    local methodToUse=$4
    local apertureUnits=$5
    local noisechisel_param=$6
    
    mycatdone=$mycatdir/done.txt
    if ! [ -d $mycatdir ]; then mkdir $mycatdir; fi
    if [ -f $mycatdone ]; then
            echo -e "\n\tSources for photometric calibration are already extracted for my image\n"
    else
        framesToUse=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToUse+=("entirecamera_$a.fits")
        done

        headerWithData=1
        printf "%s\n" "${framesToUse[@]}" | parallel -j "$num_cpus" selectStarsAndRangeForCalibrateSingleFrame {} $framesForCalibrationDir $mycatdir $headerWithData $methodToUse NO $apertureUnits $noisechisel_param
        echo done > $mycatdone
    fi
}
export -f selectStarsAndSelectionRangeOurData


matchDecalsAndSingleFrame() {
    local a=$1
    local myCatalogues=$2
    local calibrationCatalogues=$3
    local matchdir=$4
    local surveyForCalibration=$5
    local calibratingMosaic=$6
    
    base=${a%.cat}
    ourDataCatalogue=$myCatalogues/$a
    out_cat=$matchdir/match-"$base".cat
    out_fits=$matchdir/match-"$base".fits
    
    for h in $(seq 1 $num_ccd); do
        tmpCatalogue=$matchdir/match-$base-tmp_ccd"$h".cat
        out_ccd=$matchdir/match-"$base"_ccd"$h".fits
        if [[ ($surveyForCalibration = "SPECTRA") || ("$calibratingMosaic" == true) ]]; then
            calibrationCatalogue=$calibrationCatalogues/wholeFieldPhotometricCatalogue.cat
            astmatch $ourDataCatalogue --hdu=$h $calibrationCatalogue --hdu2=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 \
                --outcols=bRA,bDEC,aRA,aDEC,bMAGNITUDE,bSUM,aMAGNITUDE,aSUM -o$tmpCatalogue

        else
            calibrationCatalogue=$calibrationCatalogues/"${base%%.*}".cat
    
            astmatch $ourDataCatalogue --hdu=$h $calibrationCatalogue --hdu2=$h --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 \
                --outcols=bRA,bDEC,aRA,aDEC,bMAGNITUDE,bSUM,aMAGNITUDE,aSUM -o$tmpCatalogue
        fi
        asttable $tmpCatalogue --output=$out_ccd --colmetadata=1,RA_CALIBRATED,deg,"Right ascension DECaLs" \
                --colmetadata=2,DEC_CALIBRATED,none,"Declination DECaLs" \
                --colmetadata=3,RA_NONCALIBRATED,deg,"Right ascension data being reduced" \
                --colmetadata=4,DEC_NONCALIBRATED,none,"Declination data being reduced" \
                --colmetadata=5,MAGNITUDE_CALIBRATED,none,"Magnitude in DECaLS data" \
                --colmetadata=6,SUM_CALIBRATED,none,"Sum in DECaLS" \
                --colmetadata=7,MAGNITUDE_NONCALIBRATED,none,"Magnitude in data being reduced" \
                --colmetadata=8,SUM_NONCALIBRATED,none,"Sum in in data being reduced" 
        astfits $out_ccd --copy=1 -o$out_fits
        rm $tmpCatalogue $out_ccd
    done
    mv $out_fits $out_cat

    
}
export -f matchDecalsAndSingleFrame

matchDecalsAndOurData() {
    local myCatalogues=$1
    local calibrationCatalogues=$2
    local matchdir=$3
    local surveyForCalibration=$4 
    local calibratingMosaic=$5
    
    matchdirdone=$matchdir/done_automatic.txt
    if ! [ -d $matchdir ]; then mkdir $matchdir; fi
    if [ -f $matchdirdone ]; then
        echo -e "\n\tMatch between decals (aperture) catalog and my (aperture) catalogs already done\n"
    else
        frameNumber=()
        for a in $( ls $myCatalogues/*.cat ); do
            frameName=$( basename $a )
            frameNumber+=("$frameName")
        done
        printf "%s\n" "${frameNumber[@]}" | parallel -j "$num_cpus" matchDecalsAndSingleFrame {} $myCatalogues $calibrationCatalogues $matchdir $surveyForCalibration $calibratingMosaic
        echo done > $matchdirdone
    fi
}
export -f matchDecalsAndOurData

buildOurCatalogueOfMatchedSourcesForFrame() {
    local a=$1
    local ourDatadir=$2
    local framesForCalibrationDir=$3
    local mycatdir=$4
    local numberOfFWHMToUse=$5

    base="entirecamera_$a.fits"
    i=$framesForCalibrationDir/$base
    automaticCatalogue=$mycatdir/selected_"$base"_automatic.fits.cat

    for h in $(seq 1 $num_ccd); do

        r_myData_pix_=$(awk 'NR=='$h' {printf $1}' $mycatdir/range_"$base".txt)
        r_myData_pix=$(astarithmetic $r_myData_pix_ $numberOfFWHMToUse. x -q )

        dataHdu=$h

    # raColumnName=RA
    # decColumnName=DEC
    # photometryOnImage_noisechisel $a $ourDatadir $automaticCatalogue $i $r_myData_pix $ourDatadir/$base.cat 22.5 $dataHdu \
    #                                 $raColumnName $decColumnName

        columnWithXCoordForOutDataPx=0 # These numbers come from how the catalogue of the matches stars is built. This is not very clear right now, should be improved
        columnWithYCoordForOutDataPx=1
        columnWithXCoordForOutDataWCS=2
        columnWithYCoordForOutDataWCS=3
        photometryOnImage_photutils $a $ourDatadir $automaticCatalogue $i $r_myData_pix $ourDatadir/"$base"_ccd"$dataHdu".fits 22.5 $dataHdu \
                                $columnWithXCoordForOutDataPx $columnWithYCoordForOutDataPx $columnWithXCoordForOutDataWCS $columnWithYCoordForOutDataWCS
        astfits $ourDatadir/"$base"_ccd"$dataHdu".fits --copy=1 -o$ourDatadir/"$base".cat
        rm $ourDatadir/"$base"_ccd"$dataHdu".fits
    done
}
export -f buildOurCatalogueOfMatchedSourcesForFrame

buildOurCatalogueOfMatchedSources() {
    local ourDatadir=$1
    local framesForCalibrationDir=$2
    local mycatdir=$3
    local numberOfFWHMToUse=$4

    ourDatadone=$ourDatadir/done.txt
    if ! [ -d $ourDatadir ]; then mkdir $ourDatadir; fi
    if [ -f $ourDatadone ]; then
        echo -e "\n\tAperture catalogs in our data done\n"
    else
        framesToUse=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToUse+=("$a")
        done
        printf "%s\n" "${framesToUse[@]}" | parallel -j "$num_cpus" buildOurCatalogueOfMatchedSourcesForFrame {} $ourDatadir $framesForCalibrationDir $mycatdir $numberOfFWHMToUse
        echo done > $ourDatadone
    fi
}
export -f buildOurCatalogueOfMatchedSources

matchCalibrationStarsCatalogues() {
    local matchdir2=$1
    local ourDatadir=$2
    local decalsdir=$3
    matchdir2done=$matchdir2/done_aperture.txt

    if [ -f $matchdir2done ]; then
        echo -e "\n\tMatch between decals (aperture) catalog and our (aperture) catalogs done for extension $h\n"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            base="entirecamera_$a.fits"
            i=$ourDatadir/"$base".cat

            out_tmp=$matchdir2/"$objectName"_Decals_"$a"_tmp.cat
            out=$matchdir2/"$objectName"_Decals-"$filter"_"$a".cat

            astmatch $decalsdir/decals_"$base".cat --hdu=1 $i --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aRA,aDEC,bRA,bDEC,aMAGNITUDE,aSUM,bMAGNITUDE,bSUM -o$out_tmp
            asttable $out_tmp --output=$out --colmetadata=1,RA,deg,"Right ascension DECaLs" \
                        --colmetadata=2,DEC,none,"Declination DECaLs" \
                        --colmetadata=3,RA,deg,"Right ascension data being reduced" \
                        --colmetadata=4,DEC,none,"Declination data being reduced" \
                        --colmetadata=5,MAGNITUDE_CALIBRATED,none,"Magnitude in DECaLS data" \
                        --colmetadata=6,SUM,none,"Sum in DECaLS" \
                        --colmetadata=7,MAGNITUDE_NONCALIBRATED,none,"Magnitude in data being reduced" \
                        --colmetadata=8,SUM,none,"Sum in in data being reduced" 
            rm $out_tmp
        done
        echo done > $matchdir2done
    fi
}
export -f matchCalibrationStarsCatalogues

computeAndStoreFactors() {
    local alphatruedir=$1
    local matchdir=$2
    local brightLimit=$3
    local faintLimit=$4

    alphatruedone=$alphatruedir/done.txt
    numberOfStarsUsedToCalibrateFile=$alphatruedir/numberOfStarsUsedForCalibrate.txt

    if ! [ -d $alphatruedir ]; then mkdir $alphatruedir; fi
    if [ -f $alphatruedone ]; then
        echo -e "\n\tTrustable alphas computed for extension $h\n"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            
            f=$matchdir/match-entirecamera_$a.fits.cat

            
            for h in $(seq 1 $num_ccd); do
                alphaFile=alpha_"$a"_ccd"$h".txt
                alphatruet=$alphatruedir/"$objectName"_"$filter"_"$a"_ccd"$h".txt
                asttable $f -h$h --range=MAGNITUDE_CALIBRATED,$brightLimit,$faintLimit -o$alphatruet
                asttable $alphatruet -h1 -c1,2,'arith $6 $8 /' -o$alphatruedir/$alphaFile

                mean=$(asttable $alphatruedir/$alphaFile -c'ARITH_1' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-mean)
                std=$(asttable $alphatruedir/$alphaFile -c'ARITH_1' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
                if [ $mean == "n/a" ]; then
                    mean=$(asttable $alphatruedir/$alphaFile -c'ARITH_1' | aststatistics --sigclip-mean)
                fi
                ###This dirty if is connected to the following: our sigma clipping params are able to avoid negative values or strange values of alpha, but
                # for some reason on tables with low rows it gets "nan" (I know, stupid right?). Because of that, we decide to avoid the parameters when we get nan values and this
                # apparently solves the problem
                echo "$mean $std" >> $alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a".txt
                count=$(asttable $alphatruedir/$alphaFile -c'ARITH_1' | aststatistics --number)
                echo "Frame number $a ccd $h : $count" >> $numberOfStarsUsedToCalibrateFile
            done
        done
        echo done > $alphatruedone
    fi
}
export -f computeAndStoreFactors

combineCatalogues() {
    local outputDir=$1
    local cataloguesDir=$2
    local frame=$3
    shift 3
    bricks=("$@")

  
    catalogueName=$( echo "$frame" | awk -F'.' '{print $1}')

    # firstBrick=$(echo "$bricks" | awk '{print $1}')  
    
    # asttablePrompt="asttable $cataloguesDir/$firstBrick.cat -o$outputDir/$catalogueName.cat"

    # remainingBricks=$(echo "$bricks" | cut -d' ' -f2-)  # Get the rest of the bricks
    # for brick in $remainingBricks; do
    #     asttablePrompt+=" --catrowfile=$cataloguesDir/$brick.cat"
    # done

    firstBrick=${bricks[0]}
   
    asttablePrompt="asttable $cataloguesDir/$firstBrick.cat -o$outputDir/$catalogueName.cat"
    

    remainingBricks=$(echo "$bricks" | cut -d' ' -f2-)  # Get the rest of the bricks
    for brick in "${bricks[@]}"; do  # Iterate over remaining bricks
        asttablePrompt+=" --catrowfile=$cataloguesDir/$brick.cat"
    done

    $asttablePrompt
}
export -f combineCatalogues

combineDecalsCataloguesForSingleFrame() {
    local outputDir=$1
    local frame=$2
    local h=$3
    local bricks=$4

    catalogueName=$( echo "$frame" | awk -F'.' '{print $1}')
    firstBrick=$(echo "$bricks" | awk '{print $1}')  
    
    
    asttablePrompt="asttable $decalsCataloguesDir/$firstBrick.cat -o$outputDir/"$catalogueName"_ccd"$h".fits"

    remainingBricks=$(echo "$bricks" | cut -d' ' -f2-)  # Get the rest of the bricks
    for brick in $remainingBricks; do
        asttablePrompt+=" --catrowfile=$decalsCataloguesDir/$brick.cat"
    done
    $asttablePrompt
}

combineDecalsBricksCataloguesForEachFrame() {
    local outputDir=$1
    local frameBrickAssociationFile=$2
    local decalsCataloguesDir=$3

    combinationDone=$outputDir/done.txt
    if ! [ -d $outputDir ]; then mkdir $outputDir; fi
    if [ -f $combinationDone ]; then
        echo -e "\nCombination of the bricks catalogues for each frame already done\n"
    else
        for h in $(seq 1 $num_ccd); do
            frameBrickAssociationFile_ccd="$frameBrickAssociationFile"_ccd"$h".txt
            while IFS= read -r line; do
                currentLine=$line
                frame=$(echo "$line" | awk '{print $1}')
                frame=$( basename $frame )
                bricks=$(echo "$line" | cut -d' ' -f2-)
                combineDecalsCataloguesForSingleFrame $outputDir $frame $h "$bricks"
            done < "$frameBrickAssociationFile_ccd"
        done
        for a in $(seq 1 $totalNumberOfFrames); do
            for h in $(seq 1 $num_ccd); do
                input_cat=$outputDir/entirecamera_"$a"_ccd"$h".fits
                astfits $input_cat --copy=1 -o$outputDir/entirecamera_"$a".fits
                rm $input_cat
            done
            mv $outputDir/entirecamera_"$a".fits $outputDir/entirecamera_"$a".cat
        done
        echo "done" > $combinationDone
    fi
}
export -f combineDecalsBricksCataloguesForEachFrame

computeCommonCalibrationFactor() {
  local calibrationFactorsDir=$1
  local iteration=$2
  local objectName=$3
  local BDIR=$4 
  
  for h in $(seq 1 $num_ccd); do
    calibrationFactors=()
    for i in $( ls $calibrationFactorsDir/alpha_"$objectName"*.txt); do
        alpha=$(awk 'NR=='$h'{print $1}' $i)
        calibrationFactors+=("$alpha")
    done
    tmpTableFits=$BDIR/tableTest.fits
    printf "%s\n" "${calibrationFactors[@]}" | asttable -o "$tmpTableFits"
    commonCalibrationFactor=$( aststatistics $tmpTableFits --sigclip-median)
    calibrationFactorsStd=$( asttable $tmpTableFits | aststatistics --sclipparams=3,3 --sigclip-std)

    echo $commonCalibrationFactor $calibrationFactorsStd >> $BDIR/commonCalibrationFactor_it$iteration.txt
    rm $tmpTableFits
  done
}
export -f computeCommonCalibrationFactor

computeCalibrationFactors() {
    local surveyForCalibration=$1
    local iteration=$2
    local imagesForCalibration=$3
    local selectedDecalsStarsDir=$4
    local matchdir=$5
    local ourDataCatalogueDir=$6
    local prepareCalibrationCataloguePerFrame=$7
    local mycatdir=$8
    local rangeUsedDecalsDir=$9
    local mosaicDir=${10}
    local alphatruedir=${11}
    local brightLimit=${12}
    local faintLimit=${13}
    local apertureUnits=${14}
    local numberOfApertureUnitsForCalibration=${15}
    local calibratingMosaic=${16}
    local noisechisel_param=${17}

    

    methodToUse="sextractor"
    echo -e "\n ${GREEN} ---Selecting stars and range for our data--- ${NOCOLOUR}"
    selectStarsAndSelectionRangeOurData $iteration $imagesForCalibration $mycatdir $methodToUse $apertureUnits "'$noisechisel_param'"
     
    
    echo -e "\n ${GREEN} ---Building catalogues to our data with aperture photometry --- ${NOCOLOUR}"
    buildOurCatalogueOfMatchedSources $ourDataCatalogueDir $imagesForCalibration $mycatdir $numberOfApertureUnitsForCalibration
      
    # If we are calibrating with spectra we just have the whole catalogue of the field
    # If we are calibrating with a survey then we have a catalogue por survey's brick and we need to combine the needed bricks for build a catalogue per frame
    if ! [ -d $prepareCalibrationCataloguePerFrame ]; then mkdir $prepareCalibrationCataloguePerFrame; fi
    if [[ ("$surveyForCalibration" == "SPECTRA") || ( "$calibratingMosaic" == true) ]]; then
        cp $mosaicDir/wholeFieldPhotometricCatalogue.cat $prepareCalibrationCataloguePerFrame
    else
        echo -e "\n ${GREEN} ---Combining decals catalogues for matching each brick --- ${NOCOLOUR}"
        combineDecalsBricksCataloguesForEachFrame $prepareCalibrationCataloguePerFrame $mosaicDir/frames_bricks_association $mosaicDir/aperturePhotometryCatalogues
        
    fi
    
    echo -e "\n ${GREEN} ---Matching our aperture catalogues and Decals aperture catalogues--- ${NOCOLOUR}"
    matchDecalsAndOurData $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $matchdir $surveyForCalibration $calibratingMosaic
    
    echo -e "\n ${GREEN} ---Computing calibration factors (alpha)--- ${NOCOLOUR}"
    computeAndStoreFactors $alphatruedir $matchdir $brightLimit $faintLimit
}
export -f computeCalibrationFactors

applyCalibrationFactorsToFrame() {
    local a=$1
    local imagesForCalibration=$2
    local alphatruedir=$3
    local photCorrDir=$4
    local iteration=$5
    local applyCommon=$6
    local inverseApplication=$7

    base=entirecamera_"$a".fits
    f=$imagesForCalibration/"entirecamera_$a.fits"
    if [[ "$applyCommon" == "true" || "$applyCommon" == "True" ]]; then
        alpha_cat=$BDIR/commonCalibrationFactor_it"$iteration".txt
    else
        alpha_cat=$alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a".txt
    fi
    for h in $(seq 1 $num_ccd); do
        base_ccd=entirecamera_"$a"_ccd"$h".fits
        alpha=$(awk 'NR=='$h'{print $1}' $alpha_cat)
        if [[ "$inverseApplication" == "TRUE" ]]; then
            astarithmetic $f -h$h $alpha / float32 -o $photCorrDir/$base_ccd
        else
            astarithmetic $f -h$h $alpha x float32 -o $photCorrDir/$base_ccd
        fi
        astfits $photCorrDir/$base_ccd --copy=1 -o$photCorrDir/$base
        rm $photCorrDir/$base_ccd
    done
}
export -f applyCalibrationFactorsToFrame

applyCalibrationFactors() {
    local imagesForCalibration=$1
    local alphatruedir=$2
    local photCorrDir=$3
    local iteration=$4
    local applyCommon=$5
    local inverseApplication=${6:-false}

    muldone=$photCorrDir/done.txt
    if ! [ -d $photCorrDir ]; then mkdir $photCorrDir; fi
    if [ -f $muldone ]; then
            echo -e "\n\tMultiplication for alpha in the pointings (huge grid) is done\n"
    else
        framesToApplyFactor=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToApplyFactor+=("$a")
        done
        printf "%s\n" "${framesToApplyFactor[@]}" | parallel -j "$num_cpus" applyCalibrationFactorsToFrame {} $imagesForCalibration $alphatruedir $photCorrDir $iteration $applyCommon $inverseApplication
        echo done > $muldone
    fi
}
export -f applyCalibrationFactors

# Compute the weights o the frames based on the std of the background
# In order to perform a weighted mean
computeWeightForFrame() {
    local a=$1
    local wdir=$2
    local wonlydir=$3
    local photCorrDir=$4
    local noiseskydir=$5 
    local iteration=$6
    local minRmsFileName=$7

    
    base=entirecamera_"$a".fits
    basetmp=entirecamera_"$a"_tmp.fits

    f=$photCorrDir/$base
    rms_min=$(awk 'NR=='1'{print $1}' $BDIR/$minRmsFileName)
    for h in $(seq 1 $num_ccd); do
        rms_f=$(awk 'NR=='$h'{print $3}' $noiseskydir/entirecamera_$a.txt)

    # ****** Decision note *******
    # The weights are obtained as the quadratic ratio between the best sigma and the current sigma
    # This weights produce the optimal combinantion for a gaussian distribution 
    # Ref: https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
        weight=$(astarithmetic $rms_min 2 pow $rms_f 2 pow / --quiet) 
        echo "$weight" >> $wdir/"$objectName"_Decals-"$filter"_"$a".txt      
        
    # multiply each image for its weight
        wixi_im_tmp=$wdir/$basetmp              # frame x weight
        w_im_tmp=$wonlydir/$basetmp             # only weight
        wixi_im=$wdir/ccd"$h"_$base                     # frame x weight
        w_im=$wonlydir/ccd"$h"_$base                    # only weight

        astarithmetic $f -h$h $weight x --type=float32 -o$wixi_im_tmp 
        astarithmetic $wixi_im_tmp -h1 $f -h$h / --type=float32 -o$w_im_tmp
        astarithmetic $wixi_im_tmp float32 -g1 -o$wixi_im
        astarithmetic $w_im_tmp float32 -g1 -o$w_im
        astfits $wixi_im --copy=1 -o$wdir/$base
        astfits $w_im --copy=1 -o$wonlydir/$base
        rm -f $wixi_im_tmp $wixi_im
        rm -f $w_im_tmp $w_im
    done
    
}
export -f computeWeightForFrame

computeWeights() {
    local wdir=$1
    local wdone=$2
    local wonlydir=$3
    local wonlydone=$4
    local photCorrDir=$5
    local noiseskydir=$6 
    local iteration=$7
    local minRmsFileName=$8

    if [ -f $wdone ]; then
        echo -e "\n\tWeights computation done for extension $h\n"
    else
        framesToComputeWeight=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToComputeWeight+=("$a")
        done
        printf "%s\n" "${framesToComputeWeight[@]}" | parallel -j "$num_cpus" computeWeightForFrame {} $wdir $wonlydir $photCorrDir $noiseskydir $iteration $minRmsFileName
        echo done > $wdone
        echo done > $wonlydone
    fi
}
export -f computeWeights

# Outliers functions
buildUpperAndLowerLimitsForOutliers() {
    local clippingdir=$1
    local clippingdone=$2
    local wdir=$3
    local sigmaForStdSigclip=$4

    if ! [ -d $clippingdir ]; then mkdir $clippingdir; fi
    if [ -f $clippingdone ]; then
            echo -e "\n\tUpper and lower limits for building the masked of the weighted images already computed\n"
    else
            # Compute clipped median and std
            med_im=$clippingdir/median_image.fits
            std_im=$clippingdir/std_image.fits
            gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
            file_count=0
            output_list=""
            for file in $(ls -v $wdir/*.fits); do
                for h in $(seq 1 $num_ccd); do
                    output_list+="$file -h$h "
                    ((file_count++))
                done
            done

            if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
                
                astarithmetic $output_list $file_count $sigmaForStdSigclip 0.2 sigclip-median  --writeall -o$med_im
                astarithmetic $output_list $file_count $sigmaForStdSigclip 0.2 sigclip-std  --writeall -o$std_im
            else
                astarithmetic $output_list $file_count $sigmaForStdSigclip 0.2 sigclip-median   -o$med_im
                astarithmetic $output_list $file_count $sigmaForStdSigclip 0.2 sigclip-std   -o$std_im
            fi
            
            # Compute "borders" images
            up_lim=$clippingdir/upperlim.fits
            lo_lim=$clippingdir/lowerlim.fits
            astarithmetic 4. $std_im x -o thresh.fits
            astarithmetic $med_im thresh.fits + -g1 float32 -o $up_lim
            astarithmetic $med_im thresh.fits - -g1 float32 -o $lo_lim

            #rm -f $med_im $std_im
            rm thresh.fits
            echo done > $clippingdone
    fi
}
export -f buildUpperAndLowerLimitsForOutliers

removeOutliersFromFrame(){
    local a=$1
    local mowdir=$2
    local clippingdir=$3
    local wdir=$4

    base=$( basename $a )
    tmp_ab=$mowdir/${base%.fits}_maskabove.fits
    wom=$mowdir/$base
    
    for h in $(seq 1 $num_ccd); do
        wom_ccd=$mowdir/${base%.fits}_ccd"$h".fits
        astarithmetic $wdir/$base -h$h set-i i i $clippingdir/upperlim.fits -h1 gt nan where float32 -q -o $tmp_ab
        astarithmetic $tmp_ab -h1 set-i i i $clippingdir/lowerlim.fits -h1 lt nan where float32 -q -o$wom_ccd
        astfits $wom_ccd --copy=1 -o$wom
    # save the new mask
        #mask=$mowdir/"$objectName"_Decals-"$filter"_"$a"_ccd"$h"_mask.fits
        #astarithmetic $wom_ccd -h1 isblank float32 -o $mask
    # mask the onlyweight image
        #owom=$moonwdir/$base
        #owom_ccd=$moonwdir/ccd_$base
        #astarithmetic $wonlydir/$base -h$h $mask -h1 1 eq nan where -q float32    -o $owom_ccd
        #astfits $owom_ccd --copy=1 -o$owom

    # Remove temporary files
        rm -f $tmp_ab
        #rm -f $mask
        rm -f $wom_ccd
    done
}
export -f removeOutliersFromFrame

removeOutliersFromWeightedFrames () {
  local mowdone=$1
  local mowdir=$2
  local clippingdir=$3
  local wdir=$4

  if [ -f $mowdone ]; then
      echo -e "\n\tOutliers of the weighted images already masked\n"
  else
      framesToRemoveOutliers=()
      for a in $wdir/*.fits; do
          framesToRemoveOutliers+=("$a")
      done
      printf "%s\n" "${framesToRemoveOutliers[@]}" | parallel -j "$num_cpus" removeOutliersFromFrame {} $mowdir $clippingdir $wdir
      echo done > $mowdone 
  fi
}
export -f removeOutliersFromWeightedFrames

# Functions for applying the mask of the coadd for a second iteration
cropAndApplyMaskPerFrame() {
    local a=$1
    local dirOfFramesToMask=$2
    local dirOfFramesMasked=$3
    local wholeMask=$4
    local dirWithCropParameters=$5


    frameToMask=$dirOfFramesToMask/entirecamera_$a.fits
    fileWithCropParameters=$dirWithCropParameters/entirecamera_"$a"_cropRegion.txt
      
    
    for h in $(seq 1 $num_ccd); do
    # Parameters for identifing our frame in the full grid
        tmpMaskFile=$dirOfFramesMasked/maskFor"$a"_ccd"$h".fits
        read row_min row_max col_min col_max < <(sed -n "${h}p" "$fileWithCropParameters")
        astcrop $wholeMask --polygon=$col_min,$row_min:$col_max,$row_min:$col_max,$row_max:$col_min,$row_max --mode=img  -o $tmpMaskFile --quiet
        astarithmetic $frameToMask -h$h $tmpMaskFile -h1 1 eq nan where float32 -o $dirOfFramesMasked/entirecamera_"$a"_ccd"$h".fits -q
        astfits $dirOfFramesMasked/entirecamera_"$a"_ccd"$h".fits --copy=1 -o $dirOfFramesMasked/entirecamera_"$a".fits
        rm $tmpMaskFile $dirOfFramesMasked/entirecamera_"$a"_ccd"$h".fits
    done
}
export -f cropAndApplyMaskPerFrame

# maskPointings receives the directory with the frames in the full grid because we need it in order to know the region of the full grid
# in which the specific frame is located. That is obtained by using getRegionToCrop.py frame
maskPointings() {
    local entiredir_smallGrid=$1
    local smallPointings_maskedDir=$2
    local maskedPointingsDone=$3
    local maskName=$4
    local dirWithCropParameters=$5

    if ! [ -d $smallPointings_maskedDir ]; then mkdir $smallPointings_maskedDir; fi
    if [ -f $maskedPointingsDone ]; then
            echo -e "\nThe masks for the pointings have been already applied\n"
    else
        framesToMask=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToMask+=("$a")
        done
        printf "%s\n" "${framesToMask[@]}" | parallel -j "$num_cpus" cropAndApplyMaskPerFrame {} $entiredir_smallGrid $smallPointings_maskedDir $maskName $dirWithCropParameters
        echo done > $maskedPointingsDone 
    fi
}
export -f maskPointings

produceAstrometryCheckPlot() {
    local matchCataloguesDir=$1
    local pythonScriptsPath=$2
    local output=$3
    local pixelScale=$4
    local h=$5

    astrometryTmpDir="./astrometryDiagnosisTmp"
    if ! [ -d $astrometryTmpDir ]; then mkdir $astrometryTmpDir; fi
    python3 $pythonScriptsPath/diagnosis_deltaRAdeltaDEC.py $matchCataloguesDir $output $pixelScale $h
    rm -rf $astrometryTmpDir
}
export -f produceAstrometryCheckPlot

produceCalibrationCheckPlot() {
    local myCatalogue_nonCalibrated=$1
    local myFrames_calibrated=$2
    local aperturesForMyData_dir=$3
    local referenceCatalogueDir=$4
    local pythonScriptsPath=$5
    local output=$6
    local calibrationBrightLimit=$7
    local calibrationFaintLimit=$8
    local numberOfApertureUnitsForCalibration=$9
    local outputDir=${10}
    local survey=${11}
    local bdir=${12}
    
    calibratedCataloguesDir=$bdir/calibratedCatalogues
    
    if ! [ -d $calibratedCataloguesDir ]; then mkdir $calibratedCataloguesDir; fi
    #tmpDir="./calibrationDiagnosisTmp"
    #if ! [ -d $tmpDir ]; then mkdir $tmpDir; fi

    for i in $myCatalogue_nonCalibrated/*.cat; do
        myFrame=$i
        frameNumber=$(echo "$i" | awk -F '[/]' '{print $(NF)}' | awk -F '[.]' '{print $(1)}' | awk -F '[_]' '{print $(NF)}')

        # In the nominal resolution it takes sooo long for doing this plots. So only a set of frames are used for the
        # calibration check
        if [ $survey == "SPECTRA" ]; then
            referenceCatalogue=$referenceCatalogueDir/wholeFieldPhotometricCatalogue.cat
        else
            referenceCatalogue=$referenceCatalogueDir/*_$frameNumber.*
        fi

        myCalibratedFrame=$myFrames_calibrated/entirecamera_$frameNumber.fits
        myNonCalibratedCatalogue=$myCatalogue_nonCalibrated/entirecamera_$frameNumber.fits*
        fileWithMyApertureData=$aperturesForMyData_dir/range_entirecamera_$frameNumber*
        for h in $(seq 1 $num_ccd); do
                r_myData_pix_=$(awk 'NR=='$h' {printf $1}' $fileWithMyApertureData)
                r_myData_pix=$(astarithmetic $r_myData_pix_ $numberOfApertureUnitsForCalibration. x -q )

            # raColumnName=RA
            # decColumnName=DEC
            # photometryOnImage_noisechisel -1 $tmpDir $myNonCalibratedCatalogue $myCalibratedFrame $r_myData_pix $tmpDir/$frameNumber.cat 22.5 \
            #                                 $raColumnName $decColumnName
                dataHdu=$h
                columnWithXCoordForOutDataPx=1 # These numbers come from how the catalogue of the matches stars is built. This is not very clear right now, should be improved
                columnWithYCoordForOutDataPx=2
                columnWithXCoordForOutDataWCS=3
                columnWithYCoordForOutDataWCS=4
                photometryOnImage_photutils -1 $calibratedCataloguesDir $myNonCalibratedCatalogue $myCalibratedFrame $r_myData_pix $calibratedCataloguesDir/$frameNumber.cat 22.5 $dataHdu \
                                        $columnWithXCoordForOutDataPx $columnWithYCoordForOutDataPx $columnWithXCoordForOutDataWCS $columnWithYCoordForOutDataWCS

                astmatch $referenceCatalogue --hdu=$h $calibratedCataloguesDir/$frameNumber.cat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=1/3600 --outcols=aMAGNITUDE,bMAGNITUDE -o$calibratedCataloguesDir/"$frameNumber"_temp.fits
                asttable $calibratedCataloguesDir/"$frameNumber"_temp.fits --colmetadata=1,MAGNITUDE_CALIBRATED --colmetadata=2,MAGNITUDE_NONCALIBRATED -o$calibratedCataloguesDir/"$frameNumber"_matched_ccd"$h".fits
                astfits $calibratedCataloguesDir/"$frameNumber"_matched_ccd"$h".fits --copy=1 -o$calibratedCataloguesDir/"$frameNumber"_matched.fits
                rm $calibratedCataloguesDir/"$frameNumber"_matched_ccd"$h".fits $calibratedCataloguesDir/"$frameNumber"_temp.fits
        done
        
        mv $calibratedCataloguesDir/"$frameNumber"_matched.fits $calibratedCataloguesDir/"$frameNumber"_matched.cat
        rm $calibratedCataloguesDir/$frameNumber.cat
        
    done
    for h in $(seq 1 $num_ccd); do
        python3 $pythonScriptsPath/diagnosis_magVsDeltaMag.py $calibratedCataloguesDir $output $outputDir $calibrationBrightLimit $calibrationFaintLimit $survey $h
    done
    
}
export -f produceCalibrationCheckPlot

produceHalfMaxRadVsMagForSingleImage() {
    local image=$1 
    local outputDir=$2
    local gaiaCat=$3
    local tolerance=$4
    local pythonPath=$5
    local alternativeIdentifier=$6 # Applied when there is no number in the name
    local tileSize=$7
    local survey=$8
    a=$( echo $image | grep -oP '\d+(?=\.fits)' )
    if ! [[ -n "$a" ]]; then
        a=$alternativeIdentifier
    fi

    # header=1
    # catalogueName=$(generateCatalogueFromImage_noisechisel $image $outputDir $a $headerToUse $tileSize)
    catalogueName=$(generateCatalogueFromImage_sextractor $image $outputDir $a $survey)
    for h in $(seq 1 $num_ccd); do
        astmatch $catalogueName --hdu=$h $gaiaCat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$tolerance/3600 --outcols=aX,aY,aRA,aDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o $outputDir/match_decals_gaia_"$a"_ccd"$h".txt 
    
        plotXLowerLimit=0.5
        plotXHigherLimit=10
        plotYLowerLimit=12
        plotYHigherLimit=22
        python3 $pythonPath/diagnosis_halfMaxRadVsMag.py $catalogueName $outputDir/match_decals_gaia_"$a"_ccd"$h".txt -1 -1 -1 $outputDir/$a.png  \
            $plotXLowerLimit $plotXHigherLimit $plotYLowerLimit $plotYHigherLimit $h

        rm $catalogueName $outputDir/match_decals_gaia_"$a"_ccd"$h".txt 
    done
}
export -f produceHalfMaxRadVsMagForSingleImage


produceHalfMaxRadVsMagForOurData() {
    local imagesDir=$1
    local outputDir=$2
    local gaiaCat=$3
    local toleranceForMatching=$4
    local pythonScriptsPath=$5
    local num_cpus=$6
    local tileSize=$7

    images=()
    for i in $imagesDir/*.fits; do
        images+=("$i")
    done

    # images=("/home/sguerra/NGC598/build/photCorrSmallGrid-dir_it1/entirecamera_1.fits")
    printf "%s\n" "${images[@]}" | parallel --line-buffer -j "$num_cpus" produceHalfMaxRadVsMagForSingleImage {} $outputDir $gaiaCat $toleranceForMatching $pythonScriptsPath "-" $tileSize NO
}
export -f produceHalfMaxRadVsMagForOurData

buildCoadd() {
    local coaddir=$1
    local coaddName=$2
    local mowdir=$3
    local moonwdir=$4
    local coaddone=$5

    if ! [ -d $coaddir ]; then mkdir $coaddir; fi
    if [ -f $coaddone ]; then
            echo -e "\n\tThe first weighted (based upon std) mean of the images already done\n"
    else
            gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
            output_mow=""
            file_count_mow=0
            for file in $(ls -v $mowdir/*.fits); do
                for h in $(seq 1 $num_ccd); do
                    output_mow+="$file -h$h "
                    ((file_count_mow++))
                done
            done
            output_moonw=""
            file_count_moonw=0
            for file in $(ls -v $moonwdir/*.fits); do
                for h in $(seq 1 $num_ccd); do
                    output_moonw+="$file -h$h "
                    ((file_count_moonw++))
                done
            done
            if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
                astarithmetic $output_mow $file_count_mow sum  --writeall -o$coaddir/"$k"_wx.fits
                astarithmetic $output_moonw $file_count_moonw sum  --writeall -o$coaddir/"$k"_w.fits
            else
                astarithmetic $output_mow $file_count_mow sum   -o$coaddir/"$k"_wx.fits
                astarithmetic $output_moonw $file_count_moonw sum   -o$coaddir/"$k"_w.fits
            fi
            astarithmetic $coaddir/"$k"_wx.fits -h1 $coaddir/"$k"_w.fits -h1 / -o$coaddName
            echo done > $coaddone
    fi
}
export -f buildCoadd

subtractCoaddToFrames() {
    local dirWithFrames=$1
    local coadd=$2
    local destinationDir=$3

    for i in $dirWithFrames/*.fits; do
        name=$( basename $i )
        for h in $(seq 1 $num_ccd); do
            temp_file=$destinationDir/temp_$name
            astarithmetic $i -h$h $coadd -h1 - -o$temp_file
            astfits $temp_file --copy=1 -o $destinationDir/$name
            rm $temp_file
        done
    done
}
export -f subtractCoaddToFrames

changeNonNansOfFrameToOnes() {
  local a=$1
  local framesDir=$2
  local outputDir=$3

  frame=$framesDir/entirecamera_$a.fits
  output=$outputDir/exposure_tmp_$a.fits
  for h in $(seq 1 $num_ccd); do
    astarithmetic $frame $frame isblank not 1 where --output=$outputDir/exposure_tmp_ccd_$a.fits -g$h
    astfits $outputDir/exposure_tmp_ccd_$a.fits --copy=1 -o$output
    rm $outputDir/exposure_tmp_ccd_$a.fits
  done
}
export -f changeNonNansOfFrameToOnes

computeExposureMap() {
    local framesDir=$1
    local exposureMapDir=$2
    local exposureMapDone=$3
    local exposureMapName=$4

    if ! [ -d $exposureMapDir ]; then mkdir $exposureMapDir; fi
    if [ -f $exposureMapDone ]; then
        echo -e "\n\tThe exposure map is already done\n"
    else
      
      framesToProcess=()
      for a in $(seq 1 $totalNumberOfFrames); do
        framesToProcess+=("$a")
      done
      
      printf "%s\n" "${framesToProcess[@]}" | parallel -j "$num_cpus" changeNonNansOfFrameToOnes {} $framesDir $exposureMapDir
      gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
      exposure_frames=""
      file_count=0
      for file in $(ls -v $exposureMapDir/*.fits); do
        for h in $(seq 1 $num_ccd); do
            exposure_frames+="$file -h$h "
            ((file_count++))
        done
      done
      if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
        astarithmetic $exposure_frames $file_count sum   --writeall -o$exposureMapName
      else
        astarithmetic $exposure_frames $file_count sum   --writeall -o$exposureMapName
      fi
      rm -rf $exposureMapDir
      echo done > $exposureMapDone
    fi
}
export -f computeExposureMap






# ------------------------------------------------------------------------------------------------------
# Functions for generating the catalogue of an image
# For the moment these are in the calibration and when we need to make a half-max-radius vs magnitude plot

# In order to be used interchangeable, they have to return a catalogue with the following columns:
# 1.- X
# 2.- Y
# 3.- RA
# 4.- DEC
# 3.- Magnitude
# 4.- FWHM/2

generateCatalogueFromImage_noisechisel() { 
    local image=$1
    local outputDir=$2
    local a=$3   
    local header=$4
    local noisechisel_param=$5

    astmkprof --kernel=gaussian,1.5,3 --oversample=1 -o $outputDir/kernel_$a.fits 1>/dev/null
    astconvolve $image -h$header --kernel=$outputDir/kernel_$a.fits --domain=spatial --output=$outputDir/convolved_$a.fits 1>/dev/null
    astnoisechisel $image -h$header -o $outputDir/det_$a.fits --convolved=$outputDir/convolved_$a.fits --tilesize=$tileSize,$tileSize --numthreads=$num_cpus 1>/dev/null
    astsegment $outputDir/det_$a.fits -o $outputDir/seg_$a.fits --gthresh=-15 --objbordersn=0 1>/dev/null
    astmkcatalog $outputDir/seg_$a.fits --x --y --ra --dec --magnitude --half-max-radius --sum --clumpscat -o $outputDir/decals_$a.txt --zeropoint=22.5 1>/dev/null

    rm $outputDir/kernel_$a.fits $outputDir/convolved_$a.fits $outputDir/det_$a.fits $outputDir/seg_$a.fits  $outputDir/decals_"$a"_o.txt
    mv $outputDir/decals_"$a"_c.txt $outputDir/catalogue_$a.cat
    echo $outputDir/catalogue_$a.cat
}
export -f generateCatalogueFromImage_noisechisel


generateCatalogueFromImage_sextractor(){
    local image=$1
    local outputDir=$2
    local a=$3   
    local survey=$4
    # I specify the configuration path here because in the photometric calibration the working directoy changes. This has to be changed and use the config path given in the pipeline
    cfgPath=$ROOTDIR/"$objectName"/config
    #For panstarrs we need to use zp=25
   
    source-extractor $image -c $cfgPath/sextractor_detection.sex -CATALOG_NAME $outputDir/"$a"_tmp.cat -FILTER_NAME $cfgPath/default.conv -PARAMETERS_NAME $cfgPath/sextractor_detection.param   1>/dev/null 2>&1
    #awk '{ $6 = $6 / 2; print }' $outputDir/"$a"_tmp.cat > $outputDir/"$a".cat # I divide because SExtractor gives the FWHM and the pipeline expects half
    #If multidetector the .cat output file is a num_ccd layer .cat with the catalogues of each layer of the .fits file. We make the arith and save everything properly
    if [[ "$survey" == "YES" ]]; then
        tmp_file=$outputDir/"$a".fits
        asttable $outputDir/"$a"_tmp.cat -h1 -c1,2,3,4,5,'arith $6 2 /',7 --colmetadata=1,X --colmetadata=2,Y --colmetadata=3,RA --colmetadata=4,DEC --colmetadata=5,MAGNITUDE --colmetadata=6,FWHM,pixel --colmetadata=7,Re,pixel -o$tmp_file
    else
        for h in $(seq 1 $num_ccd); do
            tmp_file=$outputDir/"$a"_ccd"$h"_tmp.fits
            asttable $outputDir/"$a"_tmp.cat -h$h -c1,2,3,4,5,'arith $6 2 /',7 --colmetadata=1,X --colmetadata=2,Y --colmetadata=3,RA --colmetadata=4,DEC --colmetadata=5,MAGNITUDE --colmetadata=6,FWHM,pixel --colmetadata=7,Re,pixel -o$tmp_file
            astfits $tmp_file --copy=1 -o$outputDir/"$a".fits
        done
        rm $outputDir/"$a"_ccd*_tmp.fits
    fi
    
    # Headers to mimic the noisechisel format. Change between MacOS and Linux
    #With the previous solution all of this is unnecessary
    #if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS (BSD sed)
    #    sed -i '' '1s/^/# Column 1: X\
# Column 2: Y\
# Column 3: RA\
# Column 4: DEC\
# Column 5: MAGNITUDE\
# Column 6: HALF_MAX_RADIUS\
#/' "$outputDir/$a.cat"
#    else
#        sed -i '1i# Column 6: HALF_MAX_RADIUS' $outputDir/$a.cat
#        sed -i '1i# Column 5: MAGNITUDE      ' $outputDir/$a.cat
#        sed -i '1i# Column 4: DEC            ' $outputDir/$a.cat
#        sed -i '1i# Column 3: RA             ' $outputDir/$a.cat
#        sed -i '1i# Column 2: Y              ' $outputDir/$a.cat
#        sed -i '1i# Column 1: X              ' $outputDir/$a.cat
#    fi

    rm $outputDir/"$a"_tmp.cat 

    mv  $outputDir/$a.fits $outputDir/catalogue_$a.cat
    echo $outputDir/catalogue_$a.cat
}
export -f generateCatalogueFromImage_sextractor



# ------------------------------------------------------------------------------------------------------
# Functions for performing photometry in an image
# They produce an output catalogue with the format
# 1.- ID
# 2.- X
# 3.- Y
# 4.- RA
# 5.- DEC
# 6.- MAGNITUDE
# 7.- SUM

photometryOnImage_noisechisel() {
    local a=$1
    local directoryToWork=$2
    local matchedCatalogue=$3
    local imageToUse=$4
    local aperture_radius_px=$5
    local outputCatalogue=$6
    local zeropoint=$7
    local hduWithData=$8
    local raColumnName=$9
    local decColumnName=${10}

    asttable $matchedCatalogue -hSOURCE_ID -c$raColumnName,$decColumnName | awk '!/^#/{print NR, $1, $2, 5, '$aperture_radius_px', 0, 0, 1, NR, 1}' > $directoryToWork/apertures_decals_$a.txt
    astmkprof $directoryToWork/apertures_decals_$a.txt --background=$imageToUse --backhdu=$hduWithData \
            --clearcanvas --replace --type=int16 --mforflatpix \
            --mode=wcs --output=$directoryToWork/apertures_decals_$a.fits
    astmkcatalog $directoryToWork/apertures_decals_$a.fits -h1 --zeropoint=$zeropoint \
                    --valuesfile=$imageToUse --valueshdu=$hduWithData \
                    --ids --x --y --ra --dec --magnitude --sum \
                    --output=$outputCatalogue

    rm $directoryToWork/apertures_decals_$a.txt
    rm $directoryToWork/apertures_decals_$a.fits
}
export -f photometryOnImage_noisechisel

photometryOnImage_photutils() {
    local a=$1
    local directoryToWork=$2
    local matchedCatalogue=$3
    local imageToUse=$4
    local aperture_radius_px=$5
    local outputCatalogue=$6
    local zeropoint=$7
    local hduWithData=$8
    local xColumnPx=$9
    local yColumnPx=${10}
    local xColumnWCS=${11}
    local yColumnWCS=${12}

    tmpCatalogName=$directoryToWork/tmp_"$a".cat
    python3 $pythonScriptsPath/photutilsPhotometry.py $matchedCatalogue $imageToUse $aperture_radius_px $tmpCatalogName $zeropoint $hduWithData $xColumnPx $yColumnPx $xColumnWCS $yColumnWCS

    asttable $tmpCatalogName -p4 --colmetadata=2,X,px,"X" \
                            --colmetadata=3,Y,px,"Y" \
                            --colmetadata=4,RA,deg,"Right ascension" \
                            --colmetadata=5,DEC,none,"Declination" \
                            --colmetadata=6,MAGNITUDE,none,"Magnitude" \
                            --colmetadata=7,SUM,none,"sum" \
                            --output=$outputCatalogue
    rm $tmpCatalogName
}
export -f photometryOnImage_photutils

limitingSurfaceBrightness() {
    local image=$1
    local mask=$2
    local exposureMap=$3
    local directoryOfImages=$4
    local areaSB=$5
    local fracExpMap=$6
    local pixelScale=$7
    local outFile=$8

    numOfSigmasForMetric=3

    out_mask=$directoryOfImages/mask_det.fits
    astarithmetic $image -h1 $mask -h1 0 ne nan where -q --output=$out_mask >/dev/null 2>&1

    out_maskexp=$directoryOfImages/mask_exp.fits
    expMax=$(aststatistics $exposureMap --maximum -q)
    exp_fr=$(astarithmetic $expMax $fracExpMap x -q)
    astarithmetic $out_mask $exposureMap -g1 $exp_fr lt nan where --output=$out_maskexp >/dev/null 2>&1
    zp_asec=$(astarithmetic $pixelScale log10 5 x 22.5 + -q)
    sigma=$(aststatistics $out_maskexp --sigclip-std -q)
    
    sb_lim=$(astarithmetic $sigma 3 x $pixelScale x $areaSB / log10 -2.5 x $zp_asec + -q)
    echo "$sb_lim" > "$outFile"

    rm $out_mask $out_maskexp
    echo "Limiting magnitude ($numOfSigmasForMetric sigma, $areaSB x $areaSB): $sb_lim" > "$outFile"
    echo "$sb_lim" # We need to recover the value outside for adding it to the coadd header
}
export -f limitingSurfaceBrightness

prepareRingTemplate(){
    local template=$1
    local imageDir=$2
    local ringDir=$3
    x_camera=$(awk '{print $2}' $template)
    y_camera=$(awk '{print $3}' $template)
    radius=$(awk '{print $5}' $template)
    oneImage=$(ls $imageDir/*.fit* | head -1)
    
    for h in $(seq 1 $num_ccd); do
        eval $(python3 $pythonScriptsPath/singleCameratoCCD_cordtransform.py $oneImage $h $x_camera $y_camera)
        base="${template%.txt}"
        out_ring=$ringDir/"$base"_ccd"$h".txt
        echo "1 $x_det $y_det 6 $radius 1 1 1 1 1" > $out_ring
    done

}
export -f prepareRingTemplate

smallGridToFullGridSingleFrame(){
    local smallFrame=$1
    local fullDir=$2
    local fullSize=$3
    local fullRA=$4
    local fullDEC=$5

    base=$( basename $smallFrame )
    out=$fullDir/$base
    astfits $smallFrame --copy=0 --primaryimghdu -o $out
    for h in $(seq 1 $num_ccd); do
        astcrop $smallFrame -h$h --mode=wcs --center=$fullRA,$fullDEC --widthinpix --width=$fullSize,$fullSize --zeroisnotblank -o $fullDir/ccd_$base
        astfits $fullDir/ccd_$base --copy=1 -o $out
        rm $fullDir/ccd_$base
    done
}
export -f smallGridToFullGridSingleFrame

smallGridtoFullGrid(){
    local smallGridDir=$1
    local fullGridDir=$2
    local fullGridDone=$3
    local fullGridSize=$4
    local fullGridRA=$5
    local fullGridDEC=$6
    

    if [ -f $fullGridDone ]; then
        echo -e "\n\tFrames from {$smallGridDir} have been already pased into the full grid\n"
    else
        framesToGrid=()
        for frame in $smallGridDir/*.fits; do
            framesToGrid+=("$frame")
        done
        printf "%s\n" "${framesToGrid[@]}" | parallel -j "$num_cpus" smallGridToFullGridSingleFrame {} $fullGridDir $fullGridSize $fullGridRA $fullGridDEC 
        echo done > $fullGridDone
        
    fi
}
export -f smallGridtoFullGrid

computeFWHMSingleFrame(){
    local a=$1
    local framesForFWHMDir=$2
    local fwhmdir=$3
    local headerToUse=$4
    local methodToUse=$5
    local tileSize=$6           # This parameter will only be used if the catalogue is being generated with noisechisel
    local survey=$7
    
    i=$framesForFWHMDir/entirecamera_"$a"
    ##In the case of using it for Decals or Panstarrs, we need the variable survey
    
    if [[ "$methodToUse" == "sextractor" ]]; then
        outputCatalogue=$( generateCatalogueFromImage_sextractor $i $fwhmdir $a $survey )
    elif [[ "$methodToUse" == "noisechisel" ]]; then
        outputCatalogue=$( generateCatalogueFromImage_noisechisel $i $fwhmdir $a $headerToUse $tileSize  )
    else
        errorNumber=9
        echo "Error, method for selecting stars and the range in the calibration not recognised"
        echo "Exiting with error number: $erroNumber"
        exit $erroNumber
    fi
    
    if [[ "$survey" == "YES" ]]; then
        astmatch $outputCatalogue --hdu=1 $BDIR/catalogs/"$objectName"_Gaia_DR3.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aMAGNITUDE,aFWHM -o$fwhmdir/match_"$a"_my_gaia.txt
        hFWHM=$(asttable $fwhmdir/match_"$a"_my_gaia.txt -h1 -c6 --noblank=MAGNITUDE   | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
        FWHM=$(awk "BEGIN {print $hFWHM * 2}")
        echo $FWHM > $fwhmdir/fwhm_"$a".txt
        rm $fwhmdir/match_"$a"_my_gaia.txt
    else
        for h in $(seq 1 $num_ccd); do
            tmp_cat=$fwhmdir/match_"$a"_ccd"$h"_my_gaia.fits
            astmatch $outputCatalogue --hdu=$h $BDIR/catalogs/"$objectName"_Gaia_DR3.fits --hdu2=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aMAGNITUDE,a$apertureUnits -o$tmp_cat
            astfits $tmp_cat --copy=1 -o $fwhmdir/match_"$a"_my_gaia.fits
            rm $tmp_cat
        done
    
    
    # The intermediate step with awk is because I have come across an Inf value which make the std calculus fail
    # Maybe there is some beautiful way of ignoring it in gnuastro. I didn't find int, I just clean de inf fields.
        for h in $(seq 1 $num_ccd); do
            hFWHM=$(asttable $fwhmdir/match_"$a"_my_gaia.fits -h$h -c6 --noblank=MAGNITUDE   | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
            FWHM=$(awk "BEGIN {print $hFWHM * 2}")
            echo "$FWHM" >> $fwhmdir/fwhm_"$a".txt
            #mv $mycatdir/selected_"$a"_automatic.fits #$mycatdir/selected_"$a"_automatic.fits.cat
        done 
        rm $fwhmdir/match_"$a"_my_gaia.fits
    fi
    rm $outputCatalogue
}
export -f computeFWHMSingleFrame

computeGainCorrectionSingleFrame(){
    local image=$1
    local outputDir=$2
    local num_ccd=$3
    local ringDir=$4
    base=$( basename $image )
    tmpRingFits=$(echo $base | sed 's/.fits/_ring.fits/')
    tmpRingFits_single=$(echo $base | sed 's/.fits/_ring_single.fits/')
    out=$(echo $base | sed 's/.fits/.txt/')
    for h in $(seq 1 $num_ccd); do
        x_ring=$( awk ' {print $2}' $ringDir/ring_ccd"$h".txt )
        y_ring=$( awk ' {print $3}' $ringDir/ring_ccd"$h".txt )
        tmpRingDefinition=$(echo $base | sed 's/.fits/_ring_ccd.txt/')
        ringRadius=$( awk '{print $5}' $ringDir/ring_ccd"$h".txt )
        naxis1=$(fitsheader $image -e $h | grep "NAXIS1" | awk '{print $3'})
        naxis2=$(fitsheader $image -e $h | grep "NAXIS2" | awk '{print $3'})
        naxis1_r=$(fitsheader $ringDir/ring.fits -e $h | grep "NAXIS1" | awk '{print $3'})
        naxis2_r=$(fitsheader $ringDir/ring.fits -e $h | grep "NAXIS2" | awk '{print $3'})
        if [[ $naxis1 -gt $naxis2 && $naxis1_r -gt $naxis2_r ]] || [[ $naxis1 -lt $naxis2 && $naxis1_r -lt $naxis2_r ]]; then
            echo "1 $x_ring $y_ring 6 $ringRadius 1 1 1 1 1" > $ringDir/$tmpRingDefinition
        else
            #If image new is rotated, we look for the astrometrized, not rotated in order to get the correct position
            image_astro=${base#entirecamera_}
            ringCentre=$( xy2sky $BDIR/astro-ima/$image_astro,$h $x_ring $y_ring )
            ringRa=$(echo "$ringCentre" | awk '{print $1}')
            ringDec=$(echo "$ringCentre" | awk '{print $2}')
            newringCentre=$( sky2xy $image,$h $ringRa $ringDec )
            x_new=$(echo "$newringCentre" | awk '{print $5}')
            y_new=$(echo "$newringCentre" | awk '{print $6}')
            echo "1 $x_new $y_new 6 $ringRadius 1 1 1 1 1" > $ringDir/$tmpRingDefinition
                   
                        
        fi
        astmkprof --background=$image --backhdu=$h --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas --quiet -o $ringDir/$tmpRingFits_single $ringDir/$tmpRingDefinition
        rm -f $ringDir/$tmpRingDefinition
        astfits $ringDir/$tmpRingFits_single --copy=1 -o $ringDir/$tmpRingFits
    done
    
    rm -f $ringDir/$tmpRingFits_single    
    base=$( basename $image )
    for h in $(seq 1 $num_ccd); do

        astnoisechisel $image -h$h --tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2 --rawoutput --numthreads=$num_cpus -o $outputDir/tmp_$base
        astarithmetic $image -h$h $outputDir/tmp_$base -h1 0 ne nan where -q -o$outputDir/tmp2_$base
        astarithmetic $outputDir/tmp2_$base -h1 $ringDir/$tmpRingFits -h$h 0 ne nan where -q -o$outputDir/tmp3_$base
        mean_g=$(aststatistics $outputDir/tmp3_$base --sigclip-mean -q)
        echo "$mean_g" >> $outputDir/$out
        rm $outputDir/tmp*_$base
    done
    rm $ringDir/$tmpRingFits


}
export -f computeGainCorrectionSingleFrame

computeGainCorrection(){
    local inputDir=$1
    local outputDir=$2
    local num_ccd=$3
    local ringDir=$4

    if ! [ -d $outputDir ]; then mkdir $outputDir; fi
    framesToComputeGain=()
    for a in $inputDir/*.fits; do
        framesToComputeGain+=("$a")
    done
    printf "%s\n" "${framesToComputeGain[@]}" | parallel -j "$num_cpus" computeGainCorrectionSingleFrame {} $outputDir $num_ccd $ringDir
    
}
export -f computeGainCorrection


subtractStars(){
    local inputFolder_small=$1
    local starLine=$2
    local psf=$3
    local psfProfile=$4
    local outputDir_small=$5
    local starId=$6
    
    starRa=$(echo "$starLine" | awk '{print $1}')
    starDec=$(echo "$starLine" | awk '{print $2}')
    starMag=$(echo "$starLine" | awk '{print $3}')
    echo -e "Star $starId: RA=$starRa DEC=$starDec MAG=$starMag "
    echo -e "·Locating star in CCDs and computing profiles"
    profileFolder=$BDIR/profileStar_"$starId"
    profileDone=$profileFolder/done.txt
    starLocationFile=$profileFolder/starlocation.txt
    if ! [ -d $profileFolder ]; then mkdir $profileFolder; fi
    if [ -f $profileDone ]; then
        echo -e "Profiles of star {$starId} already computed"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            image=$inputFolder_small/entirecamera_"$a".fits
            #We are gonna do the following: we generate circles where PSF reaches 26.5 mag arcsec^-2, ie, ~50 times lower than typical
            # sky background (THis is literally by eye)
            circleRad=$(python3 $pythonScriptsPath/get_cropRadiusPSF.py $starMag $psfProfile 26.5 $pixelScale) 
            python3 $pythonScriptsPath/check_starLocation.py $image $starLocationFile $num_ccd $starRa $starDec $circleRad 
        done
        while IFS= read -r line; do
            #On the config file I'm gonna add some TXT where I will put if an star needs azimuth or not
            #It is true that is not the best solution, but I don't think we will need it to much
            
            imageProf=$(echo "$line" | awk '{print $1}')
            imageProf=$( basename $imageProf )
            rp_params="$line --mode=wcs --center=$starRa,$starDec --rmax=1200 --undersample=5 --measure=sigclip-mean,sigclip-std "
            starAz_file=$CDIR/star"$starId"_az.txt
            if [ -f $starAz_file ]; then
                starAz=$(awk 'NR=='1'{print $1}' $starAz_file)
                rp_params+="--azimuth=$starAz "
            fi
            astscript-radial-profile $rp_params --output=$profileFolder/RP_$imageProf
        done < $starLocationFile
        echo done > $profileDone
    fi

    echo -e "·Computing Scale Factor"
    ###NOTE: the fitting algorithm benefits from paralelization inside python, instead of parallelizing python execution
    scaleDir=$BDIR/scaleStar_$starId
    scaleDone=$scaleDir/done.txt
    if ! [ -d $scaleDir ]; then mkdir $scaleDir; fi
    if [ -f $scaleDone ]; then
        echo -e "Scale factors already computed for Star {$starId}"
    else
        ###If there is not star, we put the scale in 0
        scale_text_base=$scaleDir/scale_entirecamera
        for a in $(seq 1 $totalNumberOfFrames); do
            profile=$profileFolder/RP_entirecamera_"$a".fits
            scale_text="$scale_text_base"_"$a".txt
            if ! [ -f $profile ]; then
                echo "0.000" > $scale_text
            else
                
                ###First step: measure a rough approach of the background, using CCD2
                out_sky=$scaleDir/sky_$a.fits
                astnoisechisel $inputFolder_small/entirecamera_"$a".fits -h2  --tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2 --checksky --numthreads=$num_cpus -o $out_sky
                out_sky=${out_sky%.fits}_sky.fits
                sky_mean=$(aststatistics $out_sky -hSKY --sigclip-mean)
                ##We will range ±500
                scale=$(python3 $pythonScriptsPath/get_fitWithPSF.py $profile $psfProfile $starMag $sky_mean $scaleDir $num_cpus $pixelScale)
                echo "$scale" > $scale_text
                rm $out_sky
            fi
        done
    
        echo done > $scaleDone
    fi

    echo -e "·Subtracting star from frames"
    limMag=32 #bellow surface brightness limit

    if ! [ -d $outputDir_small ]; then mkdir $outputDir_small; fi
    #if ! [ -d $outputDir_full ]; then mkdir $outputDir_full; fi
    subtractionDone=$outputDir_small/done.txt
    if [ -f $subtractionDone ]; then
        echo -e "Subtraction of Star $starId already done"
    else
        ##First step: create the PSFfile
        psfCrop=$outputDir_small/PSF.fits
        rCrop=$(python3 $pythonScriptsPath/get_cropRadiusPSF.py $starMag $psfProfile $limMag $pixelScale )
        if (( $(echo "$rCrop == 0" | bc -l) )); then
            cp $psfFile $psfCrop
        else
            rWidth=$((2*rCrop))
            astcrop $psfFile --mode=img --center=8001,8001 --width=$rWidth,$rWidth --zeroisnotblank -o$outputDir_small/temp.fits
            echo "1 $rCrop $rCrop 5 $rCrop 0.0 0 1 1 1" | astmkprof --background=$outputDir_small/temp.fits --mforflatpix --clearcanvas -o$outputDir_small/mask.fits
            astarithmetic $outputDir_small/temp.fits $outputDir_small/mask.fits -g1 0 eq nan where -q -o$psfCrop
            rm $outputDir_small/temp.fits $outputDir_small/mask.fits
        fi
        
        
        for a in $(seq 1 $totalNumberOfFrames); do
           subtractStarFromFrame $a $psfCrop $scaleDir $inputFolder_small $outputDir_small $starRa $starDec 
           #subtractStarFromFrame $a $psfCrop $scaleDir $inputFolder_full $outputDir_full $starRa $starDec
        done
        rm $psfCrop
        echo done > $subtractionDone
    fi

}
export -f subtractStars
subtractStarFromFrame() {
    local a=$1
    local psf=$2
    local scale_folder=$3
    local input_folder=$4
    local output_folder=$5
    local star_ra=$6
    local star_dec=$7

    base=entirecamera_"$a"
    image=$input_folder/"$base".fits
    scale_text=$scale_folder/scale_"$base".txt
    output=$output_folder/"$base".fits
    scale=$(awk 'NR=='1'{print $1}' $scale_text)
    if (( $(echo "$scale == 0" | bc -l) )); then
        cp $image $output
    else
        astfits $image --copy=0 --primaryimghdu -o$output
        for h in $(seq 1 4); do
        

            output_temp=$output_folder/"$base"_temp.fits
            astscript-psf-subtract $image -h$h --scale=$scale --mode=wcs --center=$star_ra,$star_dec --psf=$psf -o $output_temp --numthreads=$num_cpus
            astfits $output_temp --copy=1 -o$output
            gain=$(astfits $image -h$h --keyvalue=GAIN -q)
            astfits $output -h$h --write=GAIN,$gain
            rm $output_temp    
        done
    fi
}
export -f subtractStarFromFrame

