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
        "·Root directory to perform the reduction:$ROOTDIR"
        "·Keyword for the airmass:$airMassKeyWord"
        "·Keyword for date:$dateHeaderKey"
        "·Keyword for RA of the pointings:$pointingRA:[$pointingRAUnits]"
        "·Keyword for DEC of the pointings:$pointingDEC:[$pointingDECUnits]"
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
        " "
        "·Indices scales for astrometrisation"
        "  Lowest index:$lowestScaleForIndex"
        "  Highest index:$highestScaleForIndex"
        " "
        "·Scale-low parameters for solve-field (astrometry.net):$solve_field_L_Param"
        "·Scale-high parameters for solve-field (astrometry.net):$solve_field_H_Param"
        "·Scale units for parameters for solve-field (astrometry.net):$solve_field_u_Param"
        "·Survey for building the catalogue introduced to solvefield:$surveyToUseInSolveField"
        " "
        "·Filter:$filter"
        "·Pixel scale:$pixelScale:[arcsec/px]"
        "·Detector width:$detectorWidth:[px]"
        "·Detector height:$detectorHeight:[px]"
        "·Is there overscan:$overscan"
        "·Keyword for illuminated section:$trimsecKey"
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
    printf "\t%-70s $ORANGE %-20s $GREEN %-10s $NOCOLOUR\n" "$text" "$value" "$unit"
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
                telescope \
                ra_gal \
                dec_gal \
                telescopeLat \
                telescopeLong \
                telescopeElevation \
                defaultNumOfCPUs \
                ROOTDIR \
                airMassKeyWord \
                dateHeaderKey \
                pointingRA \
                pointingRAUnits \
                pointingDEC \
                pointingDECUnits \
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
                filter \
                pixelScale \
                detectorWidth \
                detectorHeight \
                overscan \
                trimsecKey \
                lowestScaleForIndex \
                highestScaleForIndex \
                solve_field_L_Param \
                solve_field_H_Param \
                solve_field_u_Param \
                surveyToUseInSolveField \
                maximumBackgroundBrightness \
                maximumSeeing \
                fractionExpMap \
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

    if [[ ("$surveyForPhotometry" != "PANSTARRS") && ("$surveyForPhotometry" != "DECaLS") && ("$surveyForPhotometry" != "SDSS") && ("$surveyForPhotometry" != "SPECTRA")]]; then
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
    local filterFileNeeded=$( checkIfAllTheTransmittancesNeededAreGiven $telescope $surveyForPhotometry $folderWithTransmittances $filter )

    if [[ $filterFileNeeded == "11_1" ]]; then
        echo "Error. Calibrating with SPECTRA option but not found the transmittance of the filter to calibrate"
        exit 11
    elif [[ $filterFileNeeded == "11_2" ]]; then
        echo "Error. Calibrating with IMAGING SURVEY ($survey) option but not found the transmittance of the filter to calibrate this survey to our reference (GAIA)"
        exit 11
    else
       checkUnitsAndConvertToCommonUnitsIfNeeded $filterFileNeeded
    fi

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


# In this function we check the transmittances given. Just for having everything in the same units
# We place everything into Angstroms and transmittances from 0 to 1.
checkUnitsAndConvertToCommonUnitsIfNeeded() {
    local filterFile=$1
    errorCode=12

    # Since we may change the file (depending on the units) I create a backup in order not to loose the original file
    cp $filterFile $filterFile.original

    read units transmittanceFormat < <(head -n 1 $filterFile)
    transmittanceFormat="${transmittanceFormat//$'\r'/}" # Remove the final endline

    if [[ ("$units" != "A") && ("$units" != "nm" ) ]]; then
        echo "format ($transmittanceFormat) for transmittance not accepted. It is expected either A or nm"
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


    # If we calibrate with spectra (survey="SPECTRA") then only the transmittance of the filter to use is needed
    if [[ "$survey" == "SPECTRA" ]]; then
        filterPath=$filterFolder/"$telescope"_"$filterToUse".dat
        if ! [ -e $filterPath ]; then
            errorCode=11_1
            echo $errorCode
        else
            echo $filterPath
        fi
    else
        # If calibrating with imaging survey then the filter from the survey is needed (to calibrate it to GAIA).
        # No filter from the telescope of the data to reduce is needed because the imaging survey will be used
        filterPath=$filterFolder/"$survey"_"$filterToUse".dat
        if ! [ -e $filterPath ]; then
            errorCode=11_2
            echo $errorCode
        else
            echo $filterPath
        fi
    fi
}
export -f checkIfAllTheTransmittancesNeededAreGiven

# Functions used in Flat
maskImages() {
    local inputDirectory=$1
    local masksDirectory=$2
    local outputDirectory=$3
    local useCommonRing=$4
    local keyWordToDecideRing=$5

    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        i=$inputDirectory/$base
        out=$outputDirectory/$base
        astarithmetic $i -h1 $masksDirectory/$base -hDETECTIONS 1 eq nan where float32 -o $out -q
       
        propagateKeyword $i $dateHeaderKey $out
        propagateKeyword $i $airMassKeyWord $out
        propagateKeyword $i $dateHeaderKey $out 
        # If we are not doing a normalisation with a common ring we propagate the keyword that will be used to decide
        # which ring is to be used. This way we can check this value in a comfortable way in the normalisation section
        if [ "$useCommonRing" = false ]; then
            propagateKeyword $i $keyWordToDecideRing $out
        fi
    done
}
export -f maskImages

getInitialMidAndFinalFrameTimes() {
  local directoryWithNights=$1

  declare -a date_obs_array

  while IFS= read -r -d '' file; do
    currentDateObs=$( gethead $file DATE-OBS)
    if [[ -n "$currentDateObs" ]]; then
        unixTime=$( TZ=UTC date -d "$currentDateObs" +"%s" )
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

    valueToPropagate=$(astfits $image --keyvalue=$keyWordToPropagate --quiet)
    eval "astfits --delete=$keyWordToPropagate $out -h1 2&>/dev/null" # I redirect the error descriptor so I avoid the error message if the keyword didn't exist
    eval "astfits --write=$keyWordToPropagate,$valueToPropagate $out -h1"
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
    local noiseskydir=${10}


    if [ "$useCommonRing" = true ]; then
        # Case when we have one common normalisation ring
        tmpName=$( basename $i ) 
        maskedImageUsingRing=$noiseskydir/maskedAllButRing$tmpName
        astarithmetic $i -h1 $commonRing -h1 0 eq nan where -o $maskedImageUsingRing --quiet
        me=$(aststatistics $maskedImageUsingRing -h1 --sclipparams=3,3 --sigclip-median --quiet)
        rm $maskedImageUsingRing
    else
        # Case when we do NOT have one common normalisation ring
        # All the following logic is to decide which normalisation ring apply
        variableToDecideRingToNormalise=$(gethead $i $keyWordToDecideRing)
        firstRingLowerBound=$(echo "$keyWordValueForFirstRing - $keyWordThreshold" | bc)
        firstRingUpperBound=$(echo "$keyWordValueForFirstRing + $keyWordThreshold" | bc)
        secondRingLowerBound=$(echo "$keyWordValueForSecondRing - $keyWordThreshold" | bc)
        secondRingUpperBound=$(echo "$keyWordValueForSecondRing + $keyWordThreshold" | bc)

        tmpName=$( basename $i ) 
        maskedImageUsingRing=$noiseskydir/maskedAllButRing$tmpName

        if (( $(echo "$variableToDecideRingToNormalise >= $firstRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $firstRingUpperBound" | bc -l) )); then
            astarithmetic $i -h1 $doubleRing_first -h1 0 eq nan where -o $maskedImageUsingRing --quiet
            me=$(aststatistics $maskedImageUsingRing -h1 --sclipparams=3,3 --sigclip-median --quiet)
        elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
            astarithmetic $i -h1 $doubleRing_second -h1 0 eq nan where -o $maskedImageUsingRing --quiet
            me=$(aststatistics $maskedImageUsingRing -h1 --sclipparams=3,3 --sigclip-median --quiet)
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
    local noiseskydir=${10}

    if [ "$useCommonRing" = true ]; then
        # Case when we have one common normalisation ring
        tmpName=$( basename $i ) 
        maskedImageUsingRing=$noiseskydir/maskedAllButRing$tmpName
        astarithmetic $i -h1 $commonRing -h1 0 eq nan where -o $maskedImageUsingRing --quiet
        std=$(aststatistics $maskedImageUsingRing -h1 --sclipparams=3,3 --sigclip-std --quiet)
        rm $maskedImageUsingRing
    else
        # Case when we do NOT have one common normalisation ring
        # All the following logic is to decide which normalisation ring apply
        variableToDecideRingToNormalise=$(gethead $i $keyWordToDecideRing)
        firstRingLowerBound=$(echo "$keyWordValueForFirstRing - $keyWordThreshold" | bc)
        firstRingUpperBound=$(echo "$keyWordValueForFirstRing + $keyWordThreshold" | bc)
        secondRingLowerBound=$(echo "$keyWordValueForSecondRing - $keyWordThreshold" | bc)
        secondRingUpperBound=$(echo "$keyWordValueForSecondRing + $keyWordThreshold" | bc)

        tmpName=$( basename $i ) 
        maskedImageUsingRing=$noiseskydir/maskedAllButRing$tmpName

        if (( $(echo "$variableToDecideRingToNormalise >= $firstRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $firstRingUpperBound" | bc -l) )); then
            astarithmetic $i -h1 $doubleRing_first -h1 0 eq nan where -o $maskedImageUsingRing --quiet
            std=$(aststatistics $maskedImageUsingRing -h1 --sclipparams=3,3 --sigclip-std --quiet)
        elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
            astarithmetic $i -h1 $doubleRing_second -h1 0 eq nan where -o $maskedImageUsingRing --quiet
            std=$(aststatistics $maskedImageUsingRing -h1 --sclipparams=3,3 --sigclip-std --quiet)
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

getSkewKurtoValueFromSkyPixels(){
    local i=$1
    local constantSkyMethod=$2
    local commonRing=$3
    local doubleRing_first=$4
    local doubleRing_second=$5
    local useCommonRing=$6
    local keyWordToDecideRing=$7
    local keyWordThreshold=$8
    local keyWordValueForFirstRing=$9
    local keyWordValueForSecondRing=${10}

    if [ "$constantSkyMethod" == "ring" ]; then
        if [ "$useCommonRing" = true ]; then
                skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS $commonRing)
                kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS $commonRing)
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
                skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS $doubleRing_first)
                kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS $doubleRing_first)
                # rm ring_masked.fits
            elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
                #astarithmetic $i -h1 $doubleRing_second -h1 0 eq nan where -q -o ring_masked.fits
                skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS $doubleRing_second)
                kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS $doubleRing_second)
                # rm ring_masked.fits
            else
                errorNumber=5
                echo -e "\nMultiple normalisation ring have been tried to be used. The keyword selection value of one has not matched with the ranges provided" >&2
                echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
                exit $errorNumber
            fi
        fi
    elif [ "$constantSkyMethod" == "wholeImage" ]; then
        skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS)
        kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS)
    elif [ "$constantSkyMethod" == "noisechisel" ]; then
        skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS)
        kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS)
    else
        errorNumber=555
        echo -e "\nIn Function getSkewKurtoValueFromSkyPixels. Not identified the "constantSkyMethod" variable  value" >&2
        echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
        exit $errorNumber  
    fi
    echo "$skew $kurto"
}
export -f getSkewKurtoValueFromSkyPixels

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

    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        i=$imageDir/$base
        out=$outputDir/$base

        me=$(getMedianValueInsideRing $i $commonRing  $doubleRing_first $doubleRing_second $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $outputDir)
        astarithmetic $i -h1 $me / -o $out

        propagateKeyword $i $dateHeaderKey $out
        propagateKeyword $i $airMassKeyWord $out
        propagateKeyword $i $dateHeaderKey $out
    done
}
export -f normaliseImagesWithRing

calculateFlat() {
    local flatName="$1"
    shift
    local filesToUse="$@"
    numberOfFiles=$#

    # ****** Decision note *******
    # The rejection parameters for the construction of the flat has been chosen to be 2 sigmas
    # The running flat implies that we only have fewer frames for the flat (in our case 11 for example)
    # So we have to be a little bit aggresive in order to be able to remove the outliers
    sigmaValue=2
    iterations=10
    gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
    if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
        astarithmetic $filesToUse $numberOfFiles $sigmaValue $iterations sigclip-median -g1 --writeall -o $flatName
    else
        astarithmetic $filesToUse $numberOfFiles $sigmaValue $iterations sigclip-median -g1 -o $flatName
    fi
}
export -f calculateFlat

calculateRunningFlat() {
    local normalisedDir=$1
    local outputDir=$2
    local doneFile=$3
    local iteration=$4

    windowSize=$(( (halfWindowSize * 2) + 1 ))

    fileArray=()
    fileArray=( $(ls -v $normalisedDir/*Decals-"$filter"_n*_f*_ccd"$h".fits) )
    fileArrayLength=( $(ls -v $normalisedDir/*Decals-"$filter"_n*_f*_ccd"$h".fits | wc -l) )

    lefFlatFiles=("${fileArray[@]:0:$windowSize}")
    echo "Computing left flat - iteration $iteration"
    calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_left_ccd"$h".fits" "${lefFlatFiles[@]}"
    propagateKeyword "${lefFlatFiles[$halfWindowSize]}" $dateHeaderKey "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_left_ccd"$h".fits"


    rightFlatFiles=("${fileArray[@]:(fileArrayLength-$windowSize):fileArrayLength}")
    echo "Computing right flat - iteration $iteration"
    calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_right_ccd"$h".fits" "${rightFlatFiles[@]}"
    propagateKeyword "${rightFlatFiles[$halfWindowSize]}" $dateHeaderKey "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_right_ccd"$h".fits"


    echo "Computing non-common flats - iteration $iteration"
    for a in $(seq 1 $n_exp); do
        if [ "$a" -gt "$((halfWindowSize + 1))" ] && [ "$((a))" -lt "$(($n_exp - $halfWindowSize))" ]; then
            leftLimit=$(( a - $halfWindowSize - 1))
            calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits" "${fileArray[@]:$leftLimit:$windowSize}"

            tmpIndex=$(( a - 1 ))
            propagateKeyword "${fileArray[$tmpIndex]}" $dateHeaderKey "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits"
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

    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        i=$imageDir/$base
        out=$outputDir/$base
       
        # -- NOTE ---
        # At the beginning we were using the first commented block. Each image with its flat
        # The thing is that we use the onOff strategy sometimes the matching between the frames and the 
        # flats is not that perfect... That's why I relax the checks and if something has no counterpart the right running flat is used

        if [ "$a" -le "$((halfWindowSize + 1))" ]; then
            flatToUse=$flatDir/flat-it*_"$filter"_n"$currentNight"_left_ccd"$h".fits
        elif [ "$a" -ge "$((n_exp - halfWindowSize))" ]; then
            flatToUse=$flatDir/flat-it*_"$filter"_n"$currentNight"_right_ccd"$h".fits
        else
            flatToUse=$flatDir/flat-it*_"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        fi
        astarithmetic $i -h1 $flatToUse -h1 / -o $out

        propagateKeyword $i $dateHeaderKey $out
        propagateKeyword $i $airMassKeyWord $out
        propagateKeyword $i $pointingRA $out
        propagateKeyword $i $pointingDEC $out
    done
    echo done > $flatDone
}
export -f divideImagesByRunningFlats

divideImagesByWholeNightFlat(){
    local imageDir=$1
    local outputDir=$2
    local flatToUse=$3
    local flatDone=$4

    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        i=$imageDir/$base
        out=$outputDir/$base

        astarithmetic $i -h1 $flatToUse -h1 / -o $out
        propagateKeyword $i $dateHeaderKey $out
        propagateKeyword $i $airMassKeyWord $out
        propagateKeyword $i $pointingRA $out
        propagateKeyword $i $pointingDEC $out
    done
    echo done > $flatDone
}
export -f divideImagesByWholeNightFlat

runNoiseChiselOnFrame() {
    local baseName=$1
    local inputFileDir=$2
    local outputDir=$3
    local noiseChiselParams=$4

    imageToUse=$inputFileDir/$baseName
    output=$outputDir/$baseName
    astnoisechisel $imageToUse $noiseChiselParams --numthreads=$num_cpus -o $output
}
export -f runNoiseChiselOnFrame

# Functions for Warping the frames
getCentralCoordinate(){
    local image=$1

    NAXIS1=$(gethead $image NAXIS1)
    NAXIS2=$(gethead $image NAXIS2)

    # Calculate the center pixel coordinates
    center_x=$((NAXIS1 / 2))
    center_y=$((NAXIS2 / 2))

    # Use xy2sky to get the celestial coordinates of the center pixel
    imageCentre=$( xy2sky $image $center_x $center_y )
    echo $imageCentre
}
export -f getCentralCoordinate

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
export -f detect_swarp

warpImage() {
    local imageToSwarp=$1
    local entireDir_fullGrid=$2
    local entiredir=$3
    local ra=$4
    local dec=$5
    local coaddSizePx=$6

    # ****** Decision note *******
    # We need to regrid the frames into the final coadd grid. But if we do this right now we will be processing
    # frames too big (which are mostly Nans) and the noisechisel routine takes a looot of time
    # The approach taken is to move the frame to that grid, and then crop it to the dimension of the data itself
    # We need to store both. I have tried to store the small one and then warp it again to the big grid, it's more time consuming
    # and the nan wholes grow so we end up with less light in the final coadd.

    currentIndex=$(basename $imageToSwarp .fits)
    tmpFile1=$entiredir"/$currentIndex"_temp1.fits
    frameFullGrid=$entireDir_fullGrid/entirecamera_$currentIndex.fits

    # Resample into the final grid
    SWARP_CMD=$( detect_swarp )
    $SWARP_CMD -c $swarpcfg $imageToSwarp -CENTER $ra,$dec -IMAGE_SIZE $coaddSizePx,$coaddSizePx -IMAGEOUT_NAME $entiredir/"$currentIndex"_swarp1.fits -WEIGHTOUT_NAME $entiredir/"$currentIndex"_swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $pixelScale -PIXELSCALE_TYPE MANUAL

    # Mask bad pixels
    astarithmetic $entiredir/"$currentIndex"_swarp_w1.fits -h0 set-i i i 0 lt nan where -o$tmpFile1
    astarithmetic $entiredir/"$currentIndex"_swarp1.fits -h0 $tmpFile1 -h1 0 eq nan where -o$frameFullGrid

    regionOfDataInFullGrid=$(python3 $pythonScriptsPath/getRegionToCrop.py $frameFullGrid 1)
    read row_min row_max col_min col_max <<< "$regionOfDataInFullGrid"
    astcrop $frameFullGrid --polygon=$col_min,$row_min:$col_max,$row_min:$col_max,$row_max:$col_min,$row_max --mode=img  -o $entiredir/entirecamera_"$currentIndex".fits --quiet
    echo $row_min $row_max $col_min $col_max > $entiredir/entirecamera_"$currentIndex"_cropRegion.txt

    rm $entiredir/"$currentIndex"_swarp_w1.fits $entiredir/"$currentIndex"_swarp1.fits $tmpFile1 $frameFullGrid

    # I'm manually propagating the date because is used in some versions of the pipeline (amateur data) but  swarp for some reason propagates it incorrectly
    propagateKeyword $imageToSwarp $dateHeaderKey $entiredir/entirecamera_"$currentIndex".fits 
}
export -f warpImage


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
        echo -e "\n\tFrames from $smallGridDir have been already pased into the full grid\n"
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

removeBadFramesFromReduction() {
    local sourceToRemoveFiles=$1
    local destinationDir=$2
    local badFilesWarningDir=$3
    local badFilesWarningFile=$4
    local prefixOfFilesToRemove=$5

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
    local noisechisel_param=${15}
    local manualMaskParams=${16}


    # Masking the frames if they are not already 
    i=$entiredir/$1
    out=$(echo $base | sed 's/.fits/.txt/')

    if ! [ "$inputImagesAreMasked" = true ]; then
        tmpMask=$(echo $base | sed 's/.fits/_mask.fits/')
        tmpMaskedImage=$(echo $base | sed 's/.fits/_masked.fits/')
        astnoisechisel $i $noisechisel_param --numthreads=$num_cpus -o $noiseskydir/$tmpMask
        astarithmetic $i -h1 $noiseskydir/$tmpMask -h1 1 eq nan where float32 -o $noiseskydir/$tmpMaskedImage -quiet
        imageToUse=$noiseskydir/$tmpMaskedImage
        rm -f $noiseskydir/$tmpMask

        # manual masks defined by the user
        valueToPut=nan
        read -r -a maskArray <<< "$manualMaskParams"
        for ((i=0; i<${#maskArray[@]}; i+=5)); do
            ra_mask="${maskArray[i]}"
            dec_mask="${maskArray[i+1]}"
            r="${maskArray[i+2]}"
            axisRatio="${maskArray[i+3]}"
            pa="${maskArray[i+4]}"

            python3 $pythonScriptsPath/manualMaskRegionFromWCSArea.py $imageToUse $valueToPut $ra_mask $dec_mask $r $axisRatio $pa
        done
    else    
        imageToUse=$i
    fi


    if [ "$constantSky" = true ]; then  # Case when we subtract a constant
        # ****** Decision note *******
        # Here we have implemented two possibilities. Either the background is estimated by a constant or by a polynomial.
        # If it is a constant we store the fileName, the background value and the std. This is implemented this way because we
        # need the background to subtract but also later the std for weighing the frames
        # If it is a polynomial we only use it to subtract the background (the weighing is always with a constat) so we only store
        # the coefficients of the polynomial.
        #
        # Storing this values is also relevant for checking for potential bad frames


        # The problem is that I can't use the same ring/s as in the normalisation because here we have warped and cropped the images... So I create a new normalisation ring from the centre of the images
        # I cannot even create a common ring for all, because they are cropped based on the number of non-nan (depending on the vignetting and how the NAN are distributed), so i create a ring per image
        # For that reason the subtraction of the background using the ring is always using a ring centered in the frame
        # More logic should be implemented to use the normalisation ring(s) and recover them after the warping and cropping
        if [ "$constantSkyMethod" = "ring" ]; then
            # We generate the ring (we cannot use the normalisation ring because we have warped and cropped) and compute the background value within it
            tmpRingDefinition=$(echo $base | sed 's/.fits/_ring.txt/')
            tmpRingFits=$(echo $base | sed 's/.fits/_ring.fits/')

            naxis1=$(fitsheader $imageToUse | grep "NAXIS1" | awk '{print $3}')
            naxis2=$(fitsheader $imageToUse | grep "NAXIS2" | awk '{print $3}')
            half_naxis1=$(echo "$naxis1 / 2" | bc)
            half_naxis2=$(echo "$naxis2 / 2" | bc)

            ringRadius=$( awk '{print $5}' $ringDir/ring.txt )
            echo "1 $half_naxis1 $half_naxis2 6 $ringRadius 1 1 1 1 1" > $ringDir/$tmpRingDefinition
            astmkprof --background=$imageToUse  -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringDir/$tmpRingFits $ringDir/$tmpRingDefinition
            
            me=$(getMedianValueInsideRing $imageToUse  $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $noiseskydir)
            std=$(getStdValueInsideRing $imageToUse $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $noiseskydir)
            read skew kurto < <(getSkewKurtoValueFromSkyPixels $imageToUse $constantSkyMethod $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)

            echo "$base $me $std $skew $kurto" > $noiseskydir/$out
            rm $ringDir/$tmpRingDefinition
            rm $ringDir/$tmpRingFits

        elif [ "$constantSkyMethod" = "wholeImage" ]; then
            me=$(aststatistics $imageToUse -h1 --sclipparams=3,3 --sigclip-median --quiet)
            std=$(aststatistics $imageToUse -h1 --sclipparams=3,3 --sigclip-std --quiet)
            read skew kurto < <(getSkewKurtoValueFromSkyPixels $imageToUse $constantSkyMethod $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)
            echo "$base $me $std $skew $kurto" > $noiseskydir/$out
       
        elif [ "$constantSkyMethod" = "noisechisel" ]; then
            sky=$(echo $base | sed 's/.fits/_sky.fits/')
            astnoisechisel $imageToUse --tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2 --checksky $noisechisel_param --numthreads=$num_cpus -o $noiseskydir/$base
            mean=$(aststatistics $noiseskydir/$sky -hSKY --sigclip-mean)
            std=$(aststatistics $noiseskydir/$sky -hSTD --sigclip-mean)

            read skew kurto < <(getSkewKurtoValueFromSkyPixels $imageToUse $constantSkyMethod $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)
            echo "$base $mean $std $skew $kurto" > $noiseskydir/$out
            rm -f $noiseskydir/$sky
        
        else
            errorNumber=6
            echo -e "\nAn invalid value for the sky_estimation_method was provided" >&2
            echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
            exit $errorNumber
        fi
        
        rm $imageToUse
    else
        # Case when we fit a plane
        planeOutput=$(echo $base | sed 's/.fits/_poly.fits/')
        planeCoeffFile=$(echo $base | sed 's/.fits/.txt/')

        python3 $pythonScriptsPath/surface-fit.py -i $imageToUse -o $noiseskydir/$planeOutput -d $polyDegree -f $noiseskydir/$planeCoeffFile
        rm -f $noiseskydir/$noiseOutTmp
        rm -f $noiseskydir/$maskTmp
        rm -f $noiseskydir/$planeOutput
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
    local noisechisel_param=${15}
    local maskParams=${16}

    if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
    if [ -f $noiseskydone ]; then
        echo -e "\n\tScience images are 'noisechiseled' for constant sky substraction for extension $h\n"
    else
        framesToComputeSky=()
        for a in $( ls $framesToUseDir/*.fits ); do
            base=$( basename $a )
            framesToComputeSky+=("$base")
        done

        printf "%s\n" "${framesToComputeSky[@]}" | parallel -j "$num_cpus" computeSkyForFrame {} $framesToUseDir $noiseskydir $constantSky $constantSkyMethod $polyDegree $inputImagesAreMasked $ringDir $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth "'$noisechisel_param'" "'$maskParams'"
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
        me=$(awk 'NR=='1'{print $2}' $i)
        astarithmetic $input -h1 $me - -o$output;
    else
        i=$directoryWithSkyValues/"entirecamera_"$a"_poly.fits"
        NAXIS1_image=$(gethead $input NAXIS1); NAXIS2_image=$(gethead $input NAXIS2)
        NAXIS1_plane=$(gethead $i NAXIS1); NAXIS2_plane=$(gethead $i NAXIS2)

        if [[ "$NAXIS1_image" == "$NAXIS1_plane" ]] && [[ "$NAXIS2_image" == "$NAXIS2_plane" ]]; then
            astarithmetic $input -h1 $i -h1 - -o$output
        else
            python3 $pythonScriptsPath/moveSurfaceFitToFullGrid.py $input $i 1 $NAXIS1_image $NAXIS2_image $directoryToStoreSkySubtracted/"planeToSubtract_"$a".fits"
            astarithmetic $input -h1 $directoryToStoreSkySubtracted/"planeToSubtract_"$a".fits" -h1 - -o$output
            rm $directoryToStoreSkySubtracted/"planeToSubtract_"$a".fits"
        fi
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
    astmatch $tmpFolder/decals_c.txt --hdu=1    $BDIR/catalogs/"$objectName"_"$surveyToUseInSolveField".fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=bRA,bDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o $tmpFolder/match_decals_gaia.txt 1>/dev/null

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
        echo -e "\n\tSpectra already downloaded\n"
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


downloadCatalogue() {
    local surveyToUse=$1
    local ra=$2
    local dec=$3
    local radius=$4
    local catDir=$5
    local catName=$6

    case "$surveyToUse" in
        gaia)
            query_param="gaia --dataset=dr3 --center=$ra,$dec --radius=$radius --column=ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error"
            downloadGaiaCatalogue "$query_param" "$catDir" "$catName"
            ;;
        panstarrs)
            query_param="vizier --dataset=panstarrs1 --center=$ra,$dec --radius=$radius --column=RAJ2000,DEJ2000,gmag"
            downloadPanstarrsCatalogue "$query_param" "$catDir" "$catName"
            ;;
        *)
            echo "Unknown catalog source: $surveyToUse"
            exit 222
            ;;
    esac

}
export -f downloadCatalogue


downloadGaiaCatalogue() {
    local query=$1
    local catdir=$2
    local catName=$3

    astquery $query -o $catdir/"$objectName"_Gaia_DR3_tmp.fits
    asttable $catdir/"$objectName"_Gaia_DR3_tmp.fits -c1,2,3 -c'arith $4 abs' -c'arith $5 3 x' -c'arith $6 abs' -c'arith $7 3 x' -c'arith $8 abs' -c'arith $9 3 x' --noblank=4 -o$catdir/tmp.txt

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

    # # Here we don't demand any condition
    # asttable $catdir/tmp.txt -o $catName

    rm $catdir/test1.txt $catdir/tmp.txt $catdir/"$objectName"_Gaia_DR3_tmp.fits $catdir/test_.txt
}
export -f downloadGaiaCatalogue

downloadPanstarrsCatalogue() {
    local query=$1
    local catdir=$2
    local catName=$3

    astquery $query -o $catdir/"$objectName"_Panstarrs_S1_tmp.fits
    asttable $catdir/"$objectName"_Panstarrs_S1_tmp.fits -c1,2,3  --colmetadata=1,RA,deg --colmetadata=2,DEC,deg --colmetadata=3,phot_g_mean_mag,mag -o$catName
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
                            -E -A RA -D  DEC\
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
    local astroimadir=$8
    local sexcfg_sf=$9
    local sizeOfOurFieldDegrees=${10}
    base=$( basename $i)


    # Get the RA and Dec of the pointing. It has to be converted to deg
    LC_NUMERIC=C  # Format to get rid of scientific notation if needed

    pointingRAValue=$( astfits $i --keyvalue=$pointingRA --quiet)
    pointingRAValue=$( printf "%.8f\n" " $pointingRAValue")
    if [[ "$pointingRAUnits" == "hours" ]]; then
        pointRA=$(echo "$pointingRAValue * 15" | bc -l)
    elif [[ "$pointingRAUnits" == "deg" || "$pointingRAUnits" == "degrees" ]]; then
        pointRA="$pointingRAValue"
    else
        echo "Error: Unsupported RA units: $pointingRAUnits"
        exit 888
    fi

    pointingDecValue=$( astfits $i --keyvalue=$pointingDEC --quiet)
    pointingDecValue=$( printf "%.8f\n" " $pointingDecValue")
    if [[ "$pointingDECUnits" == "hours" ]]; then
        pointDec=$(echo "$pointingDecValue * 15" | bc -l)
    elif [[ "$pointingDECUnits" == "deg" || "$pointingDECUnits" == "degrees" ]]; then
        pointDec="$pointingDecValue"
    else
        echo "Error: Unsupported RA units: $pointingDECUnits"
        exit 888
    fi

    # The default sextractor parameter file is used.
    # I tried to use the one of the config directory (which is used in other steps), but even using the default one, it fails
    # Maybe a bug? I have not managed to make it work
    max_attempts=4
    attempt=1
    sex_path=$( which source-extractor )
    while [ $attempt -le $max_attempts ]; do
        #Sometimes the output of solve-field is not properly writen in the computer (.i.e, size of file=0). 
        #Because of that, we iterate solve-field in a maximum of 4 times until file is properly saved
        echo solve-field $i --no-plots --ra $pointRA --dec $pointDec --radius $sizeOfOurFieldDegrees
        solve-field $i --no-plots --ra $pointRA --dec $pointDec --radius $sizeOfOurFieldDegrees\
        -L $solve_field_L_Param -H $solve_field_H_Param -u $solve_field_u_Param \
        --overwrite --extension 1 --config $confFile/astrometry_$objectName.cfg --no-verify \
        --use-source-extractor --source-extractor-path=$sex_path \
        --source-extractor-config=$sexcfg_sf --x-column X_IMAGE --y-column Y_IMAGE \
        --sort-column MAG_AUTO --sort-ascending  \
        -Unone --temp-axy  -Snone -Mnone -Rnone -Bnone -N$astroimadir/$base ;
        if [ -s "$layer_temp" ]; then
            attempt=$max_attempts
        fi
            

        ((attempt++))
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
    local gain=$8 


    # Here I put the saturation threshold and the gain directly.
    # This is because it's likely that we end up forgetting about tuning the sextractor configuration file but we will be more careful with the configuration file of the reductions
    # These two values (saturation level and gain) are key for astrometrising correctly, they are used by scamp for identifying saturated sources and weighting the sources
    # I was, in fact, having frames bad astrometrised due to this parameters.
    i=$astroimadir/"$a".fits
    source-extractor $i -c $sexcfg -PARAMETERS_NAME $sexparam -FILTER_NAME $sexconv -CATALOG_NAME $sexdir/$a.cat -SATUR_LEVEL=$saturationThreshold -GAIN=$gain
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
#     SWarp -c $swarpcfg $downSampledImages -CENTER $ra,$dec -IMAGE_SIZE $mosaicSize,$mosaicSize -IMAGEOUT_NAME $swarpedImagesDir/"$a"_swarp1.fits \
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
        elif [ "$survey" = "SDSS" ]; then
            for i in $( ls $dirWithBricks/sdss_*.fits); do
                currentName=$( basename $i )
                brickList+=("$currentName")
            done
        fi

        headerWithData=0 # After decompressing the data ends up in the hdu 0
        noisechiselTileSize=50
        printf "%s\n" "${brickList[@]}" | parallel -j "$num_cpus" selectStarsAndRangeForCalibrateSingleFrame {} $dirWithBricks $cataloguedir $headerWithData $methodToUse $noisechiselTileSize $apertureUnits
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
        brickName=decompressed_decal_image_"$brick".fits
    elif [[ "$survey" = "PANSTARRS" ]]; then
        brickName=cal_"$brick".fits
    elif [[ "$survey" = "SDSS" ]]; then
        brickName=$brick.fits
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
    local numberOfApertureForRecuperateGAIA=$6

    echo -e "\n·Performing aperture photometry to the bricks of the folder: $brickDir"


    outputDone=$outputCat/done.txt
    if ! [ -d $outputCat ]; then mkdir $outputCat; fi
    if [ -f $outputDone ]; then
        echo -e "\n\tCatalogues done with aperture photometry already done\n"
    else
        brickList=()

        if [[ "$survey" = "DECaLS" ]]; then
            for a in $( ls $brickDir/decompressed*.fits); do
                brickName=$(echo "$a" | awk -F'_image_' '{print $2}' | awk -F'.fits' '{print $1}')
                brickList+=("$brickName")
            done
        elif [[ "$survey" = "PANSTARRS" ]]; then
            for a in $( ls $brickDir/cal_*.fits); do
                brickName=$(echo "$a" | awk -F'cal_' '{print $2}' | awk -F'.fits' '{print $1}')
                brickList+=("$brickName")
            done
        elif [[ "$survey" = "SDSS" ]]; then
            for a in $( ls $brickDir/sdss_*.fits); do
                brickName="sdss"$(echo "$a" | awk -F'/sdss' '{print $2}' | awk -F'.fits' '{print $1}')
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
            prepareSpectraDataForPhotometricCalibration $spectraDir $filter $ra $dec $mosaicDir $sizeOfOurFieldDegrees $surveyForSpectra $transmittanceCurveFile
            mv $mosaicDir/wholeFieldPhotometricCatalogue_$filter.cat $mosaicDir/wholeFieldPhotometricCatalogue.cat # This is not the best solution. We need here the catalogue without the filter
                                                                                                                    # but if we calibrate with pansatarrs we need the function "prepareSpectraDataForPhotometricCalibration"
                                                                                                                    # to create the wholefieldphotometricCatlaogues with the filter in the name, so I changed it here when needed

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
    local sizeOfOurFieldDegrees=$6
    local surveyForSpectra=$7
    local transmittanceCurveFile=$8

    if ! [ -d $spectraDir ]; then mkdir $spectraDir; fi
    downloadSpectra $mosaicDir $spectraDir $ra $dec $sizeOfOurFieldDegrees $surveyForSpectra

    aperturePhotDone=$mosaicDir/done_"$filter".txt
    if [ -f $aperturePhotDone ]; then
        echo -e "\n\tThe catalogue with the magnitudes of the spectra already built\n"
    else
        output_tmpCat=$mosaicDir/wholeFieldPhotometricCatalogue_tmp.cat
        outputCat=$mosaicDir/wholeFieldPhotometricCatalogue_$filter.cat
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



    sizeOfBrick_gaia=3600 #3600 is the default value because the pipeline originally worked with DECaLS and it used this brick-size. We try to keep it to 3600 when possible
    downloadSurveyData $mosaicDir $surveyImagesDir_g $bricksIdentificationFile_g "g" $ra $dec $sizeOfOurFieldDegrees $gaiaCatalogue $survey $sizeOfBrick       
    downloadSurveyData $mosaicDir $surveyImagesDirForGaiaCalibration_g $bricksIdentificationFileForGaiaCalibration_g "g" $ra $dec $sizeOfFieldForCalibratingPANSTARRStoGAIA $gaiaCatalogue $survey $sizeOfBrick_gaia
    downloadSurveyData $mosaicDir $surveyImagesDir_r $bricksIdentificationFile_r "r" $ra $dec $sizeOfOurFieldDegrees $gaiaCatalogue $survey $sizeOfBrick
    downloadSurveyData $mosaicDir $surveyImagesDirForGaiaCalibration_r $bricksIdentificationFileForGaiaCalibration_r "r" $ra $dec $sizeOfFieldForCalibratingPANSTARRStoGAIA $gaiaCatalogue $survey $sizeOfBrick_gaia

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

    # Preprocessing needed depending on what survey we are usiing
    if [ "$survey" = "DECaLS" ]; then
        # I have to decompressed them, I don't know what DECaLS does exactly but otherwise I can't run sextractor directly on the downloaded bricks
        # I can run noisechisel, but since it is quickly to decompress I prefer to simplify the logic of the code and decompress always
        decompressBricks $surveyImagesDir_g
        decompressBricks $surveyImagesDirForGaiaCalibration_g
        decompressBricks $surveyImagesDir_r
        decompressBricks $surveyImagesDirForGaiaCalibration_r
        decompressBricks $surveyImagesDir
        decompressBricks $surveyImagesDirForGaiaCalibration
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

   
    # halfMaxRad_Mag_plots=$mosaicDir/halfMaxradius_Magnitude_plots
    # if [ $survey == "DECaLS" ]; then
    #     produceHalfMaxRadiusPlotsForDecals $selectedSurveyStarsDir $halfMaxRad_Mag_plots $filter
    # elif [ "$survey" == "PANSTARRS" ]; then
    #     produceHalfMaxRadiusPlotsForPanstarrs $selectedSurveyStarsDir $halfMaxRad_Mag_plots $filter
    # fi


    # These values have been calculated for having a stable photometry (increasing aperture does not change the magnitudes anymore)
    # This should recover GAIA magnitudes (and in fact it generaly does) but sometimes an offset is observed. That's why we
    # compute and correct that offset in "calibrationToGAIA" function
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
    
    # # --- Decision note ---
    # # Now two corrections are applied. First a colour correction to match the survey filter to the filter of our telescope, and
    # # Secondly an offset to match the survey photometry to GAIA. Thus, our photometry is referenced to GAIA spectra and to our filter
    # # The order in which the corrections are applied is slightly relevant, because the offset to match the survey photometry to GAIA
    # # can vary across different filters, then you get slightly differect colours and slightly different colour corrections
    # # Since the colour corrections are obtained from g_gaia - r_gaia I do the same here for consistency

    spectraDir=$mosaicDir/gaiaSpectra
    magFromSpectraDir_g=$mosaicDir/magnitudesFromGaiaSpectra_g
    magFromSpectraDir_r=$mosaicDir/magnitudesFromGaiaSpectra_r

    
    # These two ranges (14.5-15.5 for g and 13.65-15 for r) are tested that work for calibrating panstarrs to gaia in these bands. 
    calibrationToGAIA $spectraDir $folderWithTransmittances "g" $ra $dec $mosaicDir $sizeOfFieldForCalibratingPANSTARRStoGAIA $magFromSpectraDir_g $aperturePhotDir"ForGAIACalibration_g" 14.5 15.5 $survey
    calibrationToGAIA $spectraDir $folderWithTransmittances "r" $ra $dec $mosaicDir $sizeOfFieldForCalibratingPANSTARRStoGAIA $magFromSpectraDir_r $aperturePhotDir"ForGAIACalibration_r" 14 15 $survey

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

        calibrationToGAIA $spectraDir $folderWithTransmittances $filter $ra $dec $mosaicDir $sizeOfFieldForCalibratingPANSTARRStoGAIA $magFromSpectraDir $aperturePhotDir"ForGAIACalibration" $calibrationBrightLimit $calibrationFaintLimit $survey
        read offset factorToApplyToCounts < "$mosaicDir/offsetToCorrectSurveyToGaia_"$filter".txt"
        correctOffsetFromCatalogues $aperturePhotDir $offset $factorToApplyToCounts "beforeCorrectingPanstarrsGAIAOffset"
        correctOffsetFromCatalogues $aperturePhotDir"ForGAIACalibration" $offset $factorToApplyToCounts "beforeCorrectingPanstarrsGAIAOffset"
    fi

    
    computeColoursAndAddThemToCatalogues $aperturePhotDir $aperturePhotDir"_g" $aperturePhotDir"_r" $filter
    applyColourcorrectionToAllCatalogues $aperturePhotDir "$filterCorrectionCoeff"
    
    imagesHdu=1
    brickDecalsAssociationFile=$mosaicDir/frames_bricks_association.txt
    if [[ -f $brickDecalsAssociationFile ]]; then
        rm $brickDecalsAssociationFile
    fi
    python3 $pythonScriptsPath/associateDecalsBricksToFrames.py $referenceImagesForMosaic $imagesHdu $bricksIdentificationFile $brickDecalsAssociationFile $survey
   
    # Create a catalogue of the whole field
    listOfCatalogues=()
    for i in "$aperturePhotDir"/*.cat; do
        tmpName=$( basename $i )
        listOfCatalogues+=("${tmpName%.cat}")
    done
    combineCatalogues $mosaicDir $aperturePhotDir "wholeFieldPhotometricCatalogue.cat" "${listOfCatalogues[@]}"
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

    LC_NUMERIC=C
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

            # echo $fileName
            # echo "$cataloguesToUseDir/$fileName.$filt.cat"
            # echo "$cataloguesDir_r/$fileName.r.cat"

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
    local suffix=$4

    echo -e "\n·Correcting offset between survey and GAIA for catalogues $dirWithCatalogues\n"

    if [ -f $dirWithCatalogues/"surveyAndGAIAOffsetCorrection_done.txt" ]; then
        echo -e "\n\tOffset between survey and GAIA already corrected\n"
    else
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

            mv $i "${i%.txt}"_$suffix
            mv "$i"_tmp $i
        done
    fi
    echo "done" > $dirWithCatalogues/"surveyAndGAIAOffsetCorrection_done.txt"
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
    local survey=${12}

    echo -e "\n·Computing offset between PANSTARRS data and GAIA (panstarrs catalogues used for the calculation $panstarrsCatalogueDir)"

    
    if [ -f $mosaicDir"/offsetToCorrectSurveyToGaia_"$filter".txt" ]; then
            echo -e "\n\tOffset for calibrating survey to GAIA already computed\n"
    else
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

        
        transmittanceCurveFile="$folderWithTransmittances"/"$survey"_"$filter".dat
        prepareSpectraDataForPhotometricCalibration $spectraDir $filter $ra $dec $mosaicDir $sizeOfFieldForCalibratingPANSTARRStoGAIA "GAIA" $transmittanceCurveFile
        
        offsetValues=$( python3 $pythonScriptsPath/getOffsetBetweenPANSTARRSandGAIA.py $panstarrsCatalogueDir/mergedCatalogue.cat $mosaicDir/wholeFieldPhotometricCatalogue_"$filter".cat $brightLimitToCompareGAIAandPANSTARRS $faintLimitToCompareGAIAandPANSTARRS $mosaicDir $filter )

        offset=$(echo $offsetValues | awk '{print $1}')
        factorToApplyToCounts=$(echo $offsetValues | awk '{print $2}')
        echo $offset $factorToApplyToCounts > $mosaicDir/offsetToCorrectSurveyToGaia_"$filter".txt
        rm $mosaicDir/wholeFieldPhotometricCatalogue_"$filter".cat # I delete it because later this file has to contain the catalogue of the whole field from panstarrs
    fi
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
    local methodToUse=$5
    local tileSize=$6           # This parameter will only be used if the catalogue is being generated with noisechisel
    local apertureUnits=$7

    i=$framesForCalibrationDir/$a

    if [[ "$methodToUse" == "sextractor" ]]; then
        outputCatalogue=$( generateCatalogueFromImage_sextractor $i $mycatdir $a $apertureUnits )
    elif [[ "$methodToUse" == "noisechisel" ]]; then
        outputCatalogue=$( generateCatalogueFromImage_noisechisel $i $mycatdir $a $headerToUse $tileSize $apertureUnits )
    else
        errorNumber=9
        echo "Error, method for selecting stars and the range in the calibration not recognised"
        echo "Exiting with error number: $erroNumber"
        exit $erroNumber
    fi
    
    astmatch $outputCatalogue --hdu=1 $BDIR/catalogs/"$objectName"_gaia.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o$mycatdir/match_"$a"_my_gaia.txt
    
    # The intermediate step with awk is because I have come across an Inf value which make the std calculus fail
    # Maybe there is some beautiful way of ignoring it in gnuastro. I didn't find int, I just clean de inf fields.
    s=$(asttable $mycatdir/match_"$a"_my_gaia.txt -h1 -c5 --noblank=MAGNITUDE   | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
    std=$(asttable $mycatdir/match_"$a"_my_gaia.txt -h1 -c5 --noblank=MAGNITUDE | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
    minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
    maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)


    echo $s $std $minr $maxr > $mycatdir/range_"$a".txt
    asttable $outputCatalogue --range=HALF_MAX_RADIUS,$minr,$maxr -o $mycatdir/selected_"$a"_automatic.txt
}
export -f selectStarsAndRangeForCalibrateSingleFrame

selectStarsAndSelectionRangeOurData() {
    local iteration=$1
    local framesForCalibrationDir=$2
    local mycatdir=$3
    local methodToUse=$4
    local tileSize=$5
    local apertureUnits=$6

    mycatdone=$mycatdir/done.txt
    if ! [ -d $mycatdir ]; then mkdir $mycatdir; fi
    if [ -f $mycatdone ]; then
            echo -e "\n\tSources for photometric calibration are already extracted for my image\n"
    else
        framesToUse=()
        for a in $( ls $framesForCalibrationDir/*.fits ); do
            framename=$( basename $a )
            framesToUse+=("$framename")
        done

        headerWithData=1
        printf "%s\n" "${framesToUse[@]}" | parallel -j "$num_cpus" selectStarsAndRangeForCalibrateSingleFrame {} $framesForCalibrationDir $mycatdir $headerWithData $methodToUse $tileSize $apertureUnits
        echo done > $mycatdone
    fi
}
export -f selectStarsAndSelectionRangeOurData


matchDecalsAndSingleFrame() {
    local base=$1
    local myCatalogues=$2
    local calibrationCatalogues=$3
    local matchdir=$4
    local surveyForCalibration=$5
    local calibratingMosaic=$6

    out=$matchdir/matched_"$base".cat
    ourDataCatalogue=$myCatalogues/$base

    tmpCatalogue=$matchdir/match-$base-tmp.cat
    out=$matchdir/match-"$base".cat

    # Choose between the whole field catalogue (spectra calibration or calibrating coadd prephot) or the specific catalogue (imaging calibration)
    # I have to do this independently because if it is spectra data I want to preserve the spectrograph and the objectype
    if [[ ($surveyForCalibration == "SPECTRA") || ("$calibratingMosaic" == true) ]]; then
        calibrationCat=$calibrationCatalogues/wholeFieldPhotometricCatalogue.cat
    else
        calibrationCat=$calibrationCatalogues/"${base%%.*}".cat
    fi

    if [ $surveyForCalibration == "SPECTRA" ]; then
        astmatch $ourDataCatalogue --hdu=1 $calibrationCat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 \
                --outcols=bRA,bDEC,aRA,aDEC,bMAGNITUDE,bSUM,aMAGNITUDE,aSUM -o$tmpCatalogue

        asttable $tmpCatalogue --output=$out --colmetadata=1,RA,deg,"Right ascension survey" \
            --colmetadata=2,DEC,none,"Declination survey" \
            --colmetadata=3,RA,deg,"Right ascension data being reduced" \
            --colmetadata=4,DEC,none,"Declination data being reduced" \
            --colmetadata=5,MAGNITUDE_CALIBRATED,none,"Magnitude in survey data" \
            --colmetadata=6,SUM,none,"Sum in survey" \
            --colmetadata=7,MAGNITUDE_NONCALIBRATED,none,"Magnitude in data being reduced" \
            --colmetadata=8,SUM,none,"Sum in in data being reduced"
    else
        astmatch $ourDataCatalogue --hdu=1 $calibrationCat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 \
                --outcols=bRA,bDEC,aRA,aDEC,bMAGNITUDE,bSUM,aMAGNITUDE,aSUM -o$tmpCatalogue

        asttable $tmpCatalogue --output=$out --colmetadata=1,RA,deg,"Right ascension survey" \
                    --colmetadata=2,DEC,none,"Declination survey" \
                    --colmetadata=3,RA,deg,"Right ascension data being reduced" \
                    --colmetadata=4,DEC,none,"Declination data being reduced" \
                    --colmetadata=5,MAGNITUDE_CALIBRATED,none,"Magnitude in survey data" \
                    --colmetadata=6,SUM,none,"Sum in survey" \
                    --colmetadata=7,MAGNITUDE_NONCALIBRATED,none,"Magnitude in data being reduced" \
                    --colmetadata=8,SUM,none,"Sum in in data being reduced"
    fi

    rm $tmpCatalogue
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
    local base=$1
    local ourDatadir=$2
    local framesForCalibrationDir=$3
    local mycatdir=$4
    local numberOfApertureUnitsForCalibration=$5

    i=$framesForCalibrationDir/$base
    automaticCatalogue=$mycatdir/selected_"$base"_automatic.txt

    r_myData_pix_=$(awk 'NR==1 {printf $1}' $mycatdir/range_"$base".txt)
    r_myData_pix=$(astarithmetic $r_myData_pix_ $numberOfApertureUnitsForCalibration. x -q )

    dataHdu=1

    # raColumnName=RA
    # decColumnName=DEC
    # photometryOnImage_noisechisel $a $ourDatadir $automaticCatalogue $i $r_myData_pix $ourDatadir/$base.cat 22.5 $dataHdu \
    #                                 $raColumnName $decColumnName

    columnWithXCoordForOutDataPx=0 # These numbers come from how the catalogue of the matches stars is built. This is not very clear right now, should be improved
    columnWithYCoordForOutDataPx=1
    columnWithXCoordForOutDataWCS=2
    columnWithYCoordForOutDataWCS=3
    photometryOnImage_photutils $base $ourDatadir $automaticCatalogue $i $r_myData_pix $ourDatadir/$base.cat 22.5 $dataHdu \
                                $columnWithXCoordForOutDataPx $columnWithYCoordForOutDataPx $columnWithXCoordForOutDataWCS $columnWithYCoordForOutDataWCS
}
export -f buildOurCatalogueOfMatchedSourcesForFrame

buildOurCatalogueOfMatchedSources() {
    local ourDatadir=$1
    local framesForCalibrationDir=$2
    local mycatdir=$3
    local numberOfApertureUnitsForCalibration=$4

    ourDatadone=$ourDatadir/done.txt
    if ! [ -d $ourDatadir ]; then mkdir $ourDatadir; fi
    if [ -f $ourDatadone ]; then
        echo -e "\n\tAperture catalogs in our data done\n"
    else
        framesToUse=()
        for a in $( ls $framesForCalibrationDir/*.fits); do
            frameName=$( basename $a )
            framesToUse+=("$frameName")
        done

        printf "%s\n" "${framesToUse[@]}" | parallel -j "$num_cpus" buildOurCatalogueOfMatchedSourcesForFrame {} $ourDatadir $framesForCalibrationDir $mycatdir $numberOfApertureUnitsForCalibration
        echo done > $ourDatadone
    fi
}
export -f buildOurCatalogueOfMatchedSources

matchCalibrationStarsCatalogues() {
    local matchdir2=$1
    local ourDatadir=$2
    local decalsdir=$3

    local matchdir2done=$matchdir2/done_aperture.txt


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
        for a in $( ls $matchdir/*.cat ); do                       
            baseName=$( basename $a )
            a=$( echo $baseName | awk -F'[_.]' '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/ && $(i+1)=="fits") print $i}')

            alphaFile=alpha_$a.txt
            f=$matchdir/$baseName

            alphatruet=$alphatruedir/"$objectName"_"$filter"_"$a".txt
            asttable $f -h1 --range=MAGNITUDE_CALIBRATED,$brightLimit,$faintLimit -o$alphatruet
            raCol=0
            decCol=1
            python3 $pythonScriptsPath/createDS9RegionsFromCatalogue.py "$alphatruet" "$alphatruedir/"$objectName"_"$filter"_"$a"_starsUsedForCalibrate.reg" "plain" $raCol $decCol
            asttable $alphatruet -h1 -c1,2,'arith $6 $8 /' -o$alphatruedir/$alphaFile

            # python3 /home/sguerra/pipeline/pipelineScripts/tmp_diagnosis_distributionOfCalibrationFactorsInFrame.py $alphatruedir/$alphaFile $a

            # This was done because if we use sigclipmean or median with only 3 stars it gives nan. But if we use mean always it gives wrong values when we have more stars
            if [ "$(asttable "$alphatruedir/$alphaFile" | wc -l)" -eq 3 ]; then
                mean=$(asttable "$alphatruedir/$alphaFile" -c3 | aststatistics --mean)
            else
                mean=$(asttable "$alphatruedir/$alphaFile" -c3 | aststatistics --sigclip-median)
            fi
            
            std=$(asttable $alphatruedir/$alphaFile -c3 | aststatistics --std)
            echo "$mean $std" > $alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a".txt
            count=$(asttable $alphatruedir/$alphaFile -c3 | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --number)
            echo "Frame number $a: $count" >> $numberOfStarsUsedToCalibrateFile
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
    local bricks=("$@")

    catalogueName=$( echo "$frame" | awk -F'.' '{print $1}')

    firstBrick=${bricks[0]}
    asttablePrompt="asttable $cataloguesDir/$firstBrick.cat -o$outputDir/$catalogueName.cat"

    remainingBricks=$(echo "$bricks" | cut -d' ' -f2-)  # Get the rest of the bricks
    for brick in "${bricks[@]}"; do  # Iterate over remaining bricks
        asttablePrompt+=" --catrowfile=$cataloguesDir/$brick.cat"
    done

    $asttablePrompt
    raCol=3
    decCol=4
    python3 $pythonScriptsPath/createDS9RegionsFromCatalogue.py "$outputDir/$catalogueName.cat" "$outputDir/$catalogueName.reg" "plain" $raCol $decCol
}
export -f combineCatalogues


combineDecalsBricksCataloguesForEachFrame() {
    local outputDir=$1
    local frameBrickAssociationFile=$2
    local decalsCataloguesDir=$3

    combinationDone=$outputDir/done.txt
    if ! [ -d $outputDir ]; then mkdir $outputDir; fi
    if [ -f $combinationDone ]; then
        echo -e "\nCombination of the bricks catalogues for each frame already done\n"
    else
        count=0
        while IFS= read -r line; do
            currentLine=$line
            frame=$(echo "$line" | awk '{print $1}')
            frame=$( basename $frame )
            bricks=$(echo "$line" | cut -d' ' -f2-)
            brickList=()
            for i in ${bricks[@]}; do
                brickList+=("$i")
            done
            combineCatalogues $outputDir $decalsCataloguesDir $frame "${brickList[@]}"
        done < "$frameBrickAssociationFile"
        echo "done" > $combinationDone
    fi

}
export -f combineDecalsBricksCataloguesForEachFrame

computeCommonCalibrationFactor() {
  local calibrationFactorsDir=$1
  local iteration=$2
  local objectName=$3
  local BDIR=$4 
  
  calibrationFactors=()
  for i in $( ls $calibrationFactorsDir/alpha_"$objectName"*.txt); do
    read currentCalibrationFactor currentStd < $i
    calibrationFactors+=("$currentCalibrationFactor")
  done

  tmpTableFits=$BDIR/tableTest.fits
  printf "%s\n" "${calibrationFactors[@]}" | asttable -o "$tmpTableFits"
  commonCalibrationFactor=$( aststatistics $tmpTableFits --sigclip-median)
  calibrationFactorsStd=$( asttable $tmpTableFits | aststatistics --sclipparams=3,3 --sigclip-std)

  echo $commonCalibrationFactor $calibrationFactorsStd > $BDIR/commonCalibrationFactor_it$iteration.txt
  rm $tmpTableFits
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
    local tileSize=${14}
    local apertureUnits=${15}
    local numberOfApertureUnitsForCalibration=${16}
    local calibratingMosaic=${17}
  
    methodToUse="sextractor"
    echo -e "\n ${GREEN} ---Selecting stars and range for our data--- ${NOCOLOUR}"
    selectStarsAndSelectionRangeOurData $iteration $imagesForCalibration $mycatdir $methodToUse $tileSize $apertureUnits

    
    echo -e "\n ${GREEN} ---Building catalogues for our data with aperture photometry --- ${NOCOLOUR}"
    buildOurCatalogueOfMatchedSources $ourDataCatalogueDir $imagesForCalibration $mycatdir $numberOfApertureUnitsForCalibration

    # If we are calibrating with spectra we just have the whole catalogue of the field
    # If we are calibrating with a survey then we have a catalogue por survey's brick and we need to combine the needed bricks for build a catalogue per frame
    if ! [ -d $prepareCalibrationCataloguePerFrame ]; then mkdir $prepareCalibrationCataloguePerFrame; fi
    if [[ ("$surveyForCalibration" == "SPECTRA") || ( "$calibratingMosaic" == true) ]]; then
        cp $mosaicDir/wholeFieldPhotometricCatalogue.cat $prepareCalibrationCataloguePerFrame
    else
        echo -e "\n ${GREEN} ---Combining decals catalogues for matching each brick --- ${NOCOLOUR}"
        combineDecalsBricksCataloguesForEachFrame $prepareCalibrationCataloguePerFrame $mosaicDir/frames_bricks_association.txt $mosaicDir/aperturePhotometryCatalogues
    fi
        
    echo -e "\n ${GREEN} ---Matching our aperture catalogues and Decals aperture catalogues--- ${NOCOLOUR}"
    matchDecalsAndOurData $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $matchdir $surveyForCalibration $calibratingMosaic
    
    echo -e "\n ${GREEN} ---Computing calibration factors (alpha)--- ${NOCOLOUR}"
    computeAndStoreFactors $alphatruedir $matchdir $brightLimit $faintLimit
}
export -f computeCalibrationFactors

getCalibrationFactorForIndividualFrame() {
    local a=$1
    local alphatruedir=$2

    alpha_cat=$alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a".txt
    alpha=$(awk 'NR=='1'{print $1}' $alpha_cat)
    echo $alpha
}
export -f getCalibrationFactorForIndividualFrame

getCommonCalibrationFactor() {
    local iteration=$1

    commonFactorFile=$BDIR/commonCalibrationFactor_it$iteration.txt
    alpha=$(awk 'NR=='1'{print $1}' $commonFactorFile)
    echo $alpha
}   
export -f getCommonCalibrationFactor

applyCalibrationFactorsToFrame() {
    local a=$1
    local imagesForCalibration=$2
    local alphatruedir=$3
    local photCorrDir=$4
    local iteration=$5
    local applyCommonCalibrationFactor=$6

    f=$imagesForCalibration/entirecamera_"$a".fits

    if [[ "$applyCommonCalibrationFactor" == "true" || "$applyCommonCalibrationFactor" == "True" ]]; then
        alpha=$( getCommonCalibrationFactor $iteration )
    elif  [[ "$applyCommonCalibrationFactor" == "false" || "$applyCommonCalibrationFactor" == "False" ]]; then
        alpha=$( getCalibrationFactorForIndividualFrame $a $alphatruedir )
    else
        echo "Value of variable applyCommonCalibrationFactor ($applyCommonCalibrationFactor) not recognised"
        exit 55
    fi
    astarithmetic $f -h1 $alpha x float32 -o $photCorrDir/entirecamera_"$a".fits
}
export -f applyCalibrationFactorsToFrame

applyCalibrationFactors() {
    local imagesForCalibration=$1
    local alphatruedir=$2
    local photCorrDir=$3
    local iteration=$4
    local applyCommonCalibrationFactor=$5

    muldone=$photCorrDir/done_calibration.txt
    if ! [ -d $photCorrDir ]; then mkdir $photCorrDir; fi
    if [ -f $muldone ]; then
            echo -e "\n\tMultiplication for alpha in the pointings (huge grid) is done for extension $h\n"
    else
        framesToApplyFactor=()
        for a in $(ls $imagesForCalibration/*.fits); do
            frameName=$( basename $a )
            frameNumber=$( echo $frameName | grep -oP '(?<=_)\d+(?=\.fits)' )
            framesToApplyFactor+=("$frameNumber")
        done

        printf "%s\n" "${framesToApplyFactor[@]}" | parallel -j "$num_cpus" applyCalibrationFactorsToFrame {} $imagesForCalibration $alphatruedir $photCorrDir $iteration $applyCommonCalibrationFactor
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

    h=0
    base=entirecamera_"$a".fits
    basetmp=entirecamera_"$a"_tmp.fits

    f=$photCorrDir/$base
    rms_min=$(awk 'NR=='1'{print $1}' $BDIR/$minRmsFileName)
    rms_f=$(awk 'NR=='1'{print $3}' $noiseskydir/entirecamera_$a.txt)

    # ****** Decision note *******
    # The weights are obtained as the quadratic ratio between the best sigma and the current sigma
    # This weights produce the optimal combinantion for a gaussian distribution
    # Ref: https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
    weight=$(astarithmetic $rms_min 2 pow $rms_f 2 pow / --quiet)
    echo "$weight" > $wdir/"$objectName"_Decals-"$filter"_"$a"_ccd"$h".txt

    # multiply each image for its weight
    wixi_im_tmp=$wdir/$basetmp              # frame x weight
    w_im_tmp=$wonlydir/$basetmp             # only weight
    wixi_im=$wdir/$base                     # frame x weight
    w_im=$wonlydir/$base                    # only weight

    astarithmetic $f -h1 $weight x --type=float32 -o$wixi_im_tmp
    astarithmetic $wixi_im_tmp -h1 $f -h1 / --type=float32 -o$w_im_tmp
    astarithmetic $wixi_im_tmp float32 -g1 -o$wixi_im
    astarithmetic $w_im_tmp float32 -g1 -o$w_im
    rm -f $wixi_im_tmp
    rm -f $w_im_tmp
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
            base=entirecamera_"$a".fits
            if [ -f $photCorrDir/$base ]; then
                framesToComputeWeight+=("$a")
            fi
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
    local photCorrImagesDir=$3
    local sigmaForStdSigclip=$4

    if ! [ -d $clippingdir ]; then mkdir $clippingdir; fi
    if [ -f $clippingdone ]; then
            echo -e "\n\tUpper and lower limits for building the masked of the weighted images already computed\n"
    else
            # Compute clipped median and std
            med_im=$clippingdir/median_image.fits
            std_im=$clippingdir/std_image.fits
            gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
            if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
                astarithmetic $(ls -v $photCorrImagesDir/*.fits) $(ls $photCorrImagesDir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-median -g1 --writeall -o$med_im
                astarithmetic $(ls -v $photCorrImagesDir/*.fits) $(ls $photCorrImagesDir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-std -g1 --writeall -o$std_im
            else
                astarithmetic $(ls -v $photCorrImagesDir/*.fits) $(ls $photCorrImagesDir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-median -g1  -o$med_im
                astarithmetic $(ls -v $photCorrImagesDir/*.fits) $(ls $photCorrImagesDir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-std -g1  -o$std_im
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
    local base=$1
    local imagesToMaskDir=$2
    local clippingDir=$3
    local imagesWithMaskedOutliersDir=$4

    tmp_ab=$imagesWithMaskedOutliersDir/"$objectName"-"$filter"_"$base"_ccd"$h"_maskabove.fits
    maskedImage=$imagesWithMaskedOutliersDir/$base

    astarithmetic $imagesToMaskDir/$base -h1 set-i i i $clippingDir/upperlim.fits -h1 gt nan where float32 -q -o $tmp_ab
    astarithmetic $tmp_ab -h1 set-i i i $clippingDir/lowerlim.fits -h1 lt nan where float32 -q -o$maskedImage
    rm -f $tmp_ab
}
export -f removeOutliersFromFrame

removeOutliersFromWeightedFrames () {
  local imagesToMaskOutliersDir=$1
  local clippingdir=$2
  local imagesWithMaskedOutliersDir=$3
  local imagesWithMaskedOutliersDone=$4

  if [ -f $imagesWithMaskedOutliersDone ]; then
      echo -e "\n\tOutliers of the phot corrected images already masked\n"
  else
      framesToRemoveOutliers=()
      for a in $( ls $imagesToMaskOutliersDir/*.fits ); do
        base=$( basename $a )
        framesToRemoveOutliers+=("$base")
      done
      printf "%s\n" "${framesToRemoveOutliers[@]}" | parallel -j "$num_cpus" removeOutliersFromFrame {} $imagesToMaskOutliersDir $clippingdir $imagesWithMaskedOutliersDir
      echo done > $imagesWithMaskedOutliersDone
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

    fileWithCropParameters=$dirWithCropParameters/entirecamera_"$a"_cropRegion.txt
    read row_min row_max col_min col_max < "$fileWithCropParameters"

    frameToMask=$dirOfFramesToMask/entirecamera_$a.fits
    tmpMaskFile=$dirOfFramesMasked/"maskFor"$a.fits

    astcrop $wholeMask --polygon=$col_min,$row_min:$col_max,$row_min:$col_max,$row_max:$col_min,$row_max --mode=img  -o $tmpMaskFile --quiet
    astarithmetic $frameToMask -h1 $tmpMaskFile -h1 1 eq nan where float32 -o $dirOfFramesMasked/entirecamera_$a.fits -q
    rm $tmpMaskFile
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

    astrometryTmpDir="./astrometryDiagnosisTmp"
    if ! [ -d $astrometryTmpDir ]; then mkdir $astrometryTmpDir; fi
    python3 $pythonScriptsPath/diagnosis_deltaRAdeltaDEC.py $matchCataloguesDir $output $pixelScale
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
    local calibrationBrighLimit=$7
    local calibrationFaintLimit=$8
    local numberOfApertureUnitsForCalibration=$9
    local outputDir=${10}
    local survey=${11}
    local BDIR=${12}
    local mosaicPlot=${13}
    local calibratedCataloguesDir=${14}
    local onlyPointLikeDir=${15}

    if ! [ -d $calibratedCataloguesDir ]; then mkdir $calibratedCataloguesDir; fi

    for i in $myCatalogue_nonCalibrated/*.cat; do        
        myFrame=$i
        frameNumber=$(echo "$i" | awk -F '[/]' '{print $(NF)}' | awk -F '[.]' '{print $(1)}' | awk -F '[_]' '{print $(NF)}')

        # In the nominal resolution it takes sooo long for doing this plots. So only a set of frames are used for the
        # calibration check
        if [ "$frameNumber" -gt 10 ]; then
            :
        else
            if [[ ($survey == "SPECTRA") || ("$mosaicPlot" == true) ]]; then
                referenceCatalogue=$referenceCatalogueDir/wholeFieldPhotometricCatalogue.cat
            else
                referenceCatalogue=$referenceCatalogueDir/entirecamera_$frameNumber.cat
            fi



            myCalibratedFrame=$myFrames_calibrated/entirecamera_$frameNumber.fits
            myNonCalibratedCatalogue=$myCatalogue_nonCalibrated/entirecamera_$frameNumber.fits*
            fileWithMyApertureData=$aperturesForMyData_dir/range_entirecamera_$frameNumber*

            r_myData_pix_=$(awk 'NR==1 {printf $1}' $fileWithMyApertureData)
            r_myData_pix=$(astarithmetic $r_myData_pix_ $numberOfApertureUnitsForCalibration. x -q )

     

            # raColumnName=RA
            # decColumnName=DEC
            # photometryOnImage_noisechisel -1 $calibratedCataloguesDir $myNonCalibratedCatalogue $myCalibratedFrame $r_myData_pix $calibratedCataloguesDir/$frameNumber.cat 22.5 \
            #                                 $raColumnName $decColumnName
            dataHdu=1
            columnWithXCoordForOutDataPx=1 # These numbers come from how the catalogue of the matches stars is built. This is not very clear right now, should be improved
            columnWithYCoordForOutDataPx=2
            columnWithXCoordForOutDataWCS=3
            columnWithYCoordForOutDataWCS=4

            onlyPointLikeCat=$onlyPointLikeDir/selected_entirecamera_"$frameNumber".fits_automatic.txt
            photometryOnImage_photutils -1 $calibratedCataloguesDir $myNonCalibratedCatalogue $myCalibratedFrame $r_myData_pix $calibratedCataloguesDir/$frameNumber.cat 22.5 $dataHdu \
                                        $columnWithXCoordForOutDataPx $columnWithYCoordForOutDataPx $columnWithXCoordForOutDataWCS $columnWithYCoordForOutDataWCS

            # Catalogue that contains only stars in the point-like region
            astmatch $calibratedCataloguesDir/$frameNumber.cat --hdu=1 $onlyPointLikeCat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=3/3600 --outcols=aRA,aDEC,aMAGNITUDE -o$calibratedCataloguesDir/"$frameNumber"_pointLike.cat
            astmatch $referenceCatalogue --hdu=1 $calibratedCataloguesDir/"$frameNumber"_pointLike.cat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=3/3600 --outcols=aRA,aDEC,aMAGNITUDE,bMAGNITUDE -o$calibratedCataloguesDir/"$frameNumber"_matched.cat
            rm $calibratedCataloguesDir/$frameNumber.cat $calibratedCataloguesDir/"$frameNumber"_pointLike.cat
        fi
    done


    python3 $pythonScriptsPath/diagnosis_magVsDeltaMag.py $calibratedCataloguesDir $output $outputDir $calibrationBrighLimit $calibrationFaintLimit $survey
    # rm $calibratedCataloguesDir/*.cat 
}
export -f produceCalibrationCheckPlot

produceHalfMaxRadVsMagForSingleImage() {
    local image=$1
    local outputDir=$2
    local gaiaCat=$3
    local toleranceForMatching=$4
    local pythonScriptsPath=$5
    local alternativeIdentifier=$6 # Applied when there is no number in the name
    local tileSize=$7
    local apertureUnits=$8
    local myCatDir=$9
    local brightCalibrationLimit=${10}
    local faintCalibrationLimit=${11}


    a=$( echo $image | grep -oP '\d+(?=\.fits)' )
    if ! [[ -n "$a" ]]; then
        a=$alternativeIdentifier
    fi

    # header=1
    # catalogueName=$(generateCatalogueFromImage_noisechisel $image $outputDir $a $headerToUse $tileSize $apertureUnits)
    catalogueName=$(generateCatalogueFromImage_sextractor $image $outputDir $a $apertureUnits)

    astmatch $catalogueName --hdu=1 $gaiaCat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o $outputDir/match_decals_gaia_$a.txt

    # We get the detected point-like range so we can plot it and check if everything's fine
    rangeFile=$myCatDir/"range_"$(basename $image)".txt"
    minRange=$(awk '{print $3}' $rangeFile)
    maxRange=$(awk '{print $4}' $rangeFile)

    plotXLowerLimit=0.5
    plotXHigherLimit=10
    plotYLowerLimit=14
    plotYHigherLimit=27
    python3 $pythonScriptsPath/diagnosis_halfMaxRadVsMag.py $catalogueName $outputDir/match_decals_gaia_$a.txt -1 -1 -1 $outputDir/$a.png  \
        $plotXLowerLimit $plotXHigherLimit $plotYLowerLimit $plotYHigherLimit $apertureUnits $minRange $maxRange $brightCalibrationLimit $faintCalibrationLimit

    rm $catalogueName $outputDir/match_decals_gaia_$a.txt
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
    local apertureUnits=$8
    local myCatDir=$9
    local brightestMagForCalibration=${10}
    local faintestMagForCalibration=${11}

    images=()
    for i in $imagesDir/*.fits; do
        images+=("$i")
    done

    printf "%s\n" "${images[@]}" | parallel --line-buffer -j "$num_cpus" produceHalfMaxRadVsMagForSingleImage {} $outputDir $gaiaCat $toleranceForMatching $pythonScriptsPath "-" $tileSize $apertureUnits $myCatDir $brightestMagForCalibration $faintestMagForCalibration
}
export -f produceHalfMaxRadVsMagForOurData

buildCoadd() {
    local coaddir=$1
    local coaddName=$2
    local wdir=$3
    local wonlydir=$4
    local coaddone=$5

    if ! [ -d $coaddir ]; then mkdir $coaddir; fi
    if [ -f $coaddone ]; then
            echo -e "\n\tThe first weighted (based upon std) mean of the images already done\n"
    else
            gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
            if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
                astarithmetic $(ls -v $wdir/*.fits) $(ls $wdir/*.fits | wc -l) sum -g1 --writeall -o$coaddir/"$k"_wx.fits
                astarithmetic $(ls -v $wonlydir/*.fits ) $(ls $wonlydir/*.fits | wc -l) sum -g1 --writeall -o$coaddir/"$k"_w.fits
            else
                astarithmetic $(ls -v $wdir/*.fits) $(ls $wdir/*.fits | wc -l) sum -g1  -o$coaddir/"$k"_wx.fits
                astarithmetic $(ls -v $wonlydir/*.fits ) $(ls $wonlydir/*.fits | wc -l) sum -g1  -o$coaddir/"$k"_w.fits
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
        astarithmetic $i $coadd - -o$destinationDir/$( basename $i ) -g1
    done
}
export -f subtractCoaddToFrames


computeMetricOfResiduals() {
    local dirWithFrames=$1
    local coadd=$2
    local destinationDir=$3

    fileWithMetricName="pixelsAdded.txt"
    if [ -f $destinationDir/$fileWithMetricName ]; then
        rm $destinationDir/$fileWithMetricName
    fi

    for i in $dirWithFrames/*.fits; do
        sumOfCurrentFrame=$( astarithmetic $destinationDir/$( basename $i ) -h1 sumvalue -q )
        echo $( basename $i) $sumOfCurrentFrame >> $destinationDir/$fileWithMetricName
    done
}
export -f subtractCoaddToFrames


tagIndividualResidualBasedOnPixels() {
    local individualResidual=$1
    local residualFramesTaggedDir=$2
    local individualResidualMasked=$3
    local frameNumber=$4

    numOfStdToTag=5
    median=$( astarithmetic $residualFramesTaggedDir/$individualResidualMasked medianvalue --quiet )
    std=$( astarithmetic $residualFramesTaggedDir/$individualResidualMasked stdvalue --quiet )
    upperThreshold=$( astarithmetic $median $std $numOfStdToTag x + --quiet )

    astarithmetic $individualResidual -h1 $individualResidual -h1 $upperThreshold gt $frameNumber where -o$residualFramesTaggedDir/tmp_$frameNumber.fits --quiet
    astarithmetic $residualFramesTaggedDir/tmp_$frameNumber.fits -h1 $i -h1 $upperThreshold lt 0 where -o$residualFramesTaggedDir/pxTagged_$base --quiet
    rm  $residualFramesTaggedDir/tmp_$frameNumber.fits
}
export -f tagIndividualResidualBasedOnPixels

tagIndividualResidualsBasedOnApertures() {
    local individualResidual=$1
    local residualFramesTaggedDir=$2
    local residualMask=$3
    local numberOfFWHMusedForApertures=$4
    local frameNumber=$5
    local fwhmDir=$6
    local segmentation=$7
    local residualCatalogue=$8

    fileWithFWHM=$fwhmDir/fwhm_entirecamera_$frameNumber.fits.txt
    read fwhmValue < $fileWithFWHM
    # formattedValue=$(printf "%f" "$fwhmValue") # To get rid of scientific notation
    formattedValue=$( awk -v v="$fwhmValue" 'BEGIN {printf "%f", v}' )
    aperturePx=$( echo "$formattedValue * $numberOfFWHMusedForApertures" | bc -l )

    # gnuastro code for obtaining the upper limit using the apertures
    echo "1 100 100 5 $aperturePx 0 0 1 1 1" | astmkprof --background=$individualResidual --clearcanvas --mforflatpix --type=uint8 --output=$residualFramesTaggedDir/tmp_aperture_$frameNumber.fits
    astmkcatalog $residualFramesTaggedDir/tmp_aperture_$frameNumber.fits -h1 --zeropoint=22.5 -o$residualFramesTaggedDir/sbl_$frameNumber.fits \
                    --valuesfile=$individualResidual --valueshdu=1 --upmaskfile=$residualFramesTaggedDir/$residualMask --upmaskhdu=DETECTIONS\
                    --upnsigma=5 --checkuplim=1 --upnum=1000 --ids --upperlimit-sb
    upperLimitMag=$( asttable $residualFramesTaggedDir/sbl_$frameNumber.fits -cUPPERLIMIT_SB )

    taggedResidual_tmp=$residualFramesTaggedDir/tmp_apertureTagged_$base
    taggedResidual=$residualFramesTaggedDir/apertureTagged_$base

    asttable $residualFramesTaggedDir/$residualCatalogue --range=MAGNITUDE,-999,$upperLimitMag -o$residualFramesTaggedDir/filteredCat_$base
    objectsToMask=$( asttable $residualFramesTaggedDir/filteredCat_$base --column=OBJ_ID)


    expr=""
    for label in $objectsToMask; do
        if [ -z "$expr" ]; then
            expr="o $label eq"
        else
            expr="$expr o $label eq or"
        fi
    done
    astarithmetic $individualResidual -h1 $residualFramesTaggedDir/$segmentation -h1 set-o $expr o isnotblank and $frameNumber where -o$taggedResidual_tmp 
    astarithmetic $taggedResidual_tmp -h1 $taggedResidual_tmp -h1 $frameNumber ne 0 where -o $taggedResidual

    rm $residualFramesTaggedDir/tmp_aperture_$frameNumber.fits $residualFramesTaggedDir/sbl_$frameNumber.fits \
        $residualFramesTaggedDir/filteredCat_$base $residualFramesTaggedDir/sbl_"$frameNumber"_upcheck.fits \
        $taggedResidual_tmp
}
export -f tagIndividualResidualsBasedOnApertures

prepareIndividualResidualWithTags() {
    local i=$1
    local residualFramesTaggedDir=$2
    local fwhmDir=$3
    local noisechisel_param=$4

    base=$( basename $i )
    frameNumber=$(echo "$base" | sed -E 's/.*_([0-9]+)\.fits/\1/')
    tmpMask=$(echo $base | sed 's/.fits/_mask.fits/')
    residualMasked=$(echo $base | sed 's/.fits/_masked.fits/')

    astnoisechisel $noisechisel_param $i --numthreads=$num_cpus -o $residualFramesTaggedDir/$tmpMask 
    astsegment $residualFramesTaggedDir/$tmpMask --gthresh=-15 --objbordersn=0 -o$residualFramesTaggedDir/seg_$base
    astmkcatalog $residualFramesTaggedDir/seg_$base --ids --x --y --ra --dec --magnitude --zeropoint=22.5 -o$residualFramesTaggedDir/cat_$base
    astarithmetic $i -h1 $residualFramesTaggedDir/$tmpMask -h1 1 eq nan where float32 -o $residualFramesTaggedDir/$residualMasked -quiet
    
    tagIndividualResidualBasedOnPixels $i $residualFramesTaggedDir $residualMasked $frameNumber
    numberOfFWHMusedForApertures=2
    tagIndividualResidualsBasedOnApertures $i $residualFramesTaggedDir $tmpMask $numberOfFWHMusedForApertures $frameNumber $fwhmDir seg_$base cat_$base 

    rm $residualFramesTaggedDir/$tmpMask $residualFramesTaggedDir/$residualMasked $residualFramesTaggedDir/seg_$base $residualFramesTaggedDir/cat_$base
}
export -f prepareIndividualResidualWithTags

computeSumMosaicAfterCoaddSubtractionWithTracesIndicated() {
    local residualFramesDir=$1
    local residualFramesTaggedDir=$2
    local finalResidualOutputpx=$3
    local finalResidualOutputAper=$4
    local fwhmFolder=$5
    local noisechisel_param=$6

    frames=()
    for i in $( ls $residualFramesDir/*.fits ); do
        frames+=("$i")
    done
    printf "%s\n" "${frames[@]}" | parallel -j "$num_cpus" prepareIndividualResidualWithTags {} $residualFramesTaggedDir $fwhmFolder "'$noisechisel_param'"

    astarithmetic $(ls -v $residualFramesTaggedDir/pxTagged*.fits) $(ls $residualFramesTaggedDir/pxTagged*.fits | wc -l) sum -g1 -o$finalResidualOutputpx
    astarithmetic $(ls -v $residualFramesTaggedDir/apertureTagged*.fits) $(ls $residualFramesTaggedDir/apertureTagged*.fits | wc -l) sum -g1 -o$finalResidualOutputAper
}
export -f computeSumMosaicAfterCoaddSubtractionWithTracesIndicated


changeNonNansOfFrameToOnes() {
  local a=$1
  local framesDir=$2
  local outputDir=$3

  frame=$framesDir/entirecamera_$a.fits
  output=$outputDir/exposure_tmp_$a.fits

  astarithmetic $frame $frame isnotblank 1 where --output=$output -g1
}
export -f changeNonNansOfFrameToOnes

computeExposureMap() {
    local framesDir=$1
    local exposureMapDir=$2
    local exposureMapDone=$3
    

    if ! [ -d $exposuremapDir ]; then mkdir $exposuremapDir; fi
    if [ -f $exposuremapdone ]; then
        echo -e "\n\tThe exposure map is already done\n"
    else
      
      framesToProcess=()
      for a in $(seq 1 $totalNumberOfFrames); do
        base=$framesDir/entirecamera_"$a".fits
        if [ -f $base ]; then
            framesToProcess+=("$a")
        fi
      done

      printf "%s\n" "${framesToProcess[@]}" | parallel -j "$num_cpus" changeNonNansOfFrameToOnes {} $framesDir $exposuremapDir
      gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
      if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
        astarithmetic $(ls -v $exposuremapDir/*.fits) $(ls $exposuremapDir/*.fits | wc -l) sum -g1  --writeall -o$coaddDir/exposureMap.fits
      else
        astarithmetic $(ls -v $exposuremapDir/*.fits) $(ls $exposuremapDir/*.fits | wc -l) sum -g1  --writeall -o$coaddDir/exposureMap.fits
      fi
      rm -rf $exposuremapDir
      echo done > $exposuremapdone
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
# 4.- FWHM/2 or Re


generateCatalogueFromImage_noisechisel() {
    local image=$1
    local outputDir=$2
    local a=$3

    local header=$4
    local tileSize=$5
    local apertureUnitsToCalculate=$6

    astmkprof --kernel=gaussian,1.5,3 --oversample=1 -o $outputDir/kernel_$a.fits 1>/dev/null
    astconvolve $image -h$header --kernel=$outputDir/kernel_$a.fits --domain=spatial --output=$outputDir/convolved_$a.fits 1>/dev/null
    astnoisechisel $image -h$header -o $outputDir/det_$a.fits --convolved=$outputDir/convolved_$a.fits --tilesize=$tileSize,$tileSize --numthreads=$num_cpus 1>/dev/null
    astsegment $outputDir/det_$a.fits -o $outputDir/seg_$a.fits --gthresh=-15 --objbordersn=0 1>/dev/null

    if [ $apertureUnitsToCalculate == "FWHM" ]; then
        astmkcatalog $outputDir/seg_$a.fits --x --y --ra --dec --magnitude --half-max-radius --sum --clumpscat -o $outputDir/decals_$a.txt --zeropoint=22.5 1>/dev/null
    elif [ $apertureUnitsToCalculate == "Re" ]; then
        astmkcatalog $outputDir/seg_$a.fits --x --y --ra --dec --magnitude --half-sum-radius --sum --clumpscat -o $outputDir/decals_$a.txt --zeropoint=22.5 1>/dev/null
    else
        echo "Error. Aperture Units not recognised. We should not get there never"
    fi

    rm $outputDir/kernel_$a.fits $outputDir/convolved_$a.fits $outputDir/det_$a.fits $outputDir/seg_$a.fits  $outputDir/decals_"$a"_o.txt
    mv $outputDir/decals_"$a"_c.txt $outputDir/catalogue_$a.cat
    echo $outputDir/catalogue_$a.cat
}
export -f generateCatalogueFromImage_noisechisel


generateCatalogueFromImage_sextractor(){
    local image=$1
    local outputDir=$2
    local a=$3
    local apertureUnitsToCalculate=$4

    # I specify the configuration path here because in the photometric calibration the working directory changes. This has to be changed and use the config path given in the pipeline
    cfgPath=$ROOTDIR/"$objectName"/config
    source-extractor $image -c $cfgPath/sextractor_detection.sex -CATALOG_NAME $outputDir/"$a"_tmp.cat -FILTER_NAME $cfgPath/default.conv -PARAMETERS_NAME $cfgPath/sextractor_detection.param -CATALOG_TYPE ASCII_HEAD 1>/dev/null 2>&1

    # The following code is to identify the FWHM and Re columns numbers. This is needed because it is dependant
    # on the order of the parameters in the .param file.
    headerLines=$( grep '^#' "$outputDir/"$a"_tmp.cat")
    while IFS= read -r line; do
        if [[ ("$line" == *"FWHM_IMAGE"*) ]]; then
            fwhmCol=$(echo "$line" | awk '{print $2}')
        fi

        if [[ ("$line" == *"FLUX_RADIUS"*) ]]; then
            reCol=$(echo "$line" | awk '{print $2}')
        fi
    done <<< $headerLines

    # We divide the fwhm by 2 so we have a radius
    # this is done here even if Re is chosen because then, when a column is removed, the column number changes and it's simpler to do it here
    awk -v col="$fwhmCol" '{ $col = $col / 2; print }' $outputDir/"$a"_tmp.cat > $outputDir/"$a"_tmp2.cat

    # Remove the sextractor headers to add later the noisechisel equivalents (for consistency)
    numberOfHeaders=$( grep '^#' $outputDir/"$a"_tmp2.cat | wc -l)
    sed -i "1,${numberOfHeaders}d" "$outputDir/"$a"_tmp2.cat"

    # Note. It seems like reCol and fwhmCol should go the other way around, but the awk command REMOVES the column that receives. So this is not a mistake
    if [ $apertureUnitsToCalculate == "FWHM" ]; then
        awk -v col="$reCol" '{for (i=1; i<=NF; i++) if (i != col) {printf "%s%s", $i, (i<NF || (i == NF && i != col) ? OFS : "");} print ""}' $outputDir/"$a"_tmp2.cat > $outputDir/"$a".cat
    elif [ $apertureUnitsToCalculate == "Re" ]; then
        awk -v col="$fwhmCol" '{for (i=1; i<=NF; i++) if (i != col) {printf "%s%s", $i, (i<NF || (i == NF && i != col) ? OFS : "");} print ""}' $outputDir/"$a"_tmp2.cat > $outputDir/"$a".cat
    else
        echo "Error. Aperture Units not recognised. We should not get there never"
    fi

    # Headers to mimic the noisechisel format. Change between MacOS and Linux
    if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS (BSD sed)
        sed -i '' '1s/^/# Column 1: X\
# Column 2: Y\
# Column 3: RA\
# Column 4: DEC\
# Column 5: MAGNITUDE\
# Column 6: HALF_MAX_RADIUS\
/' "$outputDir/$a.cat"
    else
        sed -i '1i# Column 6: HALF_MAX_RADIUS' $outputDir/$a.cat
        sed -i '1i# Column 5: MAGNITUDE      ' $outputDir/$a.cat
        sed -i '1i# Column 4: DEC            ' $outputDir/$a.cat
        sed -i '1i# Column 3: RA             ' $outputDir/$a.cat
        sed -i '1i# Column 2: Y              ' $outputDir/$a.cat
        sed -i '1i# Column 1: X              ' $outputDir/$a.cat
    fi

    rm $outputDir/"$a"_tmp.cat $outputDir/"$a"_tmp2.cat
    mv  $outputDir/$a.cat $outputDir/catalogue_$a.cat
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
    zp_asec=$(astarithmetic $pixelScale log10 5 x 22.5 + -q)
    expMax=$(aststatistics $exposureMap --maximum -q)
    exp_fr=$(astarithmetic $expMax $fracExpMap x -q)

    out_mask=$directoryOfImages/mask_det.fits
    astarithmetic $image -h1 $mask -hDETECTIONS 0 ne nan where -q --output=$out_mask >/dev/null 2>&1
    out_maskexp_tmp=$directoryOfImages/mask_exp_tmp.fits
    astarithmetic $out_mask $exposureMap -g1 $exp_fr lt nan where --output=$out_maskexp_tmp >/dev/null 2>&1

    # We run again noisechisel because we need the SKY_STD header. We need it because we want to use MEDSTD, it is said
    # in gnuastro documentation that it is more reliable than simply using the std of the background px of the image.
    out_maskexp=$directoryOfImages/mask_exp.fits
    astnoisechisel $out_maskexp_tmp --tilesize=20,20 -o$out_maskexp >/dev/null 2>&1
    sigma=$( astfits $out_maskexp --hdu=SKY_STD --keyvalue='MEDSTD' --quiet )

    sb_lim=$(astarithmetic $sigma 3 x $pixelScale x $areaSB / log10 -2.5 x $zp_asec + -q)

    rm $out_mask $out_maskexp_tmp $out_maskexp >/dev/null 2>&1
    echo "Limiting magnitude ($numOfSigmasForMetric sigma, $areaSB x $areaSB): $sb_lim" > "$outFile"
    echo "$sb_lim" # We need to recover the value outside for adding it to the coadd header
}
export -f limitingSurfaceBrightness


computeFWHMSingleFrame(){
    local a=$1
    local framesForFWHMDir=$2
    local fwhmdir=$3
    local headerToUse=$4
    local methodToUse=$5
    local tileSize=$6           # This parameter will only be used if the catalogue is being generated with noisechisel
    
    i=$framesForFWHMDir/$a

    # In the case of using it for Decals or Panstarrs, we need the variable survey
    if [[ "$methodToUse" == "sextractor" ]]; then
        outputCatalogue=$( generateCatalogueFromImage_sextractor $i $fwhmdir $a FWHM )
    elif [[ "$methodToUse" == "noisechisel" ]]; then
        outputCatalogue=$( generateCatalogueFromImage_noisechisel $i $fwhmdir $a $headerToUse $tileSize FWHM )
    else
        errorNumber=9
        echo "Error, method for selecting stars and the range in the calibration not recognised"
        echo "Exiting with error number: $erroNumber"
        exit $erroNumber
    fi

    astmatch $outputCatalogue --hdu=1 $BDIR/catalogs/"$objectName"_gaia.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aMAGNITUDE,aHALF_MAX_RADIUS --numthreads=$num_cpus -o$fwhmdir/match_"$a"_my_gaia.txt
    # Now we select the stars as we do for the photometry
    s=$(asttable $fwhmdir/match_"$a"_my_gaia.txt -h1 -c6 --noblank=MAGNITUDE   | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
    std=$(asttable $fwhmdir/match_"$a"_my_gaia.txt -h1 -c6 --noblank=MAGNITUDE | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
    minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
    maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)
    asttable $fwhmdir/match_"$a"_my_gaia.txt --noblank=MAGNITUDE --range=HALF_MAX_RADIUS,$minr:$maxr -c3,4,6 -c'arith $6 2 x' -o$fwhmdir/cat_fwhm_"$a".txt
    
    
    # The intermediate step with awk is because I have come across an Inf value which make the std calculus fail
    # Maybe there is some beautiful way of ignoring it in gnuastro. I didn't find int, I just clean de inf fields.
    FWHM=$(asttable $fwhmdir/cat_fwhm_"$a".txt  -c4 | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
    echo $FWHM > $fwhmdir/fwhm_"$a".txt
    rm $fwhmdir/match_"$a"_my_gaia.txt $outputCatalogue 
}
export -f computeFWHMSingleFrame
