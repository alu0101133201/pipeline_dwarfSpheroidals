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
    step=$1
    file=$2
    echo "Step: $step. Start time:  $(date +%D-%T)" >> $file
}
export -f writeTimeOfStepToFile

loadVariablesFromFile() {
  file=$1
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
        ""
        "·Keyword for the airmass:$airMassKeyWord"
        "·Keyword for date:$dateHeaderKey"
        "·Root directory to perform the reduction:$ROOTDIR"
        ""
        "·Calibration range"
        "  Bright limit:$calibrationBrightLimit:[mag]"
        "  Faint limit:$calibrationFaintLimit:[mag]"
        "·The aperture photometry for calibrating our data is in units of:$apertureUnits"
        "·The aperture photometry will be done with an aperture of:$numberOfApertureUnitsForCalibration:[$apertureUnits]"
        "·The calibration will be done using data from survey:$surveyForPhotometry"
        "·If the calibration is done with spectra, the survey to use is:$surveyForSpectra"
        "·The transmittances of the filters are in the folder:$folderWithTransmittances"
        ""
        "·Saturation threshold:$saturationThreshold:[ADU]"
        "·Gain:$gain:[e-/ADU]"
        "·Approximately size of the field:$sizeOfOurFieldDegrees:[deg]"
        "·Size of the coadd:$coaddSizePx:[px]"
        "·Vignetting threshold (mask every px with flat value below this one):$vignettingThreshold "
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
        "  If so, the window size is:$windowSize:[frames]"
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
        " "
        "·Filter:$filter"
        "·Pixel scale:$pixelScale:[arcsec/px]"
        "·Detector width:$detectorWidth:[px]"
        "·Detector height:$detectorHeight:[px]"
        " "
        "Parameters for measuring the surface brightness limit"
        "·Exp map fraction:$fractionExpMap"
        "·Area of the SB limit metric:$areaSBlimit: [arcsec]"
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
    DATEOBSValue=$1

    if [ "$DATEOBSValue" = "n/a" ]; then
        errorNumber=3
        echo -e "The file $i do not has the $dateHeaderKey, used for sorting the raw files for the pipeline"  >&2
        echo -e "Exiting with error number: $errorNumber"  >&2
        exit $errorNumber
    fi
}
export -f checkIfExist_DATEOBS

getHighestNumberFromFilesInFolder() {
    folderToCheck=$1
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
                vignettingThreshold \
                calibrationBrightLimit \
                calibrationFaintLimit \
                apertureUnits \
                numberOfApertureUnitsForCalibration \
                surveyForPhotometry \
                folderWithTransmittances \
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
                windowSize \
                halfWindowSize \
                MODEL_SKY_AS_CONSTANT \
                sky_estimation_method \
                polynomialDegree \
                filter \
                pixelScale \
                detectorWidth \
                detectorHeight \ 
                lowestScaleForIndex \
                highestScaleForIndex \ 
                solve_field_L_Param \
                solve_field_H_Param \
                solve_field_u_Param \ 
                numberOfStdForBadFrames
                fractionExpMap\
                areaSBlimit)

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
    telescope=$1
    survey=$2
    filterFolder=$3
    filterToUse=$4

    filterFileNeeded=$( checkIfAllTheTransmittancesNeededAreGiven $telescope $surveyForPhotometry $folderWithTransmittances $filter )
    checkUnitsAndConvertToCommonUnitsIfNeeded $filterFileNeeded
}
export -f checkTransmittanceFilterAndItsUnits

# In this function we check the transmittances given. Just for having everything in the same units
# We place everything into Angstroms and transmittances from 0 to 1.
checkUnitsAndConvertToCommonUnitsIfNeeded() {
    filterFile=$1
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
    telescope=$1
    survey=$2
    filterFolder=$3
    filterToUse=$4

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

# Functions used in Flat
maskImages() {
    inputDirectory=$1
    masksDirectory=$2
    outputDirectory=$3
    useCommonRing=$4
    keyWordToDecideRing=$5

    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        i=$inputDirectory/$base
        out=$outputDirectory/$base
        astarithmetic $i -h1 $masksDirectory/$base -hDETECTIONS 1 eq nan where float32 -o $out -q

        propagateKeyword $i $airMassKeyWord $out 
        # If we are not doing a normalisation with a common ring we propagate the keyword that will be used to decide
        # which ring is to be used. This way we can check this value in a comfortable way in the normalisation section
        if [ "$useCommonRing" = false ]; then
            propagateKeyword $i $keyWordToDecideRing $out
        fi
    done
}
export -f maskImages

getInitialMidAndFinalFrameTimes() {
  directoryWithNights=$1

  declare -a date_obs_array
  
  while IFS= read -r -d '' file; do
    currentDateObs=$( gethead $file DATE-OBS)

    if [[ -n "$currentDateObs" ]]; then
        unixTime=$( date -d "$currentDateObs" +"%s" )
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
    fitsFile=$1
    header=$2
    keyWord=$3
    value=$4
    comment=$5

    astfits --write=$keyWord,$value,"$comment" $fitsFile -h$header
}
export -f writeKeywordToFits

propagateKeyword() {
    image=$1
    keyWordToPropagate=$2
    out=$3

    variableToDecideRingToNormalise=$(gethead $image $keyWordToPropagate)
    eval "astfits --write=$keyWordToPropagate,$variableToDecideRingToNormalise $out -h1" 
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
    i=$1
    commonRing=$2
    doubleRing_first=$3
    doubleRing_second=$4
    useCommonRing=$5
    keyWordToDecideRing=$6
    keyWordThreshold=$7
    keyWordValueForFirstRing=$8
    keyWordValueForSecondRing=$9

    if [ "$useCommonRing" = true ]; then
            # Case when we have one common normalisation ring
            me=$(astarithmetic $i -h1 $commonRing -h1 0 eq nan where medianvalue --quiet)
    else
        # Case when we do NOT have one common normalisation ring
        # All the following logic is to decide which normalisation ring apply
        variableToDecideRingToNormalise=$(gethead $i $keyWordToDecideRing)
        firstRingLowerBound=$(echo "$keyWordValueForFirstRing - $keyWordThreshold" | bc)
        firstRingUpperBound=$(echo "$keyWordValueForFirstRing + $keyWordThreshold" | bc)
        secondRingLowerBound=$(echo "$keyWordValueForSecondRing - $keyWordThreshold" | bc)
        secondRingUpperBound=$(echo "$keyWordValueForSecondRing + $keyWordThreshold" | bc)

        if (( $(echo "$variableToDecideRingToNormalise >= $firstRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $firstRingUpperBound" | bc -l) )); then
            me=$(astarithmetic $i -h1 $doubleRing_first -h1 0 eq nan where medianvalue --quiet)
        elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
            me=$(astarithmetic $i -h1 $doubleRing_second -h1 0 eq nan where medianvalue --quiet)
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
    i=$1
    commonRing=$2
    doubleRing_first=$3
    doubleRing_second=$4
    useCommonRing=$5
    keyWordToDecideRing=$6
    keyWordThreshold=$7
    keyWordValueForFirstRing=$8
    keyWordValueForSecondRing=$9

    if [ "$useCommonRing" = true ]; then
            # Case when we have one common normalisation ring
            std=$(astarithmetic $i -h1 $commonRing -h1 0 eq nan where stdvalue --quiet)
    else
        # Case when we do NOT have one common normalisation ring
        # All the following logic is to decide which normalisation ring apply
        variableToDecideRingToNormalise=$(gethead $i $keyWordToDecideRing)
        firstRingLowerBound=$(echo "$keyWordValueForFirstRing - $keyWordThreshold" | bc)
        firstRingUpperBound=$(echo "$keyWordValueForFirstRing + $keyWordThreshold" | bc)
        secondRingLowerBound=$(echo "$keyWordValueForSecondRing - $keyWordThreshold" | bc)
        secondRingUpperBound=$(echo "$keyWordValueForSecondRing + $keyWordThreshold" | bc)

        if (( $(echo "$variableToDecideRingToNormalise >= $firstRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $firstRingUpperBound" | bc -l) )); then
            std=$(astarithmetic $i -h1 $doubleRing_first -h1 0 eq nan where stdvalue --quiet)
        elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
            std=$(astarithmetic $i -h1 $doubleRing_second -h1 0 eq nan where stdvalue --quiet)
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
    i=$1
    commonRing=$2
    doubleRing_first=$3
    doubleRing_second=$4
    useCommonRing=$5
    keyWordToDecideRing=$6
    keyWordThreshold=$7
    keyWordValueForFirstRing=$8
    keyWordValueForSecondRing=$9

    if [ "$useCommonRing" = true ]; then
            # Case when we have one common normalisation ring
            #astarithmetic $i -h1 $commonRing -h1 0 eq nan where -q -o ring_masked.fits
            skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS $commonRing)
            kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS $commonRing)
            rm ring_masked.fits
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
            rm ring_masked.fits
        elif (( $(echo "$variableToDecideRingToNormalise >= $secondRingLowerBound" | bc -l) )) && (( $(echo "$variableToDecideRingToNormalise <= $secondRingUpperBound" | bc -l) )); then
            #astarithmetic $i -h1 $doubleRing_second -h1 0 eq nan where -q -o ring_masked.fits
            skew=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i SKEWNESS $doubleRing_second)
            kurto=$(python3 $pythonScriptsPath/get_skewness_kurtosis.py $i KURTOSIS $doubleRing_second)
            rm ring_masked.fits
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
    imageDir=$1
    outputDir=$2
    useCommonRing=$3
    commonRing=$4
    doubleRing_first=$5
    doubleRing_second=$6
    keyWordToDecideRing=$7
    keyWordThreshold=$8
    keyWordValueForFirstRing=$9
    keyWordValueForSecondRing=${10}

    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        i=$imageDir/$base
        out=$outputDir/$base

        me=$(getMedianValueInsideRing $i $commonRing  $doubleRing_first $doubleRing_second $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)
        astarithmetic $i -h1 $me / -o $out
        propagateKeyword $i $airMassKeyWord $out 
    done
}
export -f normaliseImagesWithRing

calculateFlat() {
    flatName="$1"
    shift
    filesToUse="$@"
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
    normalisedDir=$1
    outputDir=$2
    doneFile=$3
    iteration=$4

    fileArray=()
    fileArray=( $(ls -v $normalisedDir/*Decals-"$filter"_n*_f*_ccd"$h".fits) )
    fileArrayLength=( $(ls -v $normalisedDir/*Decals-"$filter"_n*_f*_ccd"$h".fits | wc -l) )

    lefFlatFiles=("${fileArray[@]:0:$windowSize}")
    echo "Computing left flat - iteration $iteration"
    calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_left_ccd"$h".fits" "${lefFlatFiles[@]}"
    rightFlatFiles=("${fileArray[@]:(fileArrayLength-$windowSize):fileArrayLength}")
    echo "Computing right flat - iteration $iteration"
    calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_right_ccd"$h".fits" "${rightFlatFiles[@]}"

    echo "Computing non-common flats - iteration $iteration"
    for a in $(seq 1 $n_exp); do
        if [ "$a" -gt "$((halfWindowSize + 1))" ] && [ "$((a))" -lt "$(($n_exp - $halfWindowSize))" ]; then
            leftLimit=$(( a - $halfWindowSize - 1))
            calculateFlat "$outputDir/flat-it"$iteration"_"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits" "${fileArray[@]:$leftLimit:$windowSize}"
        fi
    done
    echo done > $doneFile
}
export -f calculateRunningFlat

divideImagesByRunningFlats(){
    imageDir=$1
    outputDir=$2
    flatDir=$3
    flatDone=$4

    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        i=$imageDir/$base
        out=$outputDir/$base

        if [ "$a" -le "$((halfWindowSize + 1))" ]; then
            flatToUse=$flatDir/flat-it*_"$filter"_n"$currentNight"_left_ccd"$h".fits
        elif [ "$a" -ge "$((n_exp - halfWindowSize))" ]; then
            flatToUse=$flatDir/flat-it*_"$filter"_n"$currentNight"_right_ccd"$h".fits
        else
            flatToUse=$flatDir/flat-it*_"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        fi
            astarithmetic $i -h1 $flatToUse -h1 / -o $out
            # This step can probably be removed
            astfits $i --copy=1 -o$out

        propagateKeyword $i $airMassKeyWord $out 
    done
    echo done > $flatDone
}
export -f divideImagesByRunningFlats

divideImagesByWholeNightFlat(){
    imageDir=$1
    outputDir=$2
    flatToUse=$3
    flatDone=$4

    for a in $(seq 1 $n_exp); do
        base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        i=$imageDir/$base
        out=$outputDir/$base

        astarithmetic $i -h1 $flatToUse -h1 / -o $out
        propagateKeyword $i $airMassKeyWord $out 
    done
    echo done > $flatDone
}
export -f divideImagesByWholeNightFlat

runNoiseChiselOnFrame() {
    baseName=$1
    inputFileDir=$2
    outputDir=$3
    noiseChiselParams=$4

    imageToUse=$inputFileDir/$baseName
    output=$outputDir/$baseName
    astnoisechisel $imageToUse $noiseChiselParams -o $output
}
export -f runNoiseChiselOnFrame

# Functions for Warping the frames
getCentralCoordinate(){
    image=$1

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

warpImage() {
    imageToSwarp=$1
    entireDir_fullGrid=$2
    entiredir=$3
    ra=$4
    dec=$5
    coaddSizePx=$6

    # ****** Decision note *******
    # We need to regrid the frames into the final coadd grid. But if we do this right now we will be processing
    # frames too big (which are mostly Nans) and the noisechisel routine takes a looot of time
    # The approach taken is to move the frame to that grid, and then crop it to the dimension of the data itself
    # We need to store both. I have tried to store the small one and then warp it again to the big grid, it's more time consuming
    # and the nan wholes grow so we end up with less light in the final coadd.

    # Parameters for identifing our frame in the full grid
    currentIndex=$(basename $imageToSwarp .fits)

    tmpFile1=$entiredir"/$currentIndex"_temp1.fits
    frameFullGrid=$entireDir_fullGrid/entirecamera_$currentIndex.fits

    # Resample into the final grid
    # Be careful with how do you have to call this package, because in the SIE sofware is "SWarp" and in the TST-ICR is "swarp"
    swarp -c $swarpcfg $imageToSwarp -CENTER $ra,$dec -IMAGE_SIZE $coaddSizePx,$coaddSizePx -IMAGEOUT_NAME $entiredir/"$currentIndex"_swarp1.fits -WEIGHTOUT_NAME $entiredir/"$currentIndex"_swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $pixelScale -PIXELSCALE_TYPE MANUAL
    
    # Mask bad pixels
    astarithmetic $entiredir/"$currentIndex"_swarp_w1.fits -h0 set-i i i 0 lt nan where -o$tmpFile1
    astarithmetic $entiredir/"$currentIndex"_swarp1.fits -h0 $tmpFile1 -h1 0 eq nan where -o$frameFullGrid

    regionOfDataInFullGrid=$(python3 $pythonScriptsPath/getRegionToCrop.py $frameFullGrid 1)
    read row_min row_max col_min col_max <<< "$regionOfDataInFullGrid"
    astcrop $frameFullGrid --polygon=$col_min,$row_min:$col_max,$row_min:$col_max,$row_max:$col_min,$row_max --mode=img  -o $entiredir/entirecamera_"$currentIndex".fits --quiet

    rm $entiredir/"$currentIndex"_swarp_w1.fits $entiredir/"$currentIndex"_swarp1.fits $tmpFile1 
}
export -f warpImage

removeBadFramesFromReduction() {
    sourceToRemoveFiles=$1
    destinationDir=$2
    badFilesWarningDir=$3
    badFilesWarningFile=$4

    filePath=$badFilesWarningDir/$badFilesWarningFile


    while IFS= read -r file_name; do
        file_name=$(basename "$file_name")
        fileName="entirecamera_${file_name%.*}".fits
        if [ -f $sourceToRemoveFiles/$fileName ]; then
            mv $sourceToRemoveFiles/$fileName $destinationDir/$fileName
        fi
    done < "$filePath"
}
export -f removeBadFramesFromReduction

# Functions for compute and subtract sky from frames
computeSkyForFrame(){
    base=$1
    entiredir=$2
    noiseskydir=$3
    constantSky=$4
    constantSkyMethod=$5
    polyDegree=$6
    inputImagesAreMasked=$7
    ringDir=$8
    useCommonRing=$9
    keyWordToDecideRing=${10}
    keyWordThreshold=${11}
    keyWordValueForFirstRing=${12}
    keyWordValueForSecondRing=${13}
    ringWidth=${14}

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
                astnoisechisel $i $noisechisel_param -o $noiseskydir/$tmpMask
                astarithmetic $i -h1 $noiseskydir/$tmpMask -h1 1 eq nan where float32 -o $noiseskydir/$tmpMaskedImage --quiet
                imageToUse=$noiseskydir/$tmpMaskedImage
                rm -f $noiseskydir/$tmpMask
            else
                imageToUse=$i
            fi

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

            me=$(getMedianValueInsideRing $imageToUse  $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)
            std=$(getStdValueInsideRing $imageToUse $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)
            read skew kurto < <(getSkewKurtoValueInsideRing $imageToUse $ringDir/$tmpRingFits "" "" true $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)

            echo "$base $me $std $skew $kurto" > $noiseskydir/$out

            rm $ringDir/$tmpRingDefinition
            rm $ringDir/$tmpRingFits

        elif [ "$constantSkyMethod" = "noisechisel" ]; then
            sky=$(echo $base | sed 's/.fits/_sky.fits/')

            # The sky substraction is done by using the --checksky option in noisechisel
            astnoisechisel $i --tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2 --checksky $noisechisel_param -o $noiseskydir/$base

            mean=$(aststatistics $noiseskydir/$sky -hSKY --sigclip-mean)
            std=$(aststatistics $noiseskydir/$sky -hSTD --sigclip-mean)
            echo "$base $mean $std" > $noiseskydir/$out
            rm -f $noiseskydir/$sky
        else
            errorNumber=6
            echo -e "\nAn invalid value for the sky_estimation_method was provided" >&2
            echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
            exit $errorNumber
        fi

    else
        # Case when we model a plane
        noiseOutTmp=$(echo $base | sed 's/.fits/_sky.fits/')
        maskTmp=$(echo $base | sed 's/.fits/_masked.fits/')
        planeOutput=$(echo $base | sed 's/.fits/_poly.fits/')
        planeCoeffFile=$(echo $base | sed 's/.fits/.txt/')

        # This conditional allows us to introduce the images already masked (masked with the mask of the coadd) in the second and next iterations
        if ! [ "$inputImagesAreMasked" = true ]; then
            astnoisechisel $i --tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2 --checksky $noisechisel_param -o $noiseskydir/$base
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
    framesToUseDir=$1
    noiseskydir=$2
    noiseskydone=$3
    constantSky=$4
    constantSkyMethod=$5
    polyDegree=$6
    inputImagesAreMasked=$7
    ringDir=$8
    useCommonRing=$9
    keyWordToDecideRing=${10}
    keyWordThreshold=${11}
    keyWordValueForFirstRing=${12}
    keyWordValueForSecondRing=${13}
    ringWidth=${14}
    
    if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
    if [ -f $noiseskydone ]; then
        echo -e "\n\tScience images are 'noisechiseled' for constant sky substraction for extension $h\n"
    else
        framesToComputeSky=()
        for a in $( ls $framesToUseDir/*.fits ); do
            base=$( basename $a )
            framesToComputeSky+=("$base")
        done
        printf "%s\n" "${framesToComputeSky[@]}" | parallel -j "$num_cpus" computeSkyForFrame {} $framesToUseDir $noiseskydir $constantSky $constantSkyMethod $polyDegree $inputImagesAreMasked $ringDir $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth
        echo done > $noiseskydone
    fi
}
export -f computeSky

subtractSkyForFrame() {
    a=$1
    directoryWithSkyValues=$2
    framesToSubtract=$3
    directoryToStoreSkySubtracted=$4
    constantSky=$5

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
    framesToSubtract=$1
    directoryToStoreSkySubtracted=$2
    directoryToStoreSkySubtracteddone=$3
    directoryWithSkyValues=$4
    constantSky=$5


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
    frame=$1
    frameBrickMapFile=$2

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
    image=$1
    gaiaCatalogue=$2
    kernel=$3
    tmpFolder=$4

    # The output of the commands are redirected to /dev/null because otherwise I cannot return the median and std.
    # Quite uncomfortable the return way of bash. Nevertheless, the error output is not modified so if an instruction fails we still get the error message.
    astconvolve $image --kernel=$kernel --domain=spatial --output=$tmpFolder/convolved.fits 1>/dev/null
    astnoisechisel $image -h1 -o $tmpFolder/det.fits --convolved=$tmpFolder/convolved.fits --tilesize=20,20 --detgrowquant=0.95 --erode=4 1>/dev/null
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
    mosaicDir=$1
    local surveyImagesDir=$2
    bricksIdentificationFile=$3
    filters=$4
    ra=$5
    dec=$6
    fieldSizeDeg=$7
    gaiaCatalogue=$8
    survey=$9
    
    echo -e "\n·Downloading ${survey} bricks"
    donwloadMosaicDone=$surveyImagesDir/done_downloads.txt

    if ! [ -d $mosaicDir ]; then mkdir $mosaicDir; fi
    if ! [ -d $surveyImagesDir ]; then mkdir $surveyImagesDir; fi
    if [ -f $donwloadMosaicDone ]; then
        echo -e "\n\tMosaic images already downloaded\n"
    else
        rm $bricksIdentificationFile # Remove the brick indentification file. This is done to avoid problems with files of previous executions        
        echo "Downloading $survey bricks for field centered at ($ra, $dec) and size $fieldSizeDeg deg; filters: " $filters
        python3 $pythonScriptsPath/downloadBricksForFrame.py $filters $surveyImagesDir $ra $dec $fieldSizeDeg $mosaicDir $bricksIdentificationFile $gaiaCatalogue $survey
        echo "done" > $donwloadMosaicDone
    fi
}
export -f downloadSurveyData

downloadSpectra() {
    mosaicDir=$1
    spectraDir=$2
    ra=$3
    dec=$4
    sizeOfOurFieldDegrees=$5
    surveyForSpectra=$6

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
    decalsImagesDir=$1
    filter1=$2
    filter2=$3

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


downloadGaiaCatalogue() {
    query=$1
    catdir=$2
    catName=$3

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

downloadIndex() {
    re=$1
    catdir=$2
    objectName=$3
    indexdir=$4

    build-astrometry-index -i $catdir/"$objectName"_Gaia_eDR3.fits -e1 \
                            -P $re \
                            -S phot_g_mean_mag \
                            -E -A RA -D  DEC\
                            -o $indexdir/index_$re.fits;
}
export -f downloadIndex

solveField() {
    i=$1
    solve_field_L_Param=$2
    solve_field_H_Param=$3
    solve_field_u_Param=$4
    ra_gal=$5
    dec_gal=$6
    confFile=$7
    astroimadir=$8

    base=$( basename $i)

    # The default sextractor parameter file is used.
    # I tried to use the one of the config directory (which is used in other steps), but even using the default one, it fails
    # Maybe a bug? I have not managed to make it work
    solve-field $i --no-plots \
    -L $solve_field_L_Param -H $solve_field_H_Param -u $solve_field_u_Param \
    --ra=$ra_gal --dec=$dec_gal --radius=3. \
    --overwrite --extension 1 --config $confFile/astrometry_$objectName.cfg --no-verify -E 1 -c 0.01 \
    --odds-to-solve 1e9 \
    --use-source-extractor --source-extractor-path=/usr/bin/source-extractor \
    -Unone --temp-axy -Snone -Mnone -Rnone -Bnone -N$astroimadir/$base ;
}
export -f solveField

runSextractorOnImage() {
    a=$1
    sexcfg=$2
    sexparam=$3
    sexconv=$4
    astroimadir=$5
    sexdir=$6
    saturationThreshold=$7
    gain=$8 

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
    frameName=$1
    dirWithbricks=$2

    funpack -O $dirWithbricks/decompressed_$frameName $dirWithbricks/$frameName 
}
export -f decompressDecalsFrame

decompressBricks() {
    dirWithBricks=$1

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
    brickName=$1
    dirWithBricks=$2

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
    dirWithBricks=$1

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
        printf "%s\n" "${brickList[@]}" | parallel -j "$num_cpus" selectStarsAndRangeForCalibrateSingleFrame {} $dirWithBricks $cataloguedir $headerWithData $methodToUse $noisechiselTileSize $apertureUnits
        echo "done" > $starSelectionDone
    fi
}
export -f selectStarsAndSelectionRangeSurvey

produceHalfMaxRadiusPlotsForDecals() {
    folderWithCatalogues=$1
    outputDir=$2
    filter=$3

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
    folderWithCatalogues=$1
    outputDir=$2
    filter=$3
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
    local numberOfApertureForRecuperateGAIA=$6
   
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
    surveyForCalibration=$1
    referenceImagesForMosaic=$2
    aperturePhotDir=$3
    filter=$4
    ra=$5
    dec=$6
    mosaicDir=$7
    selectedSurveyStarsDir=$8
    rangeUsedSurveyDir=$9
    dataPixelScale=${10}
    sizeOfOurFieldDegrees=${11} 
    gaiaCatalogue=${12}
    surveyForSpectra=${13}
    apertureUnits=${14}
    folderWithTransmittances=${15}

    if ! [ -d $mosaicDir ]; then mkdir $mosaicDir; fi

    if [[ "$surveyForCalibration" == "SPECTRA" ]]; then
        spectraDir=$mosaicDir/spectra

        writeTimeOfStepToFile "Spectra data processing" $fileForTimeStamps
        transmittanceCurveFile=$folderWithTransmittances/"$telescope"_"$filter".dat 
        prepareSpectraDataForPhotometricCalibration $spectraDir $filter $ra $dec $mosaicDir $aperturePhotDir $sizeOfOurFieldDegrees $surveyForSpectra $transmittanceCurveFile

    else
        surveyImagesDir=$mosaicDir/surveyImages
        writeTimeOfStepToFile "Survey data processing" $fileForTimeStamps

        prepareSurveyDataForPhotometricCalibration $referenceImagesForMosaic $surveyImagesDir $filter $ra $dec $mosaicDir $selectedSurveyStarsDir $rangeUsedSurveyDir \
                                            $dataPixelScale $surveyForCalibration $sizeOfOurFieldDegrees $gaiaCatalogue $aperturePhotDir $apertureUnits $folderWithTransmittances
    fi
}
export -f prepareCalibrationData

prepareSpectraDataForPhotometricCalibration() {
    spectraDir=$1
    filter=$2
    ra=$3
    dec=$4
    mosaicDir=$5
    aperturePhotDir=$6
    sizeOfOurFieldDegrees=$7
    surveyForSpectra=$8
    transmittanceCurveFile=$9

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
    referenceImagesForMosaic=$1
    surveyImagesDir=$2
    filter=$3
    ra=$4
    dec=$5
    mosaicDir=$6
    selectedSurveyStarsDir=$7
    rangeUsedSurveyDir=$8
    dataPixelScale=$9
    survey=${10}
    sizeOfOurFieldDegrees=${11} 
    gaiaCatalogue=${12}
    aperturePhotDir=${13}
    apertureUnits=${14}
    folderWithTransmittances=${15}
    

    sizeOfFieldForCalibratingPANSTARRStoGAIA=2

    echo -e "\n ${GREEN} ---Preparing ${survey} data--- ${NOCOLOUR}"
    bricksIdentificationFile=$surveyImagesDir/brickIdentification.txt

    surveyImagesDirForGaiaCalibration=$surveyImagesDir"ForGAIACalibration"


    # We work with two downloads of survey bricks. One of the field to calibrate the data to reduce and another one
    # (independent of the field to reduce) to calibrate the survey used for calibration (PANSTARRS or DECaLS) to our reference framework of GAIA

    # Could this be done more efficient and not have duplicity? yes. Would it need quite extra logic? Yes
    # Because sometimes the field for reducing the data is bigger than the field for calibrating the imaging survey, but other times
    # it is not. So it would need to do checks and stuff that we are not doing right now. Either way this is not bottle neck in the pipeline.

    if [ "$filter" = "lum" ]; then
        filters="g,r" # We download 'g' and 'r' because our images are taken with a luminance filter which response is a sort of g+r
        downloadSurveyData $mosaicDir $surveyImagesDirForGaiaCalibration $bricksIdentificationFile $filters $ra $dec $sizeOfFieldForCalibratingPANSTARRStoGAIA $gaiaCatalogue $survey  
        downloadSurveyData $mosaicDir $surveyImagesDir $bricksIdentificationFile $filters $ra $dec $sizeOfOurFieldDegrees $gaiaCatalogue $survey  
        # This step creates the images (g+r)/2. This is needed because we are using a luminance filter which is a sort of (g+r)
        # The division by 2 is because in AB system we work with Janskys, which are W Hz^-1 m^-2. So we have to give a flux per wavelenght
        # So, when we add two filters we have to take into account that we are increasing the wavelength rage. In our case, 'g' and 'r' have
        # practically the same wavelenght width, so dividing by 2 is enough
        addTwoFiltersAndDivideByTwo $surveyImagesDirForGaiaCalibration "g" "r"
        addTwoFiltersAndDivideByTwo $surveyImagesDir "g" "r"
    else 
        filters=$filter
        downloadSurveyData $mosaicDir $surveyImagesDirForGaiaCalibration $bricksIdentificationFile $filters $ra $dec $sizeOfFieldForCalibratingPANSTARRStoGAIA $gaiaCatalogue $survey
        downloadSurveyData $mosaicDir $surveyImagesDir $bricksIdentificationFile $filters $ra $dec $sizeOfOurFieldDegrees $gaiaCatalogue $survey
    fi

    # The photometric calibration is frame by frame, so we are not going to use the mosaic for calibration. 
    # This was implemented due to the LBT origin of the pipeline, but for doing it at original resolution takes time and memory so
    # it's not worth. I dont' delete it to let the option of using it
    # resampledDecalsBricks=$mosaicDir/resampled
    # buildDecalsMosaic $mosaicDir $decalsImagesDir $swarpcfg $ra $dec $filter $resampledDecalsBricks $dataPixelScale $sizeOfOurFieldDegrees   
    
    # Preprocessing needed depending on what survey we are usiing
    if [ "$survey" = "DECaLS" ]; then 
        decompressBricks $surveyImagesDirForGaiaCalibration
        decompressBricks $surveyImagesDir # I have to decompressed them, I don't know what DECaLS does exactly but otherwise I can't run sextractor directly on the downloaded bricks
                                          # I can run noisechisel, but since it is quickly to decompress I prefer to simplify the logic of the code and decompress always
    elif [ "$survey" = "PANSTARRS" ]; then
        # Panstarrs is calibrated at ZP=25. In order to be the most general possible, we transform panstarrs data so it is at zp=22.5 and we work always with the same zp
        divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDirForGaiaCalibration
        divideByExpTimeAndMoveZPForPanstarrs $surveyImagesDir
    fi
    
    methodToUse="sextractor"
    selectStarsAndSelectionRangeSurvey $surveyImagesDirForGaiaCalibration $selectedSurveyStarsDir"ForGAIACalibration" $methodToUse $survey $apertureUnits
    selectStarsAndSelectionRangeSurvey $surveyImagesDir $selectedSurveyStarsDir $methodToUse $survey $apertureUnits

    halfMaxRad_Mag_plots=$mosaicDir/halfMaxradius_Magnitude_plots
    if [ $survey == "DECaLS" ]; then
        produceHalfMaxRadiusPlotsForDecals $selectedSurveyStarsDir $halfMaxRad_Mag_plots $filter
    elif [ "$survey" == "PANSTARRS" ]; then
        produceHalfMaxRadiusPlotsForPanstarrs $selectedSurveyStarsDir $halfMaxRad_Mag_plots $filter
    fi

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
    performAperturePhotometryToBricks $surveyImagesDirForGaiaCalibration $selectedSurveyStarsDir"ForGAIACalibration" $aperturePhotDir"ForGAIACalibration" $filter $survey $numberOfApertureForRecuperateGAIA
    performAperturePhotometryToBricks $surveyImagesDir $selectedSurveyStarsDir $aperturePhotDir $filter $survey $numberOfApertureForRecuperateGAIA

    exit 0

    # As mentioned in other comments and in the README, our reference framework is gaia, so I compute any offset to the 
    # photometry of PANSTARRS and GAIA magnitudes (from spectra) and correct it
    spectraDir=$mosaicDir/gaiaSpectra
    magFromSpectraDir=$mosaicDir/magnitudesFromGaiaSpectra
    calibrationToGAIA $spectraDir $folderWithTransmittances $filter $ra $dec $mosaicDir $sizeOfFieldForCalibratingPANSTARRStoGAIA $magFromSpectraDir $aperturePhotDir"ForGAIACalibration"
    
    
    imagesHdu=1
    brickDecalsAssociationFile=$mosaicDir/frames_bricks_association.txt
    if [[ -f $brickDecalsAssociationFile ]]; then
        rm $brickDecalsAssociationFile
    fi
    python3 $pythonScriptsPath/associateDecalsBricksToFrames.py $referenceImagesForMosaic $imagesHdu $bricksIdentificationFile $brickDecalsAssociationFile $survey
}
export -f prepareSurveyDataForPhotometricCalibration

calibrationToGAIA() {
    spectraDir=$1
    folderWithTransmittances=$2
    filter=$3
    ra=$4
    dec=$5
    mosaicDir=$6
    sizeOfFieldForCalibratingPANSTARRStoGAIA=$7
    magFromSpectraDir=$8
    panstarrsCatalogueDir=$9

    brightLimitToCompareGAIAandPANSTARRS=14.5
    faintLimitToCompareGAIAandPANSTARRS=15.5

    transmittanceCurveFile="$folderWithTransmittances"/PANSTARRS_$filter.dat
    prepareSpectraDataForPhotometricCalibration $spectraDir $filter $ra $dec $mosaicDir $magFromSpectraDir $sizeOfFieldForCalibratingPANSTARRStoGAIA "GAIA" $transmittanceCurveFile
    # python3 $pythonScriptsPath/getOffsetBetweenPANSTARRSandGAIA.py $panstarrsCatalogueDir $magFromSpectraDir/wholeFieldPhotometricCatalogue.cat $brightLimitToCompareGAIAandPANSTARRS $faintLimitToCompareGAIAandPANSTARRS

    # At this point we have the catalogue from GAIA and the catalogue from PANSTARRS. It's a matter of compare them

}
export -f calibrationToGAIA

# Photometric calibration functions
# The function that is to be used (the 'public' function using OOP terminology)
# Is 'computeCalibrationFactors' and 'applyCalibrationFactors'
selectStarsAndRangeForCalibrateSingleFrame(){
    a=$1
    framesForCalibrationDir=$2
    mycatdir=$3
    headerToUse=$4
    methodToUse=$5
    tileSize=$6           # This parameter will only be used if the catalogue is being generated with noisechisel
    apertureUnits=$7


    i=$framesForCalibrationDir/$a
    ##In the case of using it for Decals or Panstarrs, we need the variable survey
    

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

    astmatch $outputCatalogue --hdu=1 $BDIR/catalogs/"$objectName"_Gaia_eDR3.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aMAGNITUDE,aHALF_MAX_RADIUS -o$mycatdir/match_"$a"_my_gaia.txt
    
    # The intermediate step with awk is because I have come across an Inf value which make the std calculus fail
    # Maybe there is some beautiful way of ignoring it in gnuastro. I didn't find int, I just clean de inf fields.
    s=$(asttable $mycatdir/match_"$a"_my_gaia.txt -h1 -c6 --noblank=MAGNITUDE   | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
    std=$(asttable $mycatdir/match_"$a"_my_gaia.txt -h1 -c6 --noblank=MAGNITUDE | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
    minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
    maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)
    echo $s $std $minr $maxr > $mycatdir/range_"$a".txt
    asttable $outputCatalogue --range=HALF_MAX_RADIUS,$minr,$maxr -o $mycatdir/selected_"$a"_automatic.txt
}
export -f selectStarsAndRangeForCalibrateSingleFrame

selectStarsAndSelectionRangeOurData() {
    iteration=$1
    framesForCalibrationDir=$2
    mycatdir=$3
    methodToUse=$4
    tileSize=$5
    apertureUnits=$6

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
        printf "%s\n" "${framesToUse[@]}" | parallel -j "$num_cpus" selectStarsAndRangeForCalibrateSingleFrame {} $framesForCalibrationDir $mycatdir $headerWithData $methodToUse $tileSize $apertureUnits
        echo done > $mycatdone
    fi
}
export -f selectStarsAndSelectionRangeOurData


matchDecalsAndSingleFrame() {
    a=$1
    myCatalogues=$2
    calibrationCatalogues=$3
    matchdir=$4
    surveyForCalibration=$5

    base="entirecamera_$a.fits"
    out=$matchdir/matched_"$base".cat
    ourDataCatalogue=$myCatalogues/entirecamera_$a.fits.cat

    tmpCatalogue=$matchdir/match-$base-tmp.cat
    out=$matchdir/match-"$base".cat

    # Choose between the whole field catalogue (spectra calibration) or the specific catalogue (imaging calibration)
    # I have to do this independently because if it is spectra data I want to preserve the spectrograph and the objectype
    if [ $surveyForCalibration == "SPECTRA" ]; then
        calibrationCat=$calibrationCatalogues/wholeFieldPhotometricCatalogue.cat
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
        calibrationCat=$calibrationCatalogues/entirecamera_$a.cat
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
    myCatalogues=$1
    calibrationCatalogues=$2
    matchdir=$3
    surveyForCalibration=$4

    matchdirdone=$matchdir/done_automatic.txt
    if ! [ -d $matchdir ]; then mkdir $matchdir; fi
    if [ -f $matchdirdone ]; then
        echo -e "\n\tMatch between decals (aperture) catalog and my (aperture) catalogs already done\n"
    else
        frameNumber=()
        for a in $(seq 1 $totalNumberOfFrames); do
            frameNumber+=("$a")
        done
        printf "%s\n" "${frameNumber[@]}" | parallel -j "$num_cpus" matchDecalsAndSingleFrame {} $myCatalogues $calibrationCatalogues $matchdir $surveyForCalibration
        echo done > $matchdirdone
    fi
}
export -f matchDecalsAndOurData

buildOurCatalogueOfMatchedSourcesForFrame() {
    a=$1
    ourDatadir=$2
    framesForCalibrationDir=$3
    mycatdir=$4
    numberOfApertureUnitsForCalibration=$5

    base="entirecamera_$a.fits"
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
    photometryOnImage_photutils $a $ourDatadir $automaticCatalogue $i $r_myData_pix $ourDatadir/$base.cat 22.5 $dataHdu \
                                $columnWithXCoordForOutDataPx $columnWithYCoordForOutDataPx $columnWithXCoordForOutDataWCS $columnWithYCoordForOutDataWCS
}
export -f buildOurCatalogueOfMatchedSourcesForFrame

buildOurCatalogueOfMatchedSources() {
    ourDatadir=$1
    framesForCalibrationDir=$2
    mycatdir=$3
    numberOfApertureUnitsForCalibration=$4

    ourDatadone=$ourDatadir/done.txt
    if ! [ -d $ourDatadir ]; then mkdir $ourDatadir; fi
    if [ -f $ourDatadone ]; then
        echo -e "\n\tAperture catalogs in our data done\n"
    else
        framesToUse=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToUse+=("$a")
        done
        printf "%s\n" "${framesToUse[@]}" | parallel -j "$num_cpus" buildOurCatalogueOfMatchedSourcesForFrame {} $ourDatadir $framesForCalibrationDir $mycatdir $numberOfApertureUnitsForCalibration
        echo done > $ourDatadone
    fi
}
export -f buildOurCatalogueOfMatchedSources

matchCalibrationStarsCatalogues() {
    matchdir2=$1
    ourDatadir=$2
    decalsdir=$3
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
    alphatruedir=$1
    matchdir=$2
    brightLimit=$3
    faintLimit=$4
    apertureCorrection=$5

    alphatruedone=$alphatruedir/done.txt
    numberOfStarsUsedToCalibrateFile=$alphatruedir/numberOfStarsUsedForCalibrate.txt

    if ! [ -d $alphatruedir ]; then mkdir $alphatruedir; fi
    if [ -f $alphatruedone ]; then
        echo -e "\n\tTrustable alphas computed for extension $h\n"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            alphaFile=alpha_$a.txt
            f=$matchdir/match-entirecamera_$a.fits.cat

            alphatruet=$alphatruedir/"$objectName"_"$filter"_"$a".txt
            asttable $f -h1 --range=MAGNITUDE_CALIBRATED,$brightLimit,$faintLimit -o$alphatruet
            if [ -z "$apertureCorrection" ]; then 
                # apertureCorrection is empty, we are calibrating with a survey
                asttable $alphatruet -h1 -c1,2,'arith $6 $8 /' -o$alphatruedir/$alphaFile
            else
                # apertureCorrection is not empty, we are calibrating with spectra
                asttable $alphatruet -h1 -c1,2,"arith \$6 \$8 $apertureCorrection / /" -o $alphatruedir/$alphaFile
            fi

            # python3 /home/sguerra/pipeline/pipelineScripts/tmp_diagnosis_distributionOfCalibrationFactorsInFrame.py $alphatruedir/$alphaFile $a
                        
            mean=$(asttable $alphatruedir/$alphaFile -c3 | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
            std=$(asttable $alphatruedir/$alphaFile -c3 | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
            echo "$mean $std" > $alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a".txt
            count=$(asttable $alphatruedir/$alphaFile -c3 | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --number)
            echo "Frame number $a: $count" >> $numberOfStarsUsedToCalibrateFile
        done
        echo done > $alphatruedone
    fi
}
export -f computeAndStoreFactors


combineDecalsCataloguesForSingleFrame() {
    outputDir=$1
    frame=$2
    bricks=$3

    catalogueName=$( echo "$frame" | awk -F'.' '{print $1}')
    firstBrick=$(echo "$bricks" | awk '{print $1}')  
    
    asttablePrompt="asttable $decalsCataloguesDir/$firstBrick.cat -o$outputDir/$catalogueName.cat"

    remainingBricks=$(echo "$bricks" | cut -d' ' -f2-)  # Get the rest of the bricks
    for brick in $remainingBricks; do
        asttablePrompt+=" --catrowfile=$decalsCataloguesDir/$brick.cat"
    done

    $asttablePrompt
}

combineDecalsBricksCataloguesForEachFrame() {
    outputDir=$1
    frameBrickAssociationFile=$2
    decalsCataloguesDir=$3

    combinationDone=$outputDir/done.txt
    if ! [ -d $outputDir ]; then mkdir $outputDir; fi
    if [ -f $combinationDone ]; then
        echo -e "\nCombination of the bricks catalogues for each frame already done\n"
    else
        while IFS= read -r line; do
            currentLine=$line
            frame=$(echo "$line" | awk '{print $1}')
            frame=$( basename $frame )
            bricks=$(echo "$line" | cut -d' ' -f2-)
            combineDecalsCataloguesForSingleFrame $outputDir $frame "$bricks"
        done < "$frameBrickAssociationFile"
        echo "done" > $combinationDone
    fi
}
export -f combineDecalsBricksCataloguesForEachFrame

computeCalibrationFactors() {
    surveyForPhotometry=$1
    iteration=$2
    imagesForCalibration=$3
    selectedDecalsStarsDir=$4
    matchdir=$5
    rangeUsedDecalsDir=$6
    mosaicDir=$7
    alphatruedir=$8
    brightLimit=$9
    faintLimit=${10}
    tileSize=${11}
    apertureUnits=${12}
    numberOfApertureUnitsForCalibration=${13}


    mycatdir=$BDIR/my-catalog-halfmaxradius_it$iteration

    methodToUse="sextractor"
    echo -e "\n ${GREEN} ---Selecting stars and range for our data--- ${NOCOLOUR}"
    selectStarsAndSelectionRangeOurData $iteration $imagesForCalibration $mycatdir $methodToUse $tileSize $apertureUnits

    ourDataCatalogueDir=$BDIR/ourData-aperture-photometry_it$iteration
    echo -e "\n ${GREEN} ---Building catalogues to our data with aperture photometry --- ${NOCOLOUR}"
    buildOurCatalogueOfMatchedSources $ourDataCatalogueDir $imagesForCalibration $mycatdir $numberOfApertureUnitsForCalibration
    
    # If we are calibrating with spectra we just have the whole catalogue of the field
    # If we are calibrating with a survey then we have a catalogue por survey's brick and we need to combine the needed bricks for build a catalogue per frame
    if [[ "$surveyForCalibration" == "SPECTRA" ]]; then
        prepareCalibrationCataloguePerFrame=$mosaicDir/aperturePhotometryCatalogues
    else
        prepareCalibrationCataloguePerFrame=$BDIR/survey-aperture-photometry_perBrick_it$iteration
        echo -e "\n ${GREEN} ---Combining decals catalogues for matching each brick --- ${NOCOLOUR}"
        combineDecalsBricksCataloguesForEachFrame $prepareCalibrationCataloguePerFrame $mosaicDir/frames_bricks_association.txt $mosaicDir/aperturePhotometryCatalogues
    fi


    echo -e "\n ${GREEN} ---Matching our aperture catalogues and Decals aperture catalogues--- ${NOCOLOUR}"
    matchDecalsAndOurData $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $matchdir $surveyForCalibration
   
    echo -e "\n ${GREEN} ---Computing calibration factors (alpha)--- ${NOCOLOUR}"
    computeAndStoreFactors $alphatruedir $matchdir $brightLimit $faintLimit $apertureCorrection
}
export -f computeCalibrationFactors

applyCalibrationFactorsToFrame() {
    a=$1
    imagesForCalibration=$2
    alphatruedir=$3
    photCorrDir=$4

    base=entirecamera_"$a".fits
    f=$imagesForCalibration/"entirecamera_$a.fits"
    alpha_cat=$alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a".txt
    alpha=$(awk 'NR=='1'{print $1}' $alpha_cat)
    astarithmetic $f -h1 $alpha x float32 -o $photCorrDir/$base
}
export -f applyCalibrationFactorsToFrame

applyCalibrationFactors() {
    imagesForCalibration=$1
    alphatruedir=$2
    photCorrDir=$3

    muldone=$photCorrDir/done.txt
    if ! [ -d $photCorrDir ]; then mkdir $photCorrDir; fi
    if [ -f $muldone ]; then
            echo -e "\n\tMultiplication for alpha in the pointings (huge grid) is done for extension $h\n"
    else
        framesToApplyFactor=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToApplyFactor+=("$a")
        done
        printf "%s\n" "${framesToApplyFactor[@]}" | parallel -j "$num_cpus" applyCalibrationFactorsToFrame {} $imagesForCalibration $alphatruedir $photCorrDir
        echo done > $muldone
    fi
}
export -f applyCalibrationFactors

# Compute the weights o the frames based on the std of the background
# In order to perform a weighted mean
computeWeightForFrame() {
    a=$1
    wdir=$2
    wonlydir=$3
    photCorrDir=$4
    noiseskydir=$5 
    iteration=$6
    minRmsFileName=$7

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
    wdir=$1
    wdone=$2
    wonlydir=$3
    wonlydone=$4
    photCorrDir=$5
    noiseskydir=$6 
    iteration=$7
    minRmsFileName=$8

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
    clippingdir=$1
    clippingdone=$2
    wdir=$3
    sigmaForStdSigclip=$4

    if ! [ -d $clippingdir ]; then mkdir $clippingdir; fi
    if [ -f $clippingdone ]; then
            echo -e "\n\tUpper and lower limits for building the masked of the weighted images already computed\n"
    else
            # Compute clipped median and std
            med_im=$clippingdir/median_image.fits
            std_im=$clippingdir/std_image.fits
            gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
            if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
                astarithmetic $(ls -v $wdir/*.fits) $(ls $wdir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-median -g1 --writeall -o$med_im
                astarithmetic $(ls -v $wdir/*.fits) $(ls $wdir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-std -g1 --writeall -o$std_im
            else
                astarithmetic $(ls -v $wdir/*.fits) $(ls $wdir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-median -g1  -o$med_im
                astarithmetic $(ls -v $wdir/*.fits) $(ls $wdir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-std -g1  -o$std_im
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
    a=$1
    mowdir=$2
    moonwdir=$3
    clippingdir=$4
    wdir=$5
    wonlydir=$6

    base=entirecamera_"$a".fits
    tmp_ab=$mowdir/"$objectName"_Decals-"$filter"_"$a"_ccd"$h"_maskabove.fits
    wom=$mowdir/$base

    astarithmetic $wdir/$base -h1 set-i i i $clippingdir/upperlim.fits -h1 gt nan where float32 -q -o $tmp_ab
    astarithmetic $tmp_ab -h1 set-i i i $clippingdir/lowerlim.fits -h1 lt nan where float32 -q -o$wom
    # save the new mask
    mask=$mowdir/"$objectName"_Decals-"$filter"_"$a"_ccd"$h"_mask.fits
    astarithmetic $wom -h1 isblank float32 -o $mask
    # mask the onlyweight image
    owom=$moonwdir/$base
    astarithmetic $wonlydir/$base $mask -g1 1 eq nan where -q float32    -o $owom

    # Remove temporary files
    rm -f $tmp_ab
    rm -f $mask
}
export -f removeOutliersFromFrame

removeOutliersFromWeightedFrames () {
  mowdone=$1
  totalNumberOfFrames=$2
  mowdir=$3
  moonwdir=$4
  clippingdir=$5
  wdir=$6
  wonlydir=$7

  if [ -f $mowdone ]; then
      echo -e "\n\tOutliers of the weighted images already masked\n"
  else
      framesToRemoveOutliers=()
      for a in $(seq 1 $totalNumberOfFrames); do
          framesToRemoveOutliers+=("$a")
      done
      printf "%s\n" "${framesToRemoveOutliers[@]}" | parallel -j "$num_cpus" removeOutliersFromFrame {} $mowdir $moonwdir $clippingdir $wdir $wonlydir
      echo done > $mowdone 
  fi
}
export -f removeOutliersFromWeightedFrames

# Functions for applying the mask of the coadd for a second iteration
cropAndApplyMaskPerFrame() {
    a=$1
    dirOfFramesToMask=$2
    dirOfFramesMasked=$3
    wholeMask=$4
    dirOfFramesFullGrid=$5


    frameToMask=$dirOfFramesToMask/entirecamera_$a.fits
    frameToObtainCropRegion=$dirOfFramesFullGrid/entirecamera_$a.fits
    tmpMaskFile=$dirOfFramesMasked/"maskFor"$a.fits

    # Parameters for identifing our frame in the full grid
    frameCentre=$( getCentralCoordinate $frameToMask )
    centralRa=$(echo "$frameCentre" | awk '{print $1}')
    centralDec=$(echo "$frameCentre" | awk '{print $2}')

    regionOfDataInFullGrid=$(python3 $pythonScriptsPath/getRegionToCrop.py $frameToObtainCropRegion 1)
    read row_min row_max col_min col_max <<< "$regionOfDataInFullGrid"
    astcrop $wholeMask --polygon=$col_min,$row_min:$col_max,$row_min:$col_max,$row_max:$col_min,$row_max --mode=img  -o $tmpMaskFile --quiet
    astarithmetic $frameToMask -h1 $tmpMaskFile -h1 1 eq nan where float32 -o $dirOfFramesMasked/entirecamera_$a.fits -q
    rm $tmpMaskFile
}
export -f cropAndApplyMaskPerFrame

# maskPointings receives the directory with the frames in the full grid because we need it in order to know the region of the full grid
# in which the specific frame is located. That is obtained by using getRegionToCrop.py frame
maskPointings() {
    entiredir_smallGrid=$1
    smallPointings_maskedDir=$2
    maskedPointingsDone=$3
    maskName=$4
    dirOfFramesFullGrid=$5

    if ! [ -d $smallPointings_maskedDir ]; then mkdir $smallPointings_maskedDir; fi
    if [ -f $maskedPointingsDone ]; then
            echo -e "\nThe masks for the pointings have been already applied\n"
    else
        framesToMask=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToMask+=("$a")
        done
        printf "%s\n" "${framesToMask[@]}" | parallel -j "$num_cpus" cropAndApplyMaskPerFrame {} $entiredir_smallGrid $smallPointings_maskedDir $maskName $dirOfFramesFullGrid
        echo done > $maskedPointingsDone 
    fi
}
export -f maskPointings

produceAstrometryCheckPlot() {
    matchCataloguesDir=$1
    pythonScriptsPath=$2
    output=$3
    pixelScale=$4

    astrometryTmpDir="./astrometryDiagnosisTmp"
    if ! [ -d $astrometryTmpDir ]; then mkdir $astrometryTmpDir; fi
    python3 $pythonScriptsPath/diagnosis_deltaRAdeltaDEC.py $matchCataloguesDir $output $pixelScale
    rm -rf $astrometryTmpDir
}
export -f produceAstrometryCheckPlot

produceCalibrationCheckPlot() {
    myCatalogue_nonCalibrated=$1
    myFrames_calibrated=$2
    aperturesForMyData_dir=$3
    referenceCatalogueDir=$4
    pythonScriptsPath=$5
    output=$6
    calibrationBrighLimit=$7
    calibrationFaintLimit=$8
    numberOfFWHMToUse=$9
    outputDir=${10}
    survey=${11}
    BDIR=${12}
    aperturecorrection=${13}

    calibratedCataloguesDir=$BDIR/calibratedCatalogues
    if ! [ -d $calibratedCataloguesDir ]; then mkdir $calibratedCataloguesDir; fi

    for i in $myCatalogue_nonCalibrated/*.cat; do
        myFrame=$i
        frameNumber=$(echo "$i" | awk -F '[/]' '{print $(NF)}' | awk -F '[.]' '{print $(1)}' | awk -F '[_]' '{print $(NF)}')

        # In the nominal resolution it takes sooo long for doing this plots. So only a set of frames are used for the
        # calibration check
        if [ "$frameNumber" -gt 5 ]; then
            :
        else
            if [ $survey == "SPECTRA" ]; then
                referenceCatalogue=$referenceCatalogueDir/wholeFieldPhotometricCatalogue.cat
            else
                referenceCatalogue=$referenceCatalogueDir/*_$frameNumber.*
            fi

            myCalibratedFrame=$myFrames_calibrated/entirecamera_$frameNumber.fits
            myNonCalibratedCatalogue=$myCatalogue_nonCalibrated/entirecamera_$frameNumber.fits*
            fileWithMyApertureData=$aperturesForMyData_dir/range_entirecamera_$frameNumber*

            r_myData_pix_=$(awk 'NR==1 {printf $1}' $fileWithMyApertureData)
            r_myData_pix=$(astarithmetic $r_myData_pix_ $numberOfFWHMToUse. x -q )

            # raColumnName=RA
            # decColumnName=DEC
            # photometryOnImage_noisechisel -1 $calibratedCataloguesDir $myNonCalibratedCatalogue $myCalibratedFrame $r_myData_pix $calibratedCataloguesDir/$frameNumber.cat 22.5 \
            #                                 $raColumnName $decColumnName
            dataHdu=1
            columnWithXCoordForOutDataPx=1 # These numbers come from how the catalogue of the matches stars is built. This is not very clear right now, should be improved
            columnWithYCoordForOutDataPx=2
            columnWithXCoordForOutDataWCS=3
            columnWithYCoordForOutDataWCS=4
            photometryOnImage_photutils -1 $calibratedCataloguesDir $myNonCalibratedCatalogue $myCalibratedFrame $r_myData_pix $calibratedCataloguesDir/"$frameNumber"_tmp.cat 22.5 $dataHdu \
                                        $columnWithXCoordForOutDataPx $columnWithYCoordForOutDataPx $columnWithXCoordForOutDataWCS $columnWithYCoordForOutDataWCS

            # To perform a fair comparison we must introduce the aperture correction as well
            if [ "$survey" == "SPECTRA" ]; then
                asttable $calibratedCataloguesDir/"$frameNumber"_tmp.cat -c1,2,3,4,5,6,7,"arith SUM $aperturecorrection /" -o$calibratedCataloguesDir/"$frameNumber"_sumCorrected.cat
                asttable $calibratedCataloguesDir/"$frameNumber"_sumCorrected.cat -c1,2,3,4,5,6,7,8,"arith ARITH_1 log10 -2.5 x 22.5 +" -o $calibratedCataloguesDir/"$frameNumber"_sumAndMagCorrected.cat
                asttable $calibratedCataloguesDir/"$frameNumber"_sumAndMagCorrected.cat -c1,2,3,4,5,"ARITH_3","ARITH_1" -o $calibratedCataloguesDir/"$frameNumber"_lackingMetaData.cat

                asttable $calibratedCataloguesDir/"$frameNumber"_lackingMetaData.cat -p4 --colmetadata=6,MAGNITUDE,mag,"MAG" \
                        --colmetadata=7,SUM,none,"sum" \
                        --output=$calibratedCataloguesDir/$frameNumber.cat
            else
                mv $calibratedCataloguesDir/"$frameNumber"_tmp.cat $calibratedCataloguesDir/$frameNumber.cat
            fi

            astmatch $referenceCatalogue --hdu=1 $calibratedCataloguesDir/$frameNumber.cat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=1/3600 --outcols=aRA,aDEC,aMAGNITUDE,bMAGNITUDE -o$calibratedCataloguesDir/"$frameNumber"_matched.cat
            rm $calibratedCataloguesDir/$frameNumber.cat $calibratedCataloguesDir/"$frameNumber"_sumCorrected.cat $calibratedCataloguesDir/"$frameNumber"_sumAndMagCorrected.cat $calibratedCataloguesDir/"$frameNumber"_lackingMetaData.cat $calibratedCataloguesDir/"$frameNumber"_tmp.cat
        fi
    done

    python3 $pythonScriptsPath/diagnosis_magVsDeltaMag.py $calibratedCataloguesDir $output $outputDir $calibrationBrighLimit $calibrationFaintLimit $survey
}
export -f produceCalibrationCheckPlot

produceHalfMaxRadVsMagForSingleImage() {
    image=$1 
    outputDir=$2
    gaiaCat=$3
    toleranceForMatching=$4
    pythonScriptsPath=$5
    alternativeIdentifier=$6 # Applied when there is no number in the name
    tileSize=$7

    a=$( echo $image | grep -oP '\d+(?=\.fits)' )
    if ! [[ -n "$a" ]]; then
        a=$alternativeIdentifier
    fi

    # header=1
    # catalogueName=$(generateCatalogueFromImage_noisechisel $image $outputDir $a $headerToUse $tileSize)
    catalogueName=$(generateCatalogueFromImage_sextractor $image $outputDir $a)

    astmatch $catalogueName --hdu=1 $gaiaCat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aX,aY,aRA,aDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o $outputDir/match_decals_gaia_$a.txt 
    
    plotXLowerLimit=0.5
    plotXHigherLimit=10
    plotYLowerLimit=12
    plotYHigherLimit=22
    python3 $pythonScriptsPath/diagnosis_halfMaxRadVsMag.py $catalogueName $outputDir/match_decals_gaia_$a.txt -1 -1 -1 $outputDir/$a.png  \
        $plotXLowerLimit $plotXHigherLimit $plotYLowerLimit $plotYHigherLimit

    rm $catalogueName $outputDir/match_decals_gaia_$a.txt 
}
export -f produceHalfMaxRadVsMagForSingleImage


produceHalfMaxRadVsMagForOurData() {
    imagesDir=$1
    outputDir=$2
    gaiaCat=$3
    toleranceForMatching=$4
    pythonScriptsPath=$5
    num_cpus=$6
    tileSize=$7

    images=()
    for i in $imagesDir/*.fits; do
        images+=("$i")
    done

    # images=("/home/sguerra/NGC598/build/photCorrSmallGrid-dir_it1/entirecamera_1.fits")
    printf "%s\n" "${images[@]}" | parallel --line-buffer -j "$num_cpus" produceHalfMaxRadVsMagForSingleImage {} $outputDir $gaiaCat $toleranceForMatching $pythonScriptsPath "-" $tileSize
}
export -f produceHalfMaxRadVsMagForOurData

buildCoadd() {
    coaddir=$1
    coaddName=$2
    mowdir=$3
    moonwdir=$4
    coaddone=$5

    if ! [ -d $coaddir ]; then mkdir $coaddir; fi
    if [ -f $coaddone ]; then
            echo -e "\n\tThe first weighted (based upon std) mean of the images already done\n"
    else
            gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
            if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
                astarithmetic $(ls -v $mowdir/*.fits) $(ls $mowdir/*.fits | wc -l) sum -g1 --writeall -o$coaddir/"$k"_wx.fits
                astarithmetic $(ls -v $moonwdir/*.fits ) $(ls $moonwdir/*.fits | wc -l) sum -g1 --writeall -o$coaddir/"$k"_w.fits
            else
                astarithmetic $(ls -v $mowdir/*.fits) $(ls $mowdir/*.fits | wc -l) sum -g1  -o$coaddir/"$k"_wx.fits
                astarithmetic $(ls -v $moonwdir/*.fits ) $(ls $moonwdir/*.fits | wc -l) sum -g1  -o$coaddir/"$k"_w.fits
            fi
            astarithmetic $coaddir/"$k"_wx.fits -h1 $coaddir/"$k"_w.fits -h1 / -o$coaddName
            echo done > $coaddone
    fi
}
export -f buildCoadd

subtractCoaddToFrames() {
    dirWithFrames=$1
    coadd=$2
    destinationDir=$3

    for i in $dirWithFrames/*.fits; do
        astarithmetic $i $coadd - -o$destinationDir/$( basename $i ) -g1
    done
}
export -f subtractCoaddToFrames

changeNonNansOfFrameToOnes() {
  a=$1
  framesDir=$2
  outputDir=$3

  frame=$framesDir/entirecamera_$a.fits
  output=$outputDir/exposure_tmp_$a.fits

  astarithmetic $frame $frame 0 gt 1 where --output=$output -g1
}
export -f changeNonNansOfFrameToOnes

computeExposureMap() {
    framesDir=$1
    exposureMapDir=$2
    exposureMapDone=$3

    if ! [ -d $exposuremapDir ]; then mkdir $exposuremapDir; fi
    if [ -f $exposuremapdone ]; then
        echo -e "\n\tThe exposure map is already done\n"
    else
      framesDir=$BDIR/pointings_fullGrid
      framesToProcess=()
      for a in $(seq 1 $totalNumberOfFrames); do
        framesToProcess+=("$a")
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
    image=$1
    outputDir=$2
    a=$3   
    header=$4
    tileSize=$5
    apertureUnitsToCalculate=$6

    astmkprof --kernel=gaussian,1.5,3 --oversample=1 -o $outputDir/kernel_$a.fits 1>/dev/null
    astconvolve $image -h$header --kernel=$outputDir/kernel_$a.fits --domain=spatial --output=$outputDir/convolved_$a.fits 1>/dev/null
    astnoisechisel $image -h$header -o $outputDir/det_$a.fits --convolved=$outputDir/convolved_$a.fits --tilesize=$tileSize,$tileSize 1>/dev/null
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
    image=$1
    outputDir=$2
    a=$3   
    apertureUnitsToCalculate=$4


    # I specify the configuration path here because in the photometric calibration the working directoy changes. This has to be changed and use the config path given in the pipeline
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
    a=$1
    directoryToWork=$2
    matchedCatalogue=$3
    imageToUse=$4
    aperture_radius_px=$5
    outputCatalogue=$6
    zeropoint=$7
    hduWithData=$8
    raColumnName=$9
    decColumnName=${10}

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
    a=$1
    directoryToWork=$2
    matchedCatalogue=$3
    imageToUse=$4
    aperture_radius_px=$5
    outputCatalogue=$6
    zeropoint=$7
    hduWithData=$8
    xColumnPx=$9
    yColumnPx=${10}
    xColumnWCS=${11}
    yColumnWCS=${12}

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
    image=$1
    mask=$2
    exposureMap=$3
    directoryOfImages=$4
    areaSB=$5
    fracExpMap=$6
    pixelScale=$7
    outFile=$8

    numOfSigmasForMetric=3

    out_mask=$directoryOfImages/mask_det.fits
    astarithmetic $image -h1 $mask -hDETECTIONS 0 ne nan where -q --output=$out_mask >/dev/null 2>&1

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
