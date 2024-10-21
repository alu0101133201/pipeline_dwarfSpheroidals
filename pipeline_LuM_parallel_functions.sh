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

measureTime() {
    time=$(date +%s)
    echo $time
}
export -f measureTime

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

writeKeywordToFits() {
    fitsFile=$1
    header=$2
    keyWord=$3
    value=$4

    astfits --write=$keyWord,$value $fitsFile -h$header
}

checkIfAllVariablesAreSet() {
    errorNumber=2
    flagToExit=""
    variablesToCheck=(objectName \
                ra_gal \
                dec_gal \
                ROOTDIR \
                airMassKeyWord \ 
                dateHeaderKey \
                coaddSizePx \
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
                numberOfStdForBadFrames)

    echo -e "\n"
    for currentVar in ${variablesToCheck[@]}; do
        [[ -z ${!currentVar} ]] && echo "${currentVar} variable not defined" && flagToExit=true
    done

    # I exit here and not when I find the variable missing because I want to show all the messages of "___ variable not defined", so the user knows all the variables that are needed
    [[ $flagToExit ]] && echo -e "Exiting with error number: $errorNumber" && exit $errorNumber
}
export -f checkIfAllVariablesAreSet




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

propagateKeyword() {
    image=$1
    keyWordToPropagate=$2
    out=$3

    variableToDecideRingToNormalise=$(gethead $image $keyWordToPropagate)
    eval "astfits --write=$keyWordToPropagate,$variableToDecideRingToNormalise $out -h1" 
}
export -f propagateKeyword


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
            errorNumber=4
            echo -e "\nMultiple normalisation ring have been tried to be used. The keyword selection value of one has not matched with the ranges provided" >&2
            echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
            exit $errorNumber 
        fi
    fi

    echo $std # This is for "returning" the value
}
export -f getStdValueInsideRing


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
    astarithmetic $filesToUse $numberOfFiles $sigmaValue $iterations sigclip-median -g1 -o $flatName
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
    SWarp -c $swarpcfg $imageToSwarp -CENTER $ra,$dec -IMAGE_SIZE $coaddSizePx,$coaddSizePx -IMAGEOUT_NAME $entiredir/"$currentIndex"_swarp1.fits -WEIGHTOUT_NAME $entiredir/"$currentIndex"_swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $pixelScale -PIXELSCALE_TYPE    MANUAL
    
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
        fileName="${file_name%.*}".fits
        mv $sourceToRemoveFiles/$fileName $destinationDir/$fileName
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
        # I cannot even create a common ring all, because they are cropped based on the number of non-nan (depending on the vignetting and how the NAN are distributed), so i create a ring per image
        # For that reason the subtraction of the background using the ring is only available if a common ring is used.
        # More logic should be implemented to use the TWO normalisation rings and recover them after the warping and cropping
        if [ "$constantSkyMethod" = "ring" ]; then
            if [ "$useCommonRing" = false ]; then
                errorNumber=5
                echo "The background estimation using multiple rings is not supported by the pipeline right now" >&2
                echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
                exit $errorNumber
            fi

            # Mask the image if needed
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
            echo "1 $half_naxis1 $half_naxis2 6 1150 1 1 1 1 1" > $ringDir/$tmpRingDefinition
            astmkprof --background=$imageToUse  -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=200 --clearcanvas -o $ringDir/$tmpRingFits $ringDir/$tmpRingDefinition

            me=$(getMedianValueInsideRing $imageToUse  $ringDir/$tmpRingFits "" "" $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)
            std=$(getStdValueInsideRing $imageToUse $ringDir/$tmpRingFits "" "" $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing)
            echo "$base $me $std" > $noiseskydir/$out

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
    

    if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
    if [ -f $noiseskydone ]; then
        echo -e "\nScience images are 'noisechiseled' for constant sky substraction for extension $h\n"
    else
        framesToComputeSky=()
        for a in $(seq 1 $totalNumberOfFrames); do
            base="entirecamera_"$a.fits
            framesToComputeSky+=("$base")
        done

        printf "%s\n" "${framesToComputeSky[@]}" | parallel -j "$num_cpus" computeSkyForFrame {} $framesToUseDir $noiseskydir $constantSky $constantSkyMethod $polyDegree $inputImagesAreMasked $ringDir $useCommonRing $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing
        echo done > $noiseskydone
    fi
}

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
        echo -e "\nSky substraction is already done for the science images for extension $h\n"
    else
    framesToSubtractSky=()
    for a in $(seq 1 $totalNumberOfFrames); do
            framesToSubtractSky+=("$a")            
    done
    printf "%s\n" "${framesToSubtractSky[@]}" | parallel -j "$num_cpus" subtractSkyForFrame {} $directoryWithSkyValues $framesToSubtract $directoryToStoreSkySubtracted $constantSky
    echo done > $directoryToStoreSkySubtracteddone
    fi
}


# Functions for decals data
# The function that is to be used (the 'public' function using OOP terminology)
# Is 'prepareDecalsDataForPhotometricCalibration'
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

downloadDecalsData() {
    referenceImagesForMosaic=$1
    mosaicDir=$2
    decalsImagesDir=$3
    frameBrickCorrespondenceFile=$4
    
    filters=$5
    ringFile=$6

    echo -e "\n-Downloading Decals bricks"

    donwloadMosaicDone=$mosaicDir/decalsImages/done_downloads.txt
    if ! [ -d $mosaicDir ]; then mkdir $mosaicDir; fi
    if ! [ -d $decalsImagesDir ]; then mkdir $decalsImagesDir; fi
    if [ -f $donwloadMosaicDone ]; then
        echo -e "\nMosaic images already downloaded\n"
    else
        rm $frameBrickCorrespondenceFile # Remove the brick map. This is done to avoid problems with files of previous executions
        for a in $(seq 1 $totalNumberOfFrames); do
            base="entirecamera_$a".fits
            echo "Downloading decals bricks for image: " $base " for filters: " $filters
            bricksOfTheFrame=$(python3 $pythonScriptsPath/downloadBricksForFrame.py $referenceImagesForMosaic/$base $ringFile $filters $decalsImagesDir)
            echo $base $bricksOfTheFrame >> $frameBrickCorrespondenceFile         # Creating the map between frames and bricks to recover it in the photometric calibration
        done
        echo done > $donwloadMosaicDone
    fi
}

add_gAndr_andDivideByTwo() {
    decalsImagesDir=$1

    addBricksDone=$mosaicDir/decalsImages/done_adding.txt
    if [ -f $addBricksDone ]; then
        echo -e "\nDecals 'g' and 'r' bricks are already added\n"
    else
        for file in "$decalsImagesDir"/*_g.fits; do
            # The following lines depend on the name of the decals images, which is defined in the python script "./decals_GetAndDownloadBricks.py"
            brickName=$(basename "$file" | cut -d '_' -f3)
            gFile="decal_image_"$brickName"_g.fits"
            rFile="decal_image_"$brickName"_r.fits"

            echo -e "Adding the files " $gFile " and the file " $rFile
            astarithmetic $decalsImagesDir/$gFile -h1 $decalsImagesDir/$rFile -h1 + 2 / -o$decalsImagesDir/"decal_image_"$brickName"_g+r_div2.fits"
        done
        echo done > $addBricksDone
    fi

}

buildDecalsMosaic() {
    # We only need the mosaic in order to download the gaia catalogue. That's why downgrade the bricks
    # Values for original decals resolution

    # decalsPxScale=0.262 # arcsec/px
    # mosaicSize=38000 # For our field, around 2.5deg

    # Values used when I was just downgrading a factor of 4
    # scaleFactor=0.25
    # decalsPxScale=1.048
    # mosaicSize=9500

    mosaicDir=$1
    decalsImagesDir=$2
    swarpcfg=$3
    ra=$4
    dec=$5
    filter=$6
    swarpedImagesDir=$7
    dataPixelScale=$8 # Pixel scale of our data. In order to do realisticly we should downgrade decals data to the same resolution as our data
    sizeOfOurFieldDegrees=$9 # Estimation of how big is our field

    originalDecalsPxScale=0.262 # arcsec/px

    buildMosaicDone=$swarpedImagesDir/done_t.xt
    if ! [ -d $swarpedImagesDir ]; then mkdir $swarpedImagesDir; fi
    if [ -f $buildMosaicDone ]; then
        echo -e "\nMosaic already built\n"
    else
        decalsPxScale=$dataPixelScale
        mosaicSize=$(echo "($sizeOfOurFieldDegrees * 3600) / $decalsPxScale" | bc)
        # scaleFactor=$(echo "scale=3; $originalDecalsPxScale / $dataPixelScale" | bc)
        scaleFactor=$(awk "BEGIN {print $originalDecalsPxScale / $dataPixelScale}")
        
        if [ "$filter" = "lum" ]; then
            bricks=$(ls -v $decalsImagesDir/*_g+r_div2.fits)
        else
            bricks=$(ls -v $decalsImagesDir/*$filter.fits)
        fi

        for a in $bricks; do
            downSampledImages="$swarpedImagesDir/originalGrid_$(basename $a)"
            astwarp $a --scale=$scaleFactor -o $downSampledImages
            SWarp -c $swarpcfg $downSampledImages -CENTER $ra,$dec -IMAGE_SIZE $mosaicSize,$mosaicSize -IMAGEOUT_NAME $swarpedImagesDir/swarp1.fits \
                                -WEIGHTOUT_NAME $swarpedImagesDir/swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $decalsPxScale -PIXELSCALE_TYPE MANUAL
            astarithmetic $swarpedImagesDir/swarp_w1.fits -h0 set-i i i 0 lt nan where -otemp1.fits
            astarithmetic $swarpedImagesDir/swarp1.fits -h0 temp1.fits -h1 0 eq nan where -o$swarpedImagesDir/commonGrid_"$(basename $a)"
            
        done
        rm $swarpedImagesDir/swarp_w1.fits $swarpedImagesDir/swarp1.fits ./temp1.fits

        sigma=2
        astarithmetic $(ls -v $swarpedImagesDir/commonGrid*.fits) $(ls -v $swarpedImagesDir/commonGrid*.fits | wc -l) -g1 $sigma 0.2 sigclip-median -o $mosaicDir/mosaic.fits
        echo done > $buildMosaicDone
    fi
}

downloadGaiaCatalogueForField() {
    mosaicDir=$1

    mos=$mosaicDir/mosaic.fits
    ref=$mos

    retcatdone=$BDIR/downloadedGaia_done_"$n".txt
    if [ -f $retcatdone ]; then
            echo -e "\ngaia dr3 catalog retrived\n"
    else
        astquery gaia --dataset=edr3 --overlapwith=$ref --column=ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error    -o$BDIR/catalogs/"$objectName"_Gaia_eDR3_.fits
        asttable $BDIR/catalogs/"$objectName"_Gaia_eDR3_.fits -c1,2,3 -c'arith $4 abs' -c'arith $5 3 x' -c'arith $6 abs' -c'arith $7 3 x' -c'arith $8 abs' -c'arith $9 3 x' --noblank=4 -otmp.txt

        # Here I demand that the gaia object fulfills simultaneously that:
        # 1.- Parallax > 3 times its error
        # 2.- Proper motion (ra) > 3 times its error
        # 3.- Proper motion (dec) > 3 times its error
        asttable tmp.txt -c1,2,3 -c'arith $4 $4 $5 gt 1000 where' -c'arith $6 $6 $7 gt 1000 where' -c'arith $8 $8 $9 gt 1000 where'    -otest_.txt
        asttable test_.txt -c1,2,3 -c'arith $4 $5 + $6 +' -otest1.txt
        asttable test1.txt -c1,2,3 --range=ARITH_2,2999,3001 -o $BDIR/catalogs/"$objectName"_Gaia_eDR3.fits

        rm test1.txt tmp.txt $BDIR/catalogs/"$objectName"_Gaia_eDR3_.fits test_.txt
        echo done > $retcatdone
    fi
}

checkIfSameBricksHaveBeenComputed() {
    currentBrick=$1
    bricksToLookFor=$2
    frameBrickCorrespondenceFile=$3

    # I build an array with the bricks to look for because the function receives an string
    bricksToLookForArray=()
    for i in $bricksToLookFor; do
        bricksToLookForArray+=("$i")
    done

    while read -r line; do
        # I read the frame number and its associatd bricks
        imageNumber=$(echo "$line" | awk -F '[_.]' '{print $2}')
        string_list=$(echo "$line" | sed -E "s/.*\[(.*)\]/\1/" | tr -d "'")
        IFS=', ' read -r -a readBricksArray <<< "$string_list"


        # Looking for coincidences
        if [[ ${#readBricksArray[@]} -eq 4 && ${#bricksToLookForArray[@]} -eq 4 ]]; then
            all_match=true

            for i in "${!bricksToLookForArray[@]}"; do
                if [[ "${bricksToLookForArray[$i]}" != "${readBricksArray[$i]}" ]]; then
                    all_match=false
                    break
                fi      
            done

            # I also check if the number is less, because I'm processing the files in numerical order so if is not less is not going to be computed yet
            if [[ "$all_match" = true && "$imageNumber" -lt "$currentBrick" ]]; then
                echo $imageNumber
                exit
            fi
        fi
    done < $frameBrickCorrespondenceFile
    echo 0
}
export -f checkIfSameBricksHaveBeenComputed

selectStarsAndSelectionRangeDecals() {
    framesForCalibrationDir=$1
    mosaicDir=$2
    decalsImagesDir=$3
    frameBrickCorrespondenceFile=$4
    selectedDecalsStarsDir=$5
    rangeUsedDecalsDir=$6
    filter=$7
    downSampleDecals=$8

    tmpFolder="./tmpFilesForPhotometricCalibration"
    selectDecalsStarsDone=$selectedDecalsStarsDir/automaticSelection_done.txt
    brickCombinationsDir=$mosaicDir/combinedBricksForImages # I save the combination of the bricks for calibrating each image because I will need it in future steps of the calibration

    if ! [ -d $rangeUsedDecalsDir ]; then mkdir $rangeUsedDecalsDir; fi
    if ! [ -d $brickCombinationsDir ]; then mkdir $brickCombinationsDir; fi
    if ! [ -d $selectedDecalsStarsDir ]; then mkdir $selectedDecalsStarsDir; fi
    if ! [ -d $tmpFolder ]; then mkdir $tmpFolder; fi

    if [ -f $selectDecalsStarsDone ]; then
            echo -e "\nDecals bricks and stars for doing the photometric calibration are already selected for each frame\n"
    else
        astmkprof --kernel=gaussian,1.5,3 --oversample=1 -o ./kernel.fits 1>/dev/null
        for a in $(seq 1 $totalNumberOfFrames); do
            base="entirecamera_"$a.fits
            bricks=$( getBricksWhichCorrespondToFrame $framesForCalibrationDir/$base $frameBrickCorrespondenceFile )
            echo "Image $a; Bricks $bricks"


            # The following conditional is not to compute twice the same decals frames
            alreadyComputedBrick=$( checkIfSameBricksHaveBeenComputed $a "$bricks" $frameBrickCorrespondenceFile )
            if [ "$alreadyComputedBrick" -gt 0 ]; then
                cp $selectedDecalsStarsDir/selected_decalsStarsForFrame_"$alreadyComputedBrick".txt $selectedDecalsStarsDir/selected_decalsStarsForFrame_"$a".txt 
                cp $rangeUsedDecalsDir/selected_rangeForFrame_"$alreadyComputedBrick".txt $rangeUsedDecalsDir/selected_rangeForFrame_"$a".txt
                cp $brickCombinationsDir/combinedBricks_"$(basename $alreadyComputedBrick)".fits $brickCombinationsDir/combinedBricks_"$a".fits
            else
                images=()
                for currentBrick in $bricks; do # I have to go through all the bricks for each frame
                    # Remark that the ' ' at the end of the strings is needed. Otherwise when calling $images in the SWarp the names are not separated
                    if [ "$filter" = "lum" ]; then
                        images+=$downSampleDecals/"originalGrid_decal_image_"$currentBrick"_g+r_div2.fits "
                    else
                        images+=$downSampleDecals/"originalGrid_decal_image_"$currentBrick"_$filter.fits "
                    fi
                done

                # Combining the downsampled decals image in one
                SWarp $images -IMAGEOUT_NAME $tmpFolder/swarp1_$a.fits -WEIGHTOUT_NAME $tmpFolder/swarp_w1_$a.fits 
                astarithmetic $tmpFolder/swarp_w1_$a.fits -h0 set-i i i 0 lt nan where -o$tmpFolder/temp1_$a.fits
                astarithmetic $tmpFolder/swarp1_$a.fits -h0 $tmpFolder/temp1_$a.fits -h1 0 eq nan where -o$brickCombinationsDir/combinedBricks_"$(basename $a)".fits

                astconvolve $brickCombinationsDir/combinedBricks_"$(basename $a)".fits --kernel=./kernel.fits --domain=spatial --output=$tmpFolder/convolved_$a.fits
                # In the following noisechisel I have reduced the tile because with higher tiles it was not finding enough neighbors
                astnoisechisel $brickCombinationsDir/combinedBricks_"$(basename $a)".fits -h1 -o $tmpFolder/det_$a.fits --convolved=$tmpFolder/convolved_$a.fits --tilesize=10,10
                astsegment $tmpFolder/det_$a.fits -o $tmpFolder/seg_$a.fits --snquant=0.1 --gthresh=-10 --objbordersn=0    --minriverlength=3
                astmkcatalog $tmpFolder/seg_$a.fits --ra --dec --magnitude --half-max-radius --sum --clumpscat -o $tmpFolder/decals_$a.txt --zeropoint=22.5
                astmatch $tmpFolder/decals_"$a"_c.txt --hdu=1 $BDIR/catalogs/"$objectName"_Gaia_eDR3.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=bRA,bDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o $tmpFolder/match_decals_gaia_$a.txt 1>/dev/null

                # The intermediate step with awk is because I have come across an Inf value which make the std calculus fail
                # Maybe there is some beautiful way of ignoring it in gnuastro. I didn't find int, I just clean de inf fields.
                s=$(asttable $tmpFolder/match_decals_gaia_$a.txt -h1 -c3 --noblank=MAGNITUDE   | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
                std=$(asttable $tmpFolder/match_decals_gaia_$a.txt -h1 -c3 --noblank=MAGNITUDE | awk '{for(i=1;i<=NF;i++) if($i!="inf") print $i}' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
                minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
                maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)

                echo $s $std $minr $maxr > $rangeUsedDecalsDir/selected_rangeForFrame_"$a".txt
                asttable $tmpFolder/decals_"$a"_c.txt --range=HALF_MAX_RADIUS,$minr,$maxr -o $selectedDecalsStarsDir/selected_decalsStarsForFrame_"$a".txt
                rm $tmpFolder/*
            fi
        done
        rmdir $tmpFolder
        rm ./kernel.fits
        echo done > $selectDecalsStarsDone
    fi
}
export -f selectStarsAndSelectionRangeDecals

prepareDecalsDataForPhotometricCalibration() {
    referenceImagesForMosaic=$1
    decalsImagesDir=$2
    filter=$3
    ra=$4
    dec=$5
    mosaicDir=$6
    selectedDecalsStarsDir=$7
    rangeUsedDecalsDir=$8
    dataPixelScale=$9

    echo -e "\n ${GREEN} ---Preparing Decals data--- ${NOCOLOUR}"

    frameBrickCorrespondenceFile=$decalsImagesDir/frameBrickMap.txt
    ringFile=$BDIR/ring/ring.txt

    # This first steps donwloads the decals frames that are needed for calibrating each of our frames
    # The file "frameBrickCorrespondenceFile" will store the correspondence between our frames and decals bricks
    # This way we can easily access to this relation in orther to do the photometric calibration

    # If the images are donwloaded but the done.txt file is no present, the images won't be donwloaded again but
    # the step takes a while even if the images are already downloaded because we have to do a query to the database
    # in order to obtain the brickname and check if it is already downloaded or no


    if [ "$filter" = "lum" ]; then
        filters="g,r" # We download 'g' and 'r' because our images are taken with a luminance filter which response is a sort of g+r
        downloadDecalsData $referenceImagesForMosaic $mosaicDir $decalsImagesDir $frameBrickCorrespondenceFile $filters $ringFile

        # This step creates the images (g+r)/2. This is needed because we are using a luminance filter which is a sort of (g+r)
        # The division by 2 is because in AB system we work with Janskys, which are W Hz^-1 m^-2. So we have to give a flux per wavelenght
        # So, when we add two filters we have to take into account that we are increasing the wavelength rage. In our case, 'g' and 'r' have
        # practically the same wavelenght width, so dividing by 2 is enough
        add_gAndr_andDivideByTwo $decalsImagesDir
    else 
        filters=$filter
        downloadDecalsData $referenceImagesForMosaic $mosaicDir $decalsImagesDir $frameBrickCorrespondenceFile $filters $ringFile
    fi

    # The photometric calibration is frame by frame, so we are not going to use the mosaic for calibration. But we build it anyway to retrieve in an easier way
    # the Gaia data of the whole field.
    sizeOfOurFieldDegrees=2.5 # Estimation of the size of the field for building the mosaic
    buildDecalsMosaic $mosaicDir $decalsImagesDir $swarpcfg $ra $dec $filter $mosaicDir/downSampled $dataPixelScale $sizeOfOurFieldDegrees

    echo -e "\n ${GREEN} ---Downloading GAIA catalogue for our field --- ${NOCOLOUR}"
    downloadGaiaCatalogueForField $mosaicDir

    # First of all remember that we need to do the photometric calibration frame by frame.
    # For each frame of the data, due to its large field, we have multiple (in our case 4 - defined in "downloadBricksForFrame.py") decals bricks
    # They may have been taken at different moments with different conditions so we have to process them one by one
    # CORRECTION. This was our initial though, but actually the bricks and the detectors of decals are not the same, so we are already mixing night conditions
    # Additionally, due to the restricted common range, one brick is not enough for obtaining a reliable calibration so now we use the four bricks
    echo -e "\n ${GREEN} --- Selecting stars and star selection range for Decals--- ${NOCOLOUR}"
    selectStarsAndSelectionRangeDecals $referenceImagesForMosaic $mosaicDir $decalsImagesDir $frameBrickCorrespondenceFile $selectedDecalsStarsDir $rangeUsedDecalsDir $filter $mosaicDir/downSampled
}
export -f prepareDecalsDataForPhotometricCalibration

# Photometric calibration functions
# The function that is to be used (the 'public' function using OOP terminology)
# Is 'computeCalibrationFactors' and 'applyCalibrationFactors'
selectStarsAndRangeForCalibrateSingleFrame(){
    a=$1
    framesForCalibrationDir=$2
    mycatdir=$3

    base="entirecamera_"$a.fits
    i=$framesForCalibrationDir/$base

    astnoisechisel $i -h1 -o det_"$a".fits
    astsegment det_"$a".fits -o seg_"$a".fits --snquant=0.1 --gthresh=-10 --objbordersn=0    --minriverlength=3
    astmkcatalog seg_"$a".fits --ra --dec --magnitude --half-max-radius --sum --clumpscat -o $mycatdir/"$base".txt
    astmatch $BDIR/catalogs/"$objectName"_Gaia_eDR3.fits --hdu=1 $mycatdir/"$base"_c.txt --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aRA,aDEC,bMAGNITUDE,bHALF_MAX_RADIUS -o$mycatdir/match_"$base"_my_gaia.txt

    s=$(asttable $mycatdir/match_"$base"_my_gaia.txt -h1 -c4 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
    std=$(asttable $mycatdir/match_"$base"_my_gaia.txt -h1 -c4 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
    minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
    maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)
    echo $s $std $minr $maxr > $mycatdir/range1_"$base".txt
    asttable $mycatdir/"$base"_c.txt    --range=HALF_MAX_RADIUS,$minr,$maxr -o $mycatdir/selected_"$base"_automatic.txt
    rm det_"$a".fits seg_"$a".fits
}
export -f selectStarsAndRangeForCalibrateSingleFrame

selectStarsAndSelectionRangeOurData() {
    iteration=$1
    framesForCalibrationDir=$2
    mycatdir=$3

    mycatdone=$mycatdir/done_ccd"$h".txt
    if ! [ -d $mycatdir ]; then mkdir $mycatdir; fi
    if [ -f $mycatdone ]; then
            echo -e "\nSources for photometric calibration are already extracted for my image\n"
    else
        framesToUse=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToUse+=("$a")
        done
        printf "%s\n" "${framesToUse[@]}" | parallel -j "$num_cpus" selectStarsAndRangeForCalibrateSingleFrame {} $framesForCalibrationDir $mycatdir
        echo done > $mycatdone
    fi
}

matchDecalsAndOurData() {
    iteration=$1
    selectedDecalsStarsDir=$2
    mycatdir=$3
    matchdir2=$4

    matchdir2done=$matchdir2/done_automatic_ccd"$h".txt
    if ! [ -d $matchdir2 ]; then mkdir $matchdir2; fi
    if [ -f $matchdir2done ]; then
        echo -e "\nMatch between decals (automatic) catalog and my (automatic) catalogs already done\n"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            base="entirecamera_$a.fits"
            out=$matchdir2/match-decals-"$base".cat

            # match the automatics catalogs THIS WAY I SELECT NON SATURATED
            out_auto=$matchdir2/match-decals-"$base"_automatic.cat
            astmatch $selectedDecalsStarsDir/selected_decalsStarsForFrame_"$a".txt --hdu=1 $mycatdir/selected_"$base"_automatic.txt --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aRA,aDEC,aMAGNITUDE,aHALF_MAX_RADIUS,bMAGNITUDE,bHALF_MAX_RADIUS -o$out_auto
        done
        echo done > $matchdir2done
    fi
}
export -f matchDecalsAndOurData

buildDecalsCatalogueOfMatchedSources() {
    decalsdir=$1
    rangeUsedDecalsDir=$2
    matchdir2=$3
    mosaicDir=$4
    decalsImagesDir=$5

    decalsdone=$decalsdir/done__ccd"$h".txt
    if ! [ -d $decalsdir ]; then mkdir $decalsdir; fi
    if [ -f $decalsdone ]; then
        echo -e "\nDecals: catalogue for the calibration stars already built\n"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            # I have to take 2 the FWHM (half-max-rad)
            # It is already saved as mean value of the point-like sources
            r_decals_pix_=$(awk 'NR==1 {printf $1}' $rangeUsedDecalsDir/"selected_rangeForFrame_"$a".txt")
            r_decals_pix=$(astarithmetic $r_decals_pix_ 2. x -q )

            base="entirecamera_$a.fits"
            out=$matchdir2/match-decals-"$base"_automatic.cat

            decalsCombinedBricks=$mosaicDir/combinedBricksForImages/combinedBricks_$a.fits
            #This seems extremely inneficient but maybe is the way idk
            asttable $out -hSOURCE_ID -cRA,DEC | awk '!/^#/{print NR, $1, $2, 5, '$r_decals_pix', 0, 0, 1, NR, 1}' > apertures_decals.txt
            astmkprof apertures_decals.txt --background=$decalsCombinedBricks --backhdu=1 \
                    --clearcanvas --replace --type=int16 --mforflatpix \
                    --mode=wcs --output=apertures_decals.fits
            astmkcatalog apertures_decals.fits -h1 --zeropoint=22.5 \
                            --valuesfile=$decalsCombinedBricks --valueshdu=1 \
                            --ids --ra --dec --magnitude --sum \
                            --output=$decalsdir/decals_"$base".cat
            rm apertures_decals.txt
            rm apertures_decals.fits
        done
        echo done > $decalsdone
    fi
}

buildOurCatalogueOfMatchedSources() {
    ourDatadir=$1
    framesForCalibrationDir=$2
    matchdir2=$3
    mycatdir=$4

    ourDatadone=$ourDatadir/done_"$filter"_ccd"$h".txt
    if ! [ -d $ourDatadir ]; then mkdir $ourDatadir; fi
    if [ -f $ourDatadone ]; then
        echo -e "\nAperture catalogs in our data done\n"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            base="entirecamera_$a.fits"
            i=$framesForCalibrationDir/$base
            out=$matchdir2/match-decals-"$base"_automatic.cat

            r_myData_pix_=$(awk 'NR==1 {printf $1}' $mycatdir/range1_"$base".txt)
            r_myData_pix=$(astarithmetic $r_myData_pix_ 2. x -q )


            asttable $out -hSOURCE_ID -cRA,DEC | awk '!/^#/{print NR, $1, $2, 5, '$r_myData_pix', 0, 0, 1, NR, 1}' > apertures.txt
            astmkprof apertures.txt --background=$i --backhdu=1 \
                    --clearcanvas --replace --type=int16 --mforflatpix \
                    --mode=wcs --output=aperture_myData.fits

            astmkcatalog aperture_myData.fits -h1 --zeropoint=0 \
                            --valuesfile=$i --valueshdu=1 \
                            --ids --ra --dec --magnitude --sum \
                            --output=$ourDatadir/$base.cat
            # asttable $lbtdir/$base_.fits -h1 --noblank=MAGNITUDE -o$lbtdir/$base.cat
        done
        echo done > $ourDatadone
    fi
}

matchCalibrationStarsCatalogues() {
    matchdir2=$1
    ourDatadir=$2
    decalsdir=$3
    matchdir2done=$matchdir2/done_aperture_ccd"$h".txt

    if [ -f $matchdir2done ]; then
        echo -e "\nMatch between decals (aperture) catalog and our (aperture) catalogs done for extension $h\n"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            base="entirecamera_$a.fits"
            i=$ourDatadir/"$base".cat
            out=$matchdir2/"$objectName"_Decals-"$filter"_"$a"_ccd"$h".cat
            astmatch $decalsdir/decals_"$base".cat --hdu=1 $i --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aRA,aDEC,aMAGNITUDE,aSUM,bMAGNITUDE,bSUM -o$out
        done
        echo done > $matchdir2done
    fi
}

computeAndStoreFactors() {
    alphatruedir=$1
    matchdir2=$2
    brightLimit=$3
    faintLimit=$4

    alphatruedone=$alphatruedir/done_ccd"$h".txt

    if ! [ -d $alphatruedir ]; then mkdir $alphatruedir; fi
    if [ -f $alphatruedone ]; then
        echo -e "\nTrustable alphas computed for extension $h\n"
    else
        for a in $(seq 1 $totalNumberOfFrames); do
            base="$a".fits
            f=$matchdir2/"$objectName"_Decals-"$filter"_"$a"_ccd"$h".cat

            alphatruet=$alphatruedir/"$objectName"_Decals-"$filter"_"$a"_ccd"$h".txt
            asttable $f -h1 --range=MAGNITUDE,$brightLimit,$faintLimit -o$alphatruet
            asttable $alphatruet -h1 -c1,2,3,'arith $4 $6 /' -o$alphatruedir/$base

            mean=$(asttable $alphatruedir/$base -c'ARITH_1' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
            std=$(asttable $alphatruedir/$base -c'ARITH_1' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
            echo "$mean $std" > $alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a"_ccd"$h".txt
            count=$(asttable $alphatruedir/$base -c'ARITH_1' | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --number)
            echo "$count" > $alphatruedir/numberOfStarsUsedForCalibrate_$a.txt
        done
        echo done > $alphatruedone
    fi
}

computeCalibrationFactors() {
    iteration=$1
    imagesForCalibration=$2
    selectedDecalsStarsDir=$3
    rangeUsedDecalsDir=$4
    mosaicDir=$5
    decalsImagesDir=$6
    alphatruedir=$7

    mycatdir=$BDIR/my-catalog-halfmaxradius_it$iteration


    # EXPLANATION AND TO DO
    # The next step performs an analog process to the one applied to decals (selection of stars and saving our star range)
    # But this step here is paralellised. This is because paralellising the step in the decals section is not straight forward
    # because I keep a record of the already studied bricks, so we are accessing a common file and two processes could work with
    # the same bricks and to paralellise it we need to give it a thought
    # Here we just have to apply the process to every single frame so we can paralellise it easily
    echo -e "\n ${GREEN} ---Selecting stars and range for our data--- ${NOCOLOUR}"
    # For fornax (around 490 frames). Deimos, 20 cores -> 40 min
    selectStarsAndSelectionRangeOurData $iteration $imagesForCalibration $mycatdir


    matchdir2=$BDIR/match-decals-myData_it$iteration

    echo -e "\n ${GREEN} ---Matching our data and Decals--- ${NOCOLOUR}"
    matchDecalsAndOurData $iteration $selectedDecalsStarsDir $mycatdir $matchdir2 >> tmp.txt

    decalsdir=$BDIR/decals-aperture-catalog_it$iteration
    echo -e "\n ${GREEN} ---Building Decals catalogue for (matched) calibration stars--- ${NOCOLOUR}"
    buildDecalsCatalogueOfMatchedSources $decalsdir $rangeUsedDecalsDir $matchdir2 $mosaicDir $decalsImagesDir

    ourDataCatalogueDir=$BDIR/ourData-catalogs-apertures_it$iteration
    echo -e "\n ${GREEN} ---Building our catalogue for calibration stars--- ${NOCOLOUR}"
    buildOurCatalogueOfMatchedSources $ourDataCatalogueDir $imagesForCalibration $matchdir2 $mycatdir

    
    echo -e "\n ${GREEN} ---Matching calibration stars catalogues--- ${NOCOLOUR}"
    matchCalibrationStarsCatalogues $matchdir2 $ourDataCatalogueDir $decalsdir

    echo -e "\n ${GREEN} ---Computing calibration factors (alpha)--- ${NOCOLOUR}"

    # ****** Decision note *******
    # The reasonable range is decideed based on the half-max-sum/mag, on the decalsmag/(decalsmag - ourmag)
    # Then we have explored this reasonable range. The bright limit is the one that it is (before saturation)
    # but the faint limit has been explored from 18 to 20. From 18.5 and fainter we do not see any difference
    # on the images or on the std-airmass plot, so we use this range.
    brightLimit=16
    faintLimit=18.5
    computeAndStoreFactors $alphatruedir $matchdir2 $brightLimit $faintLimit
}
export -f computeCalibrationFactors

applyCalibrationFactorsToFrame() {
    a=$1
    imagesForCalibration=$2
    alphatruedir=$3
    photCorrDir=$4

    h=0
    base=entirecamera_"$a".fits
    f=$imagesForCalibration/"entirecamera_$a.fits"
    alpha_cat=$alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a"_ccd"$h".txt
    alpha=$(awk 'NR=='1'{print $1}' $alpha_cat)
    echo astarithmetic $f -h1 $alpha x float32 -o $photCorrDir/$base
    astarithmetic $f -h1 $alpha x float32 -o $photCorrDir/$base
}
export -f applyCalibrationFactorsToFrame

applyCalibrationFactors() {
    imagesForCalibration=$1
    alphatruedir=$2
    photCorrDir=$3

    muldone=$photCorrDir/done_ccd"$h".txt
    if ! [ -d $photCorrDir ]; then mkdir $photCorrDir; fi
    if [ -f $muldone ]; then
            echo -e "\nMultiplication for alpha in the pointings (huge grid) is done for extension $h\n"
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
    h=0

    base=entirecamera_"$a".fits
    basetmp=entirecamera_"$a"_tmp.fits

    f=$photCorrDir/$base
    rms_min=$(awk 'NR=='1'{print $1}' $BDIR/rms_min_val-1_ccd"$h"_it$iteration.txt)
    rms_f=$(awk 'NR=='1'{print $3}' $noiseskydir/entirecamera_$a.txt)

    weight=$(astarithmetic $rms_min $rms_f / --quiet)
    echo "$weight" > $wdir/"$objectName"_Decals-"$filter"_"$a"_ccd"$h".txt        #    saving into file

    # multiply each image for its weight
    wixi_im_tmp=$wdir/$basetmp                     # frame x weight
    w_im_tmp=$wonlydir/$basetmp                     # only weight
    wixi_im=$wdir/$base                     # frame x weight
    w_im=$wonlydir/$base                     # only weight

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

    if [ -f $wdone ]; then
        echo -e "\nWeights computation done for extension $h\n"
    else
        framesToComputeWeight=()
        for a in $(seq 1 $totalNumberOfFrames); do
            framesToComputeWeight+=("$a")
        done
        printf "%s\n" "${framesToComputeWeight[@]}" | parallel -j "$num_cpus" computeWeightForFrame {} $wdir $wonlydir $photCorrDir $noiseskydir $iteration
        echo done > $wdone
        echo done > $wonlydone
    fi
}

# Outliers functions
buildUpperAndLowerLimitsForOutliers() {
    clippingdir=$1
    clippingdone=$2
    wdir=$3
    sigmaForStdSigclip=$4


    if ! [ -d $clippingdir ]; then mkdir $clippingdir; fi
    if [ -f $clippingdone ]; then
            echo -e "\nUpper and lower limits for building the masked of the weighted images already computed\n"
    else
            # Compute clipped median and std
            med_im=$clippingdir/median_image.fits
            std_im=$clippingdir/std_image.fits

            astarithmetic $(ls -v $wdir/*.fits) $(ls $wdir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-median -g1 -o$med_im
            astarithmetic $(ls -v $wdir/*.fits) $(ls $wdir/*.fits | wc -l) $sigmaForStdSigclip 0.2 sigclip-std -g1 -o$std_im
            # Compute "borders" images
            up_lim=$clippingdir/upperlim.fits
            lo_lim=$clippingdir/lowerlim.fits
            astarithmetic 4. $std_im x -o thresh.fits
            astarithmetic $med_im thresh.fits + -g1 float32 -o $up_lim
            astarithmetic $med_im thresh.fits - -g1 float32 -o $lo_lim

            #rm -f $med_im $std_im
            echo done > $clippingdone
    fi
}

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

    # # Remove temporary files
    rm -f $tmp_ab
    rm -f $mask
}
export -f removeOutliersFromFrame


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

    # # Offset for having some margin and not lose any pixel (the image will be )
    # securityOffset=200
    # detectorWidthDeg=$(echo    "(($detectorWidth + $securityOffset) * $pixelScale)" | bc )
    # detectorHeightDeg=$(echo "(($detectorHeight + $securityOffset) * $pixelScale) + $securityOffset" | bc )
    # astcrop $wholeMask --center=$centralRa,$centralDec --mode=wcs --width=$detectorHeightDeg/3600,$detectorWidthDeg/3600 -o $tmpMaskFile

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
    myCatalogue=$1
    referenceCatalogue=$2
    pythonScriptsPath=$3
    outputDir=$4
    pixelScale=$5

    astrometryTmpDir="./astrometryDiagnosisTmp"
    if ! [ -d $astrometryTmpDir ]; then mkdir $astrometryTmpDir; fi

    for i in $myCatalogue/match*.txt; do
        myFrame=$i
        frameNumber=$(echo "$i" | awk -F '[/]' '{print $(NF)}' | awk -F '[.]' '{print $(1)}' | awk -F '[_]' '{print $(NF)}')
        referenceFrame=$referenceCatalogue/*_$frameNumber.*
        astmatch $referenceFrame --hdu=1 $myFrame --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=1/3600 --outcols=aRA,aDEC,bRA,bDEC -o./$astrometryTmpDir/$frameNumber.cat
    done

    python3 $pythonScriptsPath/diagnosis_deltaRAdeltaDEC.py $astrometryTmpDir $outputDir/astrometry.png $pixelScale
    rm -r $astrometryTmpDir
}
export -f produceAstrometryCheckPlot

produceCalibrationCheckPlot() {
    myCatalogue_nonCalibrated=$1
    myFrames_calibrated=$2
    aperturesForMyData_dir=$3
    referenceCatalogueDir=$4
    pythonScriptsPath=$5
    outputDir=$6

    calibrationTmpDir="./calibrationDiagnosisTmp"
    if ! [ -d $calibrationTmpDir ]; then mkdir $calibrationTmpDir; fi

    for i in $myCatalogue_nonCalibrated/*.cat; do
        myFrame=$i
        frameNumber=$(echo "$i" | awk -F '[/]' '{print $(NF)}' | awk -F '[.]' '{print $(1)}' | awk -F '[_]' '{print $(NF)}')
        referenceCatalogue=$referenceCatalogueDir/*_$frameNumber.*

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

        astmatch $referenceCatalogue --hdu=1 $tmpDir/$frameNumber.cat --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=1/3600 --outcols=aMAGNITUDE,bMAGNITUDE -o$calibrationTmpDir/"$frameNumber"_matched.cat
        rm $tmpDir/apertures.txt
        rm $tmpDir/aperture_myData.fits
        rm $tmpDir/$frameNumber.cat
done

python3 $pythonScriptsPath/diagnosis_magVsDeltaMag.py $calibrationTmpDir $outputDir/magVsMagDiff.png
rm -rf $calibrationTmpDir

}
export -f produceCalibrationCheckPlot

# Coadd function
buildCoadd() {
    baseCoaddir=$1
    mowdir=$2
    moonwdir=$3

    if ! [ -d $baseCoaddir ]; then mkdir $baseCoaddir; fi

    coaddir=$baseCoaddir/stdWeighted
    coaddName=$coaddir/"$objectName"_coadd1_"$filter".fits
    coaddone=$coaddir/done_"$k".txt
    if ! [ -d $coaddir ]; then mkdir $coaddir; fi
    if [ -f $coaddone ]; then
            echo -e "\nThe first weighted (based upon std) mean of the images already done\n"
    else
            astarithmetic $(ls -v $mowdir/*.fits) $(ls $mowdir/*.fits | wc -l) sum -g1 -o$coaddir/"$k"_wx.fits
            astarithmetic $(ls -v $moonwdir/*.fits ) $(ls $moonwdir/*.fits    | wc -l) sum -g1 -o$coaddir/"$k"_w.fits
            astarithmetic $coaddir/"$k"_wx.fits -h1 $coaddir/"$k"_w.fits -h1 / -o$coaddName

            # Useful keywords for the final image
            numberOfFramesUsed=$(ls $mowdir/*.fits | wc -l)
            writeKeywordToFits $coaddName 1 "FRAMESCOMBINED" $numberOfFramesUsed
            echo done > $coaddone
    fi
}
