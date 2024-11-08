#!/bin/bash

# The functions needed for the pipeline are declared in another file (currently called "pipeline_LuM_parallel_functions.sh")
# This file is expected to be in the same directory as the pipeline
# The scripts that are used by the pipeline are expected to be in the directory "pipelineScripts" (value stored in pythonScriptsPath variable)

# How the pipeline expects the data -----------------------------------------------

# Two folders must exist: DATA-or and dark
# DATA-or:
#   Must contain the data. Each night placed in a folder called nightN where N is the number of the night
# dark:
#   Must contain the darks. Same format as data (i.e. nightN with N the night number)

# The pipeline also will look for a 'config' folder
# This folder must contain:
#   · Scamp conf file
#   · Swarp conf file
#   · Sextractor conf files (.conv, .param and .sex)

# The path of the ring(s) definition (.txt file) have to be provided in the configuration file of the galaxy.
# A common normalisation ring (most of the cases will be centered in the image) has to be provided (mandatory)
# Because it will be used also for selecting what decals bricks are going to be donwloaded for the photometric calibration
# The 2 rings needed for normalising with them are only requested if the normalisation is going to be done in that way (non mandatory)

# ----------------------------------------------------------------------------------------------

ORANGE='\033[0;33m'
GREEN='\033[0;32m'
RED='\033[0;31m'
NOCOLOUR='\033[0m'

export ORANGE
export GREEN
export NOCOLOUR



echo -e "\n ${GREEN} ---Loading pipeline Functions--- ${NOCOLOUR}"

# The path from which the pipeline is called is unknown. The functions for the pipeline are expected
# to be in the same folder as the pipeline, so we retrieve the path in order to run the functions file
pipelinePath=`dirname "$0"`
pipelinePath=`( cd "$pipelinePath" && pwd )`
pythonScriptsPath=$pipelinePath/pipelineScripts
export pipelinePath
export pythonScriptsPath

# Load the file with all the needed functions
source "$pipelinePath/pipeline_LuM_parallel_functions.sh"


echo -e "\n ${GREEN} ---Loading Modules--- ${NOCOLOUR}"

gnuastroModuleName="gnuastro/0.22"
load_module $gnuastroModuleName

astrometryModuleName="astrometry.net/0.94"
load_module $astrometryModuleName

# Needed if using the SIE software
pythonModuleName="python"
load_module $pythonModuleName


########## Handling options and arguments ##########
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

confFile=$1
if [[ -f $confFile ]]; then 
  source $confFile
  echo -e "\nVariables loaded from $confFile file\n"
else
  errorNumber=1
  echo -e "\nA configuration file has to be provided in order to run the pipeline"  >&2
  echo -e "Exiting with error number: $RED $errorNumber $NOCOLOUR" >&2
  exit $errorNumber
fi


startTime=$(date +%s)
echo ${GREEN} "\nStarting pipeline. Start time:  $(date +%D-%T)" ${NOCOLOUR}


# Exporting the variables from .conf file
echo -e "\n ${GREEN} ---Loading variables from conf file --- ${NOCOLOUR}"

export objectName
export ra_gal
export dec_gal

export ROOTDIR

export saturationThreshold
echo -e "\nSaturation threshold set to  " $saturationThreshold

export sizeOfOurFieldDegrees
echo -e "\nThe size of the field in degrees is " $sizeOfOurFieldDegrees

export coaddSizePx
echo -e "\nThe size in px of each side of the coadded image is " $coaddSizePx

export filter
export detectorWidth
export detectorHeight
export pixelScale


export USE_COMMON_RING
export commonRingDefinitionFile

export keyWordToDecideRing
export keyWordThreshold
export firstRingDefinitionFile
export keyWordValueForFirstRing
export secondRingDefinitionFile
export keyWordValueForSecondRing

export calibrationBrightLimit
export calibrationFaintLimit
echo -e "\nThe calibration range is from: " $ORANGE $calibrationBrightLimit $NOCOLOUR " to " $ORANGE $calibrationFaintLimit $NOCOLOUR 

echo -e "\nA common normalisation ring is going to be used?: " $ORANGE $USE_COMMON_RING $NOCOLOUR
if [ "$USE_COMMON_RING" = true ]; then
  echo -e "The file which contains the ring definition is: " $ORANGE $commonRingDefinitionFile $NOCOLOUR
else
  echo -e "The keyword for deciding the ring to use is: " $ORANGE $keyWordToDecideRing $NOCOLOUR 
  echo -e "The file containing the first ring definition is: " $ORANGE $firstRingDefinitionFile $NOCOLOUR " and the value for using this ring is " $keyWordValueForFirstRing
  echo -e "The file containing the second ring definition is: " $ORANGE  $secondRingDefinitionFile $NOCOLOUR  " and the value for using this ring is " $keyWordValueForSecondRing
  echo -e "And the treshold for detecting the value is: " $ORANGE $keyWordThreshold $NOCOLOUR
fi

echo -e "\nThe running flat is going to be used?: $ORANGE $RUNNING_FLAT $NOCOLOUR"
export RUNNING_FLAT
if [ "$RUNNING_FLAT" = true ]; then
  echo -e "The running flat will be computed with a window size of " $windowSize
  echo -e "\n"
fi
export windowSize
export halfWindowSize


export MODEL_SKY_AS_CONSTANT
export sky_estimation_method
export polynomialDegree
echo -e "\nThe background will be modelled as a constant?: $ORANGE $MODEL_SKY_AS_CONSTANT $NOCOLOUR"
echo -e "If so, the method to model the sky is: $sky_estimation_method"
echo -e "Otherwise, the polynomial degree is: $ORANGE $polynomialDegree $NOCOLOUR"


echo -e "\nThe indices that will be built for the construction of indices for astrometrisation are:"
echo -e "\tLowest index: $lowestScaleForIndex"
echo -e "\tHighest index: $highestScaleForIndex"
export lowestScaleForIndex
export highestScaleForIndex

export solve_field_L_Param
export solve_field_H_Param
export solve_field_u_Param

export numberOfStdForBadFrames

export defaultNumOfCPUs

checkIfAllVariablesAreSet
#

# The following lines are responsible of the cpu's used for paralellise
# If it is running in a system with slurm it takes the number of cpu's from the slurm job
# Otherwise it takes it from the provided configuration file
num_cpus=$SLURM_CPUS_ON_NODE
if [ -z $num_cpus ]; then
  num_cpus=$defaultNumOfCPUs
fi

echo "Number of CPUs allocated: $num_cpus"
export num_cpus

echo -e "\n-Variables defined"
echo -e "\t·Object name: ${ORANGE} ${objectName} ${NOCOLOUR}"
echo -e "\t·Filter: ${ORANGE} ${filter} ${NOCOLOUR}"
echo -e "\t·Ra: ${ORANGE} ${ra_gal} ${NOCOLOUR}"
echo -e "\t·Dec: ${ORANGE} ${dec_gal} ${NOCOLOUR}"
echo -e "\t·Detector width: ${ORANGE} ${detectorWidth} ${NOCOLOUR}"
echo -e "\t·Detector width: ${ORANGE} ${detectorHeight} ${NOCOLOUR}"

# ****** Decision note *******
# These parameters have been selected in order to obtain an aggresive mask
# Maybe they are not the best, I am not really used to use noisechisel so be careful


# I have tried with different set of parameters but the default ones work just fine with amateur data an rebinned TST data
# I think it's because with these big pixels it's easier to detect signal
# I just decreasing the erode and increasing the detgrowmaxholesize to be more conservative
noisechisel_param="--tilesize=35,35 \
                    --erode=1 \
                    --detgrowmaxholesize=5000 \
                    --rawoutput"

# These paremeters are oriented to TST data at original resolution. 
# In the nominal TST resolution the default parameters work really bad.
# I have modified them to detect fainter signal, since the pixel size is smaller I think it's harder an requires fine-tuning to detect more
# noisechisel_param="--tilesize=200,200
#                     --meanmedqdiff=0.01 \
#                     --detgrowquant=0.7 \
#                     --qthresh=0.25 \
#                     --snquant=0.98 \
#                     --erode=1 \
#                     --detgrowmaxholesize=5000
#                     --rawoutput"

export noisechisel_param

echo -e "\n-Noisechisel parameters used for masking:"
echo -e "\t" $noisechisel_param


########## Prepare data ##########

echo -e "\n ${GREEN} ---Preparing data--- ${NOCOLOUR}"

DIR=$ROOTDIR/"$objectName"
INDIRo=$ROOTDIR/"$objectName"/DATA-or
BDIR=$ROOTDIR/"$objectName"/build
INDIR=$ROOTDIR/"$objectName"/DATA
DARKDIR=$ROOTDIR/"$objectName"/dark
keyWordDirectory=$ROOTDIR/"$objectName"/keywords

export DIR
export INDIRo
export BDIR
export INDIR
export DARKDIR
export keyWordDirectory

if ! [ -d $BDIR ]; then mkdir $BDIR; fi
if ! [ -d $INDIR ]; then mkdir $INDIR; fi
if ! [ -d $filtereyWordDirectory ]; then mkdir $filtereyWordDirectory; fi

echo -e "\n-Directories defined"
echo -e "\t·Main directory (DIR): ${ORANGE} ${DIR} ${NOCOLOUR}"
echo -e "\t·Build directory (BDIR): ${ORANGE} ${BDIR} ${NOCOLOUR}"
echo -e "\t·Original data directory (INDIRo): ${ORANGE} ${INDIRo} ${NOCOLOUR}"
echo -e "\t·Data directory (INDIR): ${ORANGE} ${INDIR} ${NOCOLOUR}"
echo -e "\t·Dark Data directory (DARKDIR): ${ORANGE} ${DARKDIR} ${NOCOLOUR}"
echo -e "\t·KeyWords directory (keyWordDirectory): ${ORANGE} ${keyWordDirectory} ${NOCOLOUR}"


# Folders where the data and the results were being stored
SDIR=$ROOTDIR/"$objectName"
CDIR=$SDIR/config
export SDIR
export CDIR

echo -e "\n-Directories for results defined"
echo -e "\t·SDIR directory ${ORANGE} ${SDIR} ${NOCOLOUR}"
echo -e "\t·Config directory ${ORANGE} ${CDIR} ${NOCOLOUR}"
echo -e "INDIR and SDIR have no difference in this pipeline, we have to directories due to LBT procedence"

if ! [ -d $SDIR ]; then mkdir $SDIR; fi
if ! [ -d $CDIR ]; then mkdir $CDIR; fi

# Getting the coordinates of the galaxy
ra=$ra_gal
dec=$dec_gal
export ra
export dec

echo -e "\nCoordinates of the galaxy:"
echo -e "\t·RA: ${ORANGE} ${ra} ${NOCOLOUR}"
echo -e "\t·DEC: ${ORANGE} ${dec} ${NOCOLOUR}"

numberOfNights=$(ls -d $INDIRo/night* | wc -l)
export numberOfNights

echo -e "\nNumber of nights to reduce: ${ORANGE} $numberOfNights ${NOCOLOUR}"
echo -e "\n"

export airMassKeyWord
export dateHeaderKey

framesForCommonReductionDir=$BDIR/framesForCommonReduction
export framesForCommonReductionDir

# Function which processes a whole night
oneNightPreProcessing() {
  currentNight=$1
  framesForCommonReductionDone=$framesForCommonReductionDir/done_"$filter"_ccd"$h"_n"$currentNight".txt

  echo -e "\n\n"
  echo -e "${ORANGE} --- STARTING TO PROCESS NIGHT NUMBER $currentNight --- ${NOCOLOUR}"

  h=0

  if ! [ -d $framesForCommonReductionDir ]; then mkdir $framesForCommonReductionDir; fi
  if [ -f $framesForCommonReductionDone ]; then
    echo -e "\nScience images for night $currentNight are already processed\n"
    return 0
  fi


  # ****** Decision note *******
  # In the following, the data from "INDIRo/nightX" is placed in "INDIR/nightX". Additionally they are
  # sorted and renamed based on an objetive criteria as is the time in which the frame was taken
  # otherwise if the files are selected just by default (using glob for example) the order is OS-dependant.

  # The name of the files contains the patter nX and fY, standing for "night number X" and "frame number Y"
  currentINDIR=$INDIR/night"$currentNight"
  currentINDIRo=$INDIRo/night"$currentNight"
  renamedone=$currentINDIR/done_.txt
  if ! [ -d $currentINDIR ]; then mkdir $currentINDIR; fi
  if [ -f $renamedone ]; then
    echo -e "\nScience images for night $currentNight are already renamed\n"
  else
      for h in 0; do
          for i in $currentINDIRo/*.fits; do
            nameWithEscapedSpaces=$(escapeSpacesFromString "$i")
            DATEOBS=$(eval "astfits $nameWithEscapedSpaces -h0 --keyvalue=$dateHeaderKey --quiet")
            checkIfExist_DATEOBS $DATEOBS

            unixTimeInSeconds=$(date -d "$DATEOBS" +"%s")
            out=$currentINDIR/$unixTimeInSeconds.fits

            # HERE A CHECK IF THE DATA IS IN FLOAT32 IS NEEDED
            eval "astfits $nameWithEscapedSpaces --copy=$h -o$out"  # I run this with eval so the escaped spaces are re-parsed by bash and understood by astfits
            nameOfOriginalFile="${nameWithEscapedSpaces##*/}"
            eval "astfits --write=OriginalName,$nameOfOriginalFile $out -h0"
          done

          index=1
          for i in $(ls -v $currentINDIR/*.fits); do
            mv $i $currentINDIR/"$objectName"-Decals-"$filter"_n"$currentNight"_f"$index"_ccd"$h".fits
            index=$((index+1));
          done
      done
      echo done > $renamedone
  fi

  # -------------------------------------------------------
  # Number of exposures of the current night
  n_exp=$(ls -v $currentINDIRo/*.fits | wc -l)
  echo -e "Number of exposures ${ORANGE} ${n_exp} ${NOCOLOUR}"

  currentDARKDIR=$DARKDIR/night$currentNight
  mdadir=$BDIR/masterdark_n$currentNight

  # Loop for all the ccds
  for h in 0; do
    ########## Creating master bias ##########
    echo -e "\n ${GREEN} Creating master bias/dark-bias ${NOCOLOUR}"
    mdadone=$mdadir/mdark_"$filter"_ccd"$h".txt
    if ! [ -d $mdadir ]; then mkdir $mdadir; fi
    if [ -f $mdadone ]; then
      echo -e "\nMasterdark is already done for night $currentNight and extension $h\n"
    else
      escaped_files=""
      for file in $currentDARKDIR/*.fits; do
        escaped_files+="$(escapeSpacesFromString "$file") "
      done

      eval "astarithmetic $escaped_files $(ls -v $currentDARKDIR/* | wc -l) \
                    3 0.2 sigclip-mean -g$h \
                    -o $mdadir/mdark_"$filter"_n"$currentNight"_ccd$h.fits"
    fi
    echo done > $mdadone
  done


  ########## Save airmass ##########
  # The airmass is saved in this files on airmass-analysis_n folder but also propagated throught the steps of the pipeline until that information
  # reaches the "framesForCommonReduction", because that information needs to be used in the future for the detection of bad frames
  echo -e "\n ${GREEN} Saving airmass ${NOCOLOUR}"

  skydir=$BDIR/airmass-analysis_n$currentNight
  skydone=$skydir/done_.txt
  if ! [ -d $skydir ]; then mkdir $skydir; fi
  if [ -f $skydone ]; then
    echo -e "\nAirmass for night $currentNight already saved\n"
  else
      for i in $(ls -v $currentINDIR/*.fits ); do
        air=$(astfits $i -h1 --keyvalue=$airMassKeyWord | awk '{print $2}')
        echo $air >> $skydir/airmass.txt
    done
    echo done > $skydone
  fi

  
  ########## Subtract master bias and dark ##########
  echo -e "\n ${GREEN} Subtracting master bias/dark-bias ${NOCOLOUR}"

  # Now using indexes that could be better
  # Substracting mbias to science images.
  # and also from images for flat.
  # Also a counter variable is
  # created to rename the images. Bad and saturated pixels are masked.
  mbiascorrdir=$BDIR/bias-corrected_n$currentNight
  mbiascorrdone=$mbiascorrdir/done_"$filter"_ccd"$h".txt
  if ! [ -d $mbiascorrdir ]; then mkdir $mbiascorrdir; fi
  if [ -f $mbiascorrdone ]; then
    echo -e "\nScience images are already bias/dark corrected for night $currentNight and extension $h\n"
  else
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      i=$currentINDIR/$base
      out=$mbiascorrdir/$base
      astarithmetic $i -h1 set-i $mdadir/mdark_"$filter"_n"$currentNight"_ccd$h.fits  -h1  set-m \
                i i $saturationThreshold gt i isblank or 2 dilate nan where m -  float32  \
                -o $out

      propagateKeyword $i $airMassKeyWord $out
      # If we are not doing a normalisation with a common ring we propagate the keyword that will be used to decide
      # which ring is to be used. This way we can check this value in a comfortable way in the normalisation section
      # This is also done in the function maskImages()
      if [ "$USE_COMMON_RING" = false ]; then
        propagateKeyword $i $keyWordToDecideRing $out
      fi
    done
    echo done > $mbiascorrdone
  fi


  # until here everything is corrected of bias and renamed
  echo -e "${ORANGE} ------ FLATS ------ ${NOCOLOUR}\n"

  #now flat
  echo -e "${GREEN} --- Flat iteration 1 --- ${NOCOLOUR}"


  ########## Creating the ring mask ##########

  # Since the ring possibilities are quite flexible (using one common ring or two rings depending on the angle) I always clean and rebuild the rings (It doesn't cost almost time so is worth it)
  ringdir=$BDIR/ring
  rm -rf $ringdir
  mkdir $ringdir

  # We always need the common ring  definition always stored for photometric calibration (selection of decals bricks to download)
  cp $commonRingDefinitionFile $ringdir/ring.txt 
  # We create the .fits ring image based on how the normalisation is going to be done
  if [ "$USE_COMMON_RING" = true ]; then
    astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=200 --clearcanvas -o $ringdir/ring.fits $commonRingDefinitionFile
  else
    astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=200 --clearcanvas -o $ringdir/ring_2.fits $secondRingDefinitionFile
    astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=200 --clearcanvas -o $ringdir/ring_1.fits $firstRingDefinitionFile
  fi

  ########## Creating the it1 master flat image ##########

  # ****** Decision note *******
  # Running flat: summary
  # A running flat is a flat built with not all the frames of one night, but each frame has a local flat built from N frames based on time-near frames.
  # Caveat: This requieres to check the data of the night, this only works if the data has been taken at similar times (check DATE-OBS or the airmass)
  # The size of the window for the running flat can vambiascorrdir/$base 
  # If the running flat is activated, the whole night flat will be used to correct the running flat
  # If the running flat is not activated, the whole night flat will be used to be applied to the data

  # · Current situation with the pipeline
  # For the running flat being effective we need a great dithering pattern. The data right now has a not appropriate dithering for the running flat
  # So the whole night flat approach will be used. But In order to generalise the pipeline the option of using the running flat or not is
  # configure by the parameter "RUNNING_FLAT"


  # Creating iteration 1 flat_it1. First we need to normalise the science images.
  normit1dir=$BDIR/norm-it1-images_n$currentNight
  normit1done=$normit1dir/done_"$filter"_ccd"$h".txt
  if ! [ -d $normit1dir ]; then mkdir $normit1dir; fi
  if [ -f $normit1done ]; then
    echo -e "\nScience images are already normalized for night $currentNight and extension $h\n"
  else
    normaliseImagesWithRing $mbiascorrdir $normit1dir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
    echo done > $normit1done
  fi
 
  # Then, if the running flat is configured to be used, we combine the normalised images with a sigma clipping median
  # using the running flat strategy
  if $RUNNING_FLAT; then
    flatit1dir=$BDIR/flat-it1-Running_n$currentNight
    flatit1done=$flatit1dir/done_"$filter"_ccd"$h".txt
    iteration=1
    if ! [ -d $flatit1dir ]; then mkdir $flatit1dir; fi
    if [ -f $flatit1done ]; then
      echo -e "\nRunning flats it-1 already built for night $currentNight and extension $h\n"
    else
      calculateRunningFlat $normit1dir $flatit1dir $flatit1done $iteration
    fi
  fi

  # We compute the flat using all the frames of the night
  flatit1WholeNightdir=$BDIR/flat-it1-WholeNight_n$currentNight
  flatit1WholeNightdone=$flatit1WholeNightdir/done_"$filter"_ccd"$h".txt
  iteration=1
  if ! [ -d $flatit1WholeNightdir ]; then mkdir $flatit1WholeNightdir; fi
  if [ -f $flatit1WholeNightdone ]; then
    echo -e "\nWhole night flat it-1 already built for night $currentNight and extension $h\n"
  else
    calculateFlat $flatit1WholeNightdir/flat-it1_wholeNight_n$currentNight.fits $normit1dir/*.fits
    echo "done" >> $flatit1WholeNightdone
  fi


  # Dividing the science images for the running it1 flat
  if $RUNNING_FLAT; then
    flatit1imadir=$BDIR/flat-it1-Running-ima_n$currentNight
    flatit1imadone=$flatit1imadir/done_"$filter"_ccd"$h".txt
    if ! [ -d $flatit1imadir ]; then mkdir $flatit1imadir; fi
    if [ -f $flatit1imadone ]; then
      echo -e "\nScience images are divided by flat it1 for night $currentNight and extension $h\n"
    else
      divideImagesByRunningFlats $mbiascorrdir $flatit1imadir $flatit1dir $flatit1imadone
    fi
  fi

  # Dividing the science images for the whole night it1 flat
  flatit1WholeNightimaDir=$BDIR/flat-it1-WholeNight-ima_n$currentNight
  flatit1WholeNightimaDone=$flatit1WholeNightimaDir/done_"$filter"_ccd"$h".txt
  if ! [ -d $flatit1WholeNightimaDir ]; then mkdir $flatit1WholeNightimaDir; fi
  if [ -f $flatit1WholeNightimaDone ]; then
    echo -e "\nScience images are divided by whole night flat it1 for night $currentNight and extension $h\n"
  else
    wholeNightFlatToUse=$flatit1WholeNightdir/flat-it1_wholeNight_n$currentNight.fits
    divideImagesByWholeNightFlat $mbiascorrdir $flatit1WholeNightimaDir $wholeNightFlatToUse $flatit1WholeNightimaDone
  fi


  
  ########## Creating the it2 master flat image ##########
  echo -e "${GREEN} --- Flat iteration 2 --- ${NOCOLOUR}"


  # Obtain a mask using noisechisel on the running flat images
  if $RUNNING_FLAT; then
    noiseit2dir=$BDIR/noise-it2-Running_n$currentNight
    noiseit2done=$noiseit2dir/done_"$filter"_ccd"$h".txt
    if ! [ -d $noiseit2dir ]; then mkdir $noiseit2dir; fi
    if [ -f $noiseit2done ]; then
      echo -e "\nScience images are 'noisechiseled' for it2 running flat for night $currentNight and extension $h\n"
    else
      frameNames=()
      for a in $(seq 1 $n_exp); do
          base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
          frameNames+=("$base")
      done
      printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" runNoiseChiselOnFrame {} $flatit1imadir $noiseit2dir "'$noisechisel_param'"
      echo done > $noiseit2done
    fi
  fi


  # Obtain a mask using noisechisel on the whole night flat images
  noiseit2WholeNightDir=$BDIR/noise-it2-WholeNight_n$currentNight
  noiseit2WholeNightdone=$noiseit2WholeNightDir/done_"$filter"_ccd"$h".txt
  if ! [ -d $noiseit2WholeNightDir ]; then mkdir $noiseit2WholeNightDir; fi
  if [ -f $noiseit2WholeNightdone ]; then
    echo -e "\nScience images are 'noisechiseled' for it2 whole night flat for night $currentNight and extension $h\n"
  else
    frameNames=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      frameNames+=("$base")
    done
    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" runNoiseChiselOnFrame {} $flatit1WholeNightimaDir $noiseit2WholeNightDir "'$noisechisel_param'"
    echo done > $noiseit2WholeNightdone
  fi

  # Mask the images (running flat)
  if $RUNNING_FLAT; then
    maskedit2dir=$BDIR/masked-it2-Running_n$currentNight
    maskedit2done=$maskedit2dir/done_"$filter"_ccd"$h".txt
    if ! [ -d $maskedit2dir ]; then mkdir $maskedit2dir; fi
    if [ -f $maskedit2done ]; then
      echo -e "\nScience images are masked for running flat, night $currentNight and extension $h\n"
    else
      maskImages $mbiascorrdir $noiseit2dir $maskedit2dir $USE_COMMON_RING $keyWordToDecideRing
      echo done > $maskedit2done
    fi
  fi

  # Mask the images (whole night flat)
  maskedit2WholeNightdir=$BDIR/masked-it2-WholeNight_n$currentNight
  maskedit2WholeNightdone=$maskedit2WholeNightdir/done_"$filter"_ccd"$h".txt
  if ! [ -d $maskedit2WholeNightdir ]; then mkdir $maskedit2WholeNightdir; fi
  if [ -f $maskedit2WholeNightdone ]; then
    echo -e "\nScience images are masked for whole night flat, night $currentNight and extension $h\n"
  else
    maskImages $mbiascorrdir $noiseit2WholeNightDir $maskedit2WholeNightdir $USE_COMMON_RING $keyWordToDecideRing
    echo done > $maskedit2WholeNightdone
  fi

  
  # Normalising masked images (running flat)
  if $RUNNING_FLAT; then
    normit2dir=$BDIR/norm-it2-Running-images_n$currentNight
    normit2done=$normit2dir/done_"$filter"_ccd"$h".txt
    if ! [ -d $normit2dir ]; then mkdir $normit2dir; fi
    if [ -f $normit2done ]; then
      echo -e "\nMasked science images are normalized for running flat, night $currentNight and extension $h\n"
    else
      normaliseImagesWithRing $maskedit2dir $normit2dir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
      echo done > $normit2done
    fi
  fi

  # Normalising masked images (whole night flat)
  normit2WholeNightdir=$BDIR/norm-it2-WholeNight-images_n$currentNight
  normit2WholeNightdone=$normit2WholeNightdir/done_"$filter"_ccd"$h".txt
  if ! [ -d $normit2WholeNightdir ]; then mkdir $normit2WholeNightdir; fi
  if [ -f $normit2WholeNightdone ]; then
    echo -e "\nMasked science images are normalized for whole night flat, night $currentNight and extension $h\n"
  else
    normaliseImagesWithRing $maskedit2WholeNightdir $normit2WholeNightdir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
    echo done > $normit2WholeNightdone
  fi

  # Combining masked normalized images to make it2 running flat
  if $RUNNING_FLAT; then
    flatit2dir=$BDIR/flat-it2-Running_n$currentNight
    flatit2done=$flatit2dir/done_"$filter"_ccd"$h".txt
    iteration=2
    if ! [ -d $flatit2dir ]; then mkdir $flatit2dir; fi
    if [ -f $flatit2done ]; then
      echo -e "\nScience images are stacked for it2 running flat for night $currentNight and extension $h\n"
    else
      calculateRunningFlat $normit2dir $flatit2dir $flatit2done $iteration
    fi
  fi

  # We also compute the flat using all the frames of the night.
  flatit2WholeNightdir=$BDIR/flat-it2-WholeNight_n$currentNight
  flatit2WholeNightdone=$flatit2WholeNightdir/done_"$filter"_ccd"$h".txt
  iteration=2
  if ! [ -d $flatit2WholeNightdir ]; then mkdir $flatit2WholeNightdir; fi
  if [ -f $flatit2WholeNightdone ]; then
    echo -e "\nWhole night flat it-2 already built for night $currentNight and extension $h\n"
  else
    calculateFlat $flatit2WholeNightdir/flat-it2_wholeNight_n$currentNight.fits $normit2WholeNightdir/*.fits
    echo "done" >> $flatit2WholeNightdone
  fi

  
  # Dividing the science image by the it2 flat
  if $RUNNING_FLAT; then
    flatit2imadir=$BDIR/flat-it2-Running-ima_n$currentNight
    flatit2imadone=$flatit2imadir/done_"$filter"_ccd"$h".txt
    if ! [ -d $flatit2imadir ]; then mkdir $flatit2imadir; fi
    if [ -f $flatit2imadone ]; then
      echo -e "\nRunning flats it2-2 already built for night $currentNight and extension $h\n"
    else
      divideImagesByRunningFlats $mbiascorrdir $flatit2imadir $flatit2dir $flatit2imadone
    fi
  fi

  # Dividing the science images for the whole night it2 flat
  flatit2WholeNightimaDir=$BDIR/flat-it2-WholeNight-ima_n$currentNight
  flatit2WholeNightimaDone=$flatit2WholeNightimaDir/done_"$filter"_ccd"$h".txt
  if ! [ -d $flatit2WholeNightimaDir ]; then mkdir $flatit2WholeNightimaDir; fi
  if [ -f $flatit2WholeNightimaDone ]; then
    echo -e "\nScience images are divided by whole night flat it2 for night $currentNight and extension $h\n"
  else
    wholeNightFlatToUse=$flatit2WholeNightdir/flat-it2_wholeNight_n$currentNight.fits
    divideImagesByWholeNightFlat $mbiascorrdir $flatit2WholeNightimaDir $wholeNightFlatToUse $flatit2WholeNightimaDone
  fi




  ########## Creating the it3 master flat image ##########
  echo -e "${GREEN} --- Flat iteration 3 --- ${NOCOLOUR}"

  # Obtain a mask using noisechisel on the running flat images
  if $RUNNING_FLAT; then
    noiseit3dir=$BDIR/noise-it3-Running_n$currentNight
    noiseit3done=$noiseit3dir/done_"$filter"_ccd"$h".txt
    if ! [ -d $noiseit3dir ]; then mkdir $noiseit3dir; fi
    if [ -f $noiseit3done ]; then
      echo -e "\nScience images are 'noisechiseled' for it3 running flat for night $currentNight and extension $h\n"
    else
      frameNames=()
      for a in $(seq 1 $n_exp); do
          base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
          frameNames+=("$base")
      done
      printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" runNoiseChiselOnFrame {} $flatit2imadir $noiseit3dir "'$noisechisel_param'"
      echo done > $noiseit3done
    fi
  fi

  # Obtain a mask using noisechisel on the whole night flat images
  noiseit3WholeNightDir=$BDIR/noise-it3-WholeNight_n$currentNight
  noiseit3WholeNightdone=$noiseit3WholeNightDir/done_"$filter"_ccd"$h".txt
  if ! [ -d $noiseit3WholeNightDir ]; then mkdir $noiseit3WholeNightDir; fi
  if [ -f $noiseit3WholeNightdone ]; then
    echo -e "\nScience images are 'noisechiseled' for it3 whole night flat for night $currentNight and extension $h\n"
  else
    frameNames=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      frameNames+=("$base")

    done
    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" runNoiseChiselOnFrame {} $flatit2WholeNightimaDir $noiseit3WholeNightDir "'$noisechisel_param'"
    echo done > $noiseit3WholeNightdone
  fi

  
  # Mask the images (running flat)
  if $RUNNING_FLAT; then
    maskedit3dir=$BDIR/masked-it3-Running_n$currentNight
    maskedit3done=$maskedit3dir/done_"$filter"_ccd"$h".txt
    if ! [ -d $maskedit3dir ]; then mkdir $maskedit3dir; fi
    if [ -f $maskedit3done ]; then
      echo -e "\nScience images are masked for running flat, night $currentNight and extension $h\n"
    else
      maskImages $mbiascorrdir $noiseit3dir $maskedit3dir $USE_COMMON_RING $keyWordToDecideRing
      echo done > $maskedit3done
    fi
  fi

  # Mask the images (whole night flat)
  maskedit3WholeNightdir=$BDIR/masked-it3-WholeNight_n$currentNight
  maskedit3WholeNightdone=$maskedit3WholeNightdir/done_"$filter"_ccd"$h".txt
  if ! [ -d $maskedit3WholeNightdir ]; then mkdir $maskedit3WholeNightdir; fi
  if [ -f $maskedit3WholeNightdone ]; then
    echo -e "\nScience images are masked for whole night flat, night $currentNight and extension $h\n"
  else
    maskImages $mbiascorrdir $noiseit3WholeNightDir $maskedit3WholeNightdir $USE_COMMON_RING $keyWordToDecideRing
    echo done > $maskedit3WholeNightdone
  fi

  # Normalising masked images (running flat)
  if $RUNNING_FLAT; then
    normit3dir=$BDIR/norm-it3-Running-images_n$currentNight
    normit3done=$normit3dir/done_"$filter"_ccd"$h".txt
    if ! [ -d $normit3dir ]; then mkdir $normit3dir; fi
    if [ -f $normit3done ]; then
      echo -e "\nMasked science images are normalized for running flat, night $currentNight and extension $h\n"
    else
      normaliseImagesWithRing $maskedit3dir $normit3dir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
      echo done > $normit3done
    fi
  fi

  # Normalising masked images (whole night flat)
  normit3WholeNightdir=$BDIR/norm-it3-WholeNight-images_n$currentNight
  normit3WholeNightdone=$normit3WholeNightdir/done_"$filter"_ccd"$h".txt
  if ! [ -d $normit3WholeNightdir ]; then mkdir $normit3WholeNightdir; fi
  if [ -f $normit3WholeNightdone ]; then
    echo -e "\nMasked science images are normalized for whole night flat, night $currentNight and extension $h\n"
  else
    normaliseImagesWithRing $maskedit3WholeNightdir $normit3WholeNightdir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
    echo done > $normit3WholeNightdone
  fi


  # Combining masked normalized images to make it3 flat
  if $RUNNING_FLAT; then
    flatit3BeforeCorrectiondir=$BDIR/flat-it3-Running-BeforeCorrection_n$currentNight
    flatit3BeforeCorrectiondone=$flatit3BeforeCorrectiondir/done_"$k"_ccd"$h".txt
    iteration=3
    if ! [ -d $flatit3BeforeCorrectiondir ]; then mkdir $flatit3BeforeCorrectiondir; fi
    if [ -f $flatit3BeforeCorrectiondone ]; then
      echo -e "\nRunning flats it3 before correction are already build for night $currentNight and extension $h\n"
    else
      calculateRunningFlat $normit3dir $flatit3BeforeCorrectiondir $flatit3BeforeCorrectiondone $iteration
    fi
  fi

  # We also compute the flat using all the frames of the night.
  flatit3WholeNightdir=$BDIR/flat-it3-WholeNight_n$currentNight
  flatit3WholeNightdone=$flatit3WholeNightdir/done_"$filter"_ccd"$h".txt
  iteration=3
  if ! [ -d $flatit3WholeNightdir ]; then mkdir $flatit3WholeNightdir; fi
  if [ -f $flatit3WholeNightdone ]; then
    echo -e "\nWhole night flat it-3 already built for night $currentNight and extension $h\n"
  else
    calculateFlat $flatit3WholeNightdir/flat-it3_wholeNight_n$currentNight.fits $normit3WholeNightdir/*.fits
    echo "done" >> $flatit3WholeNightdone
  fi


  # Correct the running flats using the whole night flat
  if $RUNNING_FLAT; then
    flatit3dir=$BDIR/flat-it3-Running_n$currentNight
    flatit3done=$flatit3dir/done_"$k"_ccd"$h".txt
    if ! [ -d $flatit3dir ]; then mkdir $flatit3dir; fi
    if [ -f $flatit3done ]; then
      echo -e "\nFlats iteration 3 are corrected using the flat of the whole night for night $currentNight and extension $h\n"
    else
      for i in $flatit3BeforeCorrectiondir/*.fits; do

        tmpRatio=$flatit3dir/tmpRatio.fits
        astarithmetic $flatit3WholeNightdir/flat-it3_wholeNight_n$currentNight.fits -h1 $i -h1 / -o$tmpRatio

        # ****** Decision note *******
        # The standard deviation of the ratio between the whole flat and the running flat is around 0.03.
        # So choosing 0.85, which seems a reasonable value.
        # Chose this value based on the standard deviation of your ratios, how these vary through the night and how aggresive u want to apply the correction
        astarithmetic $i -h1 set-m  $tmpRatio -h1 set-f m f 0.85 lt nan where -o $flatit3dir/$(basename "$i")
        rm $tmpRatio
      done
      echo done > $flatit3done
    fi
  fi

  # Dividing the science image by the it3 flat
  # If running flat selected, we use it to produce the final flatted images
  # If not selcted, we applyt the whole night flat
  flatit3imadir=$BDIR/flat-it3-ima_n$currentNight
  flatit3imadone=$flatit3imadir/done_"$filter"_ccd"$h".txt
  if ! [ -d $flatit3imadir ]; then mkdir $flatit3imadir; fi
  if $RUNNING_FLAT; then
    if [ -f $flatit3imadone ]; then
      echo -e "\nScience images are divided by the it3 flat for night $currentNight and extension $h\n"
    else
      divideImagesByRunningFlats $mbiascorrdir $flatit3imadir $flatit3dir $flatit3imadone
    fi
  else
      wholeNightFlatToUse=$flatit3WholeNightdir/flat-it3_wholeNight_n$currentNight.fits
      divideImagesByWholeNightFlat $mbiascorrdir $flatit3imadir $wholeNightFlatToUse $flatit3imadone
  fi


  ########## Masking the vignetting zones ##########
  # Enmascarando las esquinas
  echo -e "${GREEN} --- Masking vignetting zones --- ${NOCOLOUR}"


  maskedcornerdir=$BDIR/masked-corner_n$currentNight
  maskedcornerdone=$maskedcornerdir/done_"$filter"_ccd"$h".txt
  if ! [ -d $maskedcornerdir ]; then mkdir $maskedcornerdir; fi
  if [ -f $maskedcornerdone ]; then
    echo -e "\nCorners are already masked for night $currentNight and extension $h\n"
  else
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits

      if $RUNNING_FLAT; then
        if [ "$a" -le "$((halfWindowSize + 1))" ]; then
          currentFlatImage=$flatit3dir/flat-it3_"$filter"_n"$currentNight"_left_ccd"$h".fits
        elif [ "$a" -ge "$((n_exp - halfWindowSize))" ]; then
          currentFlatImage=$flatit3dir/flat-it3_"$filter"_n"$currentNight"_right_ccd"$h".fits
        else
          currentFlatImage=$flatit3dir/flat-it3_"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
        fi
      else
        currentFlatImage=$flatit3WholeNightdir/flat-it3_wholeNight_n$currentNight.fits
      fi

      i=$flatit3imadir/$base
      out=$maskedcornerdir/$base
      astarithmetic $i -h1 set-m $currentFlatImage -h1 set-f m f 0.9 lt  nan where set-n n f 2. gt nan where -o $out
      propagateKeyword $i $airMassKeyWord $out 
    done
    echo done > $maskedcornerdone
  fi

  
  # At this point we can process the frames of all the nights in the same way
  # So we place all the final frames into a common folder.

  # WE SHOULD MANAGE THE RACE CONDITIONS
  if [ -f $framesForCommonReductionDone ]; then
    echo -e "\nFrames already placed in the folder for frames prepared to common reduction"
  else
    initialValue=$( getHighestNumberFromFilesInFolder $framesForCommonReductionDir )

    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      name=$(( $initialValue + $a ))
      cp $maskedcornerdir/$base $framesForCommonReductionDir/$name.fits
    done
    echo "done" > $framesForCommonReductionDone
  fi

  # Removing intermediate information to save space
  rm -rf $BDIR/bias-corrected_n$currentNight
  rm -rf $BDIR/masked-corner_n$currentNight
  rm -rf $BDIR/masterdark_n$currentNight
  rm -rf $BDIR/flat-it3-Running-BeforeCorrection_n$currentNight
  rm -rf $BDIR/flat-it3-ima_n$currentNight
  for a in $(seq 1 3); do
    # I comment the following two lines in order to keep the final flats used
    # rm -rf $BDIR/flat-it"$a"-Running_n$currentNight
    # rm -rf $BDIR/flat-it"$a"-WholeNight_n$currentNight
    rm -rf $BDIR/flat-it"$a"-Running-ima_n$currentNight
    rm -rf $BDIR/flat-it"$a"-WholeNight-ima_n$currentNight
    rm -rf $BDIR/masked-it"$a"-Running_n$currentNight
    rm -rf $BDIR/masked-it"$a"-WholeNight_n$currentNight
    rm -rf $BDIR/noise-it"$a"-Running_n$currentNight
    rm -rf $BDIR/noise-it"$a"-WholeNight_n$currentNight
    rm -rf $BDIR/norm-it"$a"-images_n$currentNight
    rm -rf $BDIR/norm-it"$a"-Running-images_n$currentNight
    rm -rf $BDIR/norm-it"$a"-WholeNight-images_n$currentNight
  done

}
export -f oneNightPreProcessing

nights=()
for currentNight in $(seq 1 $numberOfNights); do
      nights+=("$currentNight")
done
printf "%s\n" "${nights[@]}" | parallel --line-buffer -j "$num_cpus" oneNightPreProcessing {}

totalNumberOfFrames=$( ls $framesForCommonReductionDir/*.fits | wc -l)
export totalNumberOfFrames
echo $totalNumberOfFrames


# Up to this point the frame of every night has been corrected of bias-dark and flat.
# That corrections are perform night by night (because it's necessary for perform that corretions)
# Now, all the frames are "equal" so we do no distinction between nights.
# All the frames are stored together in $framesForCommonReductionDir with names 1.fits, 2.fits, 3.fits ... n.fits.

echo -e "${ORANGE} ------ ASTROMETRY AND BACKGROUND-SUBTRACTION ------ ${NOCOLOUR}\n"
for h in 0; do
  echo -e "${GREEN} --- Astrometry --- ${NOCOLOUR}"

  query_param="gaia --dataset=edr3 --center=$ra_gal,$dec_gal --radius=3.1 --column=ra,dec,phot_g_mean_mag"
  catdir=$BDIR/catalogs
  catdone=$catdir/"$objectName"_Gaia_eDR3.fits
  if ! [ -d $catdir ]; then mkdir $catdir; fi
  if [ -f $catdone ]; then
    echo -e "\nCatalog is already downloaded\n"
  else
    astquery $query_param -o $catdir/"$objectName"_Gaia_eDR3.fits
  fi


  # Making the indexes
  indexdir=$BDIR/indexes
  indexdone=$indexdir/done_"$filter".txt
  if ! [ -d $indexdir ]; then mkdir $indexdir; fi
  if [ -f $indexdone ]; then
    echo -e "\nGaia eDR3 indexes are already created\n"
  else
    # Here we build the indices for different index scales
    # The index defines the scale on which the stars are selected
    # It is recommended to build a range of scales
    indexes=()
    for re in $(seq $lowestScaleForIndex $highestScaleForIndex); do
        indexes+=("$re")
    done
    printf "%s\n" "${indexes[@]}" | parallel -j "$num_cpus" downloadIndex {} $catdir $objectName $indexdir
    echo done > $indexdone
  fi


  sexcfg=$CDIR/default.sex
  # Solving the images
  astrocfg=$CDIR/astrometry_$objectName.cfg
  
  rm $astrocfg
  echo inparallel > $astrocfg
  echo cpulimit 300 >> $astrocfg
  echo "add_path $indexdir" >> $astrocfg
  echo autoindex >> $astrocfg
  
  astroimadir=$BDIR/astro-ima
  astroimadone=$astroimadir/done_"$filter"_ccd"$h".txt
  if ! [ -d $astroimadir ]; then mkdir $astroimadir; fi
  if [ -f $astroimadone ]; then
    echo -e "\nImages are already astrometrized for extension $h\n"
  else
    frameNames=()
    for a in $(seq 1 $totalNumberOfFrames); do
        base=$a.fits
        i=$framesForCommonReductionDir/$base
        frameNames+=("$i")
    done
    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" solveField {} $solve_field_L_Param $solve_field_H_Param $solve_field_u_Param $ra_gal $dec_gal $astrocfg $astroimadir
    echo done > $astroimadone
  fi
done



# scamp swarp
CDIR=$ROOTDIR/"$objectName"/config

########## Distorsion correction ##########
echo -e "\n ${GREEN} ---Creating distorsion correction files--- ${NOCOLOUR}"


# Making sex catalogs
sexcfg=$CDIR/default.sex
sexparam=$CDIR/default.param
sexconv=$CDIR/default.conv
sexdir=$BDIR/sex-it1
sexdone=$sexdir/done_"$filter"_ccd"$h".txt
if ! [ -d $sexdir ]; then mkdir $sexdir; fi
if [ -f $sexdone ]; then
    echo -e "\nSex catalogs are already done for extension $h\n"
else
    frameNames=()
    for a in $(seq 1 $totalNumberOfFrames); do
        frameNames+=("$a")
    done
    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" runSextractorOnImage {} $sexcfg $sexparam $sexconv $astroimadir $sexdir
    echo done > $sexdone
fi



# Making scamp headers
scampcfg=$CDIR/scamp.cfg
scampdir=$BDIR/scamp-it1
scampres=$scampdir/results_Decals-"$filter"_ccd"$h"
scampdone=$scampdir/done_"$filter"_ccd"$h".txt
if ! [ -d $scampdir ]; then mkdir $scampdir; fi
if ! [ -d $scampres ]; then mkdir $scampres; fi
if [ -f $scampdone ]; then
    echo -e "\nScamp headers are already done for extension $h\n"
else
    scamp -c $scampcfg $(ls -v $sexdir/*.cat)
    mv *.pdf $scampres/
    echo done > $scampdone
fi

# We copy the files for improving the astrometry into the folder of the images that we are going to warp
cp $sexdir/*.head $astroimadir

echo -e "\n ${GREEN} ---Warping and correcting distorsion--- ${NOCOLOUR}"
# Warp the data so we can:
#     1.- Place it in a proper grid
#     2.- Improve the astrometry thanks to scamp

# I lose 4 frames here. Why?
entiredir_smallGrid=$BDIR/pointings_smallGrid
entiredir_fullGrid=$BDIR/pointings_fullGrid
entiredone=$entiredir_smallGrid/done_.txt
swarpcfg=$ROOTDIR/"$objectName"/config/swarp.cfg
export swarpcfg

if ! [ -d $entiredir_smallGrid ]; then mkdir $entiredir_smallGrid; fi
if ! [ -d $entiredir_fullGrid ]; then mkdir $entiredir_fullGrid; fi

if [ -f $entiredone ]; then
    echo -e "\nsubs_sky_it1 images already with astromety corrected using scamp-swarp and regrid to final grid (stored in pointings)\n"
else
  imagesToWarp=()
  for a in $(seq 1 $totalNumberOfFrames); do
      base="$a".fits
      imagesToWarp+=("$astroimadir/$base")
  done
  printf "%s\n" "${imagesToWarp[@]}" | parallel -j "$num_cpus" warpImage {} $entiredir_fullGrid $entiredir_smallGrid $ra $dec $coaddSizePx $pipelinePath
  echo done > $entiredone
fi

exit 0

echo -e "${GREEN} --- Compute and subtract Sky --- ${NOCOLOUR}"

noiseskydir=$BDIR/noise-sky_it1
noiseskydone=$noiseskydir/done_"$filter"_ccd"$h".txt

subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it1
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter"_ccd"$h".txt

subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it1
subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter"_ccd"$h".txt

echo -e "Modelling the background for subtracting it"
imagesAreMasked=false
ringDir=$BDIR/ring

computeSky $entiredir_smallGrid $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing


# If we have not done it already (i.e. the modelling of the background selected has been a polynomial) we estimate de background as a constant for identifying bad frames
noiseskyctedir=$BDIR/noise-sky_it1_cte
noiseskyctedone=$noiseskyctedir/done_"$filter"_ccd"$h".txt
if [ "$MODEL_SKY_AS_CONSTANT" = false ]; then
  echo -e "\nModelling the background for the bad frame detection"
  computeSky $entiredir_smallGrid $noiseskyctedir $noiseskyctedone true $sky_estimation_method -1 false $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing
fi


# Checking and removing bad frames based on the background value ------
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
badFilesWarningsFile=identifiedBadFrames_backgroundValue.txt
badFilesWarningsDone=$diagnosis_and_badFilesDir/done_backgroundValue.txt
if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi
if [ -f $badFilesWarningsDone ]; then
    echo -e "\nbadFiles warning already done\n"
else
  # We choose the right directory in order to provide the constant estimation and not the polynomial (in case that it has been selected)
  if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
    tmpDir=$noiseskydir
  else
    tmpDir=$noiseskyctedir
  fi
  python3 $pythonScriptsPath/checkForBadFrames_backgroundValueAndStd.py $tmpDir $framesForCommonReductionDir $airMassKeyWord $diagnosis_and_badFilesDir $badFilesWarningsFile $numberOfStdForBadFrames
  echo done > $badFilesWarningsDone
fi


rejectedFramesDir=$BDIR/rejectedFrames_background
if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames with backgroundValue"
removeBadFramesFromReduction $entiredir_fullGrid $rejectedFramesDir $diagnosis_and_badFilesDir $badFilesWarningsFile
removeBadFramesFromReduction $entiredir_smallGrid $rejectedFramesDir $diagnosis_and_badFilesDir $badFilesWarningsFile


echo -e "\nSubtracting background"
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT
subtractSky $entiredir_fullGrid $subskyFullGrid_dir $subskyFullGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT


#### PHOTOMETRIC CALIBRATION  ####
echo -e "${ORANGE} ------ PHOTOMETRIC CALIBRATION ------ ${NOCOLOUR}\n"

### PARAMETERS ###
toleranceForMatching=2 #arcsec
sigmaForPLRegion=1 # Parameter for deciding the selection region (half-max-rad region)
export toleranceForMatching
export sigmaForPLRegion

# Parameters for performing the sigma clipping to the different samples in aststatistics
sigmaForStdSigclip=2
iterationsForStdSigClip=3
export sigmaForStdSigclip
export iterationsForStdSigClip

### PREPARING DECALS DATA FOR CALIBRATION ###
referenceImagesForMosaic=$entiredir_smallGrid
mosaicDir=$DIR/mosaic
selectedDecalsStarsDir=$mosaicDir/automaticallySelectedStarsForCalibration
rangeUsedDecalsDir=$mosaicDir/rangesUsedForCalibration


decalsImagesDir=$mosaicDir/decalsImages
prepareDecalsDataForPhotometricCalibration $referenceImagesForMosaic $decalsImagesDir $filter $ra $dec $mosaicDir $selectedDecalsStarsDir $rangeUsedDecalsDir $pixelScale $diagnosis_and_badFilesDir $sizeOfOurFieldDegrees

iteration=1
imagesForCalibration=$subskySmallGrid_dir
alphatruedir=$BDIR/alpha-stars-true_it$iteration
computeCalibrationFactors $iteration $imagesForCalibration $selectedDecalsStarsDir $rangeUsedDecalsDir $mosaicDir $decalsImagesDir $alphatruedir $calibrationBrightLimit $calibrationFaintLimit

# Checking and removing bad frames based on the FWHM value ------
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
fwhmFolder=$BDIR/my-catalog-halfmaxradius_it1
badFilesWarningsFile=identifiedBadFrames_fwhm.txt
badFilesWarningsDone=$diagnosis_and_badFilesDir/done_fwhmValue.txt
if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi
if [ -f $badFilesWarningsDone ]; then
    echo -e "\nbadFiles warning already done\n"
else
  python3 $pythonScriptsPath/checkForBadFrames_fwhm.py $fwhmFolder $diagnosis_and_badFilesDir $badFilesWarningsFile $numberOfStdForBadFrames
  echo done > $badFilesWarningsDone
fi


rejectedFramesDir=$BDIR/rejectedFrames_FWHM
if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames with FWHM"
removeBadFramesFromReduction $subskySmallGrid_dir $rejectedFramesDir $diagnosis_and_badFilesDir $badFilesWarningsFile
removeBadFramesFromReduction $subskyFullGrid_dir $rejectedFramesDir $diagnosis_and_badFilesDir $badFilesWarningsFile


# DIAGNOSIS PLOT
# Histogram of the background values on magnitudes / arcsec²
if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
  tmpDir=$noiseskydir
else
  tmpDir=$noiseskyctedir
fi
python3 $pythonScriptsPath/diagnosis_normalisedBackgroundMagnitudes.py $tmpDir $framesForCommonReductionDir $airMassKeyWord $alphatruedir $pixelScale $diagnosis_and_badFilesDir $BDIR/rejectedFrames_background $BDIR/rejectedFrames_FWHM


echo -e "\n ${GREEN} ---Applying calibration factors--- ${NOCOLOUR}"
photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir
applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir


# DIAGNOSIS PLOT
# Astrometry
astrometryPlotDone=$diagnosis_and_badFilesDir/done_astrometryPlot.txt
if [ -f $astrometryPlotDone ]; then
    echo -e "\nAstrometry diagnosis plot already done\n"
else
  produceAstrometryCheckPlot $fwhmFolder $BDIR/decals-aperture-catalog_it1 $pythonScriptsPath $diagnosis_and_badFilesDir $pixelScale
  echo done > $astrometryPlotDone
fi

# DIAGNOSIS PLOT
# Calibration
calibrationPlotDone=$diagnosis_and_badFilesDir/done_calibrationPlot.txt
if [ -f $calibrationPlotDone ]; then
    echo -e "\nCalibration diagnosis plot already done\n"
else
    produceCalibrationCheckPlot $BDIR/ourData-catalogs-apertures_it1 $photCorrSmallGridDir $fwhmFolder $BDIR/decals-aperture-catalog_it1 $pythonScriptsPath $diagnosis_and_badFilesDir $calibrationBrightLimit $calibrationFaintLimit
    echo done > $calibrationPlotDone
fi

# DIAGNOSIS PLOT
# Half-Max-Radius vs magnitude plots of our calibrated data
halfMaxRadiusVsMagnitudeOurDataDir=$diagnosis_and_badFilesDir/halfMaxRadVsMagPlots_ourData
halfMaxRadiusVsMagnitudeOurDataDone=$halfMaxRadiusVsMagnitudeOurDataDir/done_halfMaxRadVsMagPlots.txt

if ! [ -d $halfMaxRadiusVsMagnitudeOurDataDir ]; then mkdir $halfMaxRadiusVsMagnitudeOurDataDir; fi
if [ -f $halfMaxRadiusVsMagnitudeOurDataDone ]; then
    echo -e "\nHalf max radius vs magnitude plots for our calibrated data already done"
else
    produceHalfMaxRadVsMagForOurData $photCorrSmallGridDir $halfMaxRadiusVsMagnitudeOurDataDir $catdir/"$objectName"_Gaia_eDR3.fits $toleranceForMatching $pythonScriptsPath $num_cpus
    echo done > $halfMaxRadiusVsMagnitudeOurDataDone
fi



echo -e "${ORANGE} ------ STD WEIGHT COMBINATION ------ ${NOCOLOUR}\n"
# Compute rms and of the photometrized frames
noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done_"$k"_ccd"$h".txt
# Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
computeSky $photCorrSmallGridDir $noiseskydir $noiseskydone true $sky_estimation_method -1 false $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing


# Store the minimum standard deviation of the frames in order to compute the weights
python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration

### Calculate the weights for the images based on the minimum rms ###
echo -e "\n ${GREEN} ---Computing weights for the frames--- ${NOCOLOUR}"


wdir=$BDIR/weight-dir
wdone=$wdir/done_"$k"_ccd"$h".txt
if ! [ -d $wdir ]; then mkdir $wdir; fi

wonlydir=$BDIR/only-w-dir
wonlydone=$wonlydir/done_"$k"_ccd"$h".txt
if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
# We provide the fullGrid because we are going to combine then now
computeWeights $wdir $wdone $wonlydir $wonlydone $photCorrFullGridDir $noiseskydir $iteration


echo -e "\n ${GREEN} ---Masking outliers--- ${NOCOLOUR}"
# Remove outliers before the final coadd by using sigclip-median and sigclip-std
# This is particularly important to remove cosmic rays

sigmaForStdSigclip=2
clippingdir=$BDIR/clipping-outliers
clippingdone=$clippingdir/done_"$k".txt
buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip

mowdir=$BDIR/weight-dir-no-outliers
moonwdir=$BDIR/only-weight-dir-no-outliers
mowdone=$mowdir/done_"$k"_ccd"$h".txt
if ! [ -d $mowdir ]; then mkdir $mowdir; fi
if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi

if [ -f $mowdone ]; then
    echo -e "\nOutliers of the weighted images already masked\n"
else
    framesToRemoveOutliers=()
    for a in $(seq 1 $totalNumberOfFrames); do
        framesToRemoveOutliers+=("$a")
    done
    printf "%s\n" "${framesToRemoveOutliers[@]}" | parallel -j "$num_cpus" removeOutliersFromFrame {} $mowdir $moonwdir $clippingdir $wdir $wonlydir
    echo done > $mowdone 
fi


echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"

echo -e "\nBuilding coadd"
coaddDir=$BDIR/coadds
coaddName=$coaddDir/"$objectName"_coadd_"$filter".fits
buildCoadd $coaddDir $coaddName $mowdir $moonwdir

echo -e "\nAdding keywords to the coadd"
keyWords=("FRAMES_COMBINED" "FILTER" "SATURATION_THRESHOLD" "CALIBRATION_BRIGHTLIMIT" "CALIBRATION_FAINTLIMIT" "RUNNING_FLAT" "WINDOW_SIZE" "STD_FOR_BAD_FRAMES")
numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
values=("$numberOfFramesCombined" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames")
addkeywords $coaddName keyWords values

produceHalfMaxRadVsMagForSingleImage $coaddName $halfMaxRadiusVsMagnitudeOurDataDir $catdir/"$objectName"_Gaia_eDR3.fits $toleranceForMatching $pythonScriptsPath "coadd_it1"

maskName=$coaddir/"$objectName"_coadd_"$filter"_mask.fits
if [ -f $maskName ]; then
  echo "The mask of the weighted coadd is already done"
else
  astnoisechisel $coaddName $noisechisel_param -o $maskName
fi



framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted
framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi
if [ -f $framesWithCoaddSubtractedDone ]; then
    echo -e "\nFrames with coadd subtracted already generated\n"
else
  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub.fits
  subtractCoaddToFrames $photCorrFullGridDir $coaddName $framesWithCoaddSubtractedDir
  astarithmetic $(ls -v $framesWithCoaddSubtractedDir/*.fits) $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l) sum -g1 -o$sumMosaicAfterCoaddSubtraction
  echo done > $framesWithCoaddSubtractedDone 

fi


# Subtract a plane and build the coadd. Thus we have the constant background coadd and the plane background coadd
# if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
#   planeEstimationForCoaddDir=$BDIR/planeEstimationBeforeCoadd
#   planeEstimationForCoaddDone=$planeEstimationForCoaddDir/done.txt
#   polynomialDegree=1
#   if ! [ -d $planeEstimationForCoaddDir ]; then mkdir $planeEstimationForCoaddDir; fi
#   computeSky $mowdir $planeEstimationForCoaddDir $planeEstimationForCoaddDone false $sky_estimation_method $polynomialDegree false $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing

#   planeSubtractionForCoaddDir=$BDIR/planeSubtractionBeforeCoadd
#   if ! [ -d $planeSubtractionForCoaddDir ]; then mkdir $planeSubtractionForCoaddDir; fi

#   planeSubtractionForCoaddDone=$planeSubtractionForCoaddDir/done.txt
#   subtractSky $mowdir $planeSubtractionForCoaddDir $planeSubtractionForCoaddDone $planeEstimationForCoaddDir false

#   coaddDir=$BDIR/coadds_plane
#   coaddName=$coaddDir/"$objectName"_coadd_"$filter".fits

#   if ! [ -d $coaddDir ]; then mkdir $coaddDir; fi
#   buildCoadd $coaddDir $coaddName $planeSubtractionForCoaddDir $moonwdir

  # echo -e "\nAdding keywords to the coadd"
  # addkeywords $coaddName keyWords values
# fi









# # --- Build exposure map
# exposuremapDir=$coaddDir/exposureMap
# exposuremapdone=$coaddDir/done_"$k".txt

# if ! [ -d $exposuremapDir ]; then mkdir $exposuremapDir; fi
# if [ -f $exposuremapdone ]; then
#     echo -e "\nThe first weighted (based upon std) mean of the images already done\n"
# else
#   #There should be more efficient way of doing this...
#   # Maybe in the warp of the pointings. Do a matrix of 1 and 0 and just add them
#   # Pure exposure map
#   framesDir=$BDIR/pointings_smallGrid
#   for a in $(seq 1 $totalNumberOfFrames); do
#     i=$framesDir/entirecamera_$a.fits
#     astarithmetic $i set-i i i i eq 1 where i isblank 1 where -g1 --output="./tmp.fits"
#     SWarp -c $swarpcfg "./tmp.fits" -CENTER $ra,$dec -IMAGE_SIZE $coaddSizePx,$coaddSizePx -IMAGEOUT_NAME $exposuremapDir/swarp1.fits -WEIGHTOUT_NAME $exposuremapDir/swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE 1.164 -PIXELSCALE_TYPE  MANUAL
#     astarithmetic $exposuremapDir/swarp_w1.fits -h0 set-i i i 0 lt nan where -otemp1.fits
#     astarithmetic $exposuremapDir/swarp1.fits -h0 temp1.fits -h1 0 eq nan where -o$exposuremapDir/entirecamera_"$a".fits
#   done
#   rm $exposuremapDir/swarp_w1.fits $exposuremapDir/swarp1.fits
#   astarithmetic $(ls -v $exposuremapDir/*.fits) $(ls $exposuremapDir/*.fits | wc -l) number -g1 -o$coaddDir/exposureMap_NoNans.fits
#   rm $exposuremapDir/*.fits
#   rmdir $exposuremapDir

#   echo done > $exposuremapdone
# fi


# Remove intermediate folders to save some space
find $BDIR/noise-sky_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/sub-sky-fullGrid_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/sub-sky-smallGrid_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrFullGrid-dir_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrSmallGrid-dir_it1 -type f ! -name 'done*' -exec rm {} \;

find $BDIR/my-catalog-halfmaxradius_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/match-decals-myData_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/decals-aperture-catalog_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/ourData-catalogs-apertures_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/framesWithCoaddSubtractedDir -type f ! -name 'done*' -exec rm {} \;

find $wdir -type f ! -name 'done*' -exec rm {} \;
find $wonlydir -type f ! -name 'done*' -exec rm {} \;
find $mowdir -type f ! -name 'done*' -exec rm {} \;
find $moonwdir -type f ! -name 'done*' -exec rm {} \;


####### ITERATION 2 ######

iteration=2

# We mask the pointings in order to measure (before photometric calibration) the sky accurately
entiredir_smallGrid=$BDIR/pointings_smallGrid
entiredir_fullGrid=$BDIR/pointings_fullGrid
smallPointings_maskedDir=$BDIR/pointings_smallGrid_masked_it$iteration
maskedPointingsDone=$smallPointings_maskedDir/done_.txt

maskPointings $entiredir_smallGrid  $smallPointings_maskedDir $maskedPointingsDone $maskName $entiredir_fullGrid

noiseskydir=$BDIR/noise-sky_it$iteration
noiseskydone=$noiseskydir/done_"$filter"_ccd"$h".txt

subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it$iteration
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter"_ccd"$h".txt

subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it$iteration
subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter"_ccd"$h".txt

# compute sky with frames masked with global mask
imagesAreMasked=true
computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT
subtractSky $entiredir_fullGrid $subskyFullGrid_dir $subskyFullGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT


imagesForCalibration=$subskySmallGrid_dir
alphatruedir=$BDIR/alpha-stars-true_it$iteration


computeCalibrationFactors $iteration $imagesForCalibration $selectedDecalsStarsDir $rangeUsedDecalsDir $mosaicDir $decalsImagesDir $alphatruedir $calibrationBrightLimit $calibrationFaintLimit



photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration

applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir
applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir


# We mask again the points in order to measure (after photometric calibration) the sky accurately
smallPointings_photCorr_maskedDir=$BDIR/photCorrSmallGrid_masked_it$iteration
maskedPointingsDone=$smallPointings_photCorr_maskedDir/done_.txt
maskPointings $photCorrSmallGridDir $smallPointings_photCorr_maskedDir $maskedPointingsDone $maskName $entiredir_fullGrid


noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done_"$k"_ccd"$h".txt
# Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
computeSky $smallPointings_photCorr_maskedDir $noiseskydir $noiseskydone true $sky_estimation_method -1 true $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing

python3 $pythonScriptsPath/find_rms_min.py "$filter" 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration

wdir=$BDIR/weight-dir_it$iteration
wdone=$wdir/done_"$k"_ccd"$h".txt
if ! [ -d $wdir ]; then mkdir $wdir; fi
wonlydir=$BDIR/only-w-dir_it$iteration
wonlydone=$wonlydir/done_"$k"_ccd"$h".txt
if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi

# We provide the fullGrid because we are going to combine then now
computeWeights $wdir $wdone $wonlydir $wonlydone $photCorrFullGridDir $noiseskydir $iteration

clippingdir=$BDIR/clipping-outliers_it$iteration
clippingdone=$clippingdir/done_"$k".txt
buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip


mowdir=$BDIR/weight-dir-no-outliers_it$iteration
if ! [ -d $mowdir ]; then mkdir $mowdir; fi
# only weight
moonwdir=$BDIR/only-weight-dir-no-outliers_it$iteration
if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi
mowdone=$mowdir/done_"$k"_ccd"$h".txt

if [ -f $mowdone ]; then
    echo -e "\nOutliers of the weighted images already masked\n"
else
  framesToRemoveOutliers=()
  for a in $(seq 1 $totalNumberOfFrames); do
    framesToRemoveOutliers+=("$a")
  done
  printf "%s\n" "${framesToRemoveOutliers[@]}" | parallel -j "$num_cpus" removeOutliersFromFrame {} $mowdir $moonwdir $clippingdir $wdir $wonlydir
  echo done > $mowdone 
fi


echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
echo -e "\nBuilding coadd"

coaddDir=$BDIR/coadds_it$iteration 
coaddName=$coaddDir/"$objectName"_coadd_"$filter".fits
buildCoadd $coaddDir $coaddName $mowdir $moonwdir

echo -e "\nAdding keywords to the coadd"
keyWords=("FRAMES_COMBINED" "FILTER" "SATURATION_THRESHOLD" "CALIBRATION_BRIGHTLIMIT" "CALIBRATION_FAINTLIMIT" "RUNNING_FLAT" "WINDOW_SIZE" "STD_FOR_BAD_FRAMES")
numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
values=("$numberOfFramesCombined" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames")
addkeywords $coaddName keyWords values

produceHalfMaxRadVsMagForSingleImage $coaddName $halfMaxRadiusVsMagnitudeOurDataDir $catdir/"$objectName"_Gaia_eDR3.fits $toleranceForMatching $pythonScriptsPath "coadd_it2"

framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted_it$iteration
framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi

if [ -f $framesWithCoaddSubtractedDone ]; then
    echo -e "\nFrames with coadd subtracted already generated\n"
else
  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_it$iteration.fits
  subtractCoaddToFrames $photCorrFullGridDir $coaddName $framesWithCoaddSubtractedDir
  astarithmetic $(ls -v $framesWithCoaddSubtractedDir/*.fits) $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l) sum -g1 -o$sumMosaicAfterCoaddSubtraction
fi

endTime=$(date +%D%T)
echo "Pipeline ended at : ${endTime}"


exit 0



























maskName=$coaddir/"$objectName"_coadd_"$filter"_mask.fits
if [ -f $maskName ]; then
  echo "The mask of the weighted coadd is already done"
else
  astnoisechisel $coaddName "'$noisechisel_param'" -o $maskName
fi

# Remove intermediate folders to save some space
find $BDIR/noise-sky_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/sub-sky-fullGrid_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/sub-sky-smallGrid_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrFullGrid-dir_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrSmallGrid-dir_it2 -type f ! -name 'done*' -exec rm {} \;


find $BDIR/my-catalog-halfmaxradius_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/match-decals-myData_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/decals-aperture-catalog_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/ourData-catalogs-apertures_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/framesWithCoaddSubtractedDir_it2 -type f ! -name 'done*' -exec rm {} \;

find $BDIR/pointings_smallGrid_masked_it2 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrSmallGrid_masked_it2 -type f ! -name 'done*' -exec rm {} \;

find $wdir -type f ! -name 'done*' -exec rm {} \;
find $wonlydir -type f ! -name 'done*' -exec rm {} \;
find $mowdir -type f ! -name 'done*' -exec rm {} \;
find $moonwdir -type f ! -name 'done*' -exec rm {} \;



####### ITERATION 3 ######

iteration=3
entiredir_smallGrid=$BDIR/pointings_smallGrid
entiredir_fullGrid=$BDIR/pointings_fullGrid

# We mask the pointings in order to measure (before photometric calibration) the sky accurately
entiredir_fullGrid=$BDIR/pointings_fullGrid

smallPointings_maskedDir=$BDIR/pointings_smallGrid_masked_it$iteration
maskedPointingsDone=$smallPointings_maskedDir/done_.txt
maskPointings $entiredir_smallGrid $smallPointings_maskedDir $maskedPointingsDone $maskName $entiredir_fullGrid

noiseskydir=$BDIR/noise-sky_it$iteration
noiseskydone=$noiseskydir/done_"$filter"_ccd"$h".txt

subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it$iteration
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter"_ccd"$h".txt

subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it$iteration
subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter"_ccd"$h".txt

# compute sky with frames masked with global mask
imagesAreMasked=true
computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT
subtractSky $entiredir_fullGrid $subskyFullGrid_dir $subskyFullGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT

imagesForCalibration=$subskySmallGrid_dir
alphatruedir=$BDIR/alpha-stars-true_it$iteration
computeCalibrationFactors $iteration $imagesForCalibration $selectedDecalsStarsDir $rangeUsedDecalsDir $mosaicDir $decalsImagesDir $alphatruedir $calibrationBrightLimit $calibrationFaintLimit

photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration

applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir
applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir


# We mask again the points in order to measure (after photometric calibration) the sky accurately
smallPointings_photCorr_maskedDir=$BDIR/photCorrSmallGrid_masked_it$iteration
maskedPointingsDone=$smallPointings_photCorr_maskedDir/done_.txt
maskPointings $photCorrSmallGridDir $smallPointings_photCorr_maskedDir $maskedPointingsDone $maskName $entiredir_fullGrid

noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done_"$k"_ccd"$h".txt
# Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
computeSky $smallPointings_photCorr_maskedDir $noiseskydir $noiseskydone true $sky_estimation_method -1 true $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing
python3 $pythonScriptsPath/find_rms_min.py "$filter" 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration



wdir=$BDIR/weight-dir_it$iteration
wdone=$wdir/done_"$k"_ccd"$h".txt
if ! [ -d $wdir ]; then mkdir $wdir; fi
wonlydir=$BDIR/only-w-dir_it$iteration
wonlydone=$wonlydir/done_"$k"_ccd"$h".txt
if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi

# We provide the fullGrid because we are going to combine then now
computeWeights $wdir $wdone $wonlydir $wonlydone $photCorrFullGridDir $noiseskydir $iteration

clippingdir=$BDIR/clipping-outliers_it$iteration
clippingdone=$clippingdir/done_"$k".txt
buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip



# Fornax. Around 490 frames. Deimos, 20 cores. Around 1 h and 15 min
mowdir=$BDIR/weight-dir-no-outliers_it$iteration
if ! [ -d $mowdir ]; then mkdir $mowdir; fi
# only weight
moonwdir=$BDIR/only-weight-dir-no-outliers_it$iteration
if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi
mowdone=$mowdir/done_"$k"_ccd"$h".txt

if [ -f $mowdone ]; then
    echo -e "\nOutliers of the weighted images already masked\n"
else
  framesToRemoveOutliers=()
  for a in $(seq 1 $totalNumberOfFrames); do
    framesToRemoveOutliers+=("$a")
  done
  printf "%s\n" "${framesToRemoveOutliers[@]}" | parallel -j "$num_cpus" removeOutliersFromFrame {} $mowdir $moonwdir $clippingdir $wdir $wonlydir
  echo done > $mowdone 
fi


echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
echo -e "\nBuilding coadd"

coaddDir=$BDIR/coadds_it$iteration 
coaddName=$coaddDir/"$objectName"_coadd_"$filter".fits
buildCoadd $coaddDir $coaddName $mowdir $moonwdir

echo -e "\nAdding keywords to the coadd"
keyWords=("FRAMES_COMBINED" "FILTER" "SATURATION_THRESHOLD" "CALIBRATION_BRIGHTLIMIT" "CALIBRATION_FAINTLIMIT" "RUNNING_FLAT" "WINDOW_SIZE" "STD_FOR_BAD_FRAMES")
numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
values=("$numberOfFramesCombined" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames")
addkeywords $coaddName keyWords values


if [ -f $framesWithCoaddSubtractedDone ]; then
    echo -e "\nFrames with coadd subtracted already generated\n"
else
  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_it$iteration.fits
  subtractCoaddToFrames $photCorrFullGridDir $coaddName $framesWithCoaddSubtractedDir
  astarithmetic $(ls -v $framesWithCoaddSubtractedDir/*.fits) $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l) sum -g1 -o$sumMosaicAfterCoaddSubtraction
fi