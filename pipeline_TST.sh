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
#   · Sextractor conf files (.conv, .param and .sex). Two different sets, one for the astrometry and another one for the rest of detections

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

# In this file the time at which the main steps start are stored
fileForTimeStamps="./pipelineTimeStamps.txt" 

echo -e "\n ${GREEN} ---Loading pipeline Functions--- ${NOCOLOUR}"

# The path from which the pipeline is called is unknown. The functions for the pipeline are expected
# to be in the same folder as the pipeline, so we retrieve the path in order to run the functions file
pipelinePath=`dirname "$0"`
pipelinePath=`( cd "$pipelinePath" && pwd )`
pythonScriptsPath=$pipelinePath/pipelineScripts
export pipelinePath
export pythonScriptsPath

# Load the file with all the needed functions
source "$pipelinePath/pipeline_TST_functions.sh"

echo -e "\n ${GREEN} ---Loading Modules--- ${NOCOLOUR}"

gnuastroModuleName="gnuastro/0.22"
load_module $gnuastroModuleName

astrometryModuleName="astrometry.net/0.94"
load_module $astrometryModuleName

# Needed if using the SIE software
# pythonModuleName="python"
# load_module $pythonModuleName

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

writeTimeOfStepToFile "Starting pipeline" $fileForTimeStamps
echo -e "\n ${GREEN} ---Loading variables from conf file --- ${NOCOLOUR}"

confFile=$1
loadVariablesFromFile $confFile

checkIfAllVariablesAreSet

checkIfStringVariablesHaveValidValues

checkTransmittanceFilterAndItsUnits $telescope $surveyForPhotometry $folderWithTransmittances $filter

filterCorrectionCoeff=$( checkIfNeededFilterCorrectionIsGiven $telescope $filter $surveyForPhotometry $ROOTDIR/"$objectName"/config )
if [[ $filterCorrectionCoeff == 11 ]]; then
  echo "The filter corrections for the filter $filter, telescope $telescope and survey $survey were not found"
  echo "Exiting with error code 11"
  exit 11
fi


outputConfigurationVariablesInformation


# The following lines are responsible of the cpu's used for paralellise
# If it is running in a system with slurm it takes the number of cpu's from the slurm job
# Otherwise it takes it from the provided configuration file
num_cpus=$SLURM_CPUS_ON_NODE
if [ -z $num_cpus ]; then
  num_cpus=$defaultNumOfCPUs
fi

echo -e "\nNumber of CPUs allocated: $num_cpus"
export num_cpus


# ****** Decision note *******

# Rebinned data
tileSize=35
noisechisel_param="--tilesize=$tileSize,$tileSize \
                    --detgrowmaxholesize=1000 \
                    --rawoutput"

# # These paremeters are oriented to TST data at original resolution. 

# astmkprof --kernel=gaussian,2,3 --oversample=1 -o$ROOTDIR/"$objectName"/kernel.fits 
#HIPERCAM 
#tileSize=20
#noisechisel_param="--tilesize=$tileSize,$tileSize \
#                     --rawoutput"

export noisechisel_param

echo -e "\n-Noisechisel parameters used for masking:"
echo -e "\t" $noisechisel_param


########## Prepare data ##########

echo -e "\n ${GREEN} ---Preparing data--- ${NOCOLOUR}"

DIR=$ROOTDIR/"$objectName"
CDIR=$DIR/config
INDIRo=$ROOTDIR/"$objectName"/DATA-or
BDIR=$ROOTDIR/"$objectName"/build
INDIR=$ROOTDIR/"$objectName"/DATA
DARKDIR=$ROOTDIR/"$objectName"/dark
keyWordDirectory=$ROOTDIR/"$objectName"/keywords

export DIR
export INDIRo
export BDIR
export INDIR
export CDIR
export DARKDIR
export keyWordDirectory

if ! [ -d $CDIR ]; then mkdir $CDIR; fi
if ! [ -d $BDIR ]; then mkdir $BDIR; fi
if ! [ -d $INDIR ]; then mkdir $INDIR; fi
if ! [ -d $filtereyWordDirectory ]; then mkdir $filtereyWordDirectory; fi

echo -e "\n-Directories defined"
echo -e "\t·Main directory (DIR): ${ORANGE} ${DIR} ${NOCOLOUR}"
echo -e "\t·Build directory (BDIR): ${ORANGE} ${BDIR} ${NOCOLOUR}"
echo -e "\t·Original data directory (INDIRo): ${ORANGE} ${INDIRo} ${NOCOLOUR}"
echo -e "\t·Config directory ${ORANGE} ${CDIR} ${NOCOLOUR}"
echo -e "\t·Data directory (INDIR): ${ORANGE} ${INDIR} ${NOCOLOUR}"
echo -e "\t·Dark Data directory (DARKDIR): ${ORANGE} ${DARKDIR} ${NOCOLOUR}"
echo -e "\t·KeyWords directory (keyWordDirectory): ${ORANGE} ${keyWordDirectory} ${NOCOLOUR}"

# Getting the coordinates of the galaxy
ra=$ra_gal
dec=$dec_gal
export ra
export dec

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
    echo -e "\n\tScience images for night $currentNight are already processed\n"
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
            if [[ $dateHeaderKey =~ ^MJD ]]; then
              unixTimeInSeconds=$(astarithmetic $DATEOBS 40587 - 86400 x -q)
              unixTimeInSeconds=$(printf "%.0f" "$unixTimeInSeconds")
            else
              ## MACOS does not support -d in date, so it is better to use coreutils:gdata
              if [[ $OSTYPE == 'darwin'* ]]; then
                unixTimeInSeconds=$(gdate -d "$DATEOBS" +"%s")
              else
                unixTimeInSeconds=$(date -d "$DATEOBS" +"%s")
              fi
            fi
            out=$currentINDIR/$unixTimeInSeconds.fits

            # HERE A CHECK IF THE DATA IS IN FLOAT32 IS NEEDED
            if [[ "$overscan" == "YES" ]]; then
              trsec=$(eval "astfits $nameWithEscapedSpaces -h $h --keyvalue=$trimsecKey -q" )
              trsec=${trsec//[\[\]]/}
              eval "astcrop $nameWithEscapedSpaces -h$h --mode=img --section=$trsec -o$out"
            else
              eval "astfits $nameWithEscapedSpaces --copy=$h -o$out"  # I run this with eval so the escaped spaces are re-parsed by bash and understood by astfits
            fi
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
      gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
      if awk "BEGIN {exit !($gnuastro_version > 0.22)}"; then
        eval "astarithmetic $escaped_files $(ls -v $currentDARKDIR/* | wc -l) \
                    3 0.2 sigclip-mean -g$h --writeall \
                    -o $mdadir/temp.fits"
      else
        eval "astarithmetic $escaped_files $(ls -v $currentDARKDIR/* | wc -l) \
                    3 0.2 sigclip-mean -g$h  \
                    -o $mdadir/temp.fits"
      fi
      #If there is overscan
      if [[ "$overscan" == "YES" ]]; then
        first_file=$(echo "$escaped_files" | awk '{print $1}')
        trsec=$(eval "astfits $first_file -h$h --keyvalue=$trimsecKey -q")
        trsec=${trsec//[\[\]]/}
        astcrop $mdadir/temp.fits -h1 --mode=img --section=$trsec -o$mdadir/mdark_"$filter"_n"$currentNight"_ccd$h.fits
        rm $mdadir/temp.fits
      else
        mv $mdadir/temp.fits $mdadir/mdark_"$filter"_n"$currentNight"_ccd$h.fits
      fi
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
      air=$(astfits $i -h1 --keyvalue=$airMassKeyWord 2>/dev/null | awk '{print $2}')
      if [[ $air == "n/a" ]]; then
 		    air=$(python3 $pythonScriptsPath/get_airmass_teo.py $i $dateHeaderKey $ra_gal $dec_gal)
       	astfits $i --write=$airMassKeyWord,$air,"Updated from secz"
      fi
    	
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


  echo -e "${ORANGE} ------ FLATS ------ ${NOCOLOUR}\n"
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
    astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring.fits $commonRingDefinitionFile
  else
    astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_2.fits $secondRingDefinitionFile
    astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_1.fits $firstRingDefinitionFile
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

  #  **** Decision note *****
  # We do here the check for bad frames in std and for not including them in the flat
  # This is done because the frames which have moved sections have a really bad impact in the flat
  # The background is not cleaned here because we don't want to run out of frames for the flats

  # I use the flatit2WholeNightIma because maybe the running flat has not been selected, but the whole night flat is 
  # going to be constructed always (either by user selection or for correcting the running)

  # This is a simplified version of the more thorough check that is done in the future (with normalised background, std, skewness and kurtosis)
  # but since the data here is not much processed and I'm not sure how reliable is the background for that detailed study, we use a simplified version
  # using only the std in order to remove the frames with moved sections
  tmpNoiseDir=$BDIR/noisesky_forCleaningBadFramesBeforeFlat
  tmpNoiseDone=$tmpNoiseDir/done.txt
  if ! [ -d $tmpNoiseDir ]; then mkdir $tmpNoiseDir; fi

  diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
  badFilesWarningsFile=identifiedBadFrames_preFlat_onlyStd.txt
  badFilesWarningsDone=$diagnosis_and_badFilesDir/done_badFrames_stdPreFlat.txt
  if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi
  if [ -f $badFilesWarningsDone ]; then
      echo -e "\n\tFrames with strange background value and std values already cleaned\n"
  else
    computeSky $flatit2WholeNightimaDir $tmpNoiseDir $tmpNoiseDone true $sky_estimation_method -1 false $ringdir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth
    python3 $pythonScriptsPath/checkForBadFrames_beforeFlat_std.py  $tmpNoiseDir $diagnosis_and_badFilesDir $badFilesWarningsFile $numberOfStdForBadFrames
    echo "done" > $badFilesWarningsDone
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
  
  # Remove the identified bad frames ONLY for the flat, they will still be present in following steps, but not used in the flat calculation
  diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
  badFilesWarningsFile=identifiedBadFrames_preFlat_onlyStd.txt
  rejectedFramesDir=$BDIR/rejectedFrames_std_preFlat
  if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
  removeBadFramesFromReduction $normit3dir $rejectedFramesDir $diagnosis_and_badFilesDir $badFilesWarningsFile
  removeBadFramesFromReduction $normit3WholeNightdir $rejectedFramesDir $diagnosis_and_badFilesDir $badFilesWarningsFile


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
      astarithmetic $i -h1 set-m $currentFlatImage -h1 set-f m f $vignettingThreshold lt  nan where set-n n f 2. gt nan where -o $out
      propagateKeyword $i $airMassKeyWord $out 
      propagateKeyword $i $dateHeaderKey $out
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
      astfits $framesForCommonReductionDir/$name.fits -h1 --write=ORIGINAL_FILE,$base
    done
    echo "done" > $framesForCommonReductionDone
  fi

  # # Removing intermediate information to save space
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

writeTimeOfStepToFile "Process the individual nights" $fileForTimeStamps

nights=()
for currentNight in $(seq 1 $numberOfNights); do
      nights+=("$currentNight")
done
printf "%s\n" "${nights[@]}" | parallel --line-buffer -j "$num_cpus" oneNightPreProcessing {}

totalNumberOfFrames=$( ls $framesForCommonReductionDir/*.fits | wc -l)
export totalNumberOfFrames
echo -e "* Total number of frames to combine: ${GREEN} $totalNumberOfFrames ${NOCOLOUR} *"

# Up to this point the frame of every night has been corrected of bias-dark and flat.
# That corrections are perform night by night (because it's necessary for perform that corretions)
# Now, all the frames are "equal" so we do no distinction between nights.
# All the frames are stored together in $framesForCommonReductionDir with names 1.fits, 2.fits, 3.fits ... n.fits.
echo -e "\n${GREEN} --- Astrometry --- ${NOCOLOUR}\n"

writeTimeOfStepToFile "Download Gaia catalogue" $fileForTimeStamps
echo -e "·Downloading Gaia Catalogue"

# Here I add some extra size to the field used to download the gaia catalogue
# This is because you don't want to use a field too big in order not to download bricks that you don't need
# So I expect the "sizeofOurFieldDegrees" value to be quite tight. But since the catalogue is text and it doesn't take
# long I prefer to add something and be sure that I don't lose any source because of the catalogue
radiusToDownloadCatalogue=$( echo "$sizeOfOurFieldDegrees + 0.5" | bc | awk '{printf "%.1f", $0}' ) #The awk part is to avoiod problems when R<1
query_param="gaia --dataset=edr3 --center=$ra_gal,$dec_gal --radius=$radiusToDownloadCatalogue --column=ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error"
catdir=$BDIR/catalogs
catName=$catdir/"$objectName"_Gaia_eDR3.fits
catdone=$catdir/done.txt
if ! [ -d $catdir ]; then mkdir $catdir; fi
if [ -f $catdone ]; then
  echo -e "\n\tCatalogue is already downloaded\n"
else
  downloadGaiaCatalogue "$query_param" $catdir $catName
  echo "done" > $catdone
fi


# Making the indexes
writeTimeOfStepToFile "Download Indices for astrometrisation" $fileForTimeStamps
echo -e "·Downloading Indices for astrometrisation"

indexdir=$BDIR/indexes
indexdone=$indexdir/done_"$filter".txt
if ! [ -d $indexdir ]; then mkdir $indexdir; fi
if [ -f $indexdone ]; then
  echo -e "\n\tGaia eDR3 indexes are already created\n"
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


sexcfg_sf=$CDIR/sextractor_solvefield.sex #Solving the images
writeTimeOfStepToFile "Solving fields" $fileForTimeStamps
echo -e "·Solving fields"

astrocfg=$CDIR/astrometry_$objectName.cfg

rm $astrocfg
echo inparallel > $astrocfg
echo cpulimit 300 >> $astrocfg
echo "add_path $indexdir" >> $astrocfg
echo autoindex >> $astrocfg

astroimadir=$BDIR/astro-ima
astroimadone=$astroimadir/done_"$filter".txt
if ! [ -d $astroimadir ]; then mkdir $astroimadir; fi
if [ -f $astroimadone ]; then
  echo -e "\n\tImages are already astrometrized\n"
else
  frameNames=()
  for a in $(seq 1 $totalNumberOfFrames); do
      base=$a.fits
      i=$framesForCommonReductionDir/$base
      frameNames+=("$i")
  done
  printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" solveField {} $solve_field_L_Param $solve_field_H_Param $solve_field_u_Param $ra_gal $dec_gal $CDIR $astroimadir $sexcfg_sf
  echo done > $astroimadone
fi


########## Distorsion correction ##########
echo -e "\n ${GREEN} ---Creating distorsion correction files--- ${NOCOLOUR}"


# Making sex catalogs and running scamp

# ****** Decision note *******
# In the LBT pipeline this two steps were done sequentially in two different blocks, first sextractor and then scamp
# I have put them together so we can loop them easily, because to perform an iterative astrometrisation we need to do
# sextractor + scamp, so for doing it in a confortable manner in the code both steps are in the same block
writeTimeOfStepToFile "Making sextractor catalogues and running scamp" $fileForTimeStamps
echo -e "·Creating SExtractor catalogues and running scamp"

numOfSextractorPlusScampIterations=3

sexcfg=$CDIR/sextractor_astrometry.sex
sexparam=$CDIR/sextractor_astrometry.param
sexconv=$CDIR/default.conv
sexdir=$BDIR/sex-it1

scampcfg=$CDIR/scamp.cfg
scampdir=$BDIR/scamp-it1
scampres=$scampdir/results_Decals-"$filter"
scampdone=$scampdir/done_"$filter".txt

if ! [ -d $sexdir ]; then mkdir $sexdir; fi
if ! [ -d $scampdir ]; then mkdir $scampdir; fi
if ! [ -d $scampres ]; then mkdir $scampres; fi

if [ -f $scampdone ]; then
    echo -e "\n\tSex catalogs and scamp are already done for extension $h\n"
else
  frameNames=()
  for a in $(seq 1 $totalNumberOfFrames); do
      frameNames+=("$a")
  done

  for ((i = 1; i <= numOfSextractorPlusScampIterations; i++)); do
    echo -e "\tSExtractor + scamp iteration $i"

    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" runSextractorOnImage {} $sexcfg $sexparam $sexconv $astroimadir $sexdir $saturationThreshold $gain
    scamp -c $scampcfg $(ls -v $sexdir/*.cat)
    cp $sexdir/*.head $astroimadir
    mv *.pdf $scampres/
    mv scamp.xml $scampdir
  done
  echo done > $scampdone
fi

echo -e "\n ${GREEN} ---Warping and correcting distorsion--- ${NOCOLOUR}"
writeTimeOfStepToFile "Warping frames" $fileForTimeStamps
# Warp the data so we can:
#     1.- Place it in a proper grid
#     2.- Improve the astrometry thanks to scamp

entiredir_smallGrid=$BDIR/pointings_smallGrid
entiredir_fullGrid=$BDIR/pointings_fullGrid
entiredone=$entiredir_smallGrid/done_.txt
swarpcfg=$ROOTDIR/"$objectName"/config/swarp.cfg
export swarpcfg

if ! [ -d $entiredir_smallGrid ]; then mkdir $entiredir_smallGrid; fi
if ! [ -d $entiredir_fullGrid ]; then mkdir $entiredir_fullGrid; fi

if [ -f $entiredone ]; then
    echo -e "\n\tsubs_sky_it1 images already with astromety corrected using scamp-swarp and regrid to final grid (stored in pointings)\n"
else
  imagesToWarp=()

  for a in $(seq 1 $totalNumberOfFrames); do
      base="$a".fits
      imagesToWarp+=($astroimadir/$base)
  done

  printf "%s\n" "${imagesToWarp[@]}" | parallel -j "$num_cpus" warpImage {} $entiredir_fullGrid $entiredir_smallGrid $ra $dec $coaddSizePx $pipelinePath
  echo done > $entiredone
fi


# Checking bad astrometrised frames ------
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
badFilesWarningsFile=identifiedBadFrames_astrometry.txt
badFilesWarningsDone=$diagnosis_and_badFilesDir/done_badFrames_astrometry.txt
if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi
if [ -f $badFilesWarningsDone ]; then
    echo -e "\n\tBad astrometrised frames warning already done\n"
else
  scampXMLFilePath=$scampdir/scamp.xml
  python3 $pythonScriptsPath/checkForBadFrames_badAstrometry.py $diagnosis_and_badFilesDir $scampXMLFilePath $badFilesWarningsFile $entiredir_fullGrid
  echo done > $badFilesWarningsDone
fi
#Building a .fits with the scamp contrast parameter
contrastdir_fullGrid=$BDIR/contrast_xy_maps
contrastdone_fullGrid=$contrastdir_fullGrid/done.txt
if ! [ -d $contrastdir_fullGrid ]; then mkdir $contrastdir_fullGrid; fi
if [ -f $contrastdone_fullGrid ]; then
    echo -e "\n\tMap of XY contrast parameters already done.\n"
else
  for frame in $(ls $entiredir_fullGrid/enti*.fits); do
    frame_out=$contrastdir_fullGrid/${frame#$entiredir_fullGrid/}
    xy_param=$(astfits $frame -h1 --keyvalue=XY-contrast -q)
    astarithmetic $frame -h1 0 x $xy_param +  -o $frame_out
  done
  if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
    astarithmetic $(ls $contrastdir_fullGrid/ent*.fits) $(ls $contrastdir_fullGrid/ent*.fits | wc -l) sum -g1 --writeall -o $contrastdir_fullGrid/sum_xycontrasts.fits
  else
    astarithmetic $(ls $contrastdir_fullGrid/ent*.fits) $(ls $contrastdir_fullGrid/ent*.fits | wc -l) sum -g1  -o $contrastdir_fullGrid/sum_xycontrasts.fits
  fi
  echo done > $contrastdone_fullGrid
  rm $contrastdir/fullGrid/ent*.fits
fi

echo -e "${GREEN} --- Compute and subtract Sky --- ${NOCOLOUR} \n"


noiseskydir=$BDIR/noise-sky_it1
noiseskydone=$noiseskydir/done_"$filter".txt

subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it1
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter".txt

subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it1
subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter".txt

echo -e "·Modelling the background for subtracting it"
imagesAreMasked=false
ringDir=$BDIR/ring

writeTimeOfStepToFile "Computing sky" $fileForTimeStamps
computeSky $entiredir_smallGrid $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth

# If we have not done it already (i.e. the modelling of the background selected has been a polynomial) we estimate de background as a constant for identifying bad frames
noiseskyctedir=$BDIR/noise-sky_it1_cte
noiseskyctedone=$noiseskyctedir/done_"$filter".txt
if [ "$MODEL_SKY_AS_CONSTANT" = false ]; then
  echo -e "\nModelling the background for the bad frame detection"
  computeSky $entiredir_smallGrid $noiseskyctedir $noiseskyctedone true $sky_estimation_method -1 false $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth
fi

# Checking and removing bad frames based on the background value ------
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
badFilesWarningsBackgroundValueFile=identifiedBadFrames_backgroundValue.txt
badFilesWarningsBackgroundStdFile=identifiedBadFrames_backgroundStd.txt
badFilesWarningsDone=$diagnosis_and_badFilesDir/done_badFrames_backgroundProperties.txt

if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi
if [ -f $badFilesWarningsDone ]; then
    echo -e "\n\tbadFiles warning already done\n"
else
  # We choose the right directory in order to provide the constant estimation and not the polynomial (in case that it has been selected)
  if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
    tmpDir=$noiseskydir
  else
    tmpDir=$noiseskyctedirk
  fi
  python3 $pythonScriptsPath/checkForBadFrames_backgroundValueAndStd.py $tmpDir $framesForCommonReductionDir $airMassKeyWord $diagnosis_and_badFilesDir $badFilesWarningsBackgroundValueFile $badFilesWarningsBackgroundStdFile $numberOfStdForBadFrames false
  echo done > $badFilesWarningsDone
fi

echo -e "\n·Subtracting background"
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT
subtractSky $entiredir_fullGrid $subskyFullGrid_dir $subskyFullGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT


# ### BUILD A FIRST COADD FROM SKY SUBTRACTION ####
# echo -e "${GREEN} --- Coadding before photometric calibration --- ${NOCOLOUR} \n"
# writeTimeOfStepToFile "Building coadd before photometry" $fileForTimeStamps
# iteration=1
# h=0
# minRmsFileName=min_rms_prev.txt
# noisesky_prephot=$BDIR/noise-sky_prephot
# noisesky_prephotdone=$noisesky_prephot/done_$filter.txt
# imagesAreMasked=false
# if ! [ -d $noisesky_prephot ]; then mkdir $noisesky_prephot; fi
# if [ -f $noisesky_prephot ]; then
# 	echo -e "\n Coadd pre-photometry already done\n"
# else
# 	computeSky $subskyFullGrid_dir $noisesky_prephot $noisesky_prephotdone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth

# 	python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $h $noisesky_prephot $DIR $iteration min_rms_prev.txt

# 	#Compute weights
# 	wdir=$BDIR/weight-dir-prephot
# 	wdone=$wdir/done.txt
# 	if ! [ -d $wdir ]; then mkdir $wdir; fi

# 	wonlydir=$BDIR/only-w-dir-prephot
# 	wonlydone=$wonlydir/done.txt
# 	if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi

# 	computeWeights $wdir $wdone $wonlydir $wonlydone $subskyFullGrid_dir $noisesky_prephot $iteration $minRmsFileName

# 	##Mask the outliers
# 	sigmaForStdSigclip=2
# 	clippingdir=$BDIR/clipping-outliers-prephot
# 	clippingdone=$clippingdir/done.txt
# 	buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip
	

# 	mowdir=$BDIR/weight-dir-no-outliers-prephot
# 	moonwdir=$BDIR/only-weight-dir-no-outliers-prephot
# 	mowdone=$mowdir/done.txt
# 	if ! [ -d $mowdir ]; then mkdir $mowdir; fi
# 	if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi
# 	removeOutliersFromWeightedFrames $mowdone $totalNumberOfFrames $mowdir $moonwdir $clippingdir $wdir $wonlydir

#   echo -e "\n·Removing bad frames"

#   diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
#   rejectedFramesDir=$BDIR/rejectedFrames
#   if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
#   echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"

#   rejectedByAstrometry=identifiedBadFrames_astrometry.txt
#   removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry
#   removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry

#   rejectedByBackgroundStd=identifiedBadFrames_backgroundStd.txt
#   removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd
#   removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd

#   rejectedByBackgroundValue=identifiedBadFrames_backgroundValue.txt
#   removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue
#   removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue

# 	##Make the coadd
# 	coaddDir=$BDIR/coadds-prephot
# 	coaddDone=$coaddDir/done.txt
# 	coaddName=$coaddDir/"$objectName"_coadd_"$filter"_prephot.fits
# 	buildCoadd $coaddDir $coaddName $mowdir $moonwdir $coaddDone

#   exposuremapDir=$coaddDir/"$objectName"_exposureMap
#   exposuremapdone=$coaddDir/done_exposureMap.txt
#   computeExposureMap $framesDir $exposuremapDir $exposuremapdone
# fi


#### PHOTOMETRIC CALIBRATION  ####
echo -e "${ORANGE} ------ PHOTOMETRIC CALIBRATION ------ ${NOCOLOUR}\n"

writeTimeOfStepToFile "Photometric calibration" $fileForTimeStamps

### PARAMETERS ###
toleranceForMatching=3 #arcsec
sigmaForPLRegion=3 # Parameter for deciding the selection region (half-max-rad region)
export toleranceForMatching
export sigmaForPLRegion

# Parameters for performing the sigma clipping to the different samples in aststatistics
sigmaForStdSigclip=1
iterationsForStdSigClip=3
export sigmaForStdSigclip
export iterationsForStdSigClip

### PREPARING DATA FOR CALIBRATION ###
referenceImagesForMosaic=$entiredir_smallGrid
mosaicDir=$DIR/mosaic
selectedCalibrationStarsDir=$mosaicDir/automaticallySelectedStarsForCalibration
rangeUsedCalibrationDir=$mosaicDir/rangesUsedForCalibration
aperturePhotDir=$mosaicDir/aperturePhotometryCatalogues # This is the final product that "prepareCalibrationData" produces and will be used in "computeCalibrationFactors"
mosaicDone=$mosaicDir/done_prep.txt
# ****** Decision note *******
# Since the calibration factors obtained with PANSTARRS imaging, GAIA spectra and SDDS spectra do NOT completely agree,
# we have decided to calibrate to GAIA spectra. Thus, we have estimated the aperture needed in PANSTARRS (XRe) to recover
# magnitudes obtained with GAIA spectra. When doing the tests for estimated this aperture we find that in certain fields we find and offset. 
# For solving this we compute this offset and correct it in each run of the pipeline (thus PANSTARRS always agreeing with GAIA)\\ 
# GAIA has been chosen over SDSS because we have more spectra, the calibration is more stable, and we have it in the southern hemisphere. 
#It is true that GAIA sources are quite bright (for TST is fine but would be problematic for other telescopes) but since we only need to calibrate Halpha 
# (much harder to saturate in that band) from bigger telescopes we expect to be fine.\\
# Additionally a correction between the survey filter (panstarrs, etc...) and your filter is applied. This is a offset introduced in the configuration file
prepareCalibrationData $surveyForPhotometry $referenceImagesForMosaic $aperturePhotDir $filter $ra $dec $mosaicDir $selectedCalibrationStarsDir $rangeUsedCalibrationDir \
                                            $pixelScale $sizeOfOurFieldDegrees $catName $surveyForSpectra $apertureUnits $folderWithTransmittances "$filterCorrectionCoeff" $surveyCalibrationToGaiaBrightLimit $surveyCalibrationToGaiaFaintLimit $mosaicDone


iteration=1
imagesForCalibration=$subskySmallGrid_dir
alphatruedir=$BDIR/alpha-stars-true_it$iteration
matchdir=$BDIR/match-decals-myData_it$iteration

writeTimeOfStepToFile "Computing calibration factors" $fileForTimeStamps
computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $rangeUsedCalibrationDir $mosaicDir  \
                          $alphatruedir $calibrationBrightLimit $calibrationFaintLimit $tileSize $apertureUnits $numberOfApertureUnitsForCalibration

# Creating histogram with the number of stars used for the calibratino of each frame
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi

numberOfStarsUsedInEachFramePlot=$diagnosis_and_badFilesDir/numOfStarsUsedForCalibrationHist.png
numberOfStarsUsedInEachFrameDone=$diagnosis_and_badFilesDir/done_numOfStarsUsedForCalibrate.txt
if [ -f $numberOfStarsUsedInEachFrameDone ]; then
  echo -e "\nHistogram with the number of stars used for calibrating each plot already done"
else
  python3 $pythonScriptsPath/diagnosis_numOfStarsUsedInCalibration.py $alphatruedir/numberOfStarsUsedForCalibrate.txt $numberOfStarsUsedInEachFramePlot
  echo done > $numberOfStarsUsedInEachFrameDone
fi

# Checking and removing bad frames based on the FWHM value ------
fwhmFolder=$BDIR/seeing_values
badFilesWarningsFile=identifiedBadFrames_fwhm.txt
badFilesWarningsDone=$diagnosis_and_badFilesDir/done_fwhmValue.txt
if [ -f $badFilesWarningsDone ]; then
    echo -e "\nbadFiles warning already done\n"
else
  if ! [ -d $fwhmFolder ]; then mkdir $fwhmFolder; fi
  imagesToFWHM=()
  for a in $(seq 1 $totalNumberOfFrames); do
    base="$a".fits
    imagesToFWHM+=("$base")
  done
  printf "%s\n" "${imagesToFWHM[@]}" | parallel -j "$num_cpus" computeFWHMSingleFrame {} $subskySmallGrid_dir $fwhmFolder 1 $methodToUse $tileSize 
  python3 $pythonScriptsPath/checkForBadFrames_fwhm.py $fwhmFolder $diagnosis_and_badFilesDir $badFilesWarningsFile $numberOfStdForBadFrames $framesForCommonReductionDir $pixelScale
  echo done > $badFilesWarningsDone
fi

# DIAGNOSIS PLOT
# Histogram of the background values on magnitudes / arcsec²
if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
  tmpDir=$BDIR/noise-sky_it1
else
  tmpDir=$noiseskyctedir
fi
python3 $pythonScriptsPath/diagnosis_normalisedBackgroundMagnitudes.py $tmpDir $framesForCommonReductionDir $airMassKeyWord $alphatruedir $pixelScale $diagnosis_and_badFilesDir


echo -e "\n ${GREEN} ---Applying calibration factors--- ${NOCOLOUR}"
photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir
applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir


# DIAGNOSIS PLOTs ---------------------------------------------------

# Astrometry
astrometryPlotName=$diagnosis_and_badFilesDir/astrometry.png
if [ -f $astrometryPlotName ]; then
    echo -e "\nAstrometry diagnosis plot already done\n"
else
  produceAstrometryCheckPlot $matchdir $pythonScriptsPath $astrometryPlotName $pixelScale
fi


# Calibration
aperturesFolder=$BDIR/my-catalog-halfmaxradius_it1
calibrationPlotName=$diagnosis_and_badFilesDir/calibrationPlot.png
if [ -f $calibrationPlotName ]; then
    echo -e "\nCalibration diagnosis plot already done\n"
else

    if [[ "$surveyForCalibration" == "SPECTRA" ]]; then
      dirWithReferenceCat=$mosaicDir/aperturePhotometryCatalogues
    else
      dirWithReferenceCat=$BDIR/survey-aperture-photometry_perBrick_it1
    fi

    produceCalibrationCheckPlot $BDIR/ourData-aperture-photometry_it1 $photCorrSmallGridDir $aperturesFolder $dirWithReferenceCat \
                                  $pythonScriptsPath $calibrationPlotName $calibrationBrightLimit $calibrationFaintLimit $numberOfApertureUnitsForCalibration $diagnosis_and_badFilesDir $surveyForPhotometry $BDIR  
fi


# Half-Max-Radius vs magnitude plots of our calibrated data
halfMaxRadiusVsMagnitudeOurDataDir=$diagnosis_and_badFilesDir/halfMaxRadVsMagPlots_ourData
halfMaxRadiusVsMagnitudeOurDataDone=$halfMaxRadiusVsMagnitudeOurDataDir/done_halfMaxRadVsMagPlots.txt
if ! [ -d $halfMaxRadiusVsMagnitudeOurDataDir ]; then mkdir $halfMaxRadiusVsMagnitudeOurDataDir; fi
if [ -f $halfMaxRadiusVsMagnitudeOurDataDone ]; then
    echo -e "\nHalf max radius vs magnitude plots for our calibrated data already done"
else
  
  produceHalfMaxRadVsMagForOurData $photCorrSmallGridDir $halfMaxRadiusVsMagnitudeOurDataDir $catdir/"$objectName"_Gaia_eDR3.fits $toleranceForMatching $pythonScriptsPath $num_cpus 30 $apertureUnits
  echo done > $halfMaxRadiusVsMagnitudeOurDataDone
fi

# ---------------------------------------------------

echo -e "\n${ORANGE} ------ STD WEIGHT COMBINATION ------ ${NOCOLOUR}\n"
# Compute rms and of the photometrized frames
noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done.txt
# Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
computeSky $photCorrSmallGridDir $noiseskydir $noiseskydone true $sky_estimation_method -1 false $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth


# Store the minimum standard deviation of the frames in order to compute the weights
h=0
minRmsFileName=min_rms_it1.txt
python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration min_rms_it1.txt


### Calculate the weights for the images based on the minimum rms ###
echo -e "\n ${GREEN} ---Computing weights for the frames--- ${NOCOLOUR}"
writeTimeOfStepToFile "Computing frame weights" $fileForTimeStamps

wdir=$BDIR/weight-dir
wdone=$wdir/done.txt
if ! [ -d $wdir ]; then mkdir $wdir; fi

wonlydir=$BDIR/only-w-dir
wonlydone=$wonlydir/done.txt
if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
# We provide the fullGrid because we are going to combine then now
computeWeights $wdir $wdone $wonlydir $wonlydone $photCorrFullGridDir $noiseskydir $iteration $minRmsFileName


echo -e "\n ${GREEN} ---Masking outliers--- ${NOCOLOUR}"
writeTimeOfStepToFile "Masking outliers" $fileForTimeStamps
# Remove outliers before the final coadd by using sigclip-median and sigclip-std
# This is particularly important to remove cosmic rays

sigmaForStdSigclip=2
clippingdir=$BDIR/clipping-outliers
clippingdone=$clippingdir/done.txt
buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip


mowdir=$BDIR/weight-dir-no-outliers
moonwdir=$BDIR/only-weight-dir-no-outliers
mowdone=$mowdir/done.txt
if ! [ -d $mowdir ]; then mkdir $mowdir; fi
if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi

removeOutliersFromWeightedFrames $mowdone $totalNumberOfFrames $mowdir $moonwdir $clippingdir $wdir $wonlydir


echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
writeTimeOfStepToFile "Building coadd" $fileForTimeStamps


echo -e "\n·Removing bad frames"
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
rejectedFramesDir=$BDIR/rejectedFrames
if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"

rejectedByAstrometry=identifiedBadFrames_astrometry.txt
removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry
removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry

rejectedByBackgroundStd=identifiedBadFrames_backgroundStd.txt
removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd
removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd

rejectedByBackgroundValue=identifiedBadFrames_backgroundValue.txt
#removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue
#removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue

rejectedByBackgroundFWHM=identifiedBadFrames_fwhm.txt
#removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM
#removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM


echo -e "\n·Building coadd"
coaddDir=$BDIR/coadds
coaddDone=$coaddDir/done.txt
coaddName=$coaddDir/"$objectName"_coadd_"$filter".fits
buildCoadd $coaddDir $coaddName $mowdir $moonwdir $coaddDone

maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  astnoisechisel $coaddName $noisechisel_param --numthreads=$num_cpus -o $maskName
fi

#astnoisechisel with the current parameters might fail due to long tilesize. I'm gonna make 2 checks to see if it fails, decreasing in steps of 5 in tilesize
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  #We assume that if this works for this iteration, then the next one will need at least same parameters
  tileSize=$((tileSize - 5))
  noisechisel_param="--tilesize=$tileSize,$tileSize \
                    --erode=1 \
                    --detgrowmaxholesize=5000 \
                    --rawoutput"
  astnoisechisel $coaddName $noisechisel_param --numthreads=$num_cpus -o $maskName
fi
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  #We assume that if this works for this iteration, then the next one will need at least same parameters
  tileSize=$((tileSize - 5))
  noisechisel_param="--tilesize=$tileSize,$tileSize \
                    --erode=1 \
                    --detgrowmaxholesize=5000 \
                    --rawoutput"
  astnoisechisel $coaddName $noisechisel_param --numthreads=$num_cpus  -o $maskName
fi

exposuremapDir=$coaddDir/"$objectName"_exposureMap
exposuremapdone=$coaddDir/done_exposureMap.txt
computeExposureMap $framesDir $exposuremapDir $exposuremapdone

#Compute surface brightness limit
sblimitFile=$coaddDir/"$objectName"_"$filter"_sblimit.txt
exposuremapName=$coaddDir/exposureMap.fits
if [ -f  $sblimitFile ]; then
    echo -e "\n\tSurface brightness limit for coadd already measured\n"
else
    surfaceBrightnessLimit=$( limitingSurfaceBrightness $coaddName $maskName $exposuremapName $coaddDir $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )
fi

times=($(getInitialMidAndFinalFrameTimes $INDIR))
initialTime=$( date -d @"${times[0]}" "+%Y-%m-%d_%H:%M:%S")
meanTime=$( date -d @"${times[1]}" "+%Y-%m-%d_%H:%M:%S")
finalTime=$( date -d @"${times[2]}" "+%Y-%m-%d_%H:%M:%S")

keyWords=("FRAMES_COMBINED" \
          "NUMBER_OF_DIFFERENT_NIGHTS" \
          "INITIAL_DATE_OBS" \
          "MEAN_DATE_OBS" \
          "FINAL_DATE_OBS" \
          "FILTER" \
          "VIGNETTING_THRESHOLD" \
          "SATURATION_THRESHOLD" \
          "CALIBRATED_USING" \
          "CALIBRATION_BRIGHTLIMIT" \
          "CALIBRATION_FAINTLIMIT" \
          "RUNNING_FLAT" \
          "WINDOW_SIZE" \
          "STD_FOR_BAD_FRAMES" \
          "SURFACE_BRIGHTNESS_LIMIT")

numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$vignettingThreshold" "$saturationThreshold" "$surveyForPhotometry" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames" "$surfaceBrightnessLimit")
comments=("" "" "" "" "" "" "" "" "" "" "" "Num. of tandard deviations used for rejection" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")
astfits $coaddName --write=/,"Pipeline information"
addkeywords $coaddName keyWords values comments


halfMaxRadForCoaddName=$halfMaxRadiusVsMagnitudeOurDataDir/coadd_it1.png
if [ -f $halfMaxRadForCoaddName ]; then
  echo -e "\tThe Half-Max-Rad vs Magnitude has been already generate for the coadd"
else
  produceHalfMaxRadVsMagForSingleImage $coaddName $halfMaxRadiusVsMagnitudeOurDataDir $catdir/"$objectName"_Gaia_eDR3.fits $toleranceForMatching $pythonScriptsPath "coadd_it1" 100
fi

writeTimeOfStepToFile "Producing frames with coadd subtracted" $fileForTimeStamps
framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted
framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi
if [ -f $framesWithCoaddSubtractedDone ]; then
    echo -e "\n\tFrames with coadd subtracted already generated\n"
else
  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub.fits
  subtractCoaddToFrames $photCorrFullGridDir $coaddName $framesWithCoaddSubtractedDir
  astarithmetic $(ls -v $framesWithCoaddSubtractedDir/*.fits) $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l) sum -g1 -o$sumMosaicAfterCoaddSubtraction
  echo done > $framesWithCoaddSubtractedDone 
fi


# Subtract a plane and build the coadd. Thus we have the constant background coadd and the plane background coadd
# if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
#   coaddPlaneDone=$coaddDir/done_plane.txt
#   if [ -f $coaddPlaneDone ]; then
#     echo -e "\n\tCoadd with background modelled as a plane already done"
#   else 
#     planeEstimationForCoaddDir=$BDIR/planeEstimationBeforeCoadd_it$iteration
#     planeEstimationForCoaddDone=$planeEstimationForCoaddDir/done.txt
#     polynomialDegree=1
#     if ! [ -d $planeEstimationForCoaddDir ]; then mkdir $planeEstimationForCoaddDir; fi
#     computeSky $photCorrSmallGridDir $planeEstimationForCoaddDir $planeEstimationForCoaddDone false $sky_estimation_method $polynomialDegree false $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth

#     planeSubtractionForCoaddDir=$BDIR/planeSubtractionBeforeCoadd_it$iteration
#     if ! [ -d $planeSubtractionForCoaddDir ]; then mkdir $planeSubtractionForCoaddDir; fi
#     planeSubtractionForCoaddDone=$planeSubtractionForCoaddDir/done.txt
#     subtractSky $photCorrFullGridDir $planeSubtractionForCoaddDir $planeSubtractionForCoaddDone $planeEstimationForCoaddDir false

#     noiseskyPlaneDir=$BDIR/noise-sky_forPlaneCoadd_it$iteration
#     if ! [ -d $noiseskyPlaneDir ]; then mkdir $noiseskyPlaneDir; fi
#     noiseskyPlaneDone=$noiseskyPlaneDir/done.txt
#     computeSky $planeSubtractionForCoaddDir $noiseskyPlaneDir $noiseskyPlaneDone true $sky_estimation_method -1 false $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth
    
#     minRmsFileName=min_rms_plane_it$iteration.txt
#     python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $h $noiseskyPlaneDir $DIR $iteration $minRmsFileName

#     wdir=$BDIR/weight-dir_plane_it$iteration
#     wonlydir=$BDIR/only-w-dir_plane_it$iteration
#     wdone=$wdir/done.txt
#     wonlydone=$wonlydir/done.txt
#     if ! [ -d $wdir ]; then mkdir $wdir; fi
#     if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
#     # We provide the fullGrid because we are going to combine then now
#     computeWeights $wdir $wdone $wonlydir $wonlydone $planeSubtractionForCoaddDir $noiseskyPlaneDir $iteration $minRmsFileName

#     sigmaForStdSigclip=2
#     clippingdir=$BDIR/clipping-outliers_plane_it$iteration
#     clippingdone=$clippingdir/done.txt
#     buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip

#     mowdir=$BDIR/weight-dir-no-outliers_plane_it$iteration
#     moonwdir=$BDIR/only-weight-dir-no-outliers_plane_it$iteration
#     mowdone=$mowdir/done.txt
#     if ! [ -d $mowdir ]; then mkdir $mowdir; fi
#     if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi
#     removeOutliersFromWeightedFrames $mowdone $totalNumberOfFrames $mowdir $moonwdir $clippingdir $wdir $wonlydir


    # echo -e "\n·Removing bad frames"

    # diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
    # rejectedFramesDir=$BDIR/rejectedFrames
    # if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
    # echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"

    # rejectedByAstrometry=identifiedBadFrames_astrometry.txt
    # removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry
    # removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry

    # rejectedByBackgroundStd=identifiedBadFrames_backgroundStd.txt
    # removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd
    # removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd

    # rejectedByBackgroundValue=identifiedBadFrames_backgroundValue.txt
    # removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue
    # removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue

    # rejectedByBackgroundFWHM=identifiedBadFrames_fwhm.txt
    # removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM
    # removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM


#     coaddDir=$BDIR/coadds
#     coaddDone=$coaddDir/done_plane.txt
#     coaddName=$coaddDir/"$objectName"_coadd_"$filter"_plane.fits
#     buildCoadd $coaddDir $coaddName $mowdir $moonwdir $coaddDone

    # #Compute surface brightness limit
    # sblimitFile=$coaddDir/"$objectName"_"$filter"_sblimit.txt
    # exposuremapName=$coaddDir/exposureMap.fits
    # if [ -f  $sblimitFile ]; then
    #     echo -e "\n\tSurface brightness limit for coadd already measured\n"
    # else
    #     surfaceBrightnessLimit=$( limitingSurfaceBrightness $coaddName $maskName $exposuremapName $coaddDir $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )
    # fi

    # keyWords=("FRAMES_COMBINED" \
    #           "NUMBER_OF_DIFFERENT_NIGHTS" \
    #           "INITIAL_DATE_OBS"
    #           "MEAN_DATE_OBS"
    #           "FINAL_DATE_OBS"
    #           "FILTER" \
    #           "SATURATION_THRESHOLD" \
    #           "CALIBRATION_BRIGHTLIMIT" \
    #           "CALIBRATION_FAINTLIMIT" \
    #           "RUNNING_FLAT" \
    #           "WINDOW_SIZE" \
    #           "STD_FOR_BAD_FRAMES" \
    #           "SURFACE_BRIGHTNESS_LIMIT")
    # numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)

    # values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames" "$surfaceBrightnessLimit")
    # comments=("" "" "" "" "" "" "" "" "" "" "" "Num. of tandard deviations used for rejection" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")

    # astfits $coaddName --write=/,"Pipeline information"
    # addkeywords $coaddName keyWords values comments 
#   fi
# fi



# # Remove intermediate folders to save some space
find $BDIR/noise-sky_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/noise-sky-after-photometry_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/noise-sky_forPlaneCoadd_it1 -type f ! -name 'done*' -exec rm {} \;

find $BDIR/sub-sky-fullGrid_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/sub-sky-smallGrid_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrFullGrid-dir_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrSmallGrid-dir_it1 -type f ! -name 'done*' -exec rm {} \;

find $BDIR/my-catalog-halfmaxradius_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/match-decals-myData_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/ourData-catalogs-apertures_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/framesWithCoaddSubtracted -type f ! -name 'done*' -exec rm {} \;

find $BDIR/weight-dir-no-outliers -type f ! -name 'done*' -exec rm {} \;
find $BDIR/weight-dir-no-outliers_plane_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/weight-dir -type f ! -name 'done*' -exec rm {} \;
find $BDIR/weight-dir_plane_it1 -type f ! -name 'done*' -exec rm {} \;

find $BDIR/only-weight-dir-no-outliers -type f ! -name 'done*' -exec rm {} \;
find $BDIR/only-weight-dir-no-outliers_plane_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/only-w-dir -type f ! -name 'done*' -exec rm {} \;
find $BDIR/only-w-dir_plane_it1 -type f ! -name 'done*' -exec rm {} \;


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
computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT noisechisel $polynomialDegree $imagesAreMasked $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT
subtractSky $entiredir_fullGrid $subskyFullGrid_dir $subskyFullGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT


imagesForCalibration=$subskySmallGrid_dir
alphatruedir=$BDIR/alpha-stars-true_it$iteration
matchdir=$BDIR/match-decals-myData_it$iteration
#tileSize=50
computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $rangeUsedCalibrationDir $mosaicDir  \
                          $alphatruedir $calibrationBrightLimit $calibrationFaintLimit $tileSize $apertureUnits $numberOfApertureUnitsForCalibration



photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir
applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir


# We mask again the points in order to measure (after photometric calibration) the sky accurately
smallPointings_photCorr_maskedDir=$BDIR/photCorrSmallGrid_masked_it$iteration
maskedPointingsDone=$smallPointings_photCorr_maskedDir/done_.txt
maskPointings $photCorrSmallGridDir $smallPointings_photCorr_maskedDir $maskedPointingsDone $maskName $entiredir_fullGrid


noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done.txt
# Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
computeSky $smallPointings_photCorr_maskedDir $noiseskydir $noiseskydone true noisechisel -1 true $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth

minRmsFileName="min_rms_it$iteration.txt"
python3 $pythonScriptsPath/find_rms_min.py "$filter" 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration $minRmsFileName

wdir=$BDIR/weight-dir_it$iteration
wdone=$wdir/done_"$k"_ccd"$h".txt
if ! [ -d $wdir ]; then mkdir $wdir; fi
wonlydir=$BDIR/only-w-dir_it$iteration
wonlydone=$wonlydir/done_"$k"_ccd"$h".txt
if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi

# We provide the fullGrid because we are going to combine then now
computeWeights $wdir $wdone $wonlydir $wonlydone $photCorrFullGridDir $noiseskydir $iteration $minRmsFileName

clippingdir=$BDIR/clipping-outliers_it$iteration
clippingdone=$clippingdir/done_"$k".txt
buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip

mowdir=$BDIR/weight-dir-no-outliers_it$iteration
moonwdir=$BDIR/only-weight-dir-no-outliers_it$iteration
if ! [ -d $mowdir ]; then mkdir $mowdir; fi
if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi
mowdone=$mowdir/done.txt
removeOutliersFromWeightedFrames $mowdone $totalNumberOfFrames $mowdir $moonwdir $clippingdir $wdir $wonlydir


echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
echo -e "\nBuilding coadd"

echo -e "\n·Removing bad frames"

diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
rejectedFramesDir=$BDIR/rejectedFrames_it$iteration
if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"

rejectedByAstrometry=identifiedBadFrames_astrometry.txt
removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry
removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry

rejectedByBackgroundStd=identifiedBadFrames_backgroundStd.txt
removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd
removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd

rejectedByBackgroundValue=identifiedBadFrames_backgroundValue.txt
#removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue
#removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue

rejectedByBackgroundFWHM=identifiedBadFrames_fwhm.txt
#removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM
#removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM


coaddDir=$BDIR/coadds_it$iteration 
coaddDone=$coaddDir/done.txt
coaddName=$coaddDir/"$objectName"_coadd_"$filter".fits
buildCoadd $coaddDir $coaddName $mowdir $moonwdir $coaddDone

maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
noisechisel_param="--tilesize=30,30 --detgrowmaxholesize=5000 --dthresh=0.1 --snminarea=2 --rawoutput"
export noisechisel_param
if [ -f $maskName ]; then
  echo "The mask of the weighted coadd is already done"
else
  astnoisechisel $coaddName $noisechisel_param --numthreads=$num_cpus -o $maskName
fi

exposuremapDir=$coaddDir/"$objectName"_exposureMap
exposuremapdone=$coaddDir/done_exposureMap.txt
computeExposureMap $framesDir $exposuremapDir $exposuremapdone

sblimitFile=$coaddDir/"$objectName"_"$filter"_sblimit.txt
exposuremapName=$coaddDir/exposureMap.fits
if [ -f  $sblimitFile ]; then
    echo -e "\n\tSurface brightness limit for coadd already measured\n"
else
    surfaceBrightnessLimit=$( limitingSurfaceBrightness $coaddName $maskName $exposuremapName $coaddDir $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )
fi

echo -e "\nAdding keywords to the coadd"
keyWords=("FRAMES_COMBINED" \
          "NUMBER_OF_DIFFERENT_NIGHTS" \
          "INITIAL_DATE_OBS" \
          "MEAN_DATE_OBS" \
          "FINAL_DATE_OBS" \
          "FILTER" \
          "VIGNETTING_THRESHOLD" \
          "SATURATION_THRESHOLD" \
          "CALIBRATED_USING" \
          "CALIBRATION_BRIGHTLIMIT" \
          "CALIBRATION_FAINTLIMIT" \
          "RUNNING_FLAT" \
          "WINDOW_SIZE" \
          "STD_FOR_BAD_FRAMES" \
          "SURFACE_BRIGHTNESS_LIMIT")

numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$vignettingThreshold" "$saturationThreshold" "$surveyForPhotometry" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames" "$surfaceBrightnessLimit")
comments=("" "" "" "" "" "" "" "" "" "" "" "Num. of tandard deviations used for rejection" "" "" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")
astfits $coaddName --write=/,"Pipeline information"
addkeywords $coaddName keyWords values comments

halfMaxRadForCoaddName=$halfMaxRadiusVsMagnitudeOurDataDir/coadd_it2.png
if [ -f $halfMaxRadForCoaddName ]; then
  echo -e "\tThe Half-Max-Rad vs Magnitude has been already generate for the coadd"
else
  produceHalfMaxRadVsMagForSingleImage $coaddName $halfMaxRadiusVsMagnitudeOurDataDir $catdir/"$objectName"_Gaia_eDR3.fits $toleranceForMatching $pythonScriptsPath "coadd_it2" 100
fi

framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted_it$iteration
framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi

if [ -f $framesWithCoaddSubtractedDone ]; then
    echo -e "\nFrames with coadd subtracted already generated\n"
else
  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_it$iteration.fits
  subtractCoaddToFrames $photCorrFullGridDir $coaddName $framesWithCoaddSubtractedDir
  astarithmetic $(ls -v $framesWithCoaddSubtractedDir/*.fits) $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l) sum -g1 -o$sumMosaicAfterCoaddSubtraction
  echo "done" > $framesWithCoaddSubtractedDone
fi


# Subtract a plane and build the coadd. Thus we have the constant background coadd and the plane background coadd
# if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
#   coaddPlaneDone=$coaddDir/done_plane.txt
#   if [ -f $coaddPlaneDone ]; then
#     echo -e "\n\tCoadd with background modelled as a plane already done"
#   else 
#     planeEstimationForCoaddDir=$BDIR/planeEstimationBeforeCoadd_it$iteration
#     planeEstimationForCoaddDone=$planeEstimationForCoaddDir/done.txt
#     polynomialDegree=1
#     if ! [ -d $planeEstimationForCoaddDir ]; then mkdir $planeEstimationForCoaddDir; fi
#     computeSky $photCorrSmallGridDir $planeEstimationForCoaddDir $planeEstimationForCoaddDone false $sky_estimation_method $polynomialDegree false $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth

#     planeSubtractionForCoaddDir=$BDIR/planeSubtractionBeforeCoadd_it$iteration
#     if ! [ -d $planeSubtractionForCoaddDir ]; then mkdir $planeSubtractionForCoaddDir; fi
#     planeSubtractionForCoaddDone=$planeSubtractionForCoaddDir/done.txt
#     subtractSky $photCorrFullGridDir $planeSubtractionForCoaddDir $planeSubtractionForCoaddDone $planeEstimationForCoaddDir false

#     noiseskyPlaneDir=$BDIR/noise-sky_forPlaneCoadd_it$iteration
#     if ! [ -d $noiseskyPlaneDir ]; then mkdir $noiseskyPlaneDir; fi
#     noiseskyPlaneDone=$noiseskyPlaneDir/done.txt
#     computeSky $planeSubtractionForCoaddDir $noiseskyPlaneDir $noiseskyPlaneDone true $sky_estimation_method -1 false $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth
    
#     minRmsFileName=min_rms_plane_it$iteration.txt
#     python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $h $noiseskyPlaneDir $DIR $iteration $minRmsFileName

#     wdir=$BDIR/weight-dir_plane_it$iteration
#     wonlydir=$BDIR/only-w-dir_plane_it$iteration
#     wdone=$wdir/done.txt
#     wonlydone=$wonlydir/done.txt
#     if ! [ -d $wdir ]; then mkdir $wdir; fi
#     if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
#     # We provide the fullGrid because we are going to combine then now
#     computeWeights $wdir $wdone $wonlydir $wonlydone $planeSubtractionForCoaddDir $noiseskyPlaneDir $iteration $minRmsFileName

#     sigmaForStdSigclip=2
#     clippingdir=$BDIR/clipping-outliers_plane_it$iteration
#     clippingdone=$clippingdir/done.txt
#     buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip

#     mowdir=$BDIR/weight-dir-no-outliers_plane_it$iteration
#     moonwdir=$BDIR/only-weight-dir-no-outliers_plane_it$iteration
#     mowdone=$mowdir/done.txt
#     if ! [ -d $mowdir ]; then mkdir $mowdir; fi
#     if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi
#     removeOutliersFromWeightedFrames $mowdone $totalNumberOfFrames $mowdir $moonwdir $clippingdir $wdir $wonlydir

  # echo -e "\n·Removing bad frames"

  # diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
  # rejectedFramesDir=$BDIR/rejectedFrames
  # if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
  # echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"

  # rejectedByAstrometry=identifiedBadFrames_astrometry.txt
  # removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry
  # removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry

  # rejectedByBackgroundStd=identifiedBadFrames_backgroundStd.txt
  # removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd
  # removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundStd

  # rejectedByBackgroundValue=identifiedBadFrames_backgroundValue.txt
  # removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue
  # removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue

  # rejectedByBackgroundFWHM=identifiedBadFrames_fwhm.txt
  # removeBadFramesFromReduction $mowdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM
  # removeBadFramesFromReduction $moonwdir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM

#     coaddDir=$BDIR/coadds_it$iteration
#     coaddDone=$coaddDir/done_plane.txt
#     coaddName=$coaddDir/"$objectName"_coadd_"$filter"_plane.fits
#     buildCoadd $coaddDir $coaddName $mowdir $moonwdir $coaddDone

    # #Compute surface brightness limit
    # sblimitFile=$coaddDir/"$objectName"_"$filter"_sblimit.txt
    # exposuremapName=$coaddDir/exposureMap.fits
    # if [ -f  $sblimitFile ]; then
    #     echo -e "\n\tSurface brightness limit for coadd already measured\n"
    # else
    #     surfaceBrightnessLimit=$( limitingSurfaceBrightness $coaddName $maskName $exposuremapName $coaddDir $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )
    # fi
    
      # keyWords=("FRAMES_COMBINED" \
      #           "NUMBER_OF_DIFFERENT_NIGHTS" \
      #           "INITIAL_DATE_OBS"
      #           "MEAN_DATE_OBS"
      #           "FINAL_DATE_OBS"
      #           "FILTER" \
      #           "SATURATION_THRESHOLD" \
      #           "CALIBRATION_BRIGHTLIMIT" \
      #           "CALIBRATION_FAINTLIMIT" \
      #           "RUNNING_FLAT" \
      #           "WINDOW_SIZE" \
      #           "STD_FOR_BAD_FRAMES" \
      #           "SURFACE_BRIGHTNESS_LIMIT")
      #     numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
    # values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames" "$surfaceBrightnessLimit")
    # comments=("" "" "" "" "" "" "" "" "" "" "" "Num. of tandard deviations used for rejection" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")
    # astfits $coaddName --write=/,"Pipeline information"
    # addkeywords $coaddName keyWords values comments 
#   fi
# fi

endTime=$(date +%D%T)
echo "Pipeline ended at : ${endTime}"
exit 0



# # The following code is for doing a third iteration. But from my experience 2 is enough
# Since we are not using it and we are still developing the pipeline this is not updated.
# Once the pipeline is "finished", if we want to add a third iteration then we update it.

# maskName=$coaddir/"$objectName"_coadd_"$filter"_mask.fits
# if [ -f $maskName ]; then
#   echo "The mask of the weighted coadd is already done"
# else
#   astnoisechisel $coaddName "'$noisechisel_param'" -o $maskName
# fi

# # Remove intermediate folders to save some space
# find $BDIR/noise-sky_it2 -type f ! -name 'done*' -exec rm {} \;
# find $BDIR/sub-sky-fullGrid_it2 -type f ! -name 'done*' -exec rm {} \;
# find $BDIR/sub-sky-smallGrid_it2 -type f ! -name 'done*' -exec rm {} \;
# find $BDIR/photCorrFullGrid-dir_it2 -type f ! -name 'done*' -exec rm {} \;
# find $BDIR/photCorrSmallGrid-dir_it2 -type f ! -name 'done*' -exec rm {} \;

# find $BDIR/my-catalog-halfmaxradius_it2 -type f ! -name 'done*' -exec rm {} \;
# find $BDIR/match-decals-myData_it2 -type f ! -name 'done*' -exec rm {} \;
# find $BDIR/ourData-catalogs-apertures_it2 -type f ! -name 'done*' -exec rm {} \;
# find $BDIR/framesWithCoaddSubtractedDir_it2 -type f ! -name 'done*' -exec rm {} \;

# find $BDIR/pointings_smallGrid_masked_it2 -type f ! -name 'done*' -exec rm {} \;
# find $BDIR/photCorrSmallGrid_masked_it2 -type f ! -name 'done*' -exec rm {} \;

# find $wdir -type f ! -name 'done*' -exec rm {} \;
# find $wonlydir -type f ! -name 'done*' -exec rm {} \;
# find $mowdir -type f ! -name 'done*' -exec rm {} \;
# find $moonwdir -type f ! -name 'done*' -exec rm {} \;



# ####### ITERATION 3 ######

# iteration=3
# entiredir_smallGrid=$BDIR/pointings_smallGrid
# entiredir_fullGrid=$BDIR/pointings_fullGrid

# # We mask the pointings in order to measure (before photometric calibration) the sky accurately
# entiredir_fullGrid=$BDIR/pointings_fullGrid

# smallPointings_maskedDir=$BDIR/pointings_smallGrid_masked_it$iteration
# maskedPointingsDone=$smallPointings_maskedDir/done_.txt
# maskPointings $entiredir_smallGrid $smallPointings_maskedDir $maskedPointingsDone $maskName $entiredir_fullGrid

# noiseskydir=$BDIR/noise-sky_it$iteration
# noiseskydone=$noiseskydir/done_"$filter"_ccd"$h".txt

# subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it$iteration
# subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter"_ccd"$h".txt

# subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it$iteration
# subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter"_ccd"$h".txt

# # compute sky with frames masked with global mask
# imagesAreMasked=true
# computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth
# subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT
# subtractSky $entiredir_fullGrid $subskyFullGrid_dir $subskyFullGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT

# imagesForCalibration=$subskySmallGrid_dir
# alphatruedir=$BDIR/alpha-stars-true_it$iteration
# computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $rangeUsedCalibrationDir $mosaicDir  \
#                           $alphatruedir $calibrationBrightLimit $calibrationFaintLimit $tileSize $numberOfFWHMForPhotometry 

# photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
# photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration

# applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir
# applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir


# # We mask again the points in order to measure (after photometric calibration) the sky accurately
# smallPointings_photCorr_maskedDir=$BDIR/photCorrSmallGrid_masked_it$iteration
# maskedPointingsDone=$smallPointings_photCorr_maskedDir/done_.txt
# maskPointings $photCorrSmallGridDir $smallPointings_photCorr_maskedDir $maskedPointingsDone $maskName $entiredir_fullGrid

# noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
# noiseskydone=$noiseskydir/done_"$k"_ccd"$h".txt
# # Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
# computeSky $smallPointings_photCorr_maskedDir $noiseskydir $noiseskydone true $sky_estimation_method -1 true $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth

# minRmsFileName="min_rms_it$iteration.txt"
# echo python3 $pythonScriptsPath/find_rms_min.py "$filter" 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration $minRmsFileName

# wdir=$BDIR/weight-dir_it$iteration
# wdone=$wdir/done_"$k"_ccd"$h".txt
# if ! [ -d $wdir ]; then mkdir $wdir; fi
# wonlydir=$BDIR/only-w-dir_it$iteration
# wonlydone=$wonlydir/done_"$k"_ccd"$h".txt
# if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi

# # We provide the fullGrid because we are going to combine then now
# computeWeights $wdir $wdone $wonlydir $wonlydone $photCorrFullGridDir $noiseskydir $iteration

# clippingdir=$BDIR/clipping-outliers_it$iteration
# clippingdone=$clippingdir/done_"$k".txt
# buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $wdir $sigmaForStdSigclip



# # Fornax. Around 490 frames. Deimos, 20 cores. Around 1 h and 15 min
# mowdir=$BDIR/weight-dir-no-outliers_it$iteration
# if ! [ -d $mowdir ]; then mkdir $mowdir; fi
# # only weight
# moonwdir=$BDIR/only-weight-dir-no-outliers_it$iteration
# if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi
# mowdone=$mowdir/done_"$k"_ccd"$h".txt

# if [ -f $mowdone ]; then
#     echo -e "\nOutliers of the weighted images already masked\n"
# else
#   framesToRemoveOutliers=()
#   for a in $(seq 1 $totalNumberOfFrames); do
#     framesToRemoveOutliers+=("$a")
#   done
#   printf "%s\n" "${framesToRemoveOutliers[@]}" | parallel -j "$num_cpus" removeOutliersFromFrame {} $mowdir $moonwdir $clippingdir $wdir $wonlydir
#   echo done > $mowdone 
# fi


# echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
# echo -e "\nBuilding coadd"

# coaddDir=$BDIR/coadds_it$iteration 
# coaddDone=$coaddDir/done.txt
# coaddName=$coaddDir/"$objectName"_coadd_"$filter".fits
# buildCoadd $coaddDir $coaddName $mowdir $moonwdir $coaddDone

# echo -e "\nAdding keywords to the coadd"
# keyWords=("FRAMES_COMBINED" \
#           "NUMBER_OF_DIFFERENT_NIGHTS" \
#           "INITIAL_DATE_OBS"
#           "MEAN_DATE_OBS"
#           "FINAL_DATE_OBS"
#           "FILTER" \
#           "SATURATION_THRESHOLD" \
#           "CALIBRATION_BRIGHTLIMIT" \
#           "CALIBRATION_FAINTLIMIT" \
#           "RUNNING_FLAT" \
#           "WINDOW_SIZE" \
#           "STD_FOR_BAD_FRAMES" \
#           "SURFACE_BRIGHTNESS_LIMIT")
# numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
# values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames" "$surfaceBrightnessLimit")
# comments=("" "" "" "" "" "" "" "" "" "" "" "Num. of tandard deviations used for rejection" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")
# astfits $coaddName --write=/,"Pipeline information"
# addkeywords $coaddName keyWords values comments


# if [ -f $framesWithCoaddSubtractedDone ]; then
#     echo -e "\nFrames with coadd subtracted already generated\n"
# else
#   sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_it$iteration.fits
#   subtractCoaddToFrames $photCorrFullGridDir $coaddName $framesWithCoaddSubtractedDir
#   astarithmetic $(ls -v $framesWithCoaddSubtractedDir/*.fits) $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l) sum -g1 -o$sumMosaicAfterCoaddSubtraction
# fi
