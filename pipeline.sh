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
source "$pipelinePath/pipeline_functions.sh"

echo -e "\n ${GREEN} ---Loading Modules--- ${NOCOLOUR}"

gnuastroModuleName="gnuastro/0.22"
load_module $gnuastroModuleName

astrometryModuleName="astrometry.net/0.98"
load_module $astrometryModuleName

scampModuleName="scamp/2.14.0"
load_module $scampModuleName

sexModuleName="sextractor"
load_module $sexModuleName

swarpModuleName="swarp/2.41.5"
load_module $swarpModuleName

wcsModuleName="wcstools/3.9.7"
load_module $wcsModuleName

########## Handling options and arguments ##########
RUN_NIGHT=""
RUN_COMMON=false

OPTSTRING=":hn:cp:"
while getopts ${OPTSTRING} opt; do
  case ${opt} in
    h)
      help
      exit 0;;
    n)
      RUN_NIGHT=${OPTARG};;
    c)
      RUN_COMMON=true;;
    p)
      PHASE=${OPTARG};;
    \?)
      echo "Invalid option: -${OPTARG}"
      exit 1;;
  esac
done

writeTimeOfStepToFile "Starting pipeline" $fileForTimeStamps
echo -e "\n ${GREEN} ---Loading variables from conf file --- ${NOCOLOUR}"

shift $((OPTIND -1))
confFile=$1
loadVariablesFromFile $confFile

checkIfAllVariablesAreSet

checkIfStringVariablesHaveValidValues

#checkTransmittanceFilterAndItsUnits $telescope $surveyForPhotometry $folderWithTransmittances $filter

filterCorrectionCoeff=$( checkIfNeededFilterCorrectionIsGiven $telescope $filter $surveyForPhotometry $ROOTDIR/"$objectName"/config )
if [[ $filterCorrectionCoeff == 11 ]]; then
  echo "The filter corrections for the filter $filter, telescope $telescope and survey $survey were not found"
  echo "Exiting with error code 11"
  exit 11
fi

outputConfigurationVariablesInformation


# The following lines are responsible of the cpu's used for paralellise
if [ -n "$SLURM_CPUS_ON_NODE" ]; then
  num_cpus=$SLURM_CPUS_ON_NODE
else
  num_cpus=$defaultNumOfCPUs
fi

echo -e "\nNumber of CPUs allocated: $num_cpus"
export num_cpus

# After some testing we find that it is efficient to use noisechisel with 4 threads. So we parallelise accordingly
noisechiselNumThreads=8
export noisechiselNumThreads

if (( num_cpus < noisechiselNumThreads )); then
  num_threads=$num_cpus
  num_parallel=1
else
  num_threads=$noisechiselNumThreads
  num_parallel=$(( num_cpus / noisechiselNumThreads ))
fi

export num_parallel
export num_threads

echo -e "NoiseChisel-related steps: running ${num_parallel} frames in parallel, ${num_threads} threads each"


######## Loading and transforming to needed format the user-defined masks to apply

maskParams=$(printf "%s " "${masksToApply[@]}")
echo $maskParams
export maskParams

########## Prepare data ##########

echo -e "\n ${GREEN} ---Preparing data--- ${NOCOLOUR}"

DIR=$ROOTDIR/"$objectName"
CDIR=$DIR/config
INDIRo=$ROOTDIR/"$objectName"/DATA-or
BDIR=$ROOTDIR/"$objectName"/build
INDIR=$ROOTDIR/"$objectName"/DATA
DARKDIR=$ROOTDIR/"$objectName"/dark
BIASDIR=$ROOTDIR/"$objectName"/bias
keyWordDirectory=$ROOTDIR/"$objectName"/keywords

export DIR
export INDIRo
export BDIR
export INDIR
export CDIR
export DARKDIR
export keyWordDirectory
source $ROOTDIR/env_lapalma/bin/activate
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
  local currentNight=$1
  local framesForCommonReductionDone=$framesForCommonReductionDir/done_"$filter"_ccd"$h"_n"$currentNight".txt

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
                unixTimeInSeconds=$(TZ=UTC  gdate -d "$DATEOBS" +"%s")
              else
                unixTimeInSeconds=$(TZ=UTC  date -d "$DATEOBS" +"%s")
              fi
            fi
            out=$currentINDIR/$unixTimeInSeconds.fits

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
  #rm $currentINDIRo/*.fits
  
  # -------------------------------------------------------
  # Number of exposures of the current night
  local n_exp=$(ls -v $currentINDIR/*.fits | wc -l)
  echo -e "Number of exposures ${ORANGE} ${n_exp} ${NOCOLOUR}"
  if [ -d $DARKDIR/night"$currentNight" ]; then
    currentDARKDIR=$DARKDIR/night$currentNight
  else
    currentDARKDIR=$DARKDIR
  fi
  mdadir=$BDIR/masterdark_n$currentNight
  
  ########## Creating master bias ##########
  for h in 0; do
    echo -e "\n ${GREEN} Creating master bias/dark-bias ${NOCOLOUR}"
    mdadone=$mdadir/mdark_"$filter"_ccd"$h".txt
 
    if ! [ -d $mdadir ]; then mkdir $mdadir; fi
    if [ -f $mdadone ]; then
      echo -e "\nMasterdark is already done for night $currentNight and extension $h\n"
    else
      if [ $(ls -v $currentDARKDIR/*.fits | wc -l) -eq 1 ]; then
        astfits $currentDARKDIR/*.fits --copy=0 -o $mdadir/mdark_"$filter"_n"$currentNight"_ccd$h.fits
        echo done > $mdadone
        continue
      fi
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
 		    air=$(python3 $pythonScriptsPath/get_airmass_teo.py $i $dateHeaderKey $ra_gal $dec_gal $telescopeLat $telescopeLong $telescopeElevation)
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
    framesToSubtract=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      framesToSubtract+=("$base")
    done
    dark=$mdadir/mdark_"$filter"_n"$currentNight"_ccd$h.fits 
    printf "%s\n" "${framesToSubtract[@]}" | parallel -j "$num_cpus" subtractBiasFromFrame {} $dark $saturationThreshold $currentINDIR $mbiascorrdir
    echo done > $mbiascorrdone
  fi
  rm -f $currentINDIR/*.fits
  #Since bias corrected is the only one we're not to delete, we re-compute n_exp with it
  local n_exp=$(ls $mbiascorrdir/*.fits | wc -l)
  echo -e "${ORANGE} ------ FLATS ------ ${NOCOLOUR}\n"
  echo -e "${GREEN} --- Flat iteration 1 --- ${NOCOLOUR}"

  ########## Creating the ring mask ##########
  # We always need the common ring  definition always stored for photometric calibration (selection of decals bricks to download)
  ringdir=$BDIR/ring
  mkdir -p "$ringdir"

  # NFS-Safe Lock variables
  LOCK_DIR="$ringdir/nfs_build_lock"
  DONE_FILE="$ringdir/nfs_build_done.txt"

  echo -e "\n ${GREEN} Checking/Creating normalisation rings ${NOCOLOUR}"

  # 1. Check if another night has ALREADY built the ring(s)
  if [ -f "$DONE_FILE" ]; then
      echo -e "\tRing(s) already built by another night. Reusing them for night $currentNight."
  else
      # 2. Try to acquire the lock to be the builder
      if mkdir "$LOCK_DIR" 2>/dev/null; then
          echo -e "\tNode for night $currentNight acquired the lock. Building the ring(s)..."
          
          if [[ "$USE_COMMON_RING" = true ]]; then
              cp $commonRingDefinitionFile $ringdir/ring.txt 
              astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring.fits $commonRingDefinitionFile

              if [ "$telescope" == "OSIRIS+" ]; then
                  python3 $pythonScriptsPath/cutRing.py -1 300 -1 -1
              fi
          else
              # Multiple rings scenario
              astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_2.fits $secondRingDefinitionFile
              astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_1.fits $firstRingDefinitionFile
          fi
          
          # Signal to all other nodes that the rings are fully written and ready
          touch "$DONE_FILE"
          rmdir "$LOCK_DIR"
      else
          # 3. Another node is currently building the ring. Wait for it to finish.
          echo -e "\tAnother night is currently building the ring(s). Waiting..."
          while [ ! -f "$DONE_FILE" ]; do
              sleep 2
          done
          echo -e "\tRing(s) finished by the other node. Proceeding with night $currentNight."
      fi
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
    normaliseImagesWithRing $mbiascorrdir $normit1dir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $n_exp 
    echo done > $normit1done
  fi

  # Then, if the running flat is configured to be used, we combine the normalised images with a sigma clipping median
  # using the running flat strategy
  if [[ "${RUNNING_FLAT,,}" == "true" ]]; then
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
    calculateWholeNightFlat $flatit1WholeNightdir/flat-it1_wholeNight_n$currentNight.fits $normit1dir $currentNight $flatit1WholeNightdir
    echo "done" >> $flatit1WholeNightdone
  fi

  rm -f $normit1dir/*.fits

  # Dividing the science images for the running it1 flat
  if [[ "${RUNNING_FLAT,,}" == "true" ]]; then
    flatit1imadir=$BDIR/flat-it1-Running-ima_n$currentNight
    flatit1imadone=$flatit1imadir/done_"$filter"_ccd"$h".txt
    if ! [ -d $flatit1imadir ]; then mkdir $flatit1imadir; fi
    if [ -f $flatit1imadone ]; then
      echo -e "\nScience images are divided by flat it1 for night $currentNight and extension $h\n"
    else
      divideImagesByRunningFlats $mbiascorrdir $flatit1imadir $flatit1dir $flatit1imadone $n_exp 1
    fi
  fi
  rm -f $flatit1dir/*.fits
  # Dividing the science images for the whole night it1 flat
  flatit1WholeNightimaDir=$BDIR/flat-it1-WholeNight-ima_n$currentNight
  flatit1WholeNightimaDone=$flatit1WholeNightimaDir/done_"$filter"_ccd"$h".txt
  if ! [ -d $flatit1WholeNightimaDir ]; then mkdir $flatit1WholeNightimaDir; fi
  if [ -f $flatit1WholeNightimaDone ]; then
    echo -e "\nScience images are divided by whole night flat it1 for night $currentNight and extension $h\n"
  else
    wholeNightFlatToUse=$flatit1WholeNightdir/flat-it1_wholeNight_n$currentNight.fits
    divideImagesByWholeNightFlat $mbiascorrdir $flatit1WholeNightimaDir $wholeNightFlatToUse $flatit1WholeNightimaDone $n_exp
  fi
  rm -f $flatit1WholeNightdir/*.fits



  # Iteration 2 -----

  maskAndNormaliseForFlatIteration 2 "$flatit1imadir" "$flatit1WholeNightimaDir" "$n_exp" "$currentNight"

  flatit2dir=$BDIR/flat-it2-Running_n$currentNight
  flatsForIteration 2 "$flatit2dir" "$currentNight"


  # Dividing the science image by the it2 flat
  if [[ "${RUNNING_FLAT,,}" == "true" ]]; then
    flatit2imadir=$BDIR/flat-it2-Running-ima_n$currentNight
    flatit2imadone=$flatit2imadir/done_"$filter"_ccd"$h".txt
    if ! [ -d $flatit2imadir ]; then mkdir $flatit2imadir; fi
    if [ -f $flatit2imadone ]; then
      echo -e "\nRunning flats it2-2 already built for night $currentNight and extension $h\n"
    else
      divideImagesByRunningFlats $mbiascorrdir $flatit2imadir $flatit2dir $flatit2imadone $n_exp 2
    fi
  fi
  rm -f $flatit2dir/*.fits

  # Dividing the science images for the whole night it2 flat
  flatit2WholeNightimaDir=$BDIR/flat-it2-WholeNight-ima_n$currentNight
  flatit2WholeNightimaDone=$flatit2WholeNightimaDir/done_"$filter"_ccd"$h".txt
  if ! [ -d $flatit2WholeNightimaDir ]; then mkdir $flatit2WholeNightimaDir; fi
  if [ -f $flatit2WholeNightimaDone ]; then
    echo -e "\nScience images are divided by whole night flat it2 for night $currentNight and extension $h\n"
  else
    wholeNightFlatToUse=$flatit2WholeNightdir/flat-it2_wholeNight_n$currentNight.fits
    divideImagesByWholeNightFlat $mbiascorrdir $flatit2WholeNightimaDir $wholeNightFlatToUse $flatit2WholeNightimaDone $n_exp
  fi
  rm -f $flatit2WholeNightdir/*.fits



  ########## Creating the it3 master flat image ##########
  echo -e "${GREEN} --- Flat iteration 3 --- ${NOCOLOUR}"

  maskAndNormaliseForFlatIteration 3 "$flatit2imadir" "$flatit2WholeNightimaDir" "$n_exp" "$currentNight"
 
  flatit3BeforeCorrectiondir=$BDIR/flat-it3-Running-BeforeCorrection_n$currentNight
  flatsForIteration 3 "$flatit3BeforeCorrectiondir" "$currentNight"


  # Correct the running flats using the whole night flat
  flatit3dir=$BDIR/flat-it3-Running_n$currentNight
  if [[ "${RUNNING_FLAT,,}" == "true" ]]; then
    flatit3done=$flatit3dir/done_"$filter"_ccd"$h".txt
    if ! [ -d $flatit3dir ]; then mkdir $flatit3dir; fi
    if [ -f $flatit3done ]; then
      echo -e "\nFlats iteration 3 are corrected using the flat of the whole night for night $currentNight and extension $h\n"
    else
      imagesToCorrect=()
      for i in $flatit3BeforeCorrectiondir/*.fits; do
        imagesToCorrect+=("$(basename $i)")
      done
      printf "%s\n" "${imagesToCorrect[@]}" | parallel -j "$num_cpus" correctRunningFlatWithWholeNightFlat {} $flatit3BeforeCorrectiondir $flatit3WholeNightdir/flat-it3_wholeNight_n$currentNight.fits $flatit3dir $dateHeaderKey
      echo done > $flatit3done
    fi
  fi
  rm -f $flatit3BeforeCorrectiondir/*.fits

  # Dividing the science image by the it3 flat
  # If running flat selected, we use it to produce the final flatted images
  # If not selcted, we applyt the whole night flat
  flatit3imadir=$BDIR/flat-it3-ima_n$currentNight
  flatit3imadone=$flatit3imadir/done_"$filter"_ccd"$h".txt
  if ! [ -d $flatit3imadir ]; then mkdir $flatit3imadir; fi
  if [[ "${RUNNING_FLAT,,}" == "true" ]]; then
    if [ -f $flatit3imadone ]; then
      echo -e "\nScience images are divided by the it3 flat for night $currentNight and extension $h\n"
    else
      divideImagesByRunningFlats $mbiascorrdir $flatit3imadir $flatit3dir $flatit3imadone $n_exp 3
    fi
  else
      wholeNightFlatToUse=$flatit3WholeNightdir/flat-it3_wholeNight_n$currentNight.fits
      divideImagesByWholeNightFlat $mbiascorrdir $flatit3imadir $wholeNightFlatToUse $flatit3imadone $n_exp
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
    imagesForVignetting=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      imagesForVignetting+=("$base")
    done
    printf "%s\n" "${imagesForVignetting[@]}" | parallel -j "$num_cpus" maskVignettingOnImages {} $flatit3imadir $maskedcornerdir $flatit3dir $flatit3WholeNightdir $RUNNING_FLAT $n_exp $currentNight $lowerVignettingThreshold $upperVignettingThreshold 
    echo done > $maskedcornerdone
  fi

  rm  $flatit3dir/*.fits  
  
  # At this point we can process the frames of all the nights in the same way
  # So we place all the final frames into a common folder.
  if [ -f $framesForCommonReductionDone ]; then
    echo -e "\nFrames already placed in the folder for frames prepared to common reduction"
  else
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      cp $maskedcornerdir/$base $framesForCommonReductionDir/$base
    done
    echo "done" > $framesForCommonReductionDone
    rm $maskedcornerdir/*.fits 
  fi

  # # Removing intermediate information to save space - We maintain the final flats for checking them
  # rm -rf $BDIR/masked-corner_n$currentNight
  rm -rf $BDIR/bias-corrected_n$currentNight
  rm -rf $BDIR/masterdark_n$currentNight
  rm -rf $BDIR/flat-it3-Running-BeforeCorrection_n$currentNight
  rm -rf $BDIR/flat-it3-ima_n$currentNight

  for a in $(seq 1 3); do
    rm -rf $BDIR/flat-it"$a"-Running_n$currentNight
    rm -rf $BDIR/flat-it"$a"-WholeNight_n$currentNight
    rm -rf $BDIR/noise-it"$a"-Running_n$currentNight
    rm -rf $BDIR/noise-it"$a"-WholeNight_n$currentNight
    rm -rf $BDIR/flat-it"$a"-Running-ima_n$currentNight
    rm -rf $BDIR/flat-it"$a"-WholeNight-ima_n$currentNight
    rm -rf $BDIR/masked-it"$a"-Running_n$currentNight
    rm -rf $BDIR/masked-it"$a"-WholeNight_n$currentNight
    rm -rf $BDIR/norm-it"$a"-images_n$currentNight
    rm -rf $BDIR/norm-it"$a"-Running-images_n$currentNight
    rm -rf $BDIR/norm-it"$a"-WholeNight-images_n$currentNight
  done

}
export -f oneNightPreProcessing

writeTimeOfStepToFile "Process the individual nights" $fileForTimeStamps

if [ -n "$RUN_NIGHT" ]; then
  echo -e "\n${GREEN} --- Processing night $RUN_NIGHT --- ${NOCOLOUR}\n"
  
  echo "Night $RUN_NIGHT is starting work at $(date)"
  oneNightPreProcessing "$RUN_NIGHT"
  echo "Night $RUN_NIGHT actually finished all its work naturally at $(date)"
  exit 0
fi


if [ -n "$PHASE" ]; then
  case "$PHASE" in
    astrometry-setup)
      renameFiles
      runAstrometrySetup
      exit 0
      ;;
    astrometry-solve)
      runAstrometryPhase
      exit 0
      ;;
    sextractor-scamp)
      runSextractorScampPhase
      exit 0
      ;;
    warp)
      runWarpPhase
      exit 0
      ;;
    sky)
      runMaskAndSkyPhase 1
      exit 0
      ;;
    sky-diagnostics)
      runSkyDiagnosticsPhase
      exit 0
      ;;
    prepare-calibration-data)
      runPrepareCalibrationDataPhase
      exit 0
      ;;
    compute-calibration-factors)
      runComputeCalibrationFactorsPhase 1
      exit 0
      ;;
    calibration-factors-diagnostics)
      runCalibrationFactorsDiagnosticsPhase 1
      exit 0
      ;;
    apply-calibration-and-fwhm)
      runApplyCalibrationFactorsPhase 1
      runComputeFWHMPhase
      exit 0
      ;;
    background-astrometry-calibPlot-diagnostics)
      runFwhmCheckPhase 1
      runBackgroundDiagnosisPhase 1
      runAstrometryDiagnosisPhase
      runCalibrationDiagnosisPlotPhase 1
      exit 0
      ;;
    small-to-fullGrid)
      runSmallGridToFullGridPhase 1
      exit 0
      ;;
    remove-bad-frames-and-create-blocks)
      runRemoveBadFramesAndCreateBlocksPhase 1
      exit 0
      ;;
    build-coadd-blocks)
      runBuildCoaddBlocksPhase 1
      exit 0
      ;;
    stitch-coadds)
      runStitchCoaddsPhase 1
      exit 0
      ;;
    coadd-mask-and-diagnostics)
      runCoaddMaskAndDiagnosticsPhase 1
      exit 0
      ;;
    create-it2-Mask)
      runCreateIt2MaskPhase
      exit 0
      ;;
    it2-mask-and-sky)
      runIt2MaskingPhase
      runMaskAndSkyPhase 2
      exit 0
      ;;
    it2-compute-calibration-factors)
      runComputeCalibrationFactorsPhase 2
      exit 0
      ;;
    it2-calibration-factors-diagnostics)
      rm -f $subskySmallGrid_dir/*.fits # This probably shouldn't be there
      runCalibrationFactorsDiagnosticsPhase 2
      runBackgroundDiagnosisPhase 2
      exit 0
      ;;
    it2-apply-calibration-factors)
      runApplyCalibrationFactorsPhase 2
      exit 0
      ;;
    it2-fwhm-and-calibPlot-diagnostics)
      runFwhmCheckPhase 2
      runCalibrationDiagnosisPlotPhase 2
      exit 0
      ;;
    it2-small-to-fullGrid)
      runSmallGridToFullGridPhase 2
      exit 0
      ;;
    it2-remove-bad-frames-and-create-blocks)
      runRemoveBadFramesAndCreateBlocksPhase 2
      exit 0
      ;;
    it2-build-coadd-blocks)
      runBuildCoaddBlocksPhase 2
      exit 0
      ;;
    it2-stitch-coadds)
      runStitchCoaddsPhase 2
      exit 0
      ;;
    it2-coadd-keywords)
      computeAndAddKeywordsToCoadd 2
      exit 0
      ;;
    *)
      echo "Unknown phase: $PHASE"
      exit 1
      ;;
  esac
fi






######## DO NOT REMOVE THIS COMMENTED LINES#########

# Star subtaction

####### ITERATION 2 ######
# iteration=2
# entiredir_smallGrid=$BDIR/pointings_smallGrid
# num_ccd=1
# export num_ccd

# if [[ "$subtractStarsFromRaw" == "true" ]]; then
#   echo -e "\n\t${GREEN} --- Subtract stars from frames --- ${NOCOLOUR} \n"
#   ###We will make the following:
#   ##  # Check where the star falls in a circle centered on RA, DEC and radius=RAFEC
#   ##  # If AZ exist in the catalog we will measure the profile in azimuth
#   ##  # This CCD will be used to compute scale factor between Rmin and Rmax, tunning the range with MAG
#   ##  # Finally, we will subtract from all the frames where the star is, and continue to the next one
#   ##  # For background range, we will measure a first background in the CCD, to select a range ±500ADU
#   ##
#   input_subStar_small=$entiredir_smallGrid
#   starsToSubtract=$BDIR/starsToSubtract.txt
#   psfFile=$CDIR/PSF_"$filter".fits
#   psfRadFile=$CDIR/RP_PSF_"$filter".fits
#   radiusToSearch=$(awk -v r="$sizeOfOurFieldDegrees" 'BEGIN { printf "%.6f", r/2 }')
#   query_param="gaia --dataset=dr3 --center=$ra_gal,$dec_gal --radius=$radiusToSearch --column=ra,dec,phot_g_mean_mag"
#   if ! [ -f $starsToSubtract ]; then
#     astquery $query_param -o$CDIR/starsToSubtract_temp.fits
#     asttable $CDIR/starsToSubtract_temp.fits --range=3,0:6.7 --sort=3 -o$starsToSubtract
#     rm $CDIR/starsToSubtract_temp.fits
#   fi
#   ###User may give a Saturation threshold for the stars in the config file. If not, we set it manually as the saturation theshold
#   if [ -z "$starSatThreshold" ]; then
#     starSatThreshod=$saturationThreshold
#   fi 
#   calFactor=$(getCommonCalibrationFactor 1)
  
#   starId=0
#   while IFS= read -r line; do
#     #We skip the lines that contain info about the columns
#     [[ $line =~ ^#.*$ ]] && continue

#     ((starId++))
#     outputDir_small=$BDIR/pointings_smallGrid_sub$starId
#   #  #outputDir_full=$BDIR/pointings_fullGrid_sub$starId
#     subtractStars $input_subStar_small "$line" $psfFile $psfRadFile $outputDir_small $starId $starSatThreshold $calFactor
    
#     if (( $(echo "$starId == 1" | bc -l) )); then exit 0; fi
#     #Sanity check: if something fail we insert an exit
#     for file in $outputDir_small/*.fits; do
#             nhdu=$( astfits $file --numhdus -q )
#             if [ $nhdu -lt $((num_ccd+1)) ]; then
#                     echo "Some frames have failed in the subtraction of star {$starId}"
#                     exit 23
#             fi
#     done

#     if ! (( $(echo "$starId == 1" | bc -l) )); then
#               rm $input_subStar_small/*.fits
#   #           #rm $input_subStar_full/*.fits

#     fi
#     input_subStar_small=$outputDir_small
#   #  #input_subStar_full=$outputDir_full
#   done < $starsToSubtract

#   # We mask the pointings in order to measure (before photometric calibration) the sky accurately
#   # MASK FROM THE NORMAL COADD
#   starsSub_small=$outputDir_small
# else
#   starsSub_small=$entiredir_smallGrid
# fi


# Pixel tagged residuals


#framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted_it$iteration
#framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
#if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi
#if [ -f $framesWithCoaddSubtractedDone ]; then
#    echo -e "\nFrames with coadd subtracted already generated\n"
#else
#  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_"$filter"_it$iteration.fits
#  photCorrfullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
#  
#  coaddAv=$coaddDir/"$objectName"_coadd_"$filter"_it"$iteration"_average.fits
#  astarithmetic $(ls -v $photCorrNoOutliersPxDir/*.fits) $(ls $photCorrNoOutliersPxDir/*.fits | wc -l) -g1 mean -o$coaddAv
#  subtractCoaddToFrames $photCorrfullGridDir $coaddAv $framesWithCoaddSubtractedDir
#  #wdir_tosub=$BDIR/weight-dir_outliers_it$iteration 
#  #wonlydir_tosub=$BDIR/only-w-dir_outliers_it$iteration
#  #wdone=$wdir_tosub/done.txt
#  #wonlydone=$wonlydir_tosub/done.txt
#  #if ! [ -d $wonlydir_tosub ]; then mkdir $wonlydir_tosub; fi
#  #if ! [ -d $wdir_tosub ]; then mkdir $wdir_tosub; fi
#  #computeWeights $wdir_tosub $wdone $wonlydir_tosub $wonlydone $framesWithCoaddSubtractedDir $noiseskydir $iteration $minRmsFileName
# #
#  #weighted=$coaddDir/sub_w.fits
#  #weightsonly=$coaddDir/sub_wx.fits
#  #astarithmetic $(ls -v $wdir_tosub/*.fits) $(ls $wdir_tosub/*.fits | wc -l) -g1 sum -o$weighted
#  #astarithmetic $(ls -v $wonlydir_tosub/*.fits) $(ls $wonlydir_tosub/*.fits | wc -l) -g1 sum -o$weightsonly
#
#  astarithmetic $(ls -v $framesWithCoaddSubtractedDir/*.fits) $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l) -g1 sum -o$sumMosaicAfterCoaddSubtraction
#  
#  #rm $wdir_tosub/*.fits
#  #rm $wonlydir_tosub/*.fits
#  diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
#  #computeMetricOfResiduals $photCorrfullGridDir $coaddName $framesWithCoaddSubtractedDir
#  #python3 $pythonScriptsPath/diagnosis_metricDistributionOfResiduals.py $framesWithCoaddSubtractedDir $diagnosis_and_badFilesDir
#
#  #sumMosaicAfterCoaddSubtractionPxTagged=$coaddDir/"$objectName"_sumMosaicAfterCoaddSubPxTagged_"$filter"_it$iteration.fits
#  ##sumMosaicAfterCoaddSubtractionAperTagged=$coaddDir/"$objectName"_sumMosaicAfterCoaddSubAperTagged_"$filter"_it$iteration.fits
#  #framesWithCoaddSubtractedTaggedDir=$BDIR/framesWithCoaddSubtractedTagged_it$iteration
#  #if ! [ -d $framesWithCoaddSubtractedTaggedDir ]; then mkdir $framesWithCoaddSubtractedTaggedDir; fi
#  #computeSumMosaicAfterCoaddSubtractionWithTracesIndicated $framesWithCoaddSubtractedDir $framesWithCoaddSubtractedTaggedDir $sumMosaicAfterCoaddSubtractionPxTagged $sumMosaicAfterCoaddSubtractionAperTagged $fwhmFolder "$noisechisel_param"
#
#  echo "done" > $framesWithCoaddSubtractedDone
#fi


