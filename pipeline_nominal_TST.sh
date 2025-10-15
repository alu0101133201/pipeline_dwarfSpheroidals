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
# THIS PIPELINE IS ONLY FOR NOMINAL RESOLUTION OF TST AFTER MAKING A 3x3 REBINNED REDUCTION
# In this pipeline, we will use the sub-products of the 3x3 rebinning reduction to make a nominal resolution reduction
# We will expect the following paths to exist:
#   · Folder with flat masks. There must exist a folder containing sub-folders divided by nights. I.e.: maskDir/night1 maskDir/night2 etc
#   · Folder with the re-binned sky values. Just entirecamera_N.txt
#   · Diagnosis and bad files folder. From 3by3, it must at least contain the identifiedBadFrames_*.txt files
#   · fileassociation_3by3.txt. Since more than a night might be expected, we need to have a file associating original frames with frames for common reduction.
# The pipeline workflow will be as follows:
#   1st step: oneNightPreProcessing(). It will only contain one iteration where masks from 3by3 are used.
#   2nd step: re-assign names of framesForCommonRecution according to fileassociation_3by3.txt
#   3rd step: astromety. We will benefit from original TST astrometry, avoiding solveField and only running scamp and swarp
#   4th step: sky subtraction. We will use 3by3 sky values to re-compute the sky. Prior to that, we divide the values by 9
#   5th step: photometry and coaddition as in the 3by3.

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

astrometryModuleName="astrometry.net/0.94"
load_module $astrometryModuleName



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

###For noisechisel runs: when we parallelize, we are parallelizing a process that already benefits from multi-threading
#Because of that, using $num_cpu for both parallel jobs and --numthreads when parallelizing noisechisel is not optimal at all
#We divide the num_cpus into 4 for num_threads and num_cpus/4 for jobs, if num_cpus>4
max_thread_per_process=4
if (( num_cpus < max_thread_per_process )); then
  num_threads=$num_cpus
  num_parallel=1
else
  num_threads=$max_thread_per_process
  num_parallel=$(( num_cpus / num_threads ))
fi
export num_parallel
export num_threads

# ****** Decision note *******
# Rebinned data
#tileSize=35
#noisechisel_param="--tilesize=$tileSize,$tileSize \
#                    --detgrowmaxholesize=5000 \
#                    --rawoutput"

# # These paremeters are oriented to TST data at original resolution. 
# astmkprof --kernel=gaussian,2,3 --oversample=1 -o$ROOTDIR/"$objectName"/kernel.fits 
 tileSize=90
 noisechisel_param="--tilesize=$tileSize,$tileSize \
                      --detgrowmaxholesize=5000 \
                      --rawoutput"
export noisechisel_param

echo -e "\n-Noisechisel parameters used for masking:"
echo -e "\t" $noisechisel_param

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

# Check if results from 3by3 are defined
nomToCheck=("flatMasksDir_w" "rebinSkyDir" "diagnosisPath" "fileAssociation3by3")
for currentVar in ${nomToCheck[@]}; do
    [[ -z ${!currentVar} ]] && echo "Variable $currentVar is not defined. Exiting with error code 8" && exit 8
done

echo -e "\n-Results from 3by3 rebinning reduction to be used:"
echo -e "\t·Flat masks directory (whole): ${ORANGE} ${flatMasksDir_w} ${NOCOLOUR}"
echo -e "\t·Re-binned sky directory: ${ORANGE} ${rebinSkyDir} ${NOCOLOUR}"
echo -e "\t·Diagnosis and bad files directory: ${ORANGE} ${diagnosisPath} ${NOCOLOUR}"
echo -e "\t·File association 3by3: ${ORANGE} ${fileAssociation3by3} ${NOCOLOUR}"
if $RUNNING_FLAT; then
  [[ -z ${flatMasksDir_r} ]] && echo "Variable flatMasksDir_r is not defined. Exiting with error code 8" && exit 8
  echo -e "\t·Flat masks directory (running): ${ORANGE} ${flatMasksDir_r} ${NOCOLOUR}"
fi
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


  # -------------------------------------------------------
  # Number of exposures of the current night
  n_exp=$(ls -v $currentINDIRo/*.fits | wc -l)
  echo -e "Number of exposures ${ORANGE} ${n_exp} ${NOCOLOUR}"
  
  if [ -d $DARKDIR/night"$currentNight" ]; then
    currentDARKDIR=$DARKDIR/night$currentNight
  else
    currentDARKDIR=$DARKDIR
  fi
  mdadir=$BDIR/masterdark_n$currentNight
  
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
  exit 0
  echo -e "${ORANGE} ------ FLATS ------ ${NOCOLOUR}\n"
  iteration=3
  # We always need the common ring  definition always stored for photometric calibration (selection of decals bricks to download)
  ringdir=$BDIR/ring
  if ! [ -d $ringdir ]; then mkdir $ringdir; fi
  # We create the .fits ring image based on how the normalisation is going to be done
  if [[ "$USE_COMMON_RING" = true && ! -f "$ringdir/ring.fits"  ]]; then
    cp $commonRingDefinitionFile $ringdir/ring.txt 
    astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring.fits $commonRingDefinitionFile
  else
    if [[ "$USE_COMMON_RING" = false ]]; then
      if [[ ! -f "$ringdir/ring_2.fits" || ! -f "$ringdir/ring_1.fits" ]]; then
        astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_2.fits $secondRingDefinitionFile
        astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_1.fits $firstRingDefinitionFile
      fi
    fi
  fi
  ##Define folders for masks for each night
  currentFlatMaskDir_r=$flatMasksDir_r/night"$currentNight"
  currentFlatMaskDir_w=$flatMasksDir_w/night"$currentNight"
  if $RUNNING_FLAT; then
    noiseit3dir=$BDIR/noise-it3-Running_n$currentNight
    noiseit3done=$noiseit3dir/done_"$filter"_ccd"$h".txt
    if ! [ -d $noiseit3dir ]; then mkdir $noiseit3dir; fi
    if [ -f $noiseit3done ]; then
      echo -e "\nMasks from 3by3 running flat are already warped for night $currentNight and extension $h\n"
    else
      frameNames=()
      for a in $(seq 1 $n_exp); do
          base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
          frameNames+=("$base")
      done
      printf "%s\n" "${frameNames[@]}" | parallel -j "$num_parallel" warpMaskForFrame {} $currentFlatMaskDir_r $noiseit3dir $num_threads
      echo done > $noiseit3done
    fi
  fi

  noiseit3WholeNightDir=$BDIR/noise-it3-WholeNight_n$currentNight
  noiseit3WholeNightdone=$noiseit3WholeNightDir/done_"$filter"_ccd"$h".txt
  if ! [ -d $noiseit3WholeNightDir ]; then mkdir $noiseit3WholeNightDir; fi
  if [ -f $noiseit3WholeNightdone ]; then
    echo -e "\nMasks from 3by3 whole night flat are already warped for night $currentNight and extension $h\n"
  else
    frameNames=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      frameNames+=("$base")
    done

    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_parallel" runNoiseChiselOnFrame {} $flatit2WholeNightimaDir $noiseit3WholeNightDir "'$noisechisel_param'"
    echo done > $noiseit3WholeNightdone 
  fi
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
  badFilesWarningsFile=identifiedBadFrames_preFlat_onlyStd_n$currentNight.txt
  rejectedFramesDir=$BDIR/rejectedFrames_std_preFlat_n$currentNight.txt
  if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
  removeBadFramesFromReduction $normit3dir $rejectedFramesDir $diagnosisPath $badFilesWarningsFile
  removeBadFramesFromReduction $normit3WholeNightdir $rejectedFramesDir $diagnosisPath $badFilesWarningsFile
    # Combining masked normalized images to make it3 flat
  if $RUNNING_FLAT; then
    flatit3BeforeCorrectiondir=$BDIR/flat-it3-Running-BeforeCorrection_n$currentNight
    flatit3BeforeCorrectiondone=$flatit3BeforeCorrectiondir/done_"$filter"_ccd"$h".txt
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
        propagateKeyword $i $dateHeaderKey $flatit3dir/$(basename "$i")
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
      astarithmetic $i -h1 set-m $currentFlatImage -h1 set-f m f $lowerVignettingThreshold lt nan where set-n n f $upperVignettingThreshold gt nan where -o $out
      propagateKeyword $i $airMassKeyWord $out 
      propagateKeyword $i $dateHeaderKey $out
      propagateKeyword $i $pointingRA $out
      propagateKeyword $i $pointingDEC $out
    done
    echo done > $maskedcornerdone
  fi

  
  # At this point we can process the frames of all the nights in the same way
  # So we place all the final frames into a common folder.
  if [ -f $framesForCommonReductionDone ]; then
    echo -e "\nFrames already placed in the folder for frames prepared to common reduction"
  else
   # This lockfile is created in order to handle the race conditions that could happen here
    lockfile="$framesForCommonReductionDir/lockfile"
    exec 200>$lockfile
    flock -x 200

    initialValue=$( getHighestNumberFromFilesInFolder $framesForCommonReductionDir )
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      name=$(( $initialValue + $a ))
      cp $maskedcornerdir/$base $framesForCommonReductionDir/$name.fits
      astfits $framesForCommonReductionDir/$name.fits -h1 --write=ORIGINAL_FILE,$base
    done
    echo "done" > $framesForCommonReductionDone
    
    flock -u 200 
    exec 200>&- 
  fi
}
export -f oneNightPreProcessing

writeTimeOfStepToFile "Process the individual nights" $fileForTimeStamps

nights=()
for currentNight in $(seq 1 $numberOfNights); do
      nights+=("$currentNight")
done
printf "%s\n" "${nights[@]}" | parallel --line-buffer -j "$num_cpus" oneNightPreProcessing {}
exit 0
totalNumberOfFrames=$( ls $framesForCommonReductionDir/*.fits | wc -l)
export totalNumberOfFrames
echo -e "* Total number of frames to combine: ${GREEN} $totalNumberOfFrames ${NOCOLOUR} *"

###Here we use the asignment file to rename files in the common folder
renameDone=$framesForCommonReductionDir/done_renaming.txt
if [ -f $renameDone ]; then
  echo -e "\nFiles already renamed in the folder for frames prepared to common reduction"
else
    python3 $pythonScriptsPath/rename_files_using_association.py $fileAssociation3by3 $framesForCommonReductionDir
    echo done > $renameDone
fi

echo -e "\n${GREEN} --- Astrometry --- ${NOCOLOUR}\n"

writeTimeOfStepToFile "Download Gaia catalogue" $fileForTimeStamps
echo -e "·Downloading Gaia Catalogue"

if (( sizeOfOurFieldDegrees > 1 )); then
  radiusToDownloadCatalogue=$( echo "$sizeOfOurFieldDegrees + 0.5" | bc | awk '{printf "%.1f", $0}' ) #The awk part is to avoiod problems when R<1
else
  radiusToDownloadCatalogue=$( echo "$sizeOfOurFieldDegrees + 1" | bc | awk '{printf "%.1f", $0}' ) #The awk part is to avoiod problems when R<1
fi

query_param="gaia --dataset=dr3 --center=$ra_gal,$dec_gal --radius=$radiusToDownloadCatalogue --column=ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error"
catdir=$BDIR/catalogs
catName=$catdir/"$objectName"_Gaia_DR3.fits
catRegionName=$catdir/"$objectName"_Gaia_DR3_regions.reg
catdone=$catdir/done.txt
if ! [ -d $catdir ]; then mkdir $catdir; fi
if [ -f $catdone ]; then
  echo -e "\n\tCatalogue is already downloaded\n"
else
  downloadGaiaCatalogue "$query_param" $catdir $catName
  python3 $pythonScriptsPath/createDS9RegionsFromCatalogue.py $catName $catRegionName "fits"
  echo "done" > $catdone
fi

astroimadir=$BDIR/astro-ima
astroimadone=$astroimadir/done_"$filter".txt
if ! [ -d $astroimadir ]; then mkdir $astroimadir; fi
if [ -f $astroimadone ]; then
  echo -e "\n\tImages are already astrometrized\n"
else
    #Since framesForCommonReduction is already astrometrized, we just need to copy to astroima
    cp $framesForCommonReductionDir/*.fits $astroimadir/
    echo done > $astroimadone
fi

# ########## Distorsion correction ##########
echo -e "\n ${GREEN} ---Creating distorsion correction files--- ${NOCOLOUR}"


# # Making sex catalogs and running scamp

# ****** Decision note *******
# In the LBT pipeline this two steps were done sequentially in two different blocks, first sextractor and then scamp
# I have put them together so we can loop them easily, because to perform an iterative astrometrisation we need to do
# sextractor + scamp, so for doing it in a confortable manner in the code both steps are in the same block
writeTimeOfStepToFile "Making sextractor catalogues and running scamp" $fileForTimeStamps
echo -e "·Creating SExtractor catalogues and running scamp"

numOfSextractorPlusScampIterations=2

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

# We need to warp to a huge grid (fullGrid) which covers the full field in order to combine all the images
# In order to save space and time, we crop it to the grid which contains data in each frame (smallGrid)
# and later we recover de full grid just before stacking

if ! [ -d "$BDIR/astro-ima" ]; then
  folderWithFramesToWarp=$BDIR/framesForCommonReduction
else
  folderWithFramesToWarp=$BDIR/astro-ima
fi

entiredir_smallGrid=$BDIR/pointings_smallGrid
entiredone=$entiredir_smallGrid/done_.txt
swarpcfg=$ROOTDIR/"$objectName"/config/swarp.cfg
export swarpcfg

if ! [ -d $entiredir_smallGrid ]; then mkdir $entiredir_smallGrid; fi

if [ -f $entiredone ]; then
    echo -e "\n\tImages already with astromety corrected using scamp-swarp and regrid to final grid (stored in pointings)\n"
else
  entiredir_fullGrid=$BDIR/pointings_fullGrid
  if ! [ -d $entiredir_fullGrid ]; then mkdir $entiredir_fullGrid; fi

  imagesToWarp=()
  for a in $(seq 1 $totalNumberOfFrames); do
      base="$a".fits
      imagesToWarp+=($folderWithFramesToWarp/$base)
  done

  printf "%s\n" "${imagesToWarp[@]}" | parallel -j "$num_cpus" warpImage {} $entiredir_fullGrid $entiredir_smallGrid $ra $dec $coaddSizePx $pipelinePath
  rm -rf $entiredir_fullGrid
  echo done > $entiredone
fi

echo -e "${GREEN} --- Compute and subtract Sky --- ${NOCOLOUR} \n"
noiseskydir=$BDIR/noise-sky_it1
noiseskydone=$noiseskydir/done_"$filter".txt

echo -e "·Dividing original background by 9"
imagesAreMasked=false
ringDir=$BDIR/ring
writeTimeOfStepToFile "Computing sky" $fileForTimeStamps

if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
if [ -f $noiseskydone ]; then
    echo -e "\n\tSky already computed\n"
else
  frameNames=()
  for a in $(seq 1 $totalNumberOfFrames); do
      base=entirecamera_"$a".txt
      frameNames+=("$base")
  done
  printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" divideSkyBy9 {} $rebinSkyDir $noiseskydir 
  echo done > $noiseskydone
fi

echo -e "\n·Subtracting background"
subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it1
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter".txt
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT

#### PHOTOMETRIC CALIBRATION  ####
echo -e "${ORANGE} ------ PHOTOMETRIC CALIBRATION ------ ${NOCOLOUR}\n"
writeTimeOfStepToFile "Photometric calibration" $fileForTimeStamps


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

# Calibration of individual frames
writeTimeOfStepToFile "Computing calibration factors for individual frames" $fileForTimeStamps
iteration=1
alphatruedir=$BDIR/alpha-stars-true_it$iteration
matchdir=$BDIR/match-decals-myData_it$iteration
ourDataCatalogueDir=$BDIR/ourData-aperture-photometry_it$iteration
prepareCalibrationCataloguePerFrame=$BDIR/survey-aperture-photometry_perBrick_it$iteration
mycatdir=$BDIR/my-catalog-halfmaxradius_it$iteration
imagesForCalibration=$subskySmallGrid_dir
calibratingMosaic=false

computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $mycatdir $rangeUsedCalibrationDir \
                          $mosaicDir $alphatruedir $calibrationBrightLimitIndividualFrames $calibrationFaintLimitIndividualFrames $tileSize $apertureUnits $numberOfApertureUnitsForCalibration $calibratingMosaic

echo -e "\n ${GREEN} ---Applying calibration factors--- ${NOCOLOUR}"
alphatruedir=$BDIR/alpha-stars-true_it$iteration
photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir $iteration $applyCommonCalibrationFactor

echo -e "\n${ORANGE} ------ STD WEIGHT COMBINATION ------ ${NOCOLOUR}\n"
# Compute rms and of the photometrized frames
noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done.txt
# Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
computeSky $photCorrSmallGridDir $noiseskydir $noiseskydone true $sky_estimation_method -1 false $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth "$noisechisel_param" "$maskParams"

photCorrfullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
photCorrfullGridDone=$photCorrfullGridDir/done.txt
if ! [ -d $photCorrfullGridDir ]; then mkdir $photCorrfullGridDir; fi
smallGridtoFullGrid $photCorrSmallGridDir $photCorrfullGridDir $photCorrfullGridDone $coaddSizePx $ra $dec


echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"
diagnosis_and_badFilesDir=$diagnosisPath
rejectedFramesDir=$BDIR/rejectedFrames
if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi

prefixOfTheFilesToRemove="entirecamera_"
rejectedByAstrometry=identifiedBadFrames_astrometry.txt
removeBadFramesFromReduction $photCorrfullGridDir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry $prefixOfTheFilesToRemove
removeBadFramesFromReduction $noiseskydir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry $prefixOfTheFilesToRemove
rejectedByBackgroundFWHM=identifiedBadFrames_fwhm.txt
removeBadFramesFromReduction $photCorrfullGridDir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM $prefixOfTheFilesToRemove
removeBadFramesFromReduction $noiseskydir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM $prefixOfTheFilesToRemove
rejectedByBackgroundValue=identifiedBadFrames_backgroundBrightness.txt
removeBadFramesFromReduction $photCorrfullGridDir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue $prefixOfTheFilesToRemove
removeBadFramesFromReduction $noiseskydir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundValue $prefixOfTheFilesToRemove
rejectedByCalibrationFactor=identifiedBadFrames_calibrationFactor.txt
removeBadFramesFromReduction $photCorrfullGridDir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByCalibrationFactor $prefixOfTheFilesToRemove
removeBadFramesFromReduction $noiseskydir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByCalibrationFactor $prefixOfTheFilesToRemove
rejectedByReadError=identifiedBadFrames_readError.txt
if [ -f $diagnosis_and_badFilesDir/$rejectedByReadError ]; then
  removeBadFramesFromReduction $photCorrfullGridDir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByReadError $prefixOfTheFilesToRemove
  removeBadFramesFromReduction $noiseskydir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByReadError $prefixOfTheFilesToRemove
fi

minRmsFileName="min_rms_it$iteration.txt"
python3 $pythonScriptsPath/find_rms_min.py "$filter" 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration $minRmsFileName

clippingdir=$BDIR/clipping-outliers_it$iteration
clippingdone=$clippingdir/done_"$k".txt
sigmaForStdSigclip=3
buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $photCorrfullGridDir $sigmaForStdSigclip

photCorrNoOutliersPxDir=$BDIR/photCorrFullGrid-dir_noOutliersPx_it$iteration
photCorrNoOutliersPxDone=$photCorrNoOutliersPxDir/done.txt
if ! [ -d $photCorrNoOutliersPxDir ]; then mkdir $photCorrNoOutliersPxDir; fi
removeOutliersFromWeightedFrames $photCorrfullGridDir $clippingdir $photCorrNoOutliersPxDir $photCorrNoOutliersPxDone

### Calculate the weights for the images based on the minimum rms ###
echo -e "\n ${GREEN} ---Computing weights for the frames--- ${NOCOLOUR}"
writeTimeOfStepToFile "Computing frame weights" $fileForTimeStamps

wdir=$BDIR/weight-dir_it$iteration
wonlydir=$BDIR/only-w-dir_it$iteration
wdone=$wdir/done.txt
wonlydone=$wonlydir/done.txt
if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
if ! [ -d $wdir ]; then mkdir $wdir; fi
computeWeights $wdir $wdone $wonlydir $wonlydone $photCorrNoOutliersPxDir $noiseskydir $iteration $minRmsFileName

echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
echo -e "\nBuilding coadd"
coaddDir=$BDIR/coadds_it$iteration 
coaddDone=$coaddDir/done.txt
coaddName=$coaddDir/"$objectName"_coadd_"$filter"_it$iteration.fits
buildCoadd $coaddDir $coaddName $wdir $wonlydir $coaddDone
maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
# noisechisel_param="--tilesize=30,30 --detgrowmaxholesize=5000 --dthresh=0.1 --snminarea=2 --rawoutput"
# export noisechisel_param
if [ -f $maskName ]; then
  echo "The mask of the weighted coadd is already done"
else
  astnoisechisel $coaddName $noisechisel_param --numthreads=$num_cpus -o $maskName
fi

exposuremapDir=$coaddDir/"$objectName"_exposureMap
exposuremapdone=$coaddDir/done_exposureMap.txt
computeExposureMap $wdir $exposuremapDir $exposuremapdone 


sblimitFile=$coaddDir/"$objectName"_"$filter"_sblimit.txt
exposuremapName=$coaddDir/exposureMap.fits
if [ -f  $sblimitFile ]; then
  echo -e "\n\tSurface brightness limit for coadd already measured\n"
  surfaceBrightnessLimit=$( awk '/Limiting magnitude/ { print $NF }' $sblimitFile )

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
          "LOWER_VIGNETTING_THRESHOLD" \
          "UPPER_VIGNETTING_THRESHOLD" \
          "SATURATION_THRESHOLD" \
          "CALIBRATED_USING" \
          "CALIBRATION_BRIGHTLIMIT" \
          "CALIBRATION_FAINTLIMIT" \
          "RUNNING_FLAT" \
          "WINDOW_SIZE" \
          "SURFACE_BRIGHTNESS_LIMIT")

numberOfFramesCombined=$(ls $BDIR/weight-dir_it$iteration/*.fits | wc -l)
values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$lowerVignettingThreshold" "$upperVignettingThreshold" "$saturationThreshold" "$surveyForPhotometry" "$calibrationBrightLimitIndividualFrames" "$calibrationFaintLimitIndividualFrames" "$RUNNING_FLAT" "$halfWindowSize" "$surfaceBrightnessLimit")
comments=("" "" "" "" "" "" "" "" "" "" "" "" "" "Running flat built with +-N frames" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")

astfits $coaddName --write=/,"Pipeline information"
addkeywords $coaddName keyWords values comments




framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted_it$iteration
framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi
if [ -f $framesWithCoaddSubtractedDone ]; then
    echo -e "\nFrames with coadd subtracted already generated\n"
else
  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_"$filter"_it$iteration.fits
  photCorrfullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
  
  coaddAv=$coaddDir/"$objectName"_coadd_"$filter"_it"$iteration"_average.fits
  astarithmetic $(ls -v $photCorrNoOutliersPxDir/*.fits) $(ls $photCorrNoOutliersPxDir/*.fits | wc -l) -g1 mean -o$coaddAv
  subtractCoaddToFrames $photCorrfullGridDir $coaddAv $framesWithCoaddSubtractedDir
  astarithmetic $(ls -v $framesWithCoaddSubtractedDir/*.fits) $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l) -g1 sum -o$sumMosaicAfterCoaddSubtraction
    echo done > $framesWithCoaddSubtractedDone
fi
endTime=$(TZ=UTC date +%D%T)
echo "Pipeline ended at : ${endTime}"