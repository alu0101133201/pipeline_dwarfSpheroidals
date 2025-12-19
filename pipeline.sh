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

astrometryModuleName="astrometry.net/0.94"
load_module $astrometryModuleName

scampModuleName="scamp/2.14.1"
load_module $scampModuleName

sexModuleName="sextractor/2.29.0"
load_module $sexModuleName

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



maskParams=$(printf "%s " "${maskToApply[@]}")
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
echo -e "\t·Python directory (pythonScriptsPath): ${ORANGE} ${pythonScriptsPath} ${NOCOLOUR}"

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

###Pre-processing: we generate the ring for each one of the detectors

ringtempDone=$DIR/done_templates.txt
if [ -f $ringtempDone ]; then
  echo -e "\n\tRing templates for each detector already generated"
else
  if [ "$USE_COMMON_RING" = true ]; then
    base="${commonRingDefinitionFile%.txt}"
    if [ $num_ccd -eq 1 ]; then
      cp $DIR/$commonRingDefinitionFile $DIR/"$base"_ccd"$h".txt
    else
      #We use data from night1, we in the end have data of the same camera always, and that is what matters
      prepareRingTemplate $commonRingDefinitionFile $INDIRo/night1 $DIR $overscan $trimsecKey
    fi
  else
    base_first="${firstRingDefinitionFile%.txt}"
    base_second="${secondRingDefinitionFile%.txt}"
    if [ $num_ccd -eq 1 ]; then
      cp $DIR/$firstRingDefinitionFile $DIR/"$base_first"_ccd"$h".txt
      cp $DIR/$secondRingDefinitionFile $DIR/"$base_second"_ccd"$h".txt
    else
      #We use data from night1, we in the end have data of the same camera always, and that is what matters
      prepareRingTemplate $firstRingDefinitionFile $INDIRo/night1 $DIR $overscan $trimsecKey
      prepareRingTemplate $secondRingDefinitionFile $INDIRo/night1 $DIR $overscan $trimsecKey
    fi
  fi
  echo done > $ringtempDone
fi


# Function which processes a whole night
oneNightPreProcessing() {
  currentNight=$1
  framesForCommonReductionDone=$framesForCommonReductionDir/done_"$filter"_n"$currentNight".txt

  echo -e "\n\n"
  echo -e "${ORANGE} --- STARTING TO PROCESS NIGHT NUMBER $currentNight --- ${NOCOLOUR}"

  

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
      for i in $currentINDIRo/*.fit*; do
          
        nameWithEscapedSpaces=$(escapeSpacesFromString "$i")
        DATEOBS=$(eval "astfits $nameWithEscapedSpaces -h0 --keyvalue=$dateHeaderKey --quiet")
        checkIfExist_DATEOBS $DATEOBS
        
	      if [[ $dateHeaderKey =~ ^MJD ]]; then
          unixTimeInSeconds=$(astarithmetic $DATEOBS 40587 - 86400 x -q | awk '{printf "%.0f", $1}')
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
	      for h in $(seq 0 $num_ccd); do
            	# HERE A CHECK IF THE DATA IS IN FLOAT32 IS NEEDED
	    	  if [ $h -eq 0 ]; then
                if [ $num_ccd -ne 0 ]; then
            		  eval "astfits $nameWithEscapedSpaces --copy=$h --primaryimghdu -o$out"  # I run this with eval so the escaped spaces are re-parsed by bash and understood by astfits
            		  nameOfOriginalFile="${nameWithEscapedSpaces##*/}"
            		  eval "astfits --write=OriginalName,$nameOfOriginalFile $out -h0"
                else
                  eval "astfits $nameWithEscapedSpaces --copy=$h -o$out"
                  nameOfOriginalFile="${nameWithEscapedSpaces##*/}"
                  eval "astfits --write=OriginalName,$nameOfOriginalFile $out -h0"
                fi
		      else
            
  			    if [[ "$overscan" == "YES" ]]; then
  				    trsec=$(eval "astfits $nameWithEscapedSpaces -h $h --keyvalue=$trimsecKey -q")
      				trsec=${trsec//[\[\]]/}
              temp1=$currentINDIR/"$unixTimeInSeconds"_ccd"$h"_temp1.fits
	  			    eval "astcrop $nameWithEscapedSpaces -h $h --mode=img --section=$trsec  -o $temp1"
              temp2=$currentINDIR/"$unixTimeInSeconds"_ccd"$h"_temp2.fits
              astarithmetic $temp1 -h1 float32 -o$temp2 
              rm $temp1
              astfits $temp2 --copy=1  -o$out
              rm $temp2
      				gain_h=$(eval "astfits $nameWithEscapedSpaces -h $h --keyvalue=$gain -q")
              astfits $out -h$h --write=$gain,$gain_h,"e- ADU"
      	    else
	 			      eval "astfits $nameWithEscapedSpaces --copy=$h -o$out"
     		    fi
		      fi
        done
      done

	  index=1
    for i in $(ls -v $currentINDIR/*.fits); do
      	mv $i $currentINDIR/"$objectName"-Decals-"$filter"_n"$currentNight"_f"$index".fits
      	index=$((index+1));
    done
      
    echo done > $renamedone
  fi

  # -------------------------------------------------------
  # Number of exposures of the current night
  n_exp=$(ls -v $currentINDIRo/*.fit* | wc -l)
  echo -e "Number of exposures ${ORANGE} ${n_exp} ${NOCOLOUR}"
  ###If WINDOW_SIZE >= number of exposures; then it is the same as making whole night flat
  window_size=$(( (halfWindowSize * 2) + 1 ))
  local RUNNING_FLAT_night=$RUNNING_FLAT
  if [ $n_exp -le $window_size ]; then
    RUNNING_FLAT_night=false
  fi
  # If you have only one set of darks, you can save them in dark/*
  if [ -d $DARKDIR/night"$currentNight" ]; then
    currentDARKDIR=$DARKDIR/night$currentNight
  else
    currentDARKDIR=$DARKDIR
  fi
  mdadir=$BDIR/masterdark_n$currentNight

  # Loop for all the ccds
  for h in $(seq 1 $num_ccd); do
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
      #There is OVERSCAN probably
      if [[ "$overscan" == "YES" ]]; then
          first_file=$(echo "$escaped_files" | awk '{print $1}')
          trsec=$(eval "astfits $first_file -h$h --keyvalue=$trimsecKey -q")
          trsec=${trsec//[\[\]]/}
          astcrop $mdadir/temp.fits -h1 --mode=img --section=$trsec  -o $mdadir/temp2.fits
          astfits $mdadir/temp2.fits --copy=1 -o $mdadir/mdark_"$filter"_n"$currentNight".fits
          rm $mdadir/temp2.fits
      else
          astfits $mdadir/temp.fits --copy=1 -o $mdadir/mdark_"$filter"_n"$currentNight".fits
      fi
      rm $mdadir/temp.fits
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
    if [ $num_ccd -eq 1 ]; then 
      h_air=1
    else
      h_air=0
    fi
    for i in $(ls -v $currentINDIR/*.fits ); do
      
      air=$(astfits $i -h $h_air --keyvalue=$airMassKeyWord 2>/dev/null | awk '{print $2}')
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

  if ! [ -d $mbiascorrdir ]; then mkdir $mbiascorrdir; fi
  
  mbiascorrdone=$mbiascorrdir/done_"$filter".txt
  
  if [ -f $mbiascorrdone ]; then
      echo -e "\nScience images are already bias/dark corrected for night $currentNight\n"
  else
    framesToSubtract=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
      framesToSubtract+=("$base")
    done
    dark=$mdadir/mdark_"$filter"_n"$currentNight".fits
    printf "%s\n" "${framesToSubtract[@]}" | parallel -j "$num_cpus" subtractBiasFromFrame {} $dark $saturationThreshold $currentINDIR $mbiascorrdir
    echo done > $mbiascorrdone
  fi
       



  echo -e "${ORANGE} ------ FLATS ------ ${NOCOLOUR}\n"
  echo -e "${GREEN} --- Flat iteration 1 --- ${NOCOLOUR}"

  ########## Creating the ring mask ##########
  # Since the ring possibilities are quite flexible (using one common ring or two rings depending on the angle) I always clean and rebuild the rings (It doesn't cost almost time so is worth it)
  
  
  # We always need the common ring  definition always stored for photometric calibration (selection of decals bricks to download)
  
  # We create the .fits ring image based on how the normalisation is going to be done
  ringdir=$BDIR/ring
  if ! [ -d $ringdir ]; then mkdir $ringdir; fi
  lockfile="$ringdir"/lockfile
  exec 200>$lockfile
  flock -x 200
  if ! [ -f "$ringdir/ring.fits" ] && ! [ -f "$ringdir/ring_1.fits" ] && ! [ -f "$ringdir/ring_2.fits" ]; then
    
    for h in $(seq 1 $num_ccd); do
      if [ "$USE_COMMON_RING" = true ]; then
        base="${commonRingDefinitionFile%.txt}"
        cp $DIR/"$base"_ccd"$h".txt $ringdir/ring_ccd"$h".txt
        astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1.fits --backhdu=$h --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_temp.fits $ringdir/ring_ccd"$h".txt
        astfits $ringdir/ring_temp.fits --copy=1 -o $ringdir/ring.fits
      else
        base_first="${firstRingDefinitionFile%.txt}"
        base_second="${secondRingDefinitionFile%.txt}"
        cp $DIR/"$base_first"_ccd"$h".txt $ringdir/ring_1_ccd"$h".txt
        cp $DIR/"$base_second"_ccd"$h".txt $ringdir/ring_2_ccd"$h".txt
        astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1.fits -h$h --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_2_temp.fits $ringdir/ring_2_ccd"$h".txt
        astfits $ringdir/ring_2_temp.fits --copy=1 -o $ringdir/ring_2.fits
        astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1.fits -h$h --mforflatpix --mode=img --type=uint8 --circumwidth=$ringWidth --clearcanvas -o $ringdir/ring_1_temp.fits $ringdir/ring_1_ccd"$h".txt
        astfits $ringdir/ring_1_temp.fits --copy=1 -o $ringdir/ring_1.fits
      fi
      rm $ringdir/*temp*
    done
  fi
  flock -u 200
  exec 200>&-
  
  #Parallel condition. If $ringdir/ring.fits or ring_1 or ring_2, we wait
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
  normit1done=$normit1dir/done_"$filter".txt
  if ! [ -d $normit1dir ]; then mkdir $normit1dir; fi
  if [ -f $normit1done ]; then
    echo -e "\nScience images are already normalized for night $currentNight\n"
  else
    normaliseImagesWithRing $mbiascorrdir $normit1dir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
    echo done > $normit1done
  fi

  # Then, if the running flat is configured to be used, we combine the normalised images with a sigma clipping median
  # using the running flat strategy
  if $RUNNING_FLAT_night; then
    flatit1dir=$BDIR/flat-it1-Running_n$currentNight
    flatit1done=$flatit1dir/done_"$filter".txt
    iteration=1
    if ! [ -d $flatit1dir ]; then mkdir $flatit1dir; fi
    if [ -f $flatit1done ]; then
      echo -e "\nRunning flats it-1 already built for night $currentNight\n"
    else
      calculateRunningFlat $normit1dir $flatit1dir $flatit1done $iteration
    fi
  fi

  # We compute the flat using all the frames of the night
  flatit1WholeNightdir=$BDIR/flat-it1-WholeNight_n$currentNight
  flatit1WholeNightdone=$flatit1WholeNightdir/done_"$filter".txt
  iteration=1
  if ! [ -d $flatit1WholeNightdir ]; then mkdir $flatit1WholeNightdir; fi
  if [ -f $flatit1WholeNightdone ]; then
    echo -e "\nWhole night flat it-1 already built for night $currentNight and extension $h\n"
  else
    calculateFlat $flatit1WholeNightdir/flat-it1_wholeNight_n$currentNight.fits $iteration $normit1dir/*.fits
    echo "done" >> $flatit1WholeNightdone
  fi


  # Dividing the science images for the running it1 flat
  if $RUNNING_FLAT_night; then
    flatit1imadir=$BDIR/flat-it1-Running-ima_n$currentNight
    flatit1imadone=$flatit1imadir/done_"$filter".txt
    if ! [ -d $flatit1imadir ]; then mkdir $flatit1imadir; fi
    if [ -f $flatit1imadone ]; then
      echo -e "\nScience images are divided by flat it1 for night $currentNight\n"
    else
      divideImagesByRunningFlats $mbiascorrdir $flatit1imadir $flatit1dir $flatit1imadone
    fi
  fi

  # Dividing the science images for the whole night it1 flat
  flatit1WholeNightimaDir=$BDIR/flat-it1-WholeNight-ima_n$currentNight
  flatit1WholeNightimaDone=$flatit1WholeNightimaDir/done_"$filter".txt
  if ! [ -d $flatit1WholeNightimaDir ]; then mkdir $flatit1WholeNightimaDir; fi
  if [ -f $flatit1WholeNightimaDone ]; then
    echo -e "\nScience images are divided by whole night flat it1 for night $currentNight \n"
  else
    wholeNightFlatToUse=$flatit1WholeNightdir/flat-it1_wholeNight_n$currentNight.fits
    divideImagesByWholeNightFlat $mbiascorrdir $flatit1WholeNightimaDir $wholeNightFlatToUse $flatit1WholeNightimaDone
  fi

  ########## Creating the it2 master flat image ##########
  echo -e "${GREEN} --- Flat iteration 2 --- ${NOCOLOUR}"
  # Obtain a mask using noisechisel on the running flat images
  if $RUNNING_FLAT_night; then
    noiseit2dir=$BDIR/noise-it2-Running_n$currentNight
    noiseit2done=$noiseit2dir/done_"$filter".txt
    if ! [ -d $noiseit2dir ]; then mkdir $noiseit2dir; fi
    if [ -f $noiseit2done ]; then
      echo -e "\nScience images are 'noisechiseled' for it2 running flat for night $currentNight \n"
    else
      frameNames=()
      for a in $(seq 1 $n_exp); do
          base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
          frameNames+=("$base")
      done
      printf "%s\n" "${frameNames[@]}" | parallel -j "$num_parallel" runNoiseChiselOnFrame {} $flatit1imadir $noiseit2dir $blockScale "'$noisechisel_param'"
      echo done > $noiseit2done
    fi
  fi

  # Obtain a mask using noisechisel on the whole night flat images
  noiseit2WholeNightDir=$BDIR/noise-it2-WholeNight_n$currentNight
  noiseit2WholeNightdone=$noiseit2WholeNightDir/done_"$filter".txt
  if ! [ -d $noiseit2WholeNightDir ]; then mkdir $noiseit2WholeNightDir; fi
  if [ -f $noiseit2WholeNightdone ]; then
    echo -e "\nScience images are 'noisechiseled' for it2 whole night flat for night $currentNight\n"
  else
    frameNames=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
      frameNames+=("$base")
    done
    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_parallel" runNoiseChiselOnFrame {} $flatit1WholeNightimaDir $noiseit2WholeNightDir $blockScale "'$noisechisel_param'"
    echo done > $noiseit2WholeNightdone
  fi

  # Mask the images (running flat)
  if $RUNNING_FLAT_night; then
    maskedit2dir=$BDIR/masked-it2-Running_n$currentNight
    maskedit2done=$maskedit2dir/done_"$filter".txt
    if ! [ -d $maskedit2dir ]; then mkdir $maskedit2dir; fi
    if [ -f $maskedit2done ]; then
      echo -e "\nScience images are masked for running flat, night $currentNight\n"
    else
      maskImages $mbiascorrdir $noiseit2dir $maskedit2dir $USE_COMMON_RING $keyWordToDecideRing
      echo done > $maskedit2done
    fi
  fi

  # Mask the images (whole night flat)
  maskedit2WholeNightdir=$BDIR/masked-it2-WholeNight_n$currentNight
  maskedit2WholeNightdone=$maskedit2WholeNightdir/done_"$filter".txt
  if ! [ -d $maskedit2WholeNightdir ]; then mkdir $maskedit2WholeNightdir; fi
  if [ -f $maskedit2WholeNightdone ]; then
    echo -e "\nScience images are masked for whole night flat, night $currentNight \n"
  else
    maskImages $mbiascorrdir $noiseit2WholeNightDir $maskedit2WholeNightdir $USE_COMMON_RING $keyWordToDecideRing
    echo done > $maskedit2WholeNightdone
  fi


  # Normalising masked images (running flat)
  if $RUNNING_FLAT_night; then
    normit2dir=$BDIR/norm-it2-Running-images_n$currentNight
    normit2done=$normit2dir/done_"$filter".txt
    if ! [ -d $normit2dir ]; then mkdir $normit2dir; fi
    if [ -f $normit2done ]; then
      echo -e "\nMasked science images are normalized for running flat, night $currentNight\n"
    else
      normaliseImagesWithRing $maskedit2dir $normit2dir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
      echo done > $normit2done
    fi
  fi

  # Normalising masked images (whole night flat)
  normit2WholeNightdir=$BDIR/norm-it2-WholeNight-images_n$currentNight
  normit2WholeNightdone=$normit2WholeNightdir/done_"$filter".txt
  if ! [ -d $normit2WholeNightdir ]; then mkdir $normit2WholeNightdir; fi
  if [ -f $normit2WholeNightdone ]; then
    echo -e "\nMasked science images are normalized for whole night flat, night $currentNight\n"
  else
    normaliseImagesWithRing $maskedit2WholeNightdir $normit2WholeNightdir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
    echo done > $normit2WholeNightdone
  fi

  # Combining masked normalized images to make it2 running flat
  if $RUNNING_FLAT_night; then
    flatit2dir=$BDIR/flat-it2-Running_n$currentNight
    flatit2done=$flatit2dir/done_"$filter".txt
    iteration=2
    if ! [ -d $flatit2dir ]; then mkdir $flatit2dir; fi
    if [ -f $flatit2done ]; then
      echo -e "\nScience images are stacked for it2 running flat for night $currentNight\n"
    else
      calculateRunningFlat $normit2dir $flatit2dir $flatit2done $iteration
    fi
  fi

  # We also compute the flat using all the frames of the night.
  flatit2WholeNightdir=$BDIR/flat-it2-WholeNight_n$currentNight
  flatit2WholeNightdone=$flatit2WholeNightdir/done_"$filter".txt
  iteration=2
  if ! [ -d $flatit2WholeNightdir ]; then mkdir $flatit2WholeNightdir; fi
  if [ -f $flatit2WholeNightdone ]; then
    echo -e "\nWhole night flat it-2 already built for night $currentNight \n"
  else
    calculateFlat $flatit2WholeNightdir/flat-it2_wholeNight_n$currentNight.fits $iteration $normit2WholeNightdir/*.fits
    echo "done" >> $flatit2WholeNightdone
  fi

  
  # Dividing the science image by the it2 flat
  if $RUNNING_FLAT_night; then
    flatit2imadir=$BDIR/flat-it2-Running-ima_n$currentNight
    flatit2imadone=$flatit2imadir/done_"$filter".txt
    if ! [ -d $flatit2imadir ]; then mkdir $flatit2imadir; fi
    if [ -f $flatit2imadone ]; then
      echo -e "\nRunning flats it2-2 already built for night $currentNight \n"
    else
      divideImagesByRunningFlats $mbiascorrdir $flatit2imadir $flatit2dir $flatit2imadone
    fi
  fi

  # Dividing the science images for the whole night it2 flat
  flatit2WholeNightimaDir=$BDIR/flat-it2-WholeNight-ima_n$currentNight
  flatit2WholeNightimaDone=$flatit2WholeNightimaDir/done_"$filter".txt
  if ! [ -d $flatit2WholeNightimaDir ]; then mkdir $flatit2WholeNightimaDir; fi
  if [ -f $flatit2WholeNightimaDone ]; then
    echo -e "\nScience images are divided by whole night flat it2 for night $currentNight \n"
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

  # **** Decision note about multidetector ****
  # We aim to keep the multilayer structure until the end. Because of that, it might be better to change compute sky itself and deal with the multilayer inside.
  # Since we write every value of sky into a .txt, we are gonna put them line by line for each frame
  # However, diagnosis and bad files is a problem. My approximation will be the following:
  #   1) Create (or check if exists) inside diagnosis and bad files, one directory for each ccd called CCD1, CCD2 etc. 
  #   2) Each time we call a python function, we will add the $h value and make a for loop. Might not be the finest idea, but is the one I think avoids more problems

  #Since we are not rejecting yet frames, we will comment this step for the moment
  #tmpNoiseDir=$BDIR/noisesky_forCleaningBadFramesBeforeFlat
  #tmpNoiseDone=$tmpNoiseDir/done.txt
  #if ! [ -d $tmpNoiseDir ]; then mkdir $tmpNoiseDir; fi
#
  #diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
  #badFilesWarningsFile=identifiedBadFrames_preFlat_onlyStd.txt
  #badFilesWarningsDone=$diagnosis_and_badFilesDir/done_badFrames_stdPreFlat_n"$currentNight".txt
  #if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi
  #for h in $(seq 1 $num_ccd); do
  #  diagnosis_and_badFilesDir_ccd=$diagnosis_and_badFilesDir/CCD$h
  #  if ! [ -d $diagnosis_and_badFilesDir_ccd ]; then mkdir $diagnosis_and_badFilesDir_ccd; fi
  #done
#
  #if [ -f $badFilesWarningsDone ]; then
  #    echo -e "\n\tFrames with strange background value and std values already cleaned for night $currentNight\n"
  #else
  #  computeSky $flatit2WholeNightimaDir $tmpNoiseDir $tmpNoiseDone true $sky_estimation_method -1 false $ringdir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth NO $blockScale "'$noisechisel_param'" "'$maskParams'"
  #  for h in $(seq 1 $num_ccd); do
  #    diagnosis_and_badFilesDir_ccd=$diagnosis_and_badFilesDir/CCD$h
  #    python3 $pythonScriptsPath/checkForBadFrames_beforeFlat_std.py  $tmpNoiseDir $diagnosis_and_badFilesDir_ccd $badFilesWarningsFile $numberOfStdForBadFrames $h
  #  done
  #  echo "done" > $badFilesWarningsDone
  #fi

  
  ########## Creating the it3 master flat image ##########
  echo -e "${GREEN} --- Flat iteration 3 --- ${NOCOLOUR}"

  # Obtain a mask using noisechisel on the running flat images
  if $RUNNING_FLAT_night; then
    noiseit3dir=$BDIR/noise-it3-Running_n$currentNight
    noiseit3done=$noiseit3dir/done_"$filter".txt
    if ! [ -d $noiseit3dir ]; then mkdir $noiseit3dir; fi
    if [ -f $noiseit3done ]; then
      echo -e "\nScience images are 'noisechiseled' for it3 running flat for night $currentNight\n"
    else
      frameNames=()
      for a in $(seq 1 $n_exp); do
          base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
          frameNames+=("$base")
      done
      printf "%s\n" "${frameNames[@]}" | parallel -j "$num_parallel" runNoiseChiselOnFrame {} $flatit2imadir $noiseit3dir $blockScale "'$noisechisel_param'"
      echo done > $noiseit3done
    fi
  fi


  # Obtain a mask using noisechisel on the whole night flat images
  noiseit3WholeNightDir=$BDIR/noise-it3-WholeNight_n$currentNight
  noiseit3WholeNightdone=$noiseit3WholeNightDir/done_"$filter".txt
  if ! [ -d $noiseit3WholeNightDir ]; then mkdir $noiseit3WholeNightDir; fi
  if [ -f $noiseit3WholeNightdone ]; then
    echo -e "\nScience images are 'noisechiseled' for it3 whole night flat for night $currentNight\n"
  else
    frameNames=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
      frameNames+=("$base")
    done

    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_parallel" runNoiseChiselOnFrame {} $flatit2WholeNightimaDir $noiseit3WholeNightDir $blockScale "'$noisechisel_param'"
    echo done > $noiseit3WholeNightdone 
  fi

  

  # Mask the images (running flat)
  if $RUNNING_FLAT_night; then
    maskedit3dir=$BDIR/masked-it3-Running_n$currentNight
    maskedit3done=$maskedit3dir/done_"$filter".txt
    if ! [ -d $maskedit3dir ]; then mkdir $maskedit3dir; fi
    if [ -f $maskedit3done ]; then
      echo -e "\nScience images are masked for running flat, night $currentNight\n"
    else
      maskImages $mbiascorrdir $noiseit3dir $maskedit3dir $USE_COMMON_RING $keyWordToDecideRing
      echo done > $maskedit3done
    fi
  fi

  
  # Mask the images (whole night flat)
  maskedit3WholeNightdir=$BDIR/masked-it3-WholeNight_n$currentNight
  maskedit3WholeNightdone=$maskedit3WholeNightdir/done_"$filter".txt
  if ! [ -d $maskedit3WholeNightdir ]; then mkdir $maskedit3WholeNightdir; fi
  if [ -f $maskedit3WholeNightdone ]; then
    echo -e "\nScience images are masked for whole night flat, night $currentNight\n"
  else
    maskImages $mbiascorrdir $noiseit3WholeNightDir $maskedit3WholeNightdir $USE_COMMON_RING $keyWordToDecideRing
    echo done > $maskedit3WholeNightdone
  fi

  
  # Normalising masked images (running flat)
  if $RUNNING_FLAT_night; then
    normit3dir=$BDIR/norm-it3-Running-images_n$currentNight
    normit3done=$normit3dir/done_"$filter".txt
    if ! [ -d $normit3dir ]; then mkdir $normit3dir; fi
    if [ -f $normit3done ]; then
      echo -e "\nMasked science images are normalized for running flat, night $currentNight\n"
    else
      normaliseImagesWithRing $maskedit3dir $normit3dir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
      echo done > $normit3done
    fi
  fi

  # Normalising masked images (whole night flat)
  normit3WholeNightdir=$BDIR/norm-it3-WholeNight-images_n$currentNight
  normit3WholeNightdone=$normit3WholeNightdir/done_"$filter".txt
  if ! [ -d $normit3WholeNightdir ]; then mkdir $normit3WholeNightdir; fi
  if [ -f $normit3WholeNightdone ]; then
    echo -e "\nMasked science images are normalized for whole night flat, night $currentNightCmo\n"
  else
    normaliseImagesWithRing $maskedit3WholeNightdir $normit3WholeNightdir $USE_COMMON_RING $ringdir/ring.fits $ringdir/ring_2.fits $ringdir/ring_1.fits $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing 
    echo done > $normit3WholeNightdone
  fi
  
  # Remove the identified bad frames ONLY for the flat, they will still be present in following steps, but not used in the flat calculation
  #diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
  #badFilesWarningsFile=identifiedBadFrames_preFlat_onlyStd.txt
  #rejectedFramesDir=$BDIR/rejectedFrames_std_preFlat
  #if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
  #removeBadFramesFromReduction $normit3dir $rejectedFramesDir $diagnosis_and_badFilesDir $badFilesWarningsFile
  #removeBadFramesFromReduction $normit3WholeNightdir $rejectedFramesDir $diagnosis_and_badFilesDir $badFilesWarningsFile


  # Combining masked normalized images to make it3 flat
  if $RUNNING_FLAT_night; then
    flatit3BeforeCorrectiondir=$BDIR/flat-it3-Running-BeforeCorrection_n$currentNight
    flatit3BeforeCorrectiondone=$flatit3BeforeCorrectiondir/done_"$k".txt
    iteration=3
    if ! [ -d $flatit3BeforeCorrectiondir ]; then mkdir $flatit3BeforeCorrectiondir; fi
    if [ -f $flatit3BeforeCorrectiondone ]; then
      echo -e "\nRunning flats it3 before correction are already build for night $currentNight \n"
    else
      calculateRunningFlat $normit3dir $flatit3BeforeCorrectiondir $flatit3BeforeCorrectiondone $iteration
    fi
  fi

  # We also compute the flat using all the frames of the night.
  flatit3WholeNightdir=$BDIR/flat-it3-WholeNight_n$currentNight
  flatit3WholeNightdone=$flatit3WholeNightdir/done_"$filter".txt
  iteration=3
  if ! [ -d $flatit3WholeNightdir ]; then mkdir $flatit3WholeNightdir; fi
  if [ -f $flatit3WholeNightdone ]; then
    echo -e "\nWhole night flat it-3 already built for night $currentNight \n"
  else
    calculateFlat $flatit3WholeNightdir/flat-it3_wholeNight_n$currentNight.fits $iteration $normit3WholeNightdir/*.fits
    echo "done" >> $flatit3WholeNightdone
  fi

  # Correct the running flats using the whole night flat
  flatit3dir=$BDIR/flat-it3-Running_n$currentNight
  if $RUNNING_FLAT_night; then
    
    flatit3done=$flatit3dir/done_"$k".txt
    if ! [ -d $flatit3dir ]; then mkdir $flatit3dir; fi
    if [ -f $flatit3done ]; then
      echo -e "\nFlats iteration 3 are corrected using the flat of the whole night for night $currentNight \n"
    else
      imagesToCorrect=()
      for i in $flatit3BeforeCorrectiondir/*.fits; do
        imagesToCorrect+=("$(basename $i)")
      done
      printf "%s\n" "${imagesToCorrect[@]}" | parallel -j "$num_cpus" correctRunningFlatWithWholeNightFlat {} $flatit3BeforeCorrectiondir $flatit3WholeNightdir/flat-it3_wholeNight_n$currentNight.fits $flatit3dir 
      echo done > $flatit3done
    fi
  fi

  # Dividing the science image by the it3 flat
  # If running flat selected, we use it to produce the final flatted images
  # If not selcted, we applyt the whole night flat
  flatit3imadir=$BDIR/flat-it3-ima_n$currentNight
  flatit3imadone=$flatit3imadir/done_"$filter".txt
  if ! [ -d $flatit3imadir ]; then mkdir $flatit3imadir; fi
  if $RUNNING_FLAT_night; then
    if [ -f $flatit3imadone ]; then
      echo -e "\nScience images are divided by the it3 flat for night $currentNight\n"
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
  maskedcornerdone=$maskedcornerdir/done_"$filter".txt
  if ! [ -d $maskedcornerdir ]; then mkdir $maskedcornerdir; fi
  if [ -f $maskedcornerdone ]; then
    echo -e "\nCorners are already masked for night $currentNight \n"
  else
    imagesForVignetting=()
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
      imagesForVignetting+=("$base")
    done
    printf "%s\n" "${imagesForVignetting[@]}" | parallel -j "$num_cpus" maskVignettingOnImages {} $flatit3imadir $maskedcornerdir $flatit3dir $flatit3WholeNightdir $RUNNING_FLAT_night $n_exp $currentNight $lowerVignettingThreshold $upperVignettingThreshold 
    echo done > $maskedcornerdone
  fi
  
  # At this point we can process the frames of all the nights in the same way
  # So we place all the final frames into a common folder.

  # WE SHOULD MANAGE THE RACE CONDITIONS
  if [ -f $framesForCommonReductionDone ]; then
    echo -e "\nFrames already placed in the folder for frames prepared to common reduction"
  else
    lockfile="$framesForCommonReductionDir/lockfile"
    exec 200>$lockfile
    flock -x 200

    initialValue=$( getHighestNumberFromFilesInFolder $framesForCommonReductionDir )
    
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a".fits
      name=$(( $initialValue + $a ))
      cp $maskedcornerdir/$base $framesForCommonReductionDir/$name.fits
      ###We need, so far: GAIN, DATEOBS and AIRMASS: we check if they exist and if not we propagate them
      air=$(astfits $framesForCommonReductionDir/$name.fits -h0 --keyvalue=$airMassKeyWord -q)
      dobs=$(astfits $framesForCommonReductionDir/$name.fits -h0 --keyvalue=$dateHeaderKey -q)

      if [ "$air" == "n/a" ]; then
        air=$(astfits $currentINDIR/$base -h0 --keyvalue=$airMassKeyWord -q)
        astfits $framesForCommonReductionDir/$name.fits -h0 --write=$airMassKeyWord,$air 
      fi
      if [ "$dobs" == "n/a" ]; then
        dobs=$(astfits $currentINDIR/$base -h0 --keyvalue=$dateHeaderKey -q)
        astfits $framesForCommonReductionDir/$name.fits -h0 --write=$dateHeaderKey,$dobs
      fi
      for h in $(seq 1 $num_ccd); do
        gain_h=$(astfits $framesForCommonReductionDir/$name.fits -h$h --keyvalue=$gain -q)
        if [ "$gain_h" == "n/a" ]; then
          gain_h=$(astfits $currentINDIR/$base -h$h --keyvalue=$gain -q)
          astfits $framesForCommonReductionDir/$name.fits -h$h --write=$gain,$gain_h
        fi
      done

    done
    echo "done" > $framesForCommonReductionDone

    flock -u 200
    exec 200>&-
  fi

  # # Removing intermediate information to save space
    rm -rf $BDIR/bias-corrected_n$currentNight
  rm -rf $BDIR/masterdark_n$currentNight
  rm -rf $BDIR/flat-it3-Running-BeforeCorrection_n$currentNight
  rm -rf $BDIR/flat-it3-ima_n$currentNight
  rm -rf $BDIR/flat-it1-Running_n$currentNight
  rm -rf $BDIR/flat-it1-WholeNight_n$currentNight
  rm -rf $BDIR/flat-it2-Running_n$currentNight
  rm -rf $BDIR/flat-it2-WholeNight_n$currentNight
  rm -rf $BDIR/noise-it2-Running_n$currentNight
  rm -rf $BDIR/noise-it2-WholeNight_n$currentNight

  for a in $(seq 1 3); do
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
writeTimeOfStepToFile "Download Gaia catalogue" $fileForTimeStamps
echo -e "·Downloading Gaia Catalogue"

# Here I add some extra size to the field used to download the gaia catalogue for two reasons
# 1.- This is because you don't want to use a field too big in order not to download bricks that you don't need
# So I expect the "sizeofOurFieldDegrees" value to be quite tight. But since the catalogue is text and it doesn't take
# long I prefer to add something and be sure that I don't lose any source because of the catalogue

# 2.- This gaia catalogue is used to match the survey data for calibration to gaia, calibrating thus the survey to our photometric framework
# of gaia. It has to be large enough to be to perform this calibration process

if ((  $(echo "$sizeOfOurFieldDegrees > 1.0" | bc -l) )); then
  radiusToDownloadCatalogue=$( echo "$sizeOfOurFieldDegrees + 0.25" | bc -l | awk '{printf "%.1f", $0}' ) #The awk part is to avoiod problems when R<1
  radiusToDownloadCatalogue_gaia=$radiusToDownloadCatalogue
else
  radiusToDownloadCatalogue=$( echo "$sizeOfOurFieldDegrees + 0.5" | bc -l | awk '{printf "%.1f", $0}' ) #The awk part is to avoiod problems when R<1
  radiusToDownloadCatalogue_gaia=1.5
fi

query_param_ps="vizier --dataset=panstarrs1 --center=$ra_gal,$dec_gal --radius=$radiusToDownloadCatalogue --column=RAJ2000,DEJ2000,gmag"

query_param_gaia="gaia --dataset=dr3 --center=$ra_gal,$dec_gal --radius=$radiusToDownloadCatalogue_gaia --column=ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error"

catdir=$BDIR/catalogs
catName_gaia=$catdir/"$objectName"_Gaia_DR3.fits
catName_ps=$catdir/"$objectName"_Panstarrs_S1.fits
catRegionName=$catdir/"$objectName"_Gaia_DR3_regions.reg
catdone=$catdir/done.txt
if ! [ -d $catdir ]; then mkdir $catdir; fi
if [ -f $catdone ]; then
  echo -e "\n\tCatalogue is already downloaded\n"
else
  downloadGaiaCatalogue "$query_param_gaia" $catdir $catName_gaia
  downloadPanstarrsCatalogue "$query_param_ps" $catdir $catName_ps
  python3 $pythonScriptsPath/createDS9RegionsFromCatalogue.py $catName $catRegionName "fits"
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
  printf "%s\n" "${indexes[@]}" | parallel -j "$num_cpus" downloadIndex {} $catName_ps $indexdir
  echo done > $indexdone
fi


sexcfg=$CDIR/sextractor_solvefield.sex
# Solving the images
writeTimeOfStepToFile "Solving fields" $fileForTimeStamps
echo -e "·Solving fields"

astrocfg=$CDIR/astrometry_$objectName.cfg

rm $astrocfg
echo inparallel > $astrocfg
echo cpulimit 300 >> $astrocfg
echo "add_path $indexdir" >> $astrocfg
echo autoindex >> $astrocfg

astroimadir=$BDIR/astro-ima
astroimadir_layer=$astroimadir/layers
astroimadone=$astroimadir/done_"$filter".txt
if ! [ -d $astroimadir ]; then mkdir $astroimadir; fi
if ! [ -d $astroimadir_layer ]; then mkdir $astroimadir_layer; fi
if [ -f $astroimadone ]; then
  echo -e "\n\tImages are already astrometrized\n"
else
  #LBT frames are already astrometrized, so we will take advantage of that info
  if [ "$telescope" == "LBT" ] || [ "$telescope" == "OSIRIS" ]; then
    cp $framesForCommonReductionDir/*.fits $astroimadir/
  else
    frameNames=()
    for a in $(seq 1 $totalNumberOfFrames); do
        base=$a.fits
        i=$framesForCommonReductionDir/$base
        frameNames+=("$i")
    done
    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" solveField {} $solve_field_L_Param $solve_field_H_Param $solve_field_u_Param $ra_gal $dec_gal $CDIR $astroimadir_layer $sexcfg $sizeOfOurFieldDegrees 

    for a in $(seq 1 $totalNumberOfFrames); do
        base=$a.fits
        i=$framesForCommonReductionDir/$base
        out=$astroimadir/$base
        astfits $i --copy=0 --primaryimghdu -o $out
        for h in $(seq 1 $num_ccd); do
          im_layer=$astroimadir_layer/layer"$h"_"$base"
          astfits $im_layer --copy=1 -o $out
        done
    done
  fi
  rm -rf $astroimadir_layer
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
    echo -e "\n\tSex catalogs and scamp are already done\n"
else
  frameNames=()
  for a in $(seq 1 $totalNumberOfFrames); do
      frameNames+=("$a")
  done

  for ((i = 1; i <= numOfSextractorPlusScampIterations; i++)); do
    echo -e "\tSExtractor + scamp iteration $i"

    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_cpus" runSextractorOnImage {} $sexcfg $sexparam $sexconv $astroimadir $sexdir $saturationThreshold 
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
#entiredir_fullGrid=$BDIR/pointings_fullGrid
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

  printf "%s\n" "${imagesToWarp[@]}" | parallel -j "$num_cpus" warpImage {} $entiredir_smallGrid $ra $dec $coaddSizePx $pipelinePath
  echo done > $entiredone
fi



# Checking and removing bad astrometrised frames ------
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
badFilesWarningsFile=identifiedBadFrames_astrometry.txt
badFilesWarningsDone=$diagnosis_and_badFilesDir/done_badFrames_astrometry.txt
if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi
if [ -f $badFilesWarningsDone ]; then
    echo -e "\n\tBad astrometrised frames warning already done\n"
else
  scampXMLFilePath=$scampdir/scamp.xml
  python3 $pythonScriptsPath/checkForBadFrames_badAstrometry.py $diagnosis_and_badFilesDir $scampXMLFilePath $badFilesWarningsFile
  echo done > $badFilesWarningsDone
fi


echo -e "${GREEN} --- Compute and subtract Sky --- ${NOCOLOUR} \n"

#entiredir_smallGrid=$outputDir_small
#entiredir_fullGrid=$outputDir_full
noiseskydir=$BDIR/noise-sky_it1
noiseskydone=$noiseskydir/done_"$filter".txt

subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it1
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter".txt

#subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it1
#subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter".txt

echo -e "·Modelling the background for subtracting it"
imagesAreMasked=false
ringDir=$BDIR/ring


writeTimeOfStepToFile "Computing sky" $fileForTimeStamps
computeSky $entiredir_smallGrid $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'" 

# If we have not done it already (i.e. the modelling of the background selected has been a polynomial) we estimate de background as a constant for identifying bad frames
noiseskyctedir=$BDIR/noise-sky_it1_cte
noiseskyctedone=$noiseskyctedir/done_"$filter".txt
if [ "$MODEL_SKY_AS_CONSTANT" = false ]; then
  echo -e "\nModelling the background for the bad frame detection"
  computeSky $entiredir_smallGrid $noiseskyctedir $noiseskyctedone true $sky_estimation_method -1 false $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'" 
fi

# Checking and removing bad frames based on the background value ------
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
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
  for h in $(seq 1 $num_ccd); do
    diagnosis_and_badFilesDir_ccd=$diagnosis_and_badFilesDir/CCD"$h"
    if ! [ -d $diagnosis_and_badFilesDir_ccd ]; then mkdir $diagnosis_and_badFilesDir_ccd; fi
    python3 $pythonScriptsPath/checkForBadFrames_backgroundValueAndStd.py $tmpDir $framesForCommonReductionDir $airMassKeyWord $diagnosis_and_badFilesDir_ccd $h $dateHeaderKey
  done
  echo done > $badFilesWarningsDone
fi

echo -e "\n·Subtracting background"
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT

toleranceForMatching=3 #arcsec
sigmaForPLRegion=3 # Parameter for deciding the selection region (half-max-rad region)
export toleranceForMatching
export sigmaForPLRegion

# Parameters for performing the sigma clipping to the different samples in aststatistics
sigmaForStdSigclip=3
iterationsForStdSigClip=3
export sigmaForStdSigclip
export iterationsForStdSigClip

# Checking and removing bad frames based on the FWHM value ------
fwhmFolder=$BDIR/seeing_values
fwhmDone=$fwhmFolder/done.txt
badFilesWarningsFile=identifiedBadFrames_fwhm.txt
badFilesWarningsDone=$diagnosis_and_badFilesDir/done_fwhmValue.txt
if [ -f $badFilesWarningsDone ]; then
    echo -e "\nbadFiles warning already done\n"
else
  if ! [ -d $fwhmFolder ]; then mkdir $fwhmFolder; fi
  if  [ -f $fwhmDone ]; then
    echo -e "\nFWHM already computed"
  else
    imagesToFWHM=()
    for a in $(seq 1 $totalNumberOfFrames); do
      base="$a".fits
      imagesToFWHM+=("$base")
    done
    methodToUse="sextractor"
    printf "%s\n" "${imagesToFWHM[@]}" | parallel -j "$num_cpus" computeFWHMSingleFrame {} $subskySmallGrid_dir $fwhmFolder 0 $methodToUse $tileSize "NO"
    echo done > $fwhmDone
  fi
  
  for h in $(seq 1 $num_ccd); do  
    python3 $pythonScriptsPath/checkForBadFrames_fwhm.py $fwhmFolder $diagnosis_and_badFilesDir $badFilesWarningsFile $maximumSeeing $framesForCommonReductionDir $h $airMassKeyWord $dateHeaderKey $pixelScale
    if ! [ -f "$diagnosis_and_badFilesDir/CCD"$h"/$badFilesWarningFile" ]; then
      touch $diagnosis_and_badFilesDir/CCD"$h"/$badFilesWarningFile
    fi
  done
  echo done > $badFilesWarningsDone
fi

if [[ ("$produceCoaddPrephot" = "true") || ("$produceCoaddPrephot" = "True" )]]; then
  echo -e "${GREEN} --- Coadding before photometric calibration --- ${NOCOLOUR} \n"
  writeTimeOfStepToFile "Building coadd before photometry" $fileForTimeStamps
  iteration=1
  
  coaddDir=$BDIR/coadds-prephot
  coaddDone=$coaddDir/done.txt
  minRmsFileName=min_rms_prev_it$iteration.txt
  noisesky_prephot=$BDIR/noise-sky_prephot
  noisesky_prephotdone=$noisesky_prephot/done_$filter.txt
  if ! [ -d $noisesky_prephot ]; then mkdir $noisesky_prephot; fi
  if [ -f $coaddDone ]; then
    echo -e "\n Coadd pre-photometry already done\n"
  else
    imagesAreMasked=false
    computeSky $subskySmallGrid_dir $noisesky_prephot $noisesky_prephotdone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'" 

    subskyfullGrid_dir=$BDIR/sub-sky-fullGrid_it1
    subskyfullGridDone=$subskyfullGrid_dir/done.txt
    if ! [ -d $subskyfullGrid_dir ]; then mkdir $subskyfullGrid_dir; fi
    smallGridtoFullGrid $subskySmallGrid_dir $subskyfullGrid_dir $subskyfullGridDone $coaddSizePx $ra $dec

    #rejectedFramesDir=$BDIR/rejectedFrames_prephot_it$iteration
    #echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"
    #diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
    #if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
    #
    #prefixOfTheFilesToRemove="entirecamera_"
    #rejectedByAstrometry=identifiedBadFrames_astrometry.txt
    #removeBadFramesFromReduction $subskyfullGrid_dir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry $prefixOfTheFilesToRemove
    #removeBadFramesFromReduction $noisesky_prephot $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry $prefixOfTheFilesToRemove
    #rejectedByBackgroundFWHM=identifiedBadFrames_fwhm.txt
    #removeBadFramesFromReduction $subskyfullGrid_dir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM $prefixOfTheFilesToRemove
    #removeBadFramesFromReduction $noisesky_prephot $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM $prefixOfTheFilesToRemove

    python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $noisesky_prephot $DIR $iteration $minRmsFileName
    
    echo -e "\n ${GREEN} ---Masking outliers--- ${NOCOLOUR}"
    writeTimeOfStepToFile "Masking outliers" $fileForTimeStamps
    sigmaForStdSigclip=3
    clippingdir=$BDIR/clipping-outliers-prephot
    clippingdone=$clippingdir/done.txt
    buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $subskyfullGrid_dir $sigmaForStdSigclip

    subSkyNoOutliersPxDir=$BDIR/sub-sky-fullGrid_noOutliersPx_it$iteration
    subSkyNoOutliersPxDone=$subSkyNoOutliersPxDir/done.txt
    if ! [ -d $subSkyNoOutliersPxDir ]; then mkdir $subSkyNoOutliersPxDir; fi
    removeOutliersFromWeightedFrames $subSkyNoOutliersPxDone $subSkyNoOutliersPxDir $clippingdir $subskyfullGrid_dir
    ### Calculate the weights for the images based on the minimum rms ###
    echo -e "\n ${GREEN} ---Computing weights for the frames--- ${NOCOLOUR}"
    writeTimeOfStepToFile "Computing frame weights" $fileForTimeStamps
    wdir=$BDIR/weight-dir_prephot
    wonlydir=$BDIR/only-w-dir_prephot
    wdone=$wdir/done.txt
    wonlydone=$wonlydir/done.txt
    if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
    if ! [ -d $wdir ]; then mkdir $wdir; fi
    computeWeights $wdir $wdone $wonlydir $wonlydone $subSkyNoOutliersPxDir $noisesky_prephot $iteration $minRmsFileName
  
    coaddName=$coaddDir/"$objectName"_coadd_"$filter"_prephot_it$iteration.fits
    buildCoadd $coaddDir $coaddName $wdir $wonlydir $coaddDone

    maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
    if [ -f $maskName ]; then
      echo -e "\tThe mask of the weighted coadd is already done"
    else
      #If block scale is greater than 1, we apply the block
      if [ "$blockScale" -gt 1 ]; then
        astwarp $coaddName -h1 --scale=1/$blockScale --numthreads=$num_cpus-o $coaddDir/coadd_blocked.fits
        imToMask=$coaddDir/coaddBlocked.fits
      else
        imToMask=$coaddName
      fi
      #If a kernel exists in the configuration file, we apply it
      kernelFile=$CDIR/kernel.fits
      if [ -f $kernelFile ]; then
        astconvolve $imToMask --kernel=$kernelFile --domain=spatial --numthreads=$num_cpus -o $coaddDir/coadd_convolved.fits
        imToMask=$coaddDir/coadd_convolved.fits
      fi
      astnoisechisel $imToMask $noisechisel_param --numthreads=$num_cpus -o $coaddDir/mask_warped.fits
      if [ "$blockScale" -gt 1 ]; then 
        astwarp $coaddDir/mask_warped.fits --gridfile=$coaddName --numthreads=$num_cpus -o $coaddDir/mask_unwarped.fits
        astarithmetic $coaddDir/mask_unwarped.fits -h1 set-i i i 0 gt i isnotblank and 1 where float32 -q -o $maskName
        rm $coaddDir/mask_unwarped.fits  $coaddDir/mask_warped.fits $coaddDir/coadd_blocked.fits
      else
        mv $coaddDir/mask_warped.fits $maskName
      fi
      rm $coaddDir/coadd_convolved.fits 2>/dev/null
    fi

    exposuremapDir=$coaddDir/"$objectName"_exposureMap
    exposuremapdone=$coaddDir/done_exposureMap.txt
    computeExposureMap $wdir $exposuremapDir $exposuremapdone
  fi
fi

#### PHOTOMETRIC CALIBRATION  ####
echo -e "${ORANGE} ------ PHOTOMETRIC CALIBRATION ------ ${NOCOLOUR}\n"
writeTimeOfStepToFile "Photometric calibration" $fileForTimeStamps

### PARAMETERS ###
toleranceForMatching=5 #arcsec
sigmaForPLRegion=3 # Parameter for deciding the selection region (half-max-rad region)
export toleranceForMatching
export sigmaForPLRegion

# Parameters for performing the sigma clipping to the different samples in aststatistics
sigmaForStdSigclip=1
iterationsForStdSigClip=3
export sigmaForStdSigclip
export iterationsForStdSigClip

### PREPARING DECALS DATA FOR CALIBRATION ###
referenceImagesForMosaic=$entiredir_smallGrid
mosaicDir=$DIR/mosaic
selectedCalibrationStarsDir=$mosaicDir/automaticallySelectedStarsForCalibration
rangeUsedCalibrationDir=$mosaicDir/rangesUsedForCalibration
aperturePhotDir=$mosaicDir/aperturePhotometryCatalogues
mosaicDone=$mosaicDir/done_prep.txt

if (( $(echo "$sizeOfOurFieldDegrees > 0.5" | bc -l) )); then
  sizeOfBrick=3600
else
  sizeOfBrick=1000
fi

writeTimeOfStepToFile "Survey data processing" $fileForTimeStamps
prepareCalibrationData $surveyForPhotometry $referenceImagesForMosaic $aperturePhotDir $filter $ra $dec $mosaicDir $selectedCalibrationStarsDir $rangeUsedCalibrationDir \
                                            $pixelScale $sizeOfOurFieldDegrees $catName_gaia $surveyForSpectra $apertureUnits $folderWithTransmittances "$filterCorrectionCoeff" \
                                            $surveyCalibrationToGaiaBrightLimit $surveyCalibrationToGaiaFaintLimit $mosaicDone $sizeOfBrick

#exit 0
# Calibration of coadd prephot
if [[ ("$produceCoaddPrephot" = "true") || ("$produceCoaddPrephot" = "True" )]]; then
  num_ccd_old=$num_ccd
  export num_ccd_old
  writeTimeOfStepToFile "Computing calibration factor for coadd prephot" $fileForTimeStamps
  if ! [ -d "$BDIR/coaddForCalibration_it$iteration" ]; then mkdir "$BDIR/coaddForCalibration_it$iteration"; fi
  cp $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_prephot_it$iteration.fits $BDIR/coaddForCalibration_it$iteration/entirecamera_1_tmp.fits

  # Calibrate the coadd only in the high snr section 
  expMax=$(aststatistics $BDIR/coadds-prephot/exposureMap.fits --maximum -q)
  exp_fr=$(astarithmetic $expMax 0.5 x -q)
  astarithmetic $BDIR/coaddForCalibration_it$iteration/entirecamera_1_tmp.fits $BDIR/coadds-prephot/exposureMap.fits -g1 $exp_fr lt nan where --output=$BDIR/coaddForCalibration_it$iteration/entirecamera_1.fits
  rm $BDIR/coaddForCalibration_it$iteration/entirecamera_1_tmp.fits

  iteration=1
  alphatruedir=$BDIR/alpha-stars-true_coaddPrephot_it$iteration
  matchdir=$BDIR/match-decals-myData_coaddPrephot_it$iteration
  ourDataCatalogueDir=$BDIR/ourData-aperture-photometry_coaddPrephot_it$iteration
  prepareCalibrationCataloguePerFrame=$BDIR/survey-aperture-photometry_perBrick_coaddPrephot_it$iteration
  mycatdir=$BDIR/my-catalog-halfmaxradius_coaddPrephot_it$iteration
  calibratingMosaic=true
  imagesForCalibration=$BDIR/coaddForCalibration_it$iteration
  num_ccd=1
  export num_ccd
  computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $mycatdir $rangeUsedCalibrationDir \
                            $mosaicDir $alphatruedir $calibrationBrightLimitCoaddPrephot $calibrationFaintLimitCoaddPrephot $apertureUnits $numberOfApertureUnitsForCalibration $calibratingMosaic "'$noisechisel_param'"
  num_ccd=$num_ccd_old
  export num_ccd
fi



iteration=1
alphatruedir=$BDIR/alpha-stars-true_it$iteration
matchdir=$BDIR/match-decals-myData_it$iteration
ourDataCatalogueDir=$BDIR/ourData-aperture-photometry_it$iteration
prepareCalibrationCataloguePerFrame=$BDIR/survey-aperture-photometry_perBrick_it$iteration
mycatdir=$BDIR/my-catalog-halfmaxradius_it$iteration
imagesForCalibration=$subskySmallGrid_dir
calibratingMosaic=false

writeTimeOfStepToFile "Computing calibration factors for individual frames" $fileForTimeStamps
computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $mycatdir $rangeUsedCalibrationDir \
                          $mosaicDir $alphatruedir $calibrationBrightLimitIndividualFrames $calibrationFaintLimitIndividualFrames $apertureUnits $numberOfApertureUnitsForCalibration $calibratingMosaic "'$noisechisel_param'"

#exit 0

# Creating histogram with the number of stars used for the calibratino of each frame
diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi

numberOfStarsUsedInEachFrameDone=$diagnosis_and_badFilesDir/done_numOfStarsUsedForCalibrate.txt
if [ -f $numberOfStarsUsedInEachFrameDone ]; then
  echo -e "\nHistogram with the number of stars used for calibrating each plot already done"
else
  for h in $(seq 1 $num_ccd); do
    diagnosis_and_badFilesDir_ccd=$diagnosis_and_badFilesDir/CCD$h
    if ! [ -d $diagnosis_and_badFilesDir_ccd ]; then mkdir $diagnosis_and_badFilesDir_ccd; fi
    numberOfStarsUsedInEachFramePlot=$diagnosis_and_badFilesDir_ccd/numOfStarsUsedForCalibrationHist.png
    python3 $pythonScriptsPath/diagnosis_numOfStarsUsedInCalibration.py $alphatruedir/numberOfStarsUsedForCalibrate.txt $numberOfStarsUsedInEachFramePlot $h
  done
  echo done > $numberOfStarsUsedInEachFrameDone
fi


applyCommonCalibrationFactor=true
if ! [ -f $BDIR/commonCalibrationFactor_it$iteration.txt ]; then
  if [[ ("$applyCommonCalibrationFactor" = "true") || ("$applyCommonCalibrationFactor" = "True") ]]; then
    computeCommonCalibrationFactor $alphatruedir $iteration $objectName $BDIR
  fi
fi
#exit 0
# DIAGNOSIS PLOT
# Histogram of the background values on magnitudes / arcsec²
if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
  tmpDir=$BDIR/noise-sky_it1
else
  tmpDir=$noiseskyctedir
fi

backgroundBrightnessDone=$diagnosis_and_badFilesDir/backgroundBrightness.done
if [ -f $backgroundBrightnessDone ]; then
  echo -e "\nDiagnosis based on background brightness already done"
else
  badFilesBackgroundWarningsFile=identifiedBadFrames_backgroundBrightness.txt
  badFilesCalibrationFactorFile=identifiedBadFrames_calibrationFactor.txt
  for h in $(seq 1 $num_ccd); do
    python3 $pythonScriptsPath/diagnosis_normalisedBackgroundMagnitudesAndCalibrationFactorPlots.py $tmpDir $framesForCommonReductionDir $airMassKeyWord $alphatruedir \
                                                                                                  $pixelScale $diagnosis_and_badFilesDir $maximumBackgroundBrightness $badFilesBackgroundWarningsFile \
                                                                                                  $badFilesCalibrationFactorFile $applyCommonCalibrationFactor $BDIR/commonCalibrationFactor_it$iteration.txt $h $dateHeaderKey
  done
  echo "done" > $backgroundBrightnessDone
fi 
#exit 0

echo -e "\n ${GREEN} ---Applying calibration factors--- ${NOCOLOUR}"
alphatruedir=$BDIR/alpha-stars-true_it$iteration
photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir $iteration $applyCommonCalibrationFactor
#applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir

# DIAGNOSIS PLOTs ---------------------------------------------------
##On the very this point we will need the match tables to have proper names because if not astropy will fail misearbly
for a in $(seq 1 $totalNumberOfFrames); do
  input_file=$matchdir/match-entirecamera_"$a".fits.cat
  for h in $(seq 1 $num_ccd); do
    tmp_file=$matchdir/match-entirecamera_"$a"_ccd"$h".fits
    asttable $input_file -h$h --colmetadata=1,RA_CALIBRATED,deg,"Right ascension Survey" --colmetadata=2,DEC_CALIBRATED,deg,"Declination Survey" --colmetadata=3,RA_NONCALIBRATED,deg,"Right ascension data being reduced" --colmetadata=4,DEC_NONCALIBRATED,deg,"Declination data being reduced" --colmetadata=5,MAGNITUDE_CALIBRATED,mag,"Magnitude in Survey" --colmetadata=6,SUM_CALIBRATED,none,"Sum in Survey" --colmetadata=7,MAGNITUDE_NONCALIBRATED,mag,"Magnitude in data being reduced" --colmetadata=8,SUM_NONCALIBRATED,none,"Sum in data being reduced" -o$tmp_file
    astfits $tmp_file --copy=1 -o$matchdir/match-entirecamera_"$a".fits
    rm $tmp_file
  done
  rm $input_file
  mv $matchdir/match-entirecamera_"$a".fits $input_file
done

# Astrometry
astrometryPlotName=$diagnosis_and_badFilesDir/astrometry.png
if [ -f $astrometryPlotName ]; then
    echo -e "\nAstrometry diagnosis plot already done\n"
else
  for h in $(seq 1 $num_ccd); do
    astrometryPlotName=$diagnosis_and_badFilesDir/CCD"$h"/astrometry.png
    produceAstrometryCheckPlot $matchdir $pythonScriptsPath $astrometryPlotName $pixelScale $h
  done
fi

# Calibration
apertureFolder=$BDIR/my-catalog-halfmaxradius_it1
calibrationPlotName=calibrationPlot.png
calibrationPlot_done=$diagnosis_and_badFilesDir/done_calibrationPlot.txt
if [ -f $calibrationPlot_done ]; then
    echo -e "\nCalibration diagnosis plot already done\n"
else
  if [[ "$surveyForCalibration" == "SPECTRA" ]]; then
    dirWithReferenceCat=$mosaicDir/aperturePhotometryCatalogues
  else
    dirWithReferenceCat=$BDIR/survey-aperture-photometry_perBrick_it1
  fi
  produceCalibrationCheckPlot $BDIR/ourData-aperture-photometry_it1 $photCorrSmallGridDir $apertureFolder $dirWithReferenceCat \
                                  $pythonScriptsPath $calibrationPlotName $calibrationBrightLimitIndividualFrames $calibrationFaintLimitIndividualFrames $numberOfApertureUnitsForCalibration $diagnosis_and_badFilesDir $surveyForPhotometry $BDIR
 
  echo done > $calibrationPlot_done
fi


# Half-Max-Radius vs magnitude plots of our calibrated data

halfMaxRadiusVsMagnitudeOurDataDone=$diagnosis_and_badFilesDir/done_halfMaxRadVsMagPlots.txt
if [ -f $halfMaxRadiusVsMagnitudeOurDataDone ]; then
    echo -e "\nHalf max radius vs magnitude plots for our calibrated data already done"
else
  for h in $(seq 1 $num_ccd); do
    diagnosis_and_badFilesDir_ccd=$diagnosis_and_badFilesDir/CCD"$h"
    halfMaxRadiusVsMagnitudeOurDataDir=$diagnosis_and_badFilesDir_ccd/halfMaxRadVsMagPlots_ourData
    if ! [ -d $halfMaxRadiusVsMagnitudeOurDataDir ]; then mkdir $halfMaxRadiusVsMagnitudeOurDataDir; fi

    produceHalfMaxRadVsMagForOurData $photCorrSmallGridDir $halfMaxRadiusVsMagnitudeOurDataDir $catdir/"$objectName"_Gaia_DR3.fits $toleranceForMatching $pythonScriptsPath $num_cpus 30
  done
  echo done > $halfMaxRadiusVsMagnitudeOurDataDone

fi

num_ccd_old=$num_ccd
num_ccd=1
if [[ ("$produceCoaddPrephot" = "true") || ("$produceCoaddPrephot" = "True" )]]; then
  iteration=1
  photCorrPrePhotDir=$BDIR/photCorr-coaddPrephot-dir_it$iteration
  rm $photCorrPrePhotDir/*
  rm $BDIR/coaddForCalibration_it$iteration/*

  cp $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_prephot_it$iteration.fits $BDIR/coaddForCalibration_it$iteration/entirecamera_1.fits
  alphatruedir=$BDIR/alpha-stars-true_coaddPrephot_it$iteration
  applyCalibrationFactors $BDIR/coaddForCalibration_it$iteration $alphatruedir $photCorrPrePhotDir $iteration False

  coaddPrephotDir=$BDIR/coadds-prephot
  coaddPrephotCalibratedName=$coaddPrephotDir/"$objectName"_prephot_calibrated.fits
  if [ ! -f "$coaddPrephotCalibratedName" ]; then
    cp $BDIR/photCorr-coaddPrephot-dir_it$iteration/entirecamera_1.fits $coaddPrephotCalibratedName
  fi

  # Compute surface brightness limit
  sblimitFile=$coaddPrephotDir/"$objectName"_"$filter"_sblimit.txt
  exposuremapName=$coaddPrephotDir/exposureMap.fits
  if [ -f  $sblimitFile ]; then
    echo -e "\n\tSurface brightness limit for coadd already measured\n"
    surfaceBrightnessLimit=$( awk '/Limiting magnitude/ { print $NF }' $sblimitFile )
  else
    maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
    surfaceBrightnessLimit=$( limitingSurfaceBrightness $coaddPrephotCalibratedName $maskName $exposuremapName $coaddPrephotDir $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )
  fi


  times=($(getInitialMidAndFinalFrameTimes $INDIR))
  initialTime=$( TZ=UTC  date -d @"${times[0]}" "+%Y-%m-%d_%H:%M:%S")
  meanTime=$( TZ=UTC  date -d @"${times[1]}" "+%Y-%m-%d_%H:%M:%S")
  finalTime=$( TZ=UTC  date -d @"${times[2]}" "+%Y-%m-%d_%H:%M:%S")


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

  numberOfFramesCombined=$(ls $BDIR/weight-dir_prephot/*.fits | wc -l)
  values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$lowerVignettingThreshold" "$upperVignettingThreshold" "$saturationThreshold" "$surveyForPhotometry" "$calibrationBrightLimitCoaddPrephot" "$calibrationFaintLimitCoaddPrephot" "$RUNNING_FLAT" "$halfWindowSize" "$surfaceBrightnessLimit")
  comments=("" "" "" "" "" "" "" "" "" "" "" "" "" "Running flat built with +-N frames" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")

  astfits $coaddPrephotCalibratedName --write=/,"Pipeline information"
  addkeywords $coaddPrephotCalibratedName keyWords values comments
fi # ------------------------------------------------------

num_ccd=$num_ccd_old
# ---------------------------------------------------

echo -e "\n${ORANGE} ------ STD WEIGHT COMBINATION ------ ${NOCOLOUR}\n"
# Compute rms and of the photometrized frames
noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done.txt
# Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
computeSky $photCorrSmallGridDir $noiseskydir $noiseskydone true $sky_estimation_method -1 false $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'"  

##STACKING: In a common pipeline, the process will be the following:
# 1) Put each small grid into the mosaic grid. This includes each CCD as well
# 2) Remove bad frames
# 3) Compute and mask the outliers usign the sigma-clipped median of the frames
# 4) Compute the weights of the MASKED frames
# 5) Build the coadd
#This process in the multi detector pipeline was done layer by layer, meaning that each stacking process
# was done with (N CCDs x M frames) images of size (coaddSizePx x coaddSizePx). This, of course, is not optimal
# when having a huge number of CCDs. The process is changed now to do the following:
# 1) In the process of building the fullGrid frames, we construct 2 layers per frame:
#    - The photometrized fullGrid frame were all the CCDs are placed in the common grid 
#    - Same, but the CCDs are placed multiplied by their weights (1/rms^2)
# 2) The clipping process is applied in the first layer, and the outliers are removed in the second layer
# 3) The second layer is the used to build the coadd directly
#In addition, instead of creating a weight frame to sum it, we will compute the sum of weights
minRmsFileName=min_rms_it1.txt
python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $noiseskydir $DIR $iteration min_rms_it1.txt

photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
photCorrFullGridDone=$photCorrFullGridDir/done.txt
identifiedBadDetectors=$CDIR/identifiedBadDetectors.txt #We can now provide a list of bad detectors to blank them
if ! [ -d $photCorrFullGridDir ]; then mkdir $photCorrFullGridDir; fi
smallGridtoFullGridAndWeight $photCorrSmallGridDir $photCorrFullGridDir $photCorrFullGridDone $coaddSizePx $ra $dec $noiseskydir $minRmsFileName $iteration $identifiedBadDetectors 


echo -e "\n·Removing bad frames"

diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
rejectedFramesDir=$BDIR/rejectedFrames
if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"

rejectedByAstrometry=identifiedBadFrames_astrometry.txt
removeBadFramesFromReduction $photCorrFullGridDir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry
removeBadFramesFromReduction $noiseskydir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry


# Store the minimum standard deviation of the frames in order to compute the weights



echo -e "\n ${GREEN} ---Masking outliers--- ${NOCOLOUR}"
writeTimeOfStepToFile "Masking outliers" $fileForTimeStamps
# Remove outliers before the final coadd by using sigclip-median and sigclip-std
# This is particularly important to remove cosmic rays

sigmaForStdSigclip=3
clippingdir=$BDIR/clipping-outliers
clippingdone=$clippingdir/done.txt
buildUpperAndLowerLimitsForOutliersNew $clippingdir $clippingdone $photCorrFullGridDir $sigmaForStdSigclip

mowdir=$BDIR/photCorrFullGrid-dir-no-outliers_it$iteration
mowdone=$mowdir/done.txt
if ! [ -d $mowdir ]; then mkdir $mowdir; fi
removeOutliersFromWeightedFramesNew $mowdone $mowdir $clippingdir $photCorrFullGridDir

#### Calculate the weights for the images based on the minimum rms ###
#echo -e "\n ${GREEN} ---Computing weights for the frames--- ${NOCOLOUR}"
#writeTimeOfStepToFile "Computing frame weights" $fileForTimeStamps
#
#wdir=$BDIR/weight-dir
#wdone=$wdir/done.txt
#if ! [ -d $wdir ]; then mkdir $wdir; fi
#
#wonlydir=$BDIR/only-w-dir
#wonlydone=$wonlydir/done.txt
#if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
## We provide the fullGrid because we are going to combine then now
#computeWeights $wdir $wdone $wonlydir $wonlydone $mowdir $noiseskydir $iteration $minRmsFileName





echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
writeTimeOfStepToFile "Building coadd" $fileForTimeStamps


echo -e "\n·Building coadd"
coaddDir=$BDIR/coadds
coaddDone=$coaddDir/done.txt
coaddName=$coaddDir/"$objectName"_coadd_"$filter"_it"$iteration".fits
buildCoaddNew $coaddDir $coaddName $mowdir $coaddDone
maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  #If block scale is greater than 1, we apply the block
  if [ "$blockScale" -gt 1 ]; then
    astwarp $coaddName -h1 --scale=1/$blockScale --numthreads=$num_cpus -o $coaddDir/coadd_blocked.fits
    imToMask=$coaddDir/coadd_blocked.fits
  else
    imToMask=$coaddName
  fi
  #If a kernel exists in the configuration file, we apply it
  kernelFile=$CDIR/kernel.fits
  if [ -f $kernelFile ]; then
    astconvolve $imToMask --kernel=$kernelFile --domain=spatial --numthreads=$num_cpus -o $coaddDir/coadd_convolved.fits
    imToMask=$coaddDir/coadd_convolved.fits
  fi
  astnoisechisel $imToMask $noisechisel_param --numthreads=$num_cpus -o $coaddDir/mask_warped.fits
  if [ "$blockScale" -gt 1 ]; then 
    astwarp $coaddDir/mask_warped.fits --gridfile=$coaddName --numthreads=$num_cpus -o $coaddDir/mask_unwarped.fits
    astarithmetic $coaddDir/mask_unwarped.fits -h1 set-i i i 0 gt 1 where float32 -q -o $maskName
    rm $coaddDir/mask_unwarped.fits  $coaddDir/mask_warped.fits $coaddDir/coadd_blocked.fits
  else
    mv $coaddDir/mask_warped.fits $maskName
  fi
  rm $coaddDir/coadd_convolved.fits 2>/dev/null
fi

#astnoisechisel with the current parameters might fail due to long tilesize. I'm gonna make 2 checks to see if it fails, decreasing in steps of 5 in tilesize
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  #We assume that if this works for this iteration, then the next one will need at least same parameters
  tileSize=20
  noisechisel_param="--tilesize=$tileSize,$tileSize \
                    --minskyfrac=0.9 \
                     --meanmedqdiff=0.01 \
                     --snthresh=5.2 \
                     --detgrowquant=0.7 \
                     --detgrowmaxholesize=1000 \
                     --rawoutput"
  astnoisechisel $coaddName $noisechisel_param  -o $maskName
fi
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  #We assume that if this works for this iteration, then the next one will need at least same parameters
  tileSize=$((tileSize - 5))
  noisechisel_param="--tilesize=$tileSize,$tileSize \
                    --minskyfrac=0.9 \
                     --meanmedqdiff=0.01 \
                     --snthresh=5.2 \
                     --detgrowquant=0.7 \
                     --detgrowmaxholesize=1000 \
                     --rawoutput"
  astnoisechisel $coaddName $noisechisel_param --numthreads=$num_cpus  -o $maskName
fi

if ! [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd could not be generated. Please, check manually."
  exit 47
fi
exposuremapDir=$coaddDir/"$objectName"_exposureMap
exposuremapdone=$coaddDir/done_exposureMap_eff.txt
exposureMapName=$coaddDir/"$objectName"_expMap_eff_"$filter"_it"$iteration".fits
computeExposureMapNew $mowdir $exposuremapDir $exposuremapdone $exposureMapName
#
exposuremapdone=$coaddDir/done_exposureMap.txt
exposureMapName=$coaddDir/"$objectName"_expMap_"$filter"_it"$iteration".fits
computeExposureMapNew $photCorrFullGridDir $exposuremapDir $exposuremapdone $exposureMapName

#Compute surface brightness limit
sblimitFile=$coaddDir/"$objectName"_"$filter"_sblimit.txt

if [ -f  $sblimitFile ]; then
    echo -e "\n\tSurface brightness limit for coadd already measured\n"
else
    surfaceBrightnessLimit=$( limitingSurfaceBrightness $coaddName $maskName $exposureMapName $coaddDir $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )
fi


times=($(getInitialMidAndFinalFrameTimes $INDIR $dateHeaderKey))
initialTime=$( date -d @"${times[0]}" "+%Y-%m-%d_%H:%M:%S")
meanTime=$( date -d @"${times[1]}" "+%Y-%m-%d_%H:%M:%S")
finalTime=$( date -d @"${times[2]}" "+%Y-%m-%d_%H:%M:%S")

keyWords=("FRAMES_COMBINED" \
          "NUMBER_OF_DIFFERENT_NIGHTS" \
          "INITIAL_DATE_OBS"
          "MEAN_DATE_OBS"
          "FINAL_DATE_OBS"
          "FILTER" \
          "SATURATION_THRESHOLD" \
          "CALIBRATION_BRIGHTLIMIT" \
          "CALIBRATION_FAINTLIMIT" \
          "RUNNING_FLAT" \
          "WINDOW_SIZE" \
          "STD_FOR_BAD_FRAMES" \
          "SURFACE_BRIGHTNESS_LIMIT")

numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames" "$surfaceBrightnessLimit")
comments=("" "" "" "" "" "" "" "" "" "" "" "Num. of tandard deviations used for rejection" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")
astfits $coaddName --write=/,"Pipeline information"
addkeywords $coaddName keyWords values comments

halfMaxRadForCoaddName=$diagnosis_and_badFilesDir/coadd_it1.png
if [ -f $halfMaxRadForCoaddName ]; then
  echo -e "\tThe Half-Max-Rad vs Magnitude has been already generate for the coadd"
else
  produceHalfMaxRadVsMagForSingleImage $coaddName $diagnosis_and_badFilesDir $catdir/"$objectName"_Gaia_eDR3.fits $toleranceForMatching $pythonScriptsPath "coadd_it1" 100 NO
fi
##To avoid problems, we re-name the pythonScriptsPath and toleranceForMatching because it is creating problems
#pythonScriptsPath=$pipelinePath/pipelineScripts
#toleranceForMatching=2 
#export pythonScriptsPath
#export toleranceForMatching

writeTimeOfStepToFile "Producing frames with coadd subtracted" $fileForTimeStamps
framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted
framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi
if [ -f $framesWithCoaddSubtractedDone ]; then
    echo -e "\n\tFrames with coadd subtracted already generated\n"
else
  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_"$filter"_it"$iteration".fits
  coadd_av=$coaddDir/"$objectName"_coadd_it"$iteration"_average.fits
  gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
  if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
    astarithmetic $(ls $photCorrFullGridDir/*.fits) -g1 $(ls $photCorrFullGridDir/*.fits | wc -l ) mean --writeall -o$coadd_av
  else
    astarithmetic $(ls $photCorrFullGridDir/*.fits) -g1 $(ls $photCorrFullGridDir/*.fits | wc -l ) mean -o$coadd_av
  fi    
  
  subtractCoaddToFramesNew $photCorrFullGridDir $coadd_av $framesWithCoaddSubtractedDir
  if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
    astarithmetic $framesWithCoaddSubtractedDir/*.fits -g1 $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l ) sum --writeall -o$sumMosaicAfterCoaddSubtraction
  else
    astarithmetic $framesWithCoaddSubtractedDir/*.fits -g1 $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l ) sum -o$sumMosaicAfterCoaddSubtraction
  fi
  echo done > $framesWithCoaddSubtractedDone 
fi



# # Remove intermediate folders to save some space
find $BDIR/noise-sky_it1 -type f -name '*.fits' -exec rm {} \;
find $BDIR/noise-sky-after-photometry_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/noise-sky_forPlaneCoadd_it1 -type f ! -name 'done*' -exec rm {} \;

#find $BDIR/sub-sky-fullGrid_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/sub-sky-smallGrid_it1 -type f ! -name 'done*' -exec rm {} \;
#find $BDIR/photCorrFullGrid-dir_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrSmallGrid-dir_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/photCorrFullGrid-dir-no-outliers_it1 -type f ! -name 'done*' -exec rm {} \;

find $BDIR/my-catalog-halfmaxradius_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/match-decals-myData_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/ourData-catalogs-apertures_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/framesWithCoaddSubtracted -type f ! -name 'done*' -exec rm {} \;

#find $BDIR/weight-dir-no-outliers -type f ! -name 'done*' -exec rm {} \;
#find $BDIR/weight-dir-no-outliers_plane_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/weight-dir -type f ! -name 'done*' -exec rm {} \;
find $BDIR/weight-dir_plane_it1 -type f ! -name 'done*' -exec rm {} \;

#find $BDIR/only-weight-dir-no-outliers -type f ! -name 'done*' -exec rm {} \;
#find $BDIR/only-weight-dir-no-outliers_plane_it1 -type f ! -name 'done*' -exec rm {} \;
find $BDIR/only-w-dir -type f ! -name 'done*' -exec rm {} \;
find $BDIR/only-w-dir_plane_it1 -type f ! -name 'done*' -exec rm {} \;
if [[ ("$produceCoaddPrephot" = "true") || ("$produceCoaddPrephot" = "True" )]]; then
  find $BDIR/weight-dir_prephot -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/only-w-dir_prephot -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/noise-sky_prephot -type f ! -name 'done*' -exec rm {} \;
fi
## This code is used for manually adding the user-defined masks to the mask from the coadd
# First we save the original mask that noisechisel produces
cp $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_mask.fits $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_mask_copy.fits

cp $BDIR/coadds/"$objectName"_coadd_"$filter"_mask.fits $BDIR/coadds/"$objectName"_coadd_"$filter"_mask_copy.fits

# Then we apply the user-defined masks
valueToPut=1
read -r -a maskArray <<< "$maskParams"
for ((i=0; i<${#maskArray[@]}; i+=5)); do
	ra="${maskArray[i]}"
	dec="${maskArray[i+1]}"
	r="${maskArray[i+2]}"
	axisRatio="${maskArray[i+3]}"
	pa="${maskArray[i+4]}"

	python3 $pythonScriptsPath/manualMaskRegionFromWCSArea.py $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_mask.fits $valueToPut $ra $dec $r $axisRatio $pa
	python3 $pythonScriptsPath/manualMaskRegionFromWCSArea.py $BDIR/coadds/"$objectName"_coadd_"$filter"_mask.fits $valueToPut $ra $dec $r $axisRatio $pa
done

#If a mask already exists on the configuration file, we combine both masks
if [ -f $CDIR/mask.fits ]; then
  cp $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_mask.fits $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_mask_copy.fits
  astarithmetic $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_mask_copy.fits $CDIR/mask.fits -g1 1 eq 1 where -q -o $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_mask.fits
  cp $BDIR/coadds/"$objectName"_coadd_"$filter"_mask.fits $BDIR/coadds/"$objectName"_coadd_"$filter"_mask_copy.fits
  astarithmetic $BDIR/coadds/"$objectName"_coadd_"$filter"_mask_copy.fits $CDIR/mask.fits -g1 1 eq 1 where -q -o $BDIR/coadds/"$objectName"_coadd_"$filter"_mask.fits
fi 
####### ITERATION 2 ######

iteration=2
entiredir_smallGrid=$BDIR/pointings_smallGrid

# We mask the pointings in order to measure (before photometric calibration) the sky accurately

#######################

#
#
########################
#entiredir_smallGrid=$outputDir_small

##We will conserve photCorrFullGrid-dir_it1 in order to make the crop 
smallPointings_maskedDir=$BDIR/pointings_smallGrid_masked_it$iteration
maskedPointingsDone=$smallPointings_maskedDir/done_.txt


maskPointings $entiredir_smallGrid $smallPointings_maskedDir $maskedPointingsDone $maskName $entiredir_smallGrid


noiseskydir=$BDIR/noise-sky_it$iteration
noiseskydone=$noiseskydir/done_"$filter".txt

subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it$iteration
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter".txt

##subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it$iteration
##subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter".txt


# compute sky with frames masked with global mask
imagesAreMasked=true
sky_estimation_method=wholeImage #If we trust the mask, we can use the full image
computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'"
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT

if [[ ("$produceCoaddPrephot" = "true") || ("$produceCoaddPrephot" = "True" )]]; then
  # MASK FROM THE COADD PREPHOT
  coaddDir=$BDIR/coadds-prephot
  maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
  smallPointings_maskedDir=$BDIR/pointings_smallGrid_maskedPrephot_it$iteration
  maskedPointingsDone=$smallPointings_maskedDir/done_.txt
  maskPointings $entiredir_smallGrid $smallPointings_maskedDir $maskedPointingsDone $maskName $entiredir_smallGrid

  noiseskydir=$BDIR/noise-sky_maskPrephot_it$iteration
  noiseskydone=$noiseskydir/done_"$filter"_ccd"$h".txt
  imagesAreMasked=true
  computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'"

  subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_maskPrephot_it$iteration
  subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter".txt
  subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT
  
  echo -e "${GREEN} --- Coadding before photometric calibration --- ${NOCOLOUR} \n"
  writeTimeOfStepToFile "Building coadd before photometry" $fileForTimeStamps
  iteration=2
  coaddDir=$BDIR/coadds-prephot_it"$iteration"
  coaddDone=$coaddDir/done.txt
  minRmsFileName=min_rms_prev_it$iteration.txt
  noisesky_prephot=$BDIR/noise-sky_prephot_it"$iteration"
  noisesky_prephotdone=$noisesky_prephot/done_$filter.txt
  if ! [ -d $noisesky_prephot ]; then mkdir $noisesky_prephot; fi
  if [ -f $coaddDone ]; then
    echo -e "\n Coadd pre-photometry already done\n"
  else
    maskName=$BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_mask.fits
    subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_maskPrephot_it$iteration
    subSkyPointings_maskedDir=$BDIR/sub-sky-smallGrid_maskPrephot_masked_it$iteration
    maskedPointingsDone=$subSkyPointings_maskedDir/done_.txt
    maskPointings $subskySmallGrid_dir $subSkyPointings_maskedDir $maskedPointingsDone $maskName $entiredir_smallGrid

    imagesAreMasked=true
    computeSky $subSkyPointings_maskedDir $noisesky_prephot $noisesky_prephotdone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $ringDir $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'"

    subskyfullGrid_dir=$BDIR/sub-sky-fullGrid_maskPrephot_it"$iteration"
    subskyfullGridDone=$subskyfullGrid_dir/done.txt
    if ! [ -d $subskyfullGrid_dir ]; then mkdir $subskyfullGrid_dir; fi
    smallGridtoFullGrid $subskySmallGrid_dir $subskyfullGrid_dir $subskyfullGridDone $coaddSizePx $ra $dec

    #rejectedFramesDir=$BDIR/rejectedFrames_prephot_it$iteration
    #echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"
    #diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
    #if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
    #
    #prefixOfTheFilesToRemove="entirecamera_"
    #rejectedByAstrometry=identifiedBadFrames_astrometry.txt
    #removeBadFramesFromReduction $subskyfullGrid_dir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry $prefixOfTheFilesToRemove
    #removeBadFramesFromReduction $noisesky_prephot $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry $prefixOfTheFilesToRemove
    #rejectedByBackgroundFWHM=identifiedBadFrames_fwhm.txt
    #removeBadFramesFromReduction $subskyfullGrid_dir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM $prefixOfTheFilesToRemove
    #removeBadFramesFromReduction $noisesky_prephot $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByBackgroundFWHM $prefixOfTheFilesToRemove

    python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $noisesky_prephot $DIR $iteration $minRmsFileName
    
    echo -e "\n ${GREEN} ---Masking outliers--- ${NOCOLOUR}"
    writeTimeOfStepToFile "Masking outliers" $fileForTimeStamps
    sigmaForStdSigclip=3
    clippingdir=$BDIR/clipping-outliers-prephot_it"$iteration"
    clippingdone=$clippingdir/done.txt
    buildUpperAndLowerLimitsForOutliers $clippingdir $clippingdone $subskyfullGrid_dir $sigmaForStdSigclip

    subSkyNoOutliersPxDir=$BDIR/sub-sky-fullGrid_noOutliersPx_it$iteration
    subSkyNoOutliersPxDone=$subSkyNoOutliersPxDir/done.txt
    if ! [ -d $subSkyNoOutliersPxDir ]; then mkdir $subSkyNoOutliersPxDir; fi
    removeOutliersFromWeightedFrames $subSkyNoOutliersPxDone $subSkyNoOutliersPxDir $clippingdir $subskyfullGrid_dir
    ### Calculate the weights for the images based on the minimum rms ###
    echo -e "\n ${GREEN} ---Computing weights for the frames--- ${NOCOLOUR}"
    writeTimeOfStepToFile "Computing frame weights" $fileForTimeStamps
    wdir=$BDIR/weight-dir_prephot_it"$iteration"
    wonlydir=$BDIR/only-w-dir_prephot_it"$iteration"
    wdone=$wdir/done.txt
    wonlydone=$wonlydir/done.txt
    if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
    if ! [ -d $wdir ]; then mkdir $wdir; fi
    computeWeights $wdir $wdone $wonlydir $wonlydone $subSkyNoOutliersPxDir $noisesky_prephot $iteration $minRmsFileName
  
    coaddName=$coaddDir/"$objectName"_coadd_"$filter"_prephot_it$iteration.fits
    buildCoadd $coaddDir $coaddName $wdir $wonlydir $coaddDone

    maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
    if [ -f $maskName ]; then
      echo -e "\tThe mask of the weighted coadd is already done"
    else
      #If block scale is greater than 1, we apply the block
      if [ "$blockScale" -gt 1 ]; then
        astwarp $coaddName -h1 --scale=1/$blockScale --numthreads=$num_cpus-o $coaddDir/coadd_blocked.fits
        imToMask=$coaddDir/coaddBlocked.fits
      else
        imToMask=$coaddName
      fi
      #If a kernel exists in the configuration file, we apply it
      kernelFile=$CDIR/kernel.fits
      if [ -f $kernelFile ]; then
        astconvolve $imToMask --kernel=$kernelFile --domain=spatial --numthreads=$num_cpus -o $coaddDir/coadd_convolved.fits
        imToMask=$coaddDir/coadd_convolved.fits
      fi
      astnoisechisel $imToMask $noisechisel_param --numthreads=$num_cpus -o $coaddDir/mask_warped.fits
      if [ "$blockScale" -gt 1 ]; then 
        astwarp $coaddDir/mask_warped.fits --gridfile=$coaddName --numthreads=$num_cpus -o $coaddDir/mask_unwarped.fits
        astarithmetic $coaddDir/mask_unwarped.fits -h1 set-i i i 0 gt i isnotblank and 1 where float32 -q -o $maskName
        rm $coaddDir/mask_unwarped.fits  $coaddDir/mask_warped.fits $coaddDir/coadd_blocked.fits
      else
        mv $coaddDir/mask_warped.fits $maskName
      fi
      rm $coaddDir/coadd_convolved.fits 2>/dev/null
    fi

    exposuremapDir=$coaddDir/"$objectName"_exposureMap
    exposuremapdone=$coaddDir/done_exposureMap.txt
    computeExposureMap $wdir $exposuremapDir $exposuremapdone
  fi
  num_ccd_old=$num_ccd
  export num_ccd_old
# Calibration of coadd prephot
  writeTimeOfStepToFile "Computing calibration factor for coadd prephot" $fileForTimeStamps
  if ! [ -d "$BDIR/coaddForCalibration_it$iteration" ]; then mkdir "$BDIR/coaddForCalibration_it$iteration"; fi
  cp $BDIR/coadds-prephot/"$objectName"_coadd_"$filter"_prephot_it$iteration.fits $BDIR/coaddForCalibration_it$iteration/entirecamera_1_tmp.fits

  # Calibrate the coadd only in the high snr section 
  expMax=$(aststatistics $BDIR/coadds-prephot/exposureMap.fits --maximum -q)
  exp_fr=$(astarithmetic $expMax 0.5 x -q)
  astarithmetic $BDIR/coaddForCalibration_it$iteration/entirecamera_1_tmp.fits $BDIR/coadds-prephot/exposureMap.fits -g1 $exp_fr lt nan where --output=$BDIR/coaddForCalibration_it$iteration/entirecamera_1.fits
  rm $BDIR/coaddForCalibration_it$iteration/entirecamera_1_tmp.fits

  iteration=1
  alphatruedir=$BDIR/alpha-stars-true_coaddPrephot_it$iteration
  matchdir=$BDIR/match-decals-myData_coaddPrephot_it$iteration
  ourDataCatalogueDir=$BDIR/ourData-aperture-photometry_coaddPrephot_it$iteration
  prepareCalibrationCataloguePerFrame=$BDIR/survey-aperture-photometry_perBrick_coaddPrephot_it$iteration
  mycatdir=$BDIR/my-catalog-halfmaxradius_coaddPrephot_it$iteration
  calibratingMosaic=true
  imagesForCalibration=$BDIR/coaddForCalibration_it$iteration
  num_ccd=1
  export num_ccd
  computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $mycatdir $rangeUsedCalibrationDir \
                            $mosaicDir $alphatruedir $calibrationBrightLimitCoaddPrephot $calibrationFaintLimitCoaddPrephot $apertureUnits $numberOfApertureUnitsForCalibration $calibratingMosaic "'$noisechisel_param'" "'$maskParams'"

  num_ccd=$num_ccd_old
  export num_ccd 
fi


writeTimeOfStepToFile "Computing calibration factors for individual frames" $fileForTimeStamps
iteration=2
alphatruedir=$BDIR/alpha-stars-true_it$iteration
matchdir=$BDIR/match-decals-myData_it$iteration
ourDataCatalogueDir=$BDIR/ourData-aperture-photometry_it$iteration
prepareCalibrationCataloguePerFrame=$BDIR/survey-aperture-photometry_perBrick_it$iteration
mycatdir=$BDIR/my-catalog-halfmaxradius_it$iteration
imagesForCalibration=$subskySmallGrid_dir
calibratingMosaic=false

computeCalibrationFactors $surveyForPhotometry $iteration $imagesForCalibration $selectedCalibrationStarsDir $matchdir $ourDataCatalogueDir $prepareCalibrationCataloguePerFrame $mycatdir $rangeUsedCalibrationDir \
                          $mosaicDir $alphatruedir $calibrationBrightLimitIndividualFrames $calibrationFaintLimitIndividualFrames $apertureUnits $numberOfApertureUnitsForCalibration $calibratingMosaic "'$noisechisel_param'" "'$maskParams'"




applyCommonCalibrationFactor=true
if ! [ -f $BDIR/commonCalibrationFactor_it$iteration.txt ]; then
  if [[ ("$applyCommonCalibrationFactor" = "true") || ("$applyCommonCalibrationFactor" = "True") ]]; then
    computeCommonCalibrationFactor $alphatruedir $iteration $objectName $BDIR
  fi
fi

alphatruedir=$BDIR/alpha-stars-true_it$iteration
photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir $iteration $applyCommonCalibrationFactor
echo -e "${GREEN} --- Correct difference in gain --- ${NOCOLOUR} \n"
gainCorrectionFile=$BDIR/gainCorrectionParameters_it$iteration.txt
ccd_ref=1 #Hardcoded, might change. Based on $diagnosis_and_badFilesDir/backgroundCCDcomparison.png
cfactorfile=$BDIR/commonCalibrationFactor_it$iteration.txt
backgroundCCDsDone=$diagnosis_and_badFilesDir/done_backgroundCCDs.txt
if [ -f $backgroundCCDsDone ]; then
  echo -e "\nPlot of background per CCD already done"
else
  python3 $pythonScriptsPath/diagnosis_backgroundBrightnessPerCCD.py $noiseskydir $diagnosis_and_badFilesDir $cfactorfile $pixelScale $dateHeaderKey $airMassKeyWord $framesForCommonReductionDir $num_ccd $ccd_ref $gainCorrectionFile
  
  echo "done" > $backgroundCCDsDone
fi 


  #Now that we have photometrically corrected the images, we apply a refining of the photometry based on the relative difference 
  # between background in nano-maggies: we expect that, after photometric correction, background should be approximately the same in between detectors
  # We use the sky measured before photometric correction, and the common calibration factor, to compute the gain correction

gaincordir=$BDIR/gain-corrected
gaincordone=$gaincordir/done.txt


if ! [ -d $gaincordir ]; then mkdir $gaincordir; fi
if [ -f $gaincordone ]; then
  echo -e "\nMulti-layer windows already gain corrected"
else
  frameNames=()
  for a in $entiredir_smallGrid/*.fits; do
    frameNames+=("$(basename $a)")
  done
  printf "%s\n" "${frameNames[@]}" | parallel -j "$num_parallel" gainCorrection {} $photCorrSmallGridDir $gainCorrectionFile $gaincordir $ccd_ref
  echo done > $gaincordone
fi
#smallPointings_maskedDir=$BDIR/pointings_smallGrid_masked_gain
#maskedPointingsDone=$smallPointings_maskedDir/done_.txt
#
#
#maskPointings $gaincordir $smallPointings_maskedDir $maskedPointingsDone $maskName $entiredir_smallGrid
#
#
#noiseskydir=$BDIR/noise-sky_gain
#noiseskydone=$noiseskydir/done_"$filter".txt
#
#
#
###subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it$iteration
###subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter".txt
#
#
## compute sky with frames masked with global mask
#imagesAreMasked=true
#sky_estimation_method=wholeImage #If we trust the mask, we can use the full image
#computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'"
#
#backgroundCCDsDone=$diagnosis_and_badFilesDir/done_backgroundCCDs_gaintest.txt
#if [ -f $backgroundCCDsDone ]; then
#  echo -e "\nPlot of background per CCD already done"
#else
#  python3 $pythonScriptsPath/diagnosis_backgroundBrightnessPerCCD.py $noiseskydir $diagnosis_and_badFilesDir $BDIR/commonCalibrationFactor_it$iteration.txt $pixelScale $dateHeaderKey $airMassKeyWord $framesForCommonReductionDir $num_ccd
#  
#  echo "done" > $backgroundCCDsDone
#fi 
#exit 0
#applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir
#diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles_it$iteration
#if ! [ -d $diagnosis_and_badFilesDir ]; then mkdir $diagnosis_and_badFilesDir; fi
#for h in $(seq 1 $num_ccd); do
#  if ! [ -d $diagnosis_and_badFilesDir/CCD"$h" ]; then mkdir $diagnosis_and_badFilesDir/CCD"$h"; fi
#done

#if [ "$MODEL_SKY_AS_CONSTANT" = true ]; then
#  tmpDir=$BDIR/noise-sky_it2
#else
#  tmpDir=$noiseskyctedir
#fi
#for h in $(seq 1 $num_ccd); do
#  python3 $pythonScriptsPath/diagnosis_normalisedBackgroundMagnitudes.py $tmpDir $framesForCommonReductionDir $airMassKeyWord $alphatruedir $pixelScale $diagnosis_and_badFilesDir $h $dateHeaderKey
#done
# We mask again the points in order to measure (after photometric calibration) the sky accurately
##Recover the full grid


smallPointings_photCorr_maskedDir=$BDIR/photCorrSmallGrid_masked_it$iteration
maskedPointingsDone=$smallPointings_photCorr_maskedDir/done_.txt
maskPointings $gaincordir $smallPointings_photCorr_maskedDir $maskedPointingsDone $maskName $entiredir_smallGrid

#Now we dont need the it1 ones
find $BDIR/photCorrFullGrid-dir_it1 -type f ! -name 'done*' -exec rm {} \;

noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done.txt
# Since here we compute the sky for obtaining the rms, we model it as a cte (true) and the polynomial degree is irrelevant (-1)
computeSky $smallPointings_photCorr_maskedDir $noiseskydir $noiseskydone true wholeImage -1 true $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'"


minRmsFileName=min_rms_it$iteration.txt
python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $noiseskydir $DIR $iteration min_rms_it$iteration.txt

photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
photCorrFullGridDone=$photCorrFullGridDir/done.txt
identifiedBadDetectors=$CDIR/identifiedBadDetectors.txt #We can now provide a list of bad detectors to blank them
if ! [ -d $photCorrFullGridDir ]; then mkdir $photCorrFullGridDir; fi
smallGridtoFullGridAndWeight $photCorrSmallGridDir $photCorrFullGridDir $photCorrFullGridDone $coaddSizePx $ra $dec $noiseskydir $minRmsFileName $iteration $identifiedBadDetectors 


echo -e "\n·Removing bad frames"

diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
rejectedFramesDir=$BDIR/rejectedFrames
if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"

rejectedByAstrometry=identifiedBadFrames_astrometry.txt
removeBadFramesFromReduction $photCorrFullGridDir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry
removeBadFramesFromReduction $noiseskydir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry


# Store the minimum standard deviation of the frames in order to compute the weights



echo -e "\n ${GREEN} ---Masking outliers--- ${NOCOLOUR}"
writeTimeOfStepToFile "Masking outliers" $fileForTimeStamps
# Remove outliers before the final coadd by using sigclip-median and sigclip-std
# This is particularly important to remove cosmic rays

sigmaForStdSigclip=3
clippingdir=$BDIR/clipping-outliers_it"$iteration"
clippingdone=$clippingdir/done.txt
buildUpperAndLowerLimitsForOutliersNew $clippingdir $clippingdone $photCorrFullGridDir $sigmaForStdSigclip

mowdir=$BDIR/photCorrFullGrid-dir-no-outliers_it$iteration
mowdone=$mowdir/done.txt
if ! [ -d $mowdir ]; then mkdir $mowdir; fi
removeOutliersFromWeightedFramesNew $mowdone $mowdir $clippingdir $photCorrFullGridDir





echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
writeTimeOfStepToFile "Building coadd" $fileForTimeStamps


echo -e "\n·Building coadd"
coaddDir=$BDIR/coadds_it"$iteration"
coaddDone=$coaddDir/done.txt
coaddName=$coaddDir/"$objectName"_coadd_"$filter"_it"$iteration".fits
buildCoaddNew $coaddDir $coaddName $mowdir $coaddDone

maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  #If block scale is greater than 1, we apply the block
  if [ "$blockScale" -gt 1 ]; then
    astwarp $coaddName -h1 --scale=1/$blockScale --numthreads=$num_cpus -o $coaddDir/coadd_blocked.fits
    imToMask=$coaddDir/coadd_blocked.fits
  else
    imToMask=$coaddName
  fi
  #If a kernel exists in the configuration file, we apply it
  kernelFile=$CDIR/kernel.fits
  if [ -f $kernelFile ]; then
    astconvolve $imToMask --kernel=$kernelFile --domain=spatial --numthreads=$num_cpus -o $coaddDir/coadd_convolved.fits
    imToMask=$coaddDir/coadd_convolved.fits
  fi
  astnoisechisel $imToMask $noisechisel_param --numthreads=$num_cpus -o $coaddDir/mask_warped.fits
  if [ "$blockScale" -gt 1 ]; then 
    astwarp $coaddDir/mask_warped.fits --gridfile=$coaddName --numthreads=$num_cpus -o $coaddDir/mask_unwarped.fits
    astarithmetic $coaddDir/mask_unwarped.fits -h1 set-i i i 0 gt 1 where float32 -q -o $maskName
    rm $coaddDir/mask_unwarped.fits  $coaddDir/mask_warped.fits $coaddDir/coadd_blocked.fits
  else
    mv $coaddDir/mask_warped.fits $maskName
  fi
  rm $coaddDir/coadd_convolved.fits 2>/dev/null
fi

#astnoisechisel with the current parameters might fail due to long tilesize. I'm gonna make 2 checks to see if it fails, decreasing in steps of 5 in tilesize
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  #We assume that if this works for this iteration, then the next one will need at least same parameters
  tileSize=20
  noisechisel_param="--tilesize=$tileSize,$tileSize \
                    --minskyfrac=0.9 \
                     --meanmedqdiff=0.01 \
                     --snthresh=5.2 \
                     --detgrowquant=0.7 \
                     --detgrowmaxholesize=1000 \
                     --rawoutput"
  astnoisechisel $coaddName $noisechisel_param  -o $maskName
fi
if [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd is already done"
else
  #We assume that if this works for this iteration, then the next one will need at least same parameters
  tileSize=$((tileSize - 5))
  noisechisel_param="--tilesize=$tileSize,$tileSize \
                    --minskyfrac=0.9 \
                     --meanmedqdiff=0.01 \
                     --snthresh=5.2 \
                     --detgrowquant=0.7 \
                     --detgrowmaxholesize=1000 \
                     --rawoutput"
  astnoisechisel $coaddName $noisechisel_param --numthreads=$num_cpus  -o $maskName
fi

if ! [ -f $maskName ]; then
  echo -e "\tThe mask of the weighted coadd could not be generated. Please, check manually."
  exit 47
fi
exposuremapDir=$coaddDir/"$objectName"_exposureMap
exposuremapdone=$coaddDir/done_exposureMap_eff.txt
exposureMapName=$coaddDir/"$objectName"_expMap_eff_"$filter"_it"$iteration".fits
computeExposureMapNew $mowdir $exposuremapDir $exposuremapdone $exposureMapName
#
exposuremapdone=$coaddDir/done_exposureMap.txt
exposureMapName=$coaddDir/"$objectName"_expMap_"$filter"_it"$iteration".fits
computeExposureMapNew $photCorrFullGridDir $exposuremapDir $exposuremapdone $exposureMapName

#Compute surface brightness limit
sblimitFile=$coaddDir/"$objectName"_"$filter"_sblimit.txt

if [ -f  $sblimitFile ]; then
    echo -e "\n\tSurface brightness limit for coadd already measured\n"
else
    surfaceBrightnessLimit=$( limitingSurfaceBrightness $coaddName $maskName $exposureMapName $coaddDir $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )
fi


times=($(getInitialMidAndFinalFrameTimes $INDIR $dateHeaderKey))
initialTime=$( date -d @"${times[0]}" "+%Y-%m-%d_%H:%M:%S")
meanTime=$( date -d @"${times[1]}" "+%Y-%m-%d_%H:%M:%S")
finalTime=$( date -d @"${times[2]}" "+%Y-%m-%d_%H:%M:%S")

keyWords=("FRAMES_COMBINED" \
          "NUMBER_OF_DIFFERENT_NIGHTS" \
          "INITIAL_DATE_OBS"
          "MEAN_DATE_OBS"
          "FINAL_DATE_OBS"
          "FILTER" \
          "SATURATION_THRESHOLD" \
          "CALIBRATION_BRIGHTLIMIT" \
          "CALIBRATION_FAINTLIMIT" \
          "RUNNING_FLAT" \
          "WINDOW_SIZE" \
          "STD_FOR_BAD_FRAMES" \
          "SURFACE_BRIGHTNESS_LIMIT")

numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames" "$surfaceBrightnessLimit")
comments=("" "" "" "" "" "" "" "" "" "" "" "Num. of tandard deviations used for rejection" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")
astfits $coaddName --write=/,"Pipeline information"
addkeywords $coaddName keyWords values comments

halfMaxRadForCoaddName=$diagnosis_and_badFilesDir/coadd_it1.png
if [ -f $halfMaxRadForCoaddName ]; then
  echo -e "\tThe Half-Max-Rad vs Magnitude has been already generate for the coadd"
else
  produceHalfMaxRadVsMagForSingleImage $coaddName $diagnosis_and_badFilesDir $catdir/"$objectName"_Gaia_eDR3.fits $toleranceForMatching $pythonScriptsPath "coadd_it1" 100 NO
fi
##To avoid problems, we re-name the pythonScriptsPath and toleranceForMatching because it is creating problems
#pythonScriptsPath=$pipelinePath/pipelineScripts
#toleranceForMatching=2 
#export pythonScriptsPath
#export toleranceForMatching

writeTimeOfStepToFile "Producing frames with coadd subtracted" $fileForTimeStamps
framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted_it"$iteration"
framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi
if [ -f $framesWithCoaddSubtractedDone ]; then
    echo -e "\n\tFrames with coadd subtracted already generated\n"
else
  sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_"$filter"_it"$iteration".fits
  coadd_av=$coaddDir/"$objectName"_coadd_it"$iteration"_average.fits
  gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
  if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
    astarithmetic $(ls $photCorrFullGridDir/*.fits) -g1 $(ls $photCorrFullGridDir/*.fits | wc -l ) mean --writeall -o$coadd_av
  else
    astarithmetic $(ls $photCorrFullGridDir/*.fits) -g1 $(ls $photCorrFullGridDir/*.fits | wc -l ) mean -o$coadd_av
  fi    
  
  subtractCoaddToFramesNew $photCorrFullGridDir $coadd_av $framesWithCoaddSubtractedDir
  if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
    astarithmetic $framesWithCoaddSubtractedDir/*.fits -g1 $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l ) sum --writeall -o$sumMosaicAfterCoaddSubtraction
  else
    astarithmetic $framesWithCoaddSubtractedDir/*.fits -g1 $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l ) sum -o$sumMosaicAfterCoaddSubtraction
  fi
  echo done > $framesWithCoaddSubtractedDone 
fi

if [[ "$subtractStarsFromRaw" == "true" ]]; then
  ##First we remove these files not needed as in it1
  # # Remove intermediate folders to save some space
  find $BDIR/noise-sky_it2 -type f -name '*.fits' -exec rm {} \;
  find $BDIR/noise-sky-after-photometry_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/noise-sky_forPlaneCoadd_it2 -type f ! -name 'done*' -exec rm {} \;

  #find $BDIR/sub-sky-fullGrid_it1 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/sub-sky-smallGrid_it2 -type f ! -name 'done*' -exec rm {} \;
  #find $BDIR/photCorrFullGrid-dir_it1 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/photCorrSmallGrid-dir_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/photCorrFullGrid-dir-no-outliers_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/gain-corrected -type -f ! -name 'done*' -exec rm {} \;
  find $BDIR/my-catalog-halfmaxradius_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/match-decals-myData_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/ourData-catalogs-apertures_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/framesWithCoaddSubtracted_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/pointings_smallGrid_masked_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/photCorrSmallGrid_masked_it2 -type f ! -name 'done*' -exec rm {} \;
  #find $BDIR/weight-dir-no-outliers -type f ! -name 'done*' -exec rm {} \;
  #find $BDIR/weight-dir-no-outliers_plane_it1 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/weight-dir_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/weight-dir_plane_it2 -type f ! -name 'done*' -exec rm {} \;

  #find $BDIR/only-weight-dir-no-outliers -type f ! -name 'done*' -exec rm {} \;
  #find $BDIR/only-weight-dir-no-outliers_plane_it1 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/only-w-dir_it2 -type f ! -name 'done*' -exec rm {} \;
  find $BDIR/only-w-dir_plane_it2 -type f ! -name 'done*' -exec rm {} \;
  if [[ ("$produceCoaddPrephot" = "true") || ("$produceCoaddPrephot" = "True" )]]; then
    find $BDIR/weight-dir_prephot_it2 -type f ! -name 'done*' -exec rm {} \;
    find $BDIR/only-w-dir_prephot_it2 -type f ! -name 'done*' -exec rm {} \;
    find $BDIR/noise-sky_prephot_it2 -type f ! -name 'done*' -exec rm {} \;
  fi
  cp $BDIR/coadds_it2/"$objectName"_coadd_"$filter"_mask.fits $BDIR/coadds_it2/"$objectName"_coadd_"$filter"_mask_copy.fits

# Then we apply the user-defined masks
  valueToPut=1
  read -r -a maskArray <<< "$maskParams"
  for ((i=0; i<${#maskArray[@]}; i+=5)); do
  	ra="${maskArray[i]}"
  	dec="${maskArray[i+1]}"
  	r="${maskArray[i+2]}"
  	axisRatio="${maskArray[i+3]}"
  	pa="${maskArray[i+4]}"

  	
  	python3 $pythonScriptsPath/manualMaskRegionFromWCSArea.py $BDIR/coadds_it2/"$objectName"_coadd_"$filter"_mask.fits $valueToPut $ra $dec $r $axisRatio $pa
  done

#If a mask already exists on the configuration file, we combine both masks
  if [ -f $CDIR/mask.fits ]; then
    cp $BDIR/coadds_it2/"$objectName"_coadd_"$filter"_mask.fits $BDIR/coadds_it2/"$objectName"_coadd_"$filter"_mask_copy.fits
    astarithmetic $BDIR/coadds_it2/"$objectName"_coadd_"$filter"_mask_copy.fits $CDIR/mask.fits -g1 1 eq 1 where -q -o $BDIR/coadds_it2/"$objectName"_coadd_"$filter"_mask.fits
  fi 
  ####STAR SUBTRACTION####
  iteration=2
  # We will now work in photometrized and gain corrected frames

  alphatruedir=$BDIR/alpha-stars-true_it$iteration
  photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_withSky
  applyCalibrationFactors $entiredir_smallGrid $alphatruedir $photCorrSmallGridDir $iteration $applyCommonCalibrationFactor
  gainCorrectionFile=$BDIR/gainCorrectionParameters_it$iteration.txt
  ccd_ref=1
  gaincordir=$BDIR/gain-corrected_withSky
  gaincordone=$gaincordir/done.txt


  if ! [ -d $gaincordir ]; then mkdir $gaincordir; fi
  if [ -f $gaincordone ]; then
    echo -e "\nMulti-layer windows already gain corrected"
  else
    frameNames=()
    for a in $entiredir_smallGrid/*.fits; do
      frameNames+=("$(basename $a)")
    done
    printf "%s\n" "${frameNames[@]}" | parallel -j "$num_parallel" gainCorrection {} $photCorrSmallGridDir $gainCorrectionFile $gaincordir $ccd_ref
    echo done > $gaincordone
  fi
  echo -e "\n\t${GREEN} --- Subtract stars from frames --- ${NOCOLOUR} \n"

  psfFile=$CDIR/PSF_"$filter".fits
  psfRadFile=$CDIR/RP_PSF_"$filter".fits
####So far I will do in this way, I will for sure change it
###Catalog will be: RA DEC MAG
###We will make the following:
##  # Check where the star falls in a circle centered on RA, DEC and radius=RAFEC
##  # If AZ exist in the catalog we will measure the profile in azimuth
##  # This CCD will be used to compute scale factor between Rmin and Rmax, tunning the range with MAG
##  # Finally, we will subtract from all the frames where the star is, and continue to the next one
##  # For background range, we will measure a first background in the CCD, to select a range ±500ADU
##
  input_subStar_small=$gaincordir
##input_subStar_full=$normalised_fullGrid
  starsToSubtract=$BDIR/starsToSubtract.txt
  radiusToSearch=$(awk -v r="$sizeOfOurFieldDegrees" 'BEGIN { printf "%.6f", r/2 }')
  query_param="gaia --dataset=dr3 --center=$ra_gal,$dec_gal --radius=$radiusToSearch --column=ra,dec,phot_g_mean_mag"
  if ! [ -f $starsToSubtract ]; then
    astquery $query_param -o$BDIR/starsToSubtract_temp.fits
    asttable $BDIR/starsToSubtract_temp.fits --range=3,0:13.5 --sort=3 -o$starsToSubtract
    rm $BDIR/starsToSubtract_temp.fits
  fi
  if [ -z "$starSatThreshold" ]; then
    starSatThreshod=$saturationThreshold
  fi 
#
  starId=0
  while IFS= read -r line; do
    #We skip the lines that contain info about the columns
    [[ $line =~ ^#.*$ ]] && continue
#
    ((starId++))
    outputDir_small=$BDIR/pointings_smallGrid_sub$starId
    subtractStars $input_subStar_small "$line" $psfFile $psfRadFile $outputDir_small $starId $starSatThreshold $BDIR/commonCalibrationFactor_it2.txt $gainCorrectionFile
#  
    #if (( $(echo "$starId == 2" | bc -l) )); then exit; fi
    if ! (( $(echo "$starId == 1" | bc -l) )); then
	    rm $input_subStar_small/*.fits
    fi
    #Sanity check: if something failed, we stop the process
    for file in $outputDir_small/*.fits; do
      nhdu=$( astfits $file --numhdus -q )
      if [ $nhdu -lt $((num_ccd+1)) ]; then
        echo "Some frames have failed in the subtraction of star $starId"
        exit 23
      fi
    done
    input_subStar_small=$outputDir_small

  done < $starsToSubtract
  starsSub_small=$outputDir_small

  iteration=3
  smallPointings_maskedDir=$BDIR/pointings_smallGrid_masked_it$iteration
  maskedPointingsDone=$smallPointings_maskedDir/done_.txt
  maskName=$BDIR/coadds_it2/"$objectName"_coadd_"$filter"_mask.fits

  maskPointings $starsSub_small $smallPointings_maskedDir $maskedPointingsDone $maskName $entiredir_smallGrid


  noiseskydir=$BDIR/noise-sky_it$iteration
  noiseskydone=$noiseskydir/done_"$filter".txt

  subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it$iteration
  subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter".txt

  imagesAreMasked=true
  sky_estimation_method=wholeImage #If we trust the mask, we can use the full image
  computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone $MODEL_SKY_AS_CONSTANT $sky_estimation_method $polynomialDegree $imagesAreMasked $BDIR/ring $USE_COMMON_RING $keyWordToDecideRing $keyWordThreshold $keyWordValueForFirstRing $keyWordValueForSecondRing $ringWidth YES $blockScale "'$noisechisel_param'" "'$maskParams'"
  subtractSky $starsSub_small $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir $MODEL_SKY_AS_CONSTANT
  minRmsFileName=min_rms_it$iteration.txt
  python3 $pythonScriptsPath/find_rms_min.py $filter 1 $totalNumberOfFrames $noiseskydir $DIR $iteration min_rms_it$iteration.txt
  
  photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
  photCorrFullGridDone=$photCorrFullGridDir/done.txt
  identifiedBadDetectors=$CDIR/identifiedBadDetectors.txt #We can now provide a list of bad detectors to blank them
  if ! [ -d $photCorrFullGridDir ]; then mkdir $photCorrFullGridDir; fi
  smallGridtoFullGridAndWeight $subskySmallGrid_dir $photCorrFullGridDir $photCorrFullGridDone $coaddSizePx $ra $dec $noiseskydir $minRmsFileName $iteration $identifiedBadDetectors 


  echo -e "\n·Removing bad frames"

  diagnosis_and_badFilesDir=$BDIR/diagnosis_and_badFiles
  rejectedFramesDir=$BDIR/rejectedFrames
  if ! [ -d $rejectedFramesDir ]; then mkdir $rejectedFramesDir; fi
  echo -e "\nRemoving (moving to $rejectedFramesDir) the frames that have been identified as bad frames"

  rejectedByAstrometry=identifiedBadFrames_astrometry.txt
  removeBadFramesFromReduction $photCorrFullGridDir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry
  removeBadFramesFromReduction $noiseskydir $rejectedFramesDir $diagnosis_and_badFilesDir $rejectedByAstrometry


  # Store the minimum standard deviation of the frames in order to compute the weights



  echo -e "\n ${GREEN} ---Masking outliers--- ${NOCOLOUR}"
  writeTimeOfStepToFile "Masking outliers" $fileForTimeStamps
  # Remove outliers before the final coadd by using sigclip-median and sigclip-std
  # This is particularly important to remove cosmic rays

  sigmaForStdSigclip=3
  clippingdir=$BDIR/clipping-outliers_it"$iteration"
  clippingdone=$clippingdir/done.txt
  buildUpperAndLowerLimitsForOutliersNew $clippingdir $clippingdone $photCorrFullGridDir $sigmaForStdSigclip
  mowdir=$BDIR/photCorrFullGrid-dir-no-outliers_it$iteration
  mowdone=$mowdir/done.txt
  if ! [ -d $mowdir ]; then mkdir $mowdir; fi
  removeOutliersFromWeightedFramesNew $mowdone $mowdir $clippingdir $photCorrFullGridDir




  echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"
  writeTimeOfStepToFile "Building coadd" $fileForTimeStamps


  echo -e "\n·Building coadd"
  coaddDir=$BDIR/coadds_it"$iteration"
  coaddDone=$coaddDir/done.txt
  coaddName=$coaddDir/"$objectName"_coadd_"$filter"_it"$iteration".fits
  buildCoaddNew $coaddDir $coaddName $mowdir $coaddDone
  
  maskName=$coaddDir/"$objectName"_coadd_"$filter"_mask.fits
  if [ -f $maskName ]; then
    echo -e "\tThe mask of the weighted coadd is already done"
  else
    #If block scale is greater than 1, we apply the block
    if [ "$blockScale" -gt 1 ]; then
      astwarp $coaddName -h1 --scale=1/$blockScale --numthreads=$num_cpus -o $coaddDir/coadd_blocked.fits
      imToMask=$coaddDir/coadd_blocked.fits
    else
      imToMask=$coaddName
    fi
    #If a kernel exists in the configuration file, we apply it
    kernelFile=$CDIR/kernel.fits
    if [ -f $kernelFile ]; then
      astconvolve $imToMask --kernel=$kernelFile --domain=spatial --numthreads=$num_cpus -o $coaddDir/coadd_convolved.fits
      imToMask=$coaddDir/coadd_convolved.fits
    fi
    astnoisechisel $imToMask $noisechisel_param --numthreads=$num_cpus -o $coaddDir/mask_warped.fits
    if [ "$blockScale" -gt 1 ]; then 
      astwarp $coaddDir/mask_warped.fits --gridfile=$coaddName --numthreads=$num_cpus -o $coaddDir/mask_unwarped.fits
      astarithmetic $coaddDir/mask_unwarped.fits -h1 set-i i i 0 gt 1 where float32 -q -o $maskName
      rm $coaddDir/mask_unwarped.fits  $coaddDir/mask_warped.fits $coaddDir/coadd_blocked.fits
    else
      mv $coaddDir/mask_warped.fits $maskName
    fi
    rm $coaddDir/coadd_convolved.fits 2>/dev/null
  fi

  #astnoisechisel with the current parameters might fail due to long tilesize. I'm gonna make 2 checks to see if it fails, decreasing in steps of 5 in tilesize
  if [ -f $maskName ]; then
    echo -e "\tThe mask of the weighted coadd is already done"
  else
    #We assume that if this works for this iteration, then the next one will need at least same parameters
    tileSize=20
    noisechisel_param="--tilesize=$tileSize,$tileSize \
                      --minskyfrac=0.9 \
                       --meanmedqdiff=0.01 \
                       --snthresh=5.2 \
                       --detgrowquant=0.7 \
                       --detgrowmaxholesize=1000 \
                       --rawoutput"
    astnoisechisel $coaddName $noisechisel_param  -o $maskName
  fi
  if [ -f $maskName ]; then
    echo -e "\tThe mask of the weighted coadd is already done"
  else
    #We assume that if this works for this iteration, then the next one will need at least same parameters
    tileSize=$((tileSize - 5))
    noisechisel_param="--tilesize=$tileSize,$tileSize \
                      --minskyfrac=0.9 \
                       --meanmedqdiff=0.01 \
                       --snthresh=5.2 \
                       --detgrowquant=0.7 \
                       --detgrowmaxholesize=1000 \
                       --rawoutput"
    astnoisechisel $coaddName $noisechisel_param --numthreads=$num_cpus  -o $maskName
  fi

  if ! [ -f $maskName ]; then
    echo -e "\tThe mask of the weighted coadd could not be generated. Please, check manually."
    exit 47
  fi
  exposuremapDir=$coaddDir/"$objectName"_exposureMap
  exposuremapdone=$coaddDir/done_exposureMap_eff.txt
  exposureMapName=$coaddDir/"$objectName"_expMap_eff_"$filter"_it"$iteration".fits
  computeExposureMapNew $mowdir $exposuremapDir $exposuremapdone $exposureMapName
  #
  exposuremapdone=$coaddDir/done_exposureMap.txt
  exposureMapName=$coaddDir/"$objectName"_expMap_"$filter"_it"$iteration".fits
  computeExposureMapNew $photCorrFullGridDir $exposuremapDir $exposuremapdone $exposureMapName

  #Compute surface brightness limit
  sblimitFile=$coaddDir/"$objectName"_"$filter"_sblimit.txt

  if [ -f  $sblimitFile ]; then
      echo -e "\n\tSurface brightness limit for coadd already measured\n"
  else
      surfaceBrightnessLimit=$( limitingSurfaceBrightness $coaddName $maskName $exposureMapName $coaddDir $areaSBlimit $fractionExpMap $pixelScale $sblimitFile )
  fi


  times=($(getInitialMidAndFinalFrameTimes $INDIR $dateHeaderKey))
  initialTime=$( date -d @"${times[0]}" "+%Y-%m-%d_%H:%M:%S")
  meanTime=$( date -d @"${times[1]}" "+%Y-%m-%d_%H:%M:%S")
  finalTime=$( date -d @"${times[2]}" "+%Y-%m-%d_%H:%M:%S")

  keyWords=("FRAMES_COMBINED" \
            "NUMBER_OF_DIFFERENT_NIGHTS" \
            "INITIAL_DATE_OBS"
            "MEAN_DATE_OBS"
            "FINAL_DATE_OBS"
            "FILTER" \
            "SATURATION_THRESHOLD" \
            "CALIBRATION_BRIGHTLIMIT" \
            "CALIBRATION_FAINTLIMIT" \
            "RUNNING_FLAT" \
            "WINDOW_SIZE" \
            "STD_FOR_BAD_FRAMES" \
            "SURFACE_BRIGHTNESS_LIMIT")

  numberOfFramesCombined=$(ls $mowdir/*.fits | wc -l)
  values=("$numberOfFramesCombined" "$numberOfNights" "$initialTime" "$meanTime" "$finalTime" "$filter" "$saturationThreshold" "$calibrationBrightLimit" "$calibrationFaintLimit" "$RUNNING_FLAT" "$windowSize" "$numberOfStdForBadFrames" "$surfaceBrightnessLimit")
  comments=("" "" "" "" "" "" "" "" "" "" "" "Num. of tandard deviations used for rejection" "[mag/arcsec^2](3sig;"$areaSBlimit"x"$areaSBlimit" arcsec)")
  astfits $coaddName --write=/,"Pipeline information"
  addkeywords $coaddName keyWords values comments

  
  ##To avoid problems, we re-name the pythonScriptsPath and toleranceForMatching because it is creating problems
  #pythonScriptsPath=$pipelinePath/pipelineScripts
  #toleranceForMatching=2 
  #export pythonScriptsPath
  #export toleranceForMatching

  writeTimeOfStepToFile "Producing frames with coadd subtracted" $fileForTimeStamps
  framesWithCoaddSubtractedDir=$BDIR/framesWithCoaddSubtracted_it"$iteration"
  framesWithCoaddSubtractedDone=$framesWithCoaddSubtractedDir/done_framesWithCoaddSubtracted.txt
  if ! [ -d $framesWithCoaddSubtractedDir ]; then mkdir $framesWithCoaddSubtractedDir; fi
  if [ -f $framesWithCoaddSubtractedDone ]; then
      echo -e "\n\tFrames with coadd subtracted already generated\n"
  else
    sumMosaicAfterCoaddSubtraction=$coaddDir/"$objectName"_sumMosaicAfterCoaddSub_"$filter"_it"$iteration".fits
    coadd_av=$coaddDir/"$objectName"_coadd_it"$iteration"_average.fits
    gnuastro_version=$(astarithmetic --version | head -n1 | awk '{print $NF}')
    if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
      astarithmetic $(ls $photCorrFullGridDir/*.fits) -g1 $(ls $photCorrFullGridDir/*.fits | wc -l ) mean --writeall -o$coadd_av
    else
      astarithmetic $(ls $photCorrFullGridDir/*.fits) -g1 $(ls $photCorrFullGridDir/*.fits | wc -l ) mean -o$coadd_av
    fi    

    subtractCoaddToFramesNew $photCorrFullGridDir $coadd_av $framesWithCoaddSubtractedDir
    if [ "$(echo "$gnuastro_version > 0.22" | bc)" -eq 1 ]; then
      astarithmetic $framesWithCoaddSubtractedDir/*.fits -g1 $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l ) sum --writeall -o$sumMosaicAfterCoaddSubtraction
    else
      astarithmetic $framesWithCoaddSubtractedDir/*.fits -g1 $(ls $framesWithCoaddSubtractedDir/*.fits | wc -l ) sum -o$sumMosaicAfterCoaddSubtraction
    fi
    echo done > $framesWithCoaddSubtractedDone 
  fi
fi
endTime=$(date +%D%T)
echo "Pipeline ended at : ${endTime}"
exit 0



