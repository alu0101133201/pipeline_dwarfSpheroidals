#!/bin/bash

# Information of the instruments used for adquiring the data

# Montura Astro-physics, modelo mach one
# Telescopio f5 , apertura 130 mm y largo focal 655 mm de la marca
# Stellarvue SVS 130.
# Camara fotografica ASI ZWO 1600GT monocromática con filtros RGB+L (bit depth of 12bits, seems to have been multiplied to 2⁴ to obtain a typical 0-2^16 range)
# PIXEL 3,8UM
# Resolucion 1,2"/pixel (after astrometrising it seems to be 0.164"/pixel)
# SQM 20,9 - 21,0
# Bortle 3
# Filtro de luminancia, marca Astrodon



ORANGE='\033[0;33m'
GREEN='\033[0;32m'
NOCOLOUR='\033[0m'

export ORANGE
export GREEN
export NOCOLOUR

########## Loading modules ##########
echo -e "\n ${GREEN} ---Loading Modules--- ${NOCOLOUR}"

module load gnuastro/0.22
module load astronomy.net/0.93

########## Functions ##########

# Functions of general sue
escapeSpacesFromString() {
  local input="$1"
  escaped_string="${input// /\\ }"
  echo $escaped_string
}
export -f escapeSpacesFromString

checkIfExist_DATEOBS() {
  DATEOBSValue=$1

  if [ "$DATEOBSValue" = "n/a" ]; then
    errorNumber=2
    echo -e "The file $i do not has the $dateHeaderKey, used for sorting the raw files for the pipeline"
    echo -e "Exiting with error number: $errorNumber"
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

retrieveChosenDecalsBrickForFrame() {
  frame=$1
  mapFile=$2

  chosenBrick=$( grep "^$frame:" $mapFile | cut -d ':' -f 2 )
  echo $chosenBrick
}
export -f retrieveChosenDecalsBrickForFrame

# Functions used in Flat
maskImages() {
  inputDirectory=$1
  masksDirectory=$2
  outputDirectory=$3

  for a in $(seq 1 $n_exp); do
    base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
    i=$inputDirectory/$base
    astarithmetic $i -h1 $masksDirectory/$base -h1 1 eq nan where float32 -o $outputDirectory/$base -q
  done
}
export -f maskImages

normaliseImagesWithRing() {
  imageDir=$1
  outputDir=$2
  ringFile=$3

  for a in $(seq 1 $n_exp); do
    base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
    i=$imageDir/$base
    out=$outputDir/$base

    me=$(astarithmetic $i -h1 $ringFile -h1 0 eq nan where medianvalue --quiet)
    astarithmetic $i -h1 $me / -o $out
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
  done
  echo done > $flatDone
}
export -f divideImagesByWholeNightFlat


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
  frameCentre=$( getCentralCoordinate $imageToSwarp )
  centralRa=$(echo "$frameCentre" | awk '{print $1}')
  centralDec=$(echo "$frameCentre" | awk '{print $2}')

  tmpFile1=$entiredir"/$currentIndex"_temp1.fits
  frameFullGrid=$entireDir_fullGrid/entirecamera_$currentIndex.fits

  # Resample into the final grid
  SWarp -c $swarpcfg $imageToSwarp -CENTER $ra,$dec -IMAGE_SIZE $coaddSizePx,$coaddSizePx -IMAGEOUT_NAME $entiredir/"$currentIndex"_swarp1.fits -WEIGHTOUT_NAME $entiredir/"$currentIndex"_swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $pixelScale -PIXELSCALE_TYPE  MANUAL
  # Mask bad pixels
  astarithmetic $entiredir/"$currentIndex"_swarp_w1.fits -h0 set-i i i 0 lt nan where -o$tmpFile1
  astarithmetic $entiredir/"$currentIndex"_swarp1.fits -h0 $tmpFile1 -h1 0 eq nan where -o$frameFullGrid

  # Offset for having some margin and not lose any pixel
  securityOffset=50
  detectorWidthDeg=$(echo  "(($detectorWidth + $securityOffset) * $pixelScale)" | bc )
  detectorHeightDeg=$(echo "(($detectorHeight + $securityOffset) * $pixelScale) + $securityOffset" | bc )
  astcrop $frameFullGrid --center=$centralRa,$centralDec --mode=wcs --width=$detectorHeightDeg/3600,$detectorWidthDeg/3600  -o $entiredir/entirecamera_"$currentIndex".fits

  # SWarp -c $swarpcfg $entiredir/entirecamera_"$currentIndex".fits -CENTER $ra,$dec -IMAGE_SIZE $coaddSizePx,$coaddSizePx -IMAGEOUT_NAME $entiredir/"$currentIndex"_back.fits -WEIGHTOUT_NAME $entiredir/"$currentIndex"_swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $pixelScale -PIXELSCALE_TYPE  MANUAL
  rm $entiredir/"$currentIndex"_swarp_w1.fits $entiredir/"$currentIndex"_swarp1.fits $tmpFile1 # $tmpFile2
}
export -f warpImage


# Functions for compute and subtract sky from frames
computeSkyForFrame(){
  base=$1
  entiredir=$2
  noiseskydir=$3

  i=$entiredir/$1
  sky=$(echo $base | sed 's/.fits/_sky.fits/')
  out=$(echo $base | sed 's/.fits/.txt/')

  # The sky substraction is done by using the --checksky option in noisechisel.
  astnoisechisel $i --tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2 --checksky $noisechisel_param -o $noiseskydir/$base

  # save mean and std
  m=$(aststatistics $noiseskydir/$sky -hSKY --sigclip-mean)
  s=$(aststatistics $noiseskydir/$sky -hSTD --sigclip-mean)
  echo "$base $m $s" > $noiseskydir/$out

  # Remove fist file for saving storage
  rm -f $noiseskydir/$sky
}
export -f computeSkyForFrame

computeSky() {
  framesToUseDir=$1
  noiseskydir=$2
  noiseskydone=$3


  if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
  if [ -f $noiseskydone ]; then
    echo -e "\nScience images are 'noisechiseled' for sky substraction for extension $h\n"
  else
    framesToComputeSky=()
    for a in $(seq 1 $totalNumberOfFrames); do
      base="entirecamera_"$a.fits
      framesToComputeSky+=("$base")
    done
    printf "%s\n" "${framesToComputeSky[@]}" | parallel -j "$num_cpus" computeSkyForFrame {} $framesToUseDir $noiseskydir
    echo done > $noiseskydone
  fi
}

subtractSkyForFrame() {
  a=$1
  directoryWithSkyValues=$2
  framesToSubtract=$3
  directoryToStoreSkySubtracted=$4

  base="entirecamera_"$a.fits
  i=$directoryWithSkyValues/"entirecamera_"$a.txt
  me=$(awk 'NR=='1'{print $2}' $i)
  input=$framesToSubtract/$base
  output=$directoryToStoreSkySubtracted/$base
  astarithmetic $input -h1 $me - -o$output;
}
export -f subtractSkyForFrame

subtractSky() {
  framesToSubtract=$1
  directoryToStoreSkySubtracted=$2
  directoryToStoreSkySubtracteddone=$3
  directoryWithSkyValues=$4
 
  if ! [ -d $directoryToStoreSkySubtracted ]; then mkdir $directoryToStoreSkySubtracted; fi
  if [ -f $directoryToStoreSkySubtracteddone ]; then
    echo -e "\nSky substraction is already done for the science images for extension $h\n"
  else
  framesToSubtractSky=()
  for a in $(seq 1 $totalNumberOfFrames); do
      framesToSubtractSky+=("$a")      
  done
  printf "%s\n" "${framesToSubtractSky[@]}" | parallel -j "$num_cpus" subtractSkyForFrame {} $directoryWithSkyValues $framesToSubtract $directoryToStoreSkySubtracted
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
  astsegment $tmpFolder/det.fits -o $tmpFolder/seg.fits --snquant=0.1 --gthresh=-10 --objbordersn=0  --minriverlength=3 1>/dev/null
  astmkcatalog $tmpFolder/seg.fits --ra --dec --magnitude --half-max-radius --sum --clumpscat -o $tmpFolder/decals.txt --zeropoint=22.5 1>/dev/null
  astmatch $tmpFolder/decals_c.txt --hdu=1  $BDIR/Gaia_eDR3.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=bRA,bDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o $tmpFolder/match_decals_gaia.txt 1>/dev/null

  numOfStars=$( cat $tmpFolder/match_decals_gaia.txt | wc -l )
  median=$( asttable $tmpFolder/match_decals_gaia.txt -h1 -c3 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median )
  std=$( asttable $tmpFolder/match_decals_gaia.txt -h1 -c3 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std )
  rm $tmpFolder/*
  echo $median $std $numOfStars
}

chooseBrickToUse() {
  bricks=$1
  brickDir=$2
  kernel=$3
  tmpFolder=$4
  bestBrickRecord=$5

  minStd=9999
  bestBrick=$( grep "^$bricks:" $bestBrickRecord | cut -d ':' -f 2 )

  # The following conditional is not to compute the same decision twice, so we recover previous brick-decision from "$bestBrickRecord" file
  if [ -n "$bestBrick" ]; then
    echo $bestBrick
  else
    for i in $bricks; do
      brickFullName_g=$brickDir/"decal_image_"$i"_g.fits"
      brickFullName_r=$brickDir/"decal_image_"$i"_r.fits"
      parameters_g=$(getParametersFromHalfMaxRadius $brickFullName_g "$BDIR/Gaia_eDR3.fits" $kernel $tmpFolder)
      parameters_r=$(getParametersFromHalfMaxRadius $brickFullName_r "$BDIR/Gaia_eDR3.fits" $kernel $tmpFolder)

      # The following code converts the parameters obtained in the previous lines to the number format needed for numerical comparison
      # The sed - bc - awk combination is needed because the values are returned in scientific notation and I have not found an easy way of
      # coverting them. In all likelihood this can be simplified

      read mean_g std_g numStars_g <<< "$parameters_g"
      mean_g=$( echo $mean_g | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l | awk '{printf "%.8f\n", $0}')
      std_g=$( echo $std_g | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l | awk '{printf "%.8f\n", $0}')

      read mean_r std_r numStars_r <<< "$parameters_r"
      mean_r=$( echo $mean_r | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l | awk '{printf "%.8f\n", $0}')
      std_r=$( echo $std_r | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l | awk '{printf "%.8f\n", $0}')

      # ****** Decision note *******
      # The criteria for selecting the best brick is having at least 50 stars and, from those that fulfill this condition, the lower std
      # where std is std_r + std_g
      std=$(echo "$std_g + $std_r" | bc -l )
      if (( $(echo "$std < $minStd" | bc -l) )) && (( $(echo "$numStars_g > 50" | bc -l) )) && (( $(echo "$numStars_r > 50" | bc -l) )); then
        minStd=$std
        bestBrick=$i
      fi
    done
    echo "$bricks:$bestBrick" >> $bestBrickRecord
    echo $bestBrick
  fi
}

downloadDecalsData() {
  referenceImagesForMosaic=$1
  mosaicDir=$2
  decalsImagesDir=$3
  frameBrickCorrespondenceFile=$4
  filters=$filters
  ringFile=$ringFile

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
      bricksOfTheFrame=$(python3 downloadBricksForFrame.py $referenceImagesForMosaic/$base $ringFile $filters $decalsImagesDir)
      echo $base $bricksOfTheFrame >> $frameBrickCorrespondenceFile     # Creating the map between frames and bricks to recover it in the photometric calibration
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
  # decalsPxScale=0.262 #arcsec/px
  # mosaicSize=38000 # For our field

  mosaicDir=$1
  decalsImagesDir=$2
  swarpcfg=$3
  ra=$4
  dec=$5

  buildMosaicDone=$mosaicDir/swarped/done_t.xt
  swarpedImagesDir=$mosaicDir/swarped
  if ! [ -d $swarpedImagesDir ]; then mkdir $swarpedImagesDir; fi
  if [ -f $buildMosaicDone ]; then
    echo -e "\nMosaic already built\n"
  else
    scaleFactor=0.25
    decalsPxScale=1.048
    mosaicSize=9500

    for a in $(ls -v $decalsImagesDir/*_g+r_div2.fits); do
      downSampledImages="$swarpedImagesDir/downSampled_$(basename $a)"
      astwarp $a --scale=$scaleFactor -o $downSampledImages
      SWarp -c $swarpcfg $downSampledImages -CENTER $ra,$dec -IMAGE_SIZE $mosaicSize,$mosaicSize -IMAGEOUT_NAME $swarpedImagesDir/swarp1.fits \
                -WEIGHTOUT_NAME $swarpedImagesDir/swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE $decalsPxScale -PIXELSCALE_TYPE  MANUAL
      astarithmetic $swarpedImagesDir/swarp_w1.fits -h0 set-i i i 0 lt nan where -otemp1.fits
      astarithmetic $swarpedImagesDir/swarp1.fits -h0 temp1.fits -h1 0 eq nan where -o$swarpedImagesDir/swarped_"$(basename $a)"
      rm $downSampledImages
    done
    rm $swarpedImagesDir/swarp_w1.fits $swarpedImagesDir/swarp1.fits ./temp1.fits

    sigma=2
    astarithmetic $(ls -v $swarpedImagesDir/*.fits) $(ls -v $swarpedImagesDir/*.fits | wc -l) -g1 $sigma 0.2 sigclip-median -o $mosaicDir/mosaic.fits
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
    astquery gaia --dataset=edr3 --overlapwith=$ref --column=ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error  -o$BDIR/Gaia_eDR3_.fits
    asttable $BDIR/Gaia_eDR3_.fits -c1,2,3 -c'arith $4 abs' -c'arith $5 3 x' -c'arith $6 abs' -c'arith $7 3 x' -c'arith $8 abs' -c'arith $9 3 x' --noblank=4 -otmp.txt

    # Here I demand that the gaia object fulfills simultaneously that:
    # 1.- Parallax > 3 times its error
    # 2.- Proper motion (ra) > 3 times its error
    # 3.- Proper motion (dec) > 3 times its error
    asttable tmp.txt -c1,2,3 -c'arith $4 $4 $5 gt 1000 where' -c'arith $6 $6 $7 gt 1000 where' -c'arith $8 $8 $9 gt 1000 where'  -otest_.txt
    asttable test_.txt -c1,2,3 -c'arith $4 $5 + $6 +' -otest1.txt
    asttable test1.txt -c1,2,3 --range=ARITH_2,2999,3001 -o $BDIR/Gaia_eDR3.fits

    rm test1.txt tmp.txt $BDIR/Gaia_eDR3_.fits test_.txt
    echo done > $retcatdone
  fi
}

selectStarsAndSelectionRangeDecals() {
  framesForCalibrationDir=$1
  mosaicDir=$2
  decalsImagesDir=$3
  frameBrickCorrespondenceFile=$4
  selectedDecalsStarsDir=$5
  rangeUsedDecalsDir=$6
  bestBrickRecord=$7
  frameChosenBrickMap=$8

  tmpFolder="./tmpFilesForPhotometricCalibration"
  selectDecalsStarsDone=$selectedDecalsStarsDir/automaticSelection_done.txt

  # Code to store clumps catalogue. So I can check the range for the photometric calibration
  # decalsClumpsCat=$mosaicDir/clumpsCats
  # if ! [ -d $decalsClumpsCat ]; then mkdir $decalsClumpsCat; fi

  if ! [ -d $rangeUsedDecalsDir ]; then mkdir $rangeUsedDecalsDir; fi
  if ! [ -d $selectedDecalsStarsDir ]; then mkdir $selectedDecalsStarsDir; fi
  if ! [ -d $tmpFolder ]; then mkdir $tmpFolder; fi

  if [ -f $selectDecalsStarsDone ]; then
      echo -e "\nDecals bricks and stars for doing the photometric calibration are already selected for each frame\n"
  else
    # Remove the record and map files. This is done to avoid problems with files from previous executions
    rm $bestBrickRecord
    rm $frameChosenBrickMap

    astmkprof --kernel=gaussian,1.5,3 --oversample=1 -o ./kernel.fits 1>/dev/null
    for a in $(seq 1 $totalNumberOfFrames); do
      base="entirecamera_"$a.fits
      bricks=$( getBricksWhichCorrespondToFrame $framesForCalibrationDir/$base $frameBrickCorrespondenceFile )
      chosenBrick=$( chooseBrickToUse "${bricks[@]}" $decalsImagesDir ./kernel.fits $tmpFolder $bestBrickRecord)
      echo "$base:$chosenBrick" >> $frameChosenBrickMap

      echo "Image $a; Bricks $bricks; Chosen $chosenBrick"

      # The following conditional is not to compute twice the same decals (g+r)/2 frame
      matchingSelectedStars=$(find $selectedDecalsStarsDir -maxdepth 1 -type f -name "*${chosenBrick}*" | head -n 1)
      matchingRanges=$(find $rangeUsedDecalsDir -maxdepth 1 -type f -name "*${chosenBrick}*" | head -n 1)
      if [ -n "$matchingSelectedStars" ]; then
        cp $matchingSelectedStars $selectedDecalsStarsDir/selected_decalsStarsForFrame_"$a"_andBrick_"$chosenBrick".txt
        cp $matchingRanges $rangeUsedDecalsDir/selected_rangeForFrame_"$a"_andBrick_"$chosenBrick".txt

        #
        # tmp=$(( $a - 1 ))
        # cp $decalsClumpsCat/decals_"$tmp".fits_$chosenBrick.txt $decalsClumpsCat/decals_"$base"_$chosenBrick.txt
        #

      else
        decalsImageToUse=$decalsImagesDir/"decal_image_"$chosenBrick"_g+r_div2.fits"
        astconvolve $decalsImageToUse --kernel=./kernel.fits --domain=spatial --output=$tmpFolder/convolved.fits
        astnoisechisel $decalsImageToUse -h1 -o $tmpFolder/det.fits --convolved=$tmpFolder/convolved.fits --tilesize=30,30 --detgrowquant=0.95 --erode=4
        astsegment $tmpFolder/det.fits -o $tmpFolder/seg.fits --snquant=0.1 --gthresh=-10 --objbordersn=0  --minriverlength=3
        astmkcatalog $tmpFolder/seg.fits --ra --dec --magnitude --half-max-radius --sum --clumpscat -o $tmpFolder/decals.txt --zeropoint=22.5

        #
        # echo cp $tmpFolder/decals_c.txt $decalsClumpsCat/decals_"$base"_$chosenBrick.txt
        # cp $tmpFolder/decals_c.txt $decalsClumpsCat/decals_"$base"_$chosenBrick.txt
        #

        astmatch $tmpFolder/decals_c.txt --hdu=1  $BDIR/Gaia_eDR3.fits --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=bRA,bDEC,aHALF_MAX_RADIUS,aMAGNITUDE -o $tmpFolder/match_decals_gaia.txt 1>/dev/null

        s=$(asttable $tmpFolder/match_decals_gaia.txt -h1 -c3 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
        std=$(asttable $tmpFolder/match_decals_gaia.txt -h1 -c3 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
        minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
        maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)

        echo $s $std $minr $maxr > $rangeUsedDecalsDir/selected_rangeForFrame_"$a"_andBrick_"$chosenBrick".txt
        asttable $tmpFolder/decals_c.txt --range=HALF_MAX_RADIUS,$minr,$maxr -o $selectedDecalsStarsDir/selected_decalsStarsForFrame_"$a"_andBrick_"$chosenBrick".txt
        rm $tmpFolder/*
      fi
    done
    rmdir $tmpFolder
    rm ./kernel.fits
    echo done > $selectDecalsStarsDone
  fi


}

prepareDecalsDataForPhotometricCalibration() {
  referenceImagesForMosaic=$1
  decalsImagesDir=$2
  ra=$3
  dec=$4
  mosaicDir=$5
  selectedDecalsStarsDir=$6
  rangeUsedDecalsDir=$7
  bestBrickRecord=$8
  frameChosenBrickMap=$9


  echo -e "\n ${GREEN} ---Preparing Decals data--- ${NOCOLOUR}"

  frameBrickCorrespondenceFile=$decalsImagesDir/frameBrickMap.txt
  ringFile=$CDIR/flat_ring_ccd"$h".txt

  # This first steps donwloads the decals frames that are needed for calibrating each of our frames
  # The file "frameBrickCorrespondenceFile" will store the correspondence between our frames and decals bricks
  # This way we can easily access to this relation in orther to do the photometric calibration

  # If the images are donwloaded but the done.txt file is no present, the images won't be donwloaded again but
  # the step takes a while even if the images are already downloaded because we have to do a query to the database
  # in order to obtain the brickname and check if it is already downloaded or no

  # We download 'g' and 'r' because our images are taken with a luminance filter which response is a sort of g+r
  filters="g,r"
  downloadDecalsData $referenceImagesForMosaic $mosaicDir $decalsImagesDir $frameBrickCorrespondenceFile $filters $ringFile

  # This step creates the images (g+r)/2. This is needed because we are using a luminance filter which is a sort of (g+r)
  # The division by 2 is because in AB system we work with Janskys, which are W Hz^-1 m^-2. So we have to give a flux per wavelenght
  # So, when we add two filters we have to take into account that we are increasing the wavelength rage. In our case, 'g' and 'r' have
  # practically the same wavelenght width, so dividing by 2 is enough
  add_gAndr_andDivideByTwo $decalsImagesDir

  # The photometric calibration is frame by frame, so we are not going to use the mosaic for calibration. But we build it anyway to retrieve in an easier way
  # the Gaia data of the whole field.
  buildDecalsMosaic $mosaicDir $decalsImagesDir $swarpcfg $ra $dec

  echo -e "\n ${GREEN} ---Downloading GAIA catalogue for our field --- ${NOCOLOUR}"
  downloadGaiaCatalogueForField $mosaicDir


  # First of all remember that we need to do the photometric calibration frame by frame.
  # For each frame of the data, due to its large field, we have multiple (in our case 4 - defined in "downloadBricksForFrame.py") decals bricks
  # They may have been taken at different moments with different conditions so we have to process them one by one

  # With only one field is enough, so in order to avoid extra work having to treat the four bricks independently, we will
  # choose one of them (based on the std and defining a minimum number of stars needed) and do the calibration of the frame with that brick
  echo -e "\n ${GREEN} --- Selecting stars and star selection range for Decals--- ${NOCOLOUR}"
  selectStarsAndSelectionRangeDecals $referenceImagesForMosaic $mosaicDir $decalsImagesDir $frameBrickCorrespondenceFile $selectedDecalsStarsDir $rangeUsedDecalsDir $bestBrickRecord $frameChosenBrickMap
}

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
  astsegment det_"$a".fits -o seg_"$a".fits --snquant=0.1 --gthresh=-10 --objbordersn=0  --minriverlength=3
  astmkcatalog seg_"$a".fits --ra --dec --magnitude --half-max-radius --sum --clumpscat -o $mycatdir/"$base".txt
  astmatch $BDIR/Gaia_eDR3.fits --hdu=1 $mycatdir/"$base"_c.txt --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aRA,aDEC,bMAGNITUDE,bHALF_MAX_RADIUS -o$mycatdir/match_"$base"_my_gaia.txt

  s=$(asttable $mycatdir/match_"$base"_my_gaia.txt -h1 -c4 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
  std=$(asttable $mycatdir/match_"$base"_my_gaia.txt -h1 -c4 --noblank=MAGNITUDE | aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
  minr=$(astarithmetic $s $sigmaForPLRegion $std x - -q)
  maxr=$(astarithmetic $s $sigmaForPLRegion $std x + -q)
  echo $s $std $minr $maxr > $mycatdir/range1_"$base".txt
  asttable $mycatdir/"$base"_c.txt  --range=HALF_MAX_RADIUS,$minr,$maxr -o $mycatdir/selected_"$base"_automatic.txt
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

      # match the aoutomatics catalogs THIS WAY I SELECT NON SATURATED
      out_auto=$matchdir2/match-decals-"$base"_automatic.cat
      astmatch $selectedDecalsStarsDir/selected_decalsStarsForFrame_"$a"_*.txt --hdu=1 $mycatdir/selected_"$base"_automatic.txt --hdu=1 --ccol1=RA,DEC --ccol2=RA,DEC --aperture=$toleranceForMatching/3600 --outcols=aRA,aDEC,aMAGNITUDE,aHALF_MAX_RADIUS,bMAGNITUDE,bHALF_MAX_RADIUS -o$out_auto
    done
    echo done > $matchdir2done
  fi
}

buildDecalsCatalogueOfMatchedSources() {
  decalsdir=$1
  rangeUsedDecalsDir=$2
  matchdir2=$3
  frameChosenBrickMap=$4
  decalsImagesDir=$5

  decalsdone=$decalsdir/done__ccd"$h".txt
  if ! [ -d $decalsdir ]; then mkdir $decalsdir; fi
  if [ -f $decalsdone ]; then
    echo -e "\nDecals: catalogue for the calibration stars already built\n"
  else
    for a in $(seq 1 $totalNumberOfFrames); do
      # I have to take 2 the FWHM (half-max-rad)
      # It is already saved as mean value of the point-like sources
      r_decals_pix_=$(awk 'NR==1 {printf $1}' $rangeUsedDecalsDir/"selected_rangeForFrame_"$a"_andBrick_"*".txt")
      r_decals_pix=$(astarithmetic $r_decals_pix_ 2. x -q )

      base="entirecamera_$a.fits"
      out=$matchdir2/match-decals-"$base"_automatic.cat
      decalsBrick=$( retrieveChosenDecalsBrickForFrame $base $frameChosenBrickMap )
      decalsBrick=$decalsImagesDir"/decal_image_"$decalsBrick"_g+r_div2.fits"

      #This seems extremely inneficient but maybe is the way idk
      asttable $out -hSOURCE_ID -cRA,DEC | awk '!/^#/{print NR, $1, $2, 5, '$r_decals_pix', 0, 0, 1, NR, 1}' > apertures_decals.txt
      astmkprof apertures_decals.txt --background=$decalsBrick --backhdu=1 \
          --clearcanvas --replace --type=int16 --mforflatpix \
          --mode=wcs --output=apertures_decals.fits
      astmkcatalog apertures_decals.fits -h1 --zeropoint=22.5 \
              --valuesfile=$decalsBrick --valueshdu=1 \
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

      mean=$(asttable $alphatruedir/$base -c'ARITH_1'| aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-median)
      std=$(asttable $alphatruedir/$base -c'ARITH_1'| aststatistics --sclipparams=$sigmaForStdSigclip,$iterationsForStdSigClip --sigclip-std)
      echo "$mean $std" > $alphatruedir/alpha_"$objectName"_Decals-"$filter"_"$a"_ccd"$h".txt

    done
    echo done > $alphatruedone
  fi
}

computeCalibrationFactors() {
  iteration=$1
  imagesForCalibration=$2
  selectedDecalsStarsDir=$3
  rangeUsedDecalsDir=$4
  frameChosenBrickMap=$5
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
  matchDecalsAndOurData $iteration $selectedDecalsStarsDir $mycatdir $matchdir2

  decalsdir=$BDIR/decals-aperture-catalog_it$iteration
  echo -e "\n ${GREEN} ---Building Decals catalogue for (matched) calibration stars--- ${NOCOLOUR}"
  buildDecalsCatalogueOfMatchedSources $decalsdir $rangeUsedDecalsDir $matchdir2 $frameChosenBrickMap $decalsImagesDir

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
  brightLimit=16.3
  faintLimit=18.5
  computeAndStoreFactors $alphatruedir $matchdir2 $brightLimit $faintLimit
}

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
  f=$photCorrDir/$base
  rms_min=$(awk 'NR=='1'{print $1}' $BDIR/rms_min_val-1_ccd"$h"_it$iteration.txt)
  rms_f=$(awk 'NR=='1'{print $3}' $noiseskydir/entirecamera_$a.txt)

  weight=$(astarithmetic $rms_min $rms_f / --quiet)
  echo "$weight" > $wdir/"$objectName"_Decals-"$filter"_"$a"_ccd"$h".txt    #  saving into file

  # multiply each image for its weight
  wixi_im=$wdir/$base           # frame x weight
  w_im=$wonlydir/$base           # only weight
  astarithmetic $f -h1 $weight x -o$wixi_im
  astarithmetic $wixi_im -h1 $f -h1 / -o$w_im
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
  astarithmetic $wonlydir/$base $mask -g1 1 eq nan where -q float32  -o $owom

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

  frameToMask=$dirOfFramesToMask/entirecamera_$a.fits
  tmpMaskFile=$dirOfFramesMasked/"maskFor"$a.fits

  # Parameters for identifing our frame in the full grid
  frameCentre=$( getCentralCoordinate $frameToMask )
  centralRa=$(echo "$frameCentre" | awk '{print $1}')
  centralDec=$(echo "$frameCentre" | awk '{print $2}')

  # Offset for having some margin and not lose any pixel
  securityOffset=50
  detectorWidthDeg=$(echo  "(($detectorWidth + $securityOffset) * $pixelScale)" | bc )
  detectorHeightDeg=$(echo "(($detectorHeight + $securityOffset) * $pixelScale) + $securityOffset" | bc )

  astcrop $wholeMask --center=$centralRa,$centralDec --mode=wcs --width=$detectorHeightDeg/3600,$detectorWidthDeg/3600  -o $tmpMaskFile
  astarithmetic $frameToMask -h1 $tmpMaskFile -h1 1 eq nan where float32 -o $dirOfFramesMasked/entirecamera_$a.fits -q
  rm $tmpMaskFile
}
export -f cropAndApplyMaskPerFrame

maskPointings() {
  entiredir_smallGrid=$1
  smallPointings_maskedDir=$2
  maskedPointingsDone=$3
  maskName=$4

  if ! [ -d $smallPointings_maskedDir ]; then mkdir $smallPointings_maskedDir; fi
  if [ -f $maskedPointingsDone ]; then
      echo -e "\nThe masks for the pointings have been already applied\n"
  else
    framesToMask=()
    for a in $(seq 1 $totalNumberOfFrames); do
      framesToMask+=("$a")
    done
    printf "%s\n" "${framesToMask[@]}" | parallel -j "$num_cpus" cropAndApplyMaskPerFrame {} $entiredir_smallGrid $smallPointings_maskedDir $maskName
    echo done > $maskedPointingsDone 
  fi
}

# Coadd function
buildCoadd() {
  baseCoaddir=$1
  mowdir=$2
  moonwdir=$3

  if ! [ -d $baseCoaddir ]; then mkdir $baseCoaddir; fi

  coaddir=$baseCoaddir/stdWeighted1
  coaddName=$coaddir/"$objectName"_coadd1_"$filter".fits
  coaddone=$coaddir/done_"$k".txt
  if ! [ -d $coaddir ]; then mkdir $coaddir; fi
  if [ -f $coaddone ]; then
      echo -e "\nThe first weighted (based upon std) mean of the images already done\n"
  else
      astarithmetic $(ls -v $mowdir/*.fits) $(ls $mowdir/*.fits | wc -l) sum -g1 -o$coaddir/"$k"_wx.fits
      astarithmetic $(ls -v $moonwdir/*.fits ) $(ls $moonwdir/*.fits  | wc -l) sum -g1 -o$coaddir/"$k"_w.fits
      astarithmetic $coaddir/"$k"_wx.fits -h1 $coaddir/"$k"_w.fits -h1 / -o$coaddName
      echo done > $coaddone
  fi
}

########## Variables ##########

objectName=Fornax
filter=lum

ra_gal=39.9515
dec_gal=-34.5110

ROOTDIR=/scratch/sguerra

defaultNumOfCPUs=12
num_cpus=$SLURM_CPUS_ON_NODE
if [ -z $num_cpus ]; then
  num_cpus=$defaultNumOfCPUs
fi

echo "Number of CPUs allocated: $num_cpus"

export objectName
export filter
export ra_gal
export dec_gal

pixelScale=1.164 # arcsec / px
detectorWidth=4656
detectorHeight=3520

export detectorWidth
export detectorHeight
export pixelScale

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
noisechisel_param="--tilesize=25,25  \
                    --meanmedqdiff=0.01 \
                    --detgrowquant=0.7 \
                    --detgrowmaxholesize=1000 \
                    --qthresh=0.25 \
                    --snquant=0.98 \
                    --rawoutput"
export noisechisel_param

echo -e "\n-Noisechisel parameters used for masking:"
echo -e "\t" $noisechisel_param


coaddSizePx=9001
echo -e "The size in px of each side of the coadded image is " $coaddSizePx
export coaddSizePx

# This variables account for the falt that is used in the reduction
# Depending on your dithering patter you may want to perform running flat or the typical whole night flat (RUNNING_FLAT variable)
# The windowSize determines how many flats are you going to use for the running flat.
RUNNING_FLAT=false
windowSize=11
halfWindowSize=5

export RUNNING_FLAT
export windowSize
export halfWindowSize

echo -e "\nThe running flat is going to be used?: $ORANGE $RUNNING_FLAT $NOCOLOUR"
echo -e "If so, the running flat will be computed with a window size of " $windowSize
echo -e "\n"


########## Prepare data ##########

echo -e "\n ${GREEN} ---Preparing data--- ${NOCOLOUR}"

DIR=$ROOTDIR/"$objectName"
INDIRo=$ROOTDIR/"$objectName"/DATA-or
BDIR=$ROOTDIR/"$objectName"/build
INDIR=$ROOTDIR/"$objectName"/DATA
DARKDIR=$ROOTDIR/"$objectName"/dark
keyWordDirectory=$ROOTDIR/"$objectName"/keywords

export ROOTDIR
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

dateHeaderKey="DATE-OBS"
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
  # sorted and renamed based in an objetive criteria es the time in which the frame was taken
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
      astarithmetic $(ls -v $currentDARKDIR/* ) \
                    $(ls -v $currentDARKDIR/* | wc -l) \
                    3 0.2 sigclip-mean -g$h \
                    -o $mdadir/mdark_"$filter"_n"$currentNight"_ccd$h.fits
    fi
    echo done > $mdadone
  done


  ########## Save airmass ##########
  echo -e "\n ${GREEN} Saving airmass ${NOCOLOUR}"

  skydir=$BDIR/airmass-analysis_n$currentNight
  skydone=$skydir/done_.txt
  if ! [ -d $skydir ]; then mkdir $skydir; fi
  if [ -f $skydone ]; then
    echo -e "\nAirmass for night $currentNight already saved\n"
  else
      for i in $(ls -v $currentINDIR/*.fits ); do
        air=$(astfits $i -h1 --keyvalue=AIRMASS | awk '{print $2}')
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
                i i 55000 gt i isblank or 2 dilate nan where m -  float32  \
                -o $mbiascorrdir/$base

    done
    echo done > $mbiascorrdone
  fi

  # until here everything is corrected of bias and renamed
  echo -e "${ORANGE} ------ FLATS ------ ${NOCOLOUR}\n"

  #now flat
  echo -e "${GREEN} --- Flat iteration 1 --- ${NOCOLOUR}"


  ########## Creating the ring mask ##########

  ####### THIS HAS TO BE CHECKED, MAYBE USING A NORMAL MEDIAN IS BETTER TO NORMALIZE
  ringdir=$BDIR/ring
  ringFile=$CDIR/flat_ring_ccd"$h".txt
  ringdone=$ringdir/ring_ccd"$h".txt
  if ! [ -d $ringdir ]; then mkdir $ringdir; fi
  if [ -f $ringdone ]; then
    echo -e "\nRing mask is already created for extension $h\n"
  else
    astmkprof --background=$mbiascorrdir/"$objectName"-Decals-"$filter"_n"$currentNight"_f1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=200 --clearcanvas -o $ringdir/ring_ccd"$h".fits $ringFile
    echo "done" >> $ringdone
  fi


  ########## Creating the it1 master flat image ##########

  # ****** Decision note *******
  # Running flat: summary
  # A running flat is a flat built with not all the frames of one night, but each frame has a local flat built from N frames based on time-near frames.
  # Caveat: This requieres to check the data of the night, this only works if the data has been taken at similar times (check DATE-OBS or the airmass)
  # The size of the window for the running flat can vary, here a window of 11 is used because I have been told that that's the minimum number of frames
  # that can be successfully used for building the flat. But in general is a compromise between having a lot of frames and characterising the local sky conditions

  # This pipeline also computes a whole night flat.
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
    normaliseImagesWithRing $mbiascorrdir $normit1dir $ringdir/ring_ccd"$h".fits
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
      for a in $(seq 1 $n_exp); do
          base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
          i=$flatit1imadir/$base
          astnoisechisel $i $noisechisel_param -o $noiseit2dir/$base
      done
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
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      i=$flatit1WholeNightimaDir/$base
      astnoisechisel $i $noisechisel_param -o $noiseit2WholeNightDir/$base
    done
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
      maskImages $mbiascorrdir $noiseit2dir $maskedit2dir
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
    maskImages $mbiascorrdir $noiseit2WholeNightDir $maskedit2WholeNightdir
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
      normaliseImagesWithRing $maskedit2dir $normit2dir $ringdir/ring_ccd"$h".fits
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
    normaliseImagesWithRing $maskedit2WholeNightdir $normit2WholeNightdir $ringdir/ring_ccd"$h".fits
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
        for a in $(seq 1 $n_exp); do
            base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
            i=$flatit2imadir/$base
            astnoisechisel $i $noisechisel_param -o $noiseit3dir/$base
      done
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
    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      i=$flatit2WholeNightimaDir/$base
      astnoisechisel $i $noisechisel_param -o $noiseit3WholeNightDir/$base
    done
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
      maskImages $mbiascorrdir $noiseit3dir $maskedit3dir
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
    maskImages $mbiascorrdir $noiseit3WholeNightDir $maskedit3WholeNightdir
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
      normaliseImagesWithRing $maskedit3dir $normit3dir $ringdir/ring_ccd"$h".fits
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
    normaliseImagesWithRing $maskedit3WholeNightdir $normit3WholeNightdir $ringdir/ring_ccd"$h".fits
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
      astfits $i --copy=1 -o$out
    done
    echo done > $maskedcornerdone
  fi


  # At this point we can process the frames of all the nights in the same way
  # So we place all the final frames into a common folder.

  # WE SHOULD LOCK THE FILE IN ORDER TO AVOID RACE CONDITIONS
  if [ -f $framesForCommonReductionDone ]; then
    echo -e "\nFrames already placed in the folder for frames prepared to common reduction"
  else
    initialValue=$( getHighestNumberFromFilesInFolder $framesForCommonReductionDir )

    for a in $(seq 1 $n_exp); do
      base="$objectName"-Decals-"$filter"_n"$currentNight"_f"$a"_ccd"$h".fits
      name=$(( $initialValue + $a ))
      mv $maskedcornerdir/$base $framesForCommonReductionDir/$name.fits
    done
    echo "done" > $framesForCommonReductionDone
  fi


  # Removing everything but the final frames
  rm -rf $BDIR/bias-corrected_n$currentNight
  rm -rf $BDIR/masked-corner_n$currentNight
  rm -rf $BDIR/masterdark_n$currentNight
  rm -rf $BDIR/flat-it3-Running-BeforeCorrection_n$currentNight
  rm -rf $BDIR/flat-it3-ima_n$currentNight
  for a in $(seq 1 3); do
    rm -rf $BDIR/flat-it"$a"-Running_n$currentNight
    rm -rf $BDIR/flat-it"$a"-Running-ima_n$currentNight
    rm -rf $BDIR/flat-it"$a"-WholeNight_n$currentNight
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
printf "%s\n" "${nights[@]}" | parallel -j "$num_cpus" oneNightPreProcessing {}

totalNumberOfFrames=$( ls $framesForCommonReductionDir/*.fits | wc -l)
export totalNumberOfFrames
echo $totalNumberOfFrames





# Up to this point the frame of every night has been corrected of bias-dark and flat.
# That corrections are perform night by night (because it's necessary for perform that corretions)
# Now, all the frames are "equall" so we do no distinction between nights.
# All the frames are stored together in $framesForCommonReductionDir with names 1.fits, 2.fits, 3.fits ... n.fits.

echo -e "${ORANGE} ------ ASTROMETRY AND BACKGROUND-SUBTRACTION ------ ${NOCOLOUR}\n"
for h in 0; do
  ########## Astrometry ##########
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
    for re in $(seq -5 6); do
      build-astrometry-index -i $catdir/"$objectName"_Gaia_eDR3.fits -e1 \
                              -P $re \
                              -S phot_g_mean_mag \
                              -E -A RA -D  DEC\
                              -o $indexdir/index_$re.fits;
    done
    echo done > $indexdone
  fi

  sexcfg=$CDIR/default.sex
  # Solving the images
  astrocfg=$CDIR/astrometry_$objectName.cfg
  if [ -f $astrocfg ]; then
    echo -e "\nAstrometry config file are already created\n"
  else
    echo inparallel > $astrocfg
    echo cpulimit 300 >> $astrocfg
    echo "add_path $indexdir" >> $astrocfg
    echo autoindex >> $astrocfg
  fi
  astroimadir=$BDIR/astro-ima
  astroimadone=$astroimadir/done_"$filter"_ccd"$h".txt
  if ! [ -d $astroimadir ]; then mkdir $astroimadir; fi
  if [ -f $astroimadone ]; then
    echo -e "\nImages are already astrometrized for extension $h\n"
  else
    for a in $(seq 1 $totalNumberOfFrames); do
        base=$a.fits
        i=$framesForCommonReductionDir/$base
        solve-field $i --no-plots \
        -L 1.15 -H 1.18 -u arcsecperpix \
        --ra=$ra_gal --dec=$dec_gal --radius=3. \
        --overwrite --extension 1 --config $astrocfg --no-verify -E 3 -c 0.03 \
        --odds-to-solve 1e7 \
        --use-source-extractor --source-extractor-path=/usr/bin/source-extractor \
        -Unone --temp-axy -Snone -Mnone -Rnone -Bnone -N$astroimadir/$base ;

        # --source-extractor-config=./config/default.sex \
        # AQUÍ HAY ALGO RARO - cuando uso la configuración por defecto de sextractor funciona. cuando uso la configuración del ./config no...
    done
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
    for a in $(seq 1 $totalNumberOfFrames); do
        base="$a".fits
        i=$astroimadir/$base
        source-extractor $i -c $sexcfg -PARAMETERS_NAME $sexparam -FILTER_NAME $sexconv -CATALOG_NAME $sexdir/$a.cat
    done
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

# For fornax (around 490 frames). Deimos, 20 cores -> 1h and 20 min

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
  printf "%s\n" "${imagesToWarp[@]}" | parallel -j "$num_cpus" warpImage {} $entiredir_fullGrid $entiredir_smallGrid $ra $dec $coaddSizePx 
  echo done > $entiredone
fi


# For fornax (around 490 frames). Deimos, 20 cores -> 25 min
echo -e "${GREEN} --- Compute and subtract Sky --- ${NOCOLOUR}"


# totalNumberOfFrames=25
# export totalNumberOfFrames
echo $totalNumberOfFrames

noiseskydir=$BDIR/noise-sky_it1
noiseskydone=$noiseskydir/done_"$filter"_ccd"$h".txt
computeSky $entiredir_smallGrid $noiseskydir $noiseskydone 

subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it1
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter"_ccd"$h".txt
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir

subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it1
subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter"_ccd"$h".txt
subtractSky $entiredir_fullGrid $subskyFullGrid_dir $subskyFullGrid_done $noiseskydir



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

# The file "bestBrickRecord" will store the sets of decals bricks already evaluated and the decision made
# This way we avoid computing the same more than once and saves a considerable amount of time
# -
# The file "frameChosenBrickMap" will store the correspondence between our frames and the chosen decals brick
# It's true that "frameChosenBrickMap" does not store information not present in "bestBrickRecord" together with "frameBrickMap"
# But it makes our life easier and the amount of space is negligible (a plain text file with some hundreds of lines)
bestBrickRecord=$decalsImagesDir/bestBrickRecord.txt
frameChosenBrickMap=$decalsImagesDir/frameChosenBrickMap.txt

prepareDecalsDataForPhotometricCalibration $referenceImagesForMosaic $decalsImagesDir $ra $dec $mosaicDir $selectedDecalsStarsDir $rangeUsedDecalsDir $bestBrickRecord $frameChosenBrickMap



iteration=1
imagesForCalibration=$subskySmallGrid_dir
alphatruedir=$BDIR/alpha-stars-true_it$iteration
computeCalibrationFactors $iteration $imagesForCalibration $selectedDecalsStarsDir $rangeUsedDecalsDir $frameChosenBrickMap $decalsImagesDir $alphatruedir


echo -e "\n ${GREEN} ---Applying calibration factors--- ${NOCOLOUR}"
photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration
applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir
applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir


echo -e "${ORANGE} ------ STD WEIGHT COMBINATION ------ ${NOCOLOUR}\n"

# Compute rms and of the photometrized frames
noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done_"$k"_ccd"$h".txt
computeSky $photCorrSmallGridDir $noiseskydir $noiseskydone

# Store the minimum standard deviation of the frames in order to compute the weights
python3 find_rms_min.py $filter 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration

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


# Fornax. Around 490 frames. Deimos, 20 cores. Around 1 h and 15 min
mowdir=$BDIR/weight-dir-no-outliers
if ! [ -d $mowdir ]; then mkdir $mowdir; fi
# only weight
moonwdir=$BDIR/only-weight-dir-no-outliers
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



### Coadd ###
echo -e "\n ${GREEN} ---Coadding--- ${NOCOLOUR}"

baseCoaddir=$BDIR/coadds
buildCoadd $baseCoaddir $mowdir $moonwdir

maskName=$coaddir/"$objectName"_coadd1_"$filter"_mask.fits
if [ -f $maskName ]; then
  echo "The mask of the weighted coadd is already done"
else
  astnoisechisel $coaddName $noisechisel_param -o $maskName
fi



# # --- Build exposure map
exposuremapDir=$baseCoaddir/exposureMap
exposuremapdone=$baseCoaddir/done_"$k".txt

if ! [ -d $exposuremapDir ]; then mkdir $exposuremapDir; fi
if [ -f $exposuremapdone ]; then
    echo -e "\nThe first weighted (based upon std) mean of the images already done\n"
else
  #There should be more efficient way of doing this...
  # Pure exposure map
  framesDir=$BDIR/pointings_smallGrid
  for a in $(seq 1 $totalNumberOfFrames); do
    i=$framesDir/entirecamera_$a.fits
    astarithmetic $i set-i i i i eq 1 where i isblank 1 where -g1 --output="./tmp.fits"
    SWarp -c $swarpcfg "./tmp.fits" -CENTER $ra,$dec -IMAGE_SIZE $coaddSizePx,$coaddSizePx -IMAGEOUT_NAME $exposuremapDir/swarp1.fits -WEIGHTOUT_NAME $exposuremapDir/swarp_w1.fits -SUBTRACT_BACK N -PIXEL_SCALE 1.164 -PIXELSCALE_TYPE  MANUAL
    astarithmetic $exposuremapDir/swarp_w1.fits -h0 set-i i i 0 lt nan where -otemp1.fits
    astarithmetic $exposuremapDir/swarp1.fits -h0 temp1.fits -h1 0 eq nan where -o$exposuremapDir/entirecamera_"$a".fits
  done
  rm $exposuremapDir/swarp_w1.fits $exposuremapDir/swarp1.fits
  astarithmetic $(ls -v $exposuremapDir/*.fits) $(ls $exposuremapDir/*.fits | wc -l) number -g1 -o$baseCoaddir/exposureMap_NoNans.fits
  rm $exposuremapDir/*.fits
  rmdir $exposuremapDir

  echo done > $exposuremapdone
fi



# Remove intermediate folders to save some space
rm -rf $BDIR/sub-sky-fullGrid_it1
rm -rf $BDIR/sub-sky-smallGrid_it1
rm -rf $BDIR/photCorrFullGrid-dir_it1
rm -rf $BDIR/photCorrSmallGrid-dir_it1

rm -rf $wdir
rm -rf $wonlydir
rm -rf $mowdir
rm -rf $moonwdir

####### ITERATION 2 ######

iteration=2

# We mask the pointings in order to measure (before photometric calibration) the sky accurately
smallPointings_maskedDir=$BDIR/pointings_smallGrid_masked
maskedPointingsDone=$smallPointings_maskedDir/done_.txt
maskPointings $entiredir_smallGrid $smallPointings_maskedDir $maskedPointingsDone $maskName

noiseskydir=$BDIR/noise-sky_it$iteration
noiseskydone=$noiseskydir/done_"$filter"_ccd"$h".txt
computeSky $smallPointings_maskedDir $noiseskydir $noiseskydone # compute sky with frames masked with global mask

subskySmallGrid_dir=$BDIR/sub-sky-smallGrid_it$iteration
subskySmallGrid_done=$subskySmallGrid_dir/done_"$filter"_ccd"$h".txt
subtractSky $entiredir_smallGrid $subskySmallGrid_dir $subskySmallGrid_done $noiseskydir

subskyFullGrid_dir=$BDIR/sub-sky-fullGrid_it$iteration
subskyFullGrid_done=$subskyFullGrid_dir/done_"$filter"_ccd"$h".txt
subtractSky $entiredir_fullGrid $subskyFullGrid_dir $subskyFullGrid_done $noiseskydir

imagesForCalibration=$subskySmallGrid_dir
alphatruedir=$BDIR/alpha-stars-true_it$iteration
computeCalibrationFactors $iteration $imagesForCalibration $selectedDecalsStarsDir $rangeUsedDecalsDir $frameChosenBrickMap $decalsImagesDir $alphatruedir

photCorrSmallGridDir=$BDIR/photCorrSmallGrid-dir_it$iteration
photCorrFullGridDir=$BDIR/photCorrFullGrid-dir_it$iteration

applyCalibrationFactors $subskySmallGrid_dir $alphatruedir $photCorrSmallGridDir
applyCalibrationFactors $subskyFullGrid_dir $alphatruedir $photCorrFullGridDir


# We mask again the points in order to measure (after photometric calibration) the sky accurately
smallPointings_photCorr_maskedDir=$BDIR/photCorrSmallGrid_masked
maskedPointingsDone=$smallPointings_photCorr_maskedDir/done_.txt
maskPointings $photCorrSmallGridDir $smallPointings_photCorr_maskedDir $maskedPointingsDone $maskName


noiseskydir=$BDIR/noise-sky-after-photometry_it$iteration
noiseskydone=$noiseskydir/done_"$k"_ccd"$h".txt
computeSky $smallPointings_photCorr_maskedDir $noiseskydir $noiseskydone

python3 find_rms_min.py "$filter" 1 $totalNumberOfFrames $h $noiseskydir $DIR $iteration



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
baseCoaddir=$BDIR/coadds_it$iteration 
buildCoadd $baseCoaddir $mowdir $moonwdir



exit 0


























