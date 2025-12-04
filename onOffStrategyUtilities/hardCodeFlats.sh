#!/bin/bash

BDIR=$1
filter=$2
buildDirWithFlatFolders=$3
numberOfNights=$4

createFolderAndDone () {
    foldName=$1
    if ! [ -d $foldName ]; then mkdir $foldName; fi
    echo "done" > $foldName/done_"$filter"_ccd0.txt
}
export -f createFolderAndDone

foldersToCreate=(
    "norm-it1-images_n"
    "flat-it1-Running_n"
    "flat-it1-WholeNight_n"
    "flat-it1-Running-ima_n"
    "flat-it1-WholeNight-ima_n"
    "noise-it2-Running_n"
    "noise-it2-WholeNight_n"
    "masked-it2-Running_n"
    "masked-it2-WholeNight_n"
    "norm-it2-Running-images_n"
    "norm-it2-WholeNight-images_n"
    "flat-it2-Running_n"
    "flat-it2-WholeNight_n"
    "flat-it2-Running-ima_n"
    "flat-it2-WholeNight-ima_n"
    "noisesky_forCleaningBadFramesBeforeFlat_n"
    "noise-it3-Running_n"
    "noise-it3-WholeNight_n"
    "masked-it3-Running_n"
    "masked-it3-WholeNight_n"
    "norm-it3-Running-images_n"
    "norm-it3-WholeNight-images_n"
    "flat-it3-Running-BeforeCorrection_n"
    "flat-it3-WholeNight_n"
    "flat-it3-Running_n"
)

if ! [ -d $BDIR ]; then mkdir $BDIR; fi


for (( currentNight=1; currentNight<=numberOfNights; currentNight++ )); do

    # Create all folders and its done files in order to be skipped by the pipeline
    for folder in "${foldersToCreate[@]}"; do
        createFolderAndDone "$BDIR/${folder}${currentNight}"
    done
    createFolderAndDone $BDIR/diagnosis_and_badFiles
    echo "done" > $BDIR/diagnosis_and_badFiles/done_badFrames_stdPreFlat_n$currentNight.txt

    # Copy the flatfields that are going to be used
    cp -r $buildDirWithFlatFolders/flat-it3-Running_n$currentNight $BDIR
done


