#!/bin/bash

escapeSpacesFromString() {
    local input="$1"
    escaped_string="${input// /\\ }"
    echo $escaped_string
}
export -f escapeSpacesFromString

getNameAndDateObs(){
    folderWithRaws=$1
    dateHeaderKey=$2
    header=$3
    local -n outputArray=$4

    for i in $folderWithRaws/*.fits; do
       
        nameWithEscapedSpaces=$(escapeSpacesFromString "$i")

        # If lightbridges couldn't reduce a frame for whatever reason, it will look like "noRed...", the DATA-OBS
        # is still present, but in the first hdu

        if [[ "${nameWithEscapedSpaces:0:5}" == "noRed" ]]; then
            DATEOBS=$(eval "astfits $nameWithEscapedSpaces -h0 --keyvalue=$dateHeaderKey --quiet")
        else
            DATEOBS=$(eval "astfits $nameWithEscapedSpaces -h$header --keyvalue=$dateHeaderKey --quiet")

        fi
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
        nameOfOriginalFile="${nameWithEscapedSpaces##*/}"

        outputArray+=("$nameOfOriginalFile" "$unixTimeInSeconds")
    done
}
export -f getNameAndDateObs

orderFramesByDate() {
    local -n array=$1

    lines=()
    for ((i=0; i<${#array[@]}; i+=2)); do
        lines+=("${array[i]} ${array[i+1]}")
    done

    IFS=$'\n' sorted_lines=($(sort -k2 <<<"${lines[*]}"))
    unset IFS

    array=()
    for line in "${sorted_lines[@]}"; do
        file=$(awk '{print $1}' <<< "$line")
        keyword=$(awk '{print $2}' <<< "$line")
        array+=("$file" "$keyword")
    done
}
export -f orderFramesByDate

getRunningFlatSizeUsed() {
    flatInfo=$1

    min=9999999
    for ((i=0; i<${#flatInfo[@]}; i+=2)); do
        filename="${flatInfo[i]}"
        if [[ $filename =~ f([0-9]+)_ccd0\.fits ]]; then
            num=${BASH_REMATCH[1]}
            (( num < min )) && min=$num
        fi
    done

    # We need to subtract 2. 
    # The running flat is applied so you have [-n frames, current frame, +n frames]. So if the running flat has n=12, then the 
    # so called "left flat" is built on frame 13. So the lowest value that we find in the flatnames (the one that we match in the previous for loop)
    # is the next to the left flat. i.e. 14.
    echo $(( min - 2))
}
export -f getRunningFlatSizeUsed

mergeArrays() {
    local -n arr1=$1
    local -n arr2=$2
    local -n out=$3

    local i=0 j=0
    out=()

    while (( i < ${#arr1[@]} && j < ${#arr2[@]} )); do
        key1="${arr1[i+1]}"
        key2="${arr2[j+1]}"

        if [[ "$key1" < "$key2" ]]; then
            out+=("${arr1[i]}" "$key1")
            ((i+=2))
        else
            out+=("${arr2[j]}" "$key2")
            ((j+=2))
        fi
    done

    while (( i < ${#arr1[@]} )); do
        out+=("${arr1[i]}" "${arr1[i+1]}")
        ((i+=2))
    done

    while (( j < ${#arr2[@]} )); do
        out+=("${arr2[j]}" "${arr2[j+1]}")
        ((j+=2))
    done
}
export -f mergeArrays

correctObjectsAndFlatsMissing() {
    local -n objAndFlatArray=$1
    local -n newArray=$2

    flag=""
    for ((i=0; i<${#objectsAndFlatsCombined[@]}; i+=2)); do
        currentFile="${objectsAndFlatsCombined[i]}"
        currentKey="${objectsAndFlatsCombined[i+1]}"
        currentFlag=${objectsAndFlatsCombined[i]:0:4}

        if [[ "$flag" == "$currentFlag" ]]; then
            if [[ "$flag" == "flat" ]]; then
                unset 'newArray[-1]'  
                unset 'newArray[-1]'
            else
                if (( i+3 < ${#objectsAndFlatsCombined[@]} )); then
                    nextFile="${objectsAndFlatsCombined[i+2]}"
                    nextKey="${objectsAndFlatsCombined[i+3]}"
                    nextkeyMod=$(( nextKey - 140 )) 
                    newArray+=("$nextFile" "$nextkeyMod")
                fi
            fi
        fi

        newArray+=("$currentFile" "$currentKey")
        flag=$currentFlag
    done
}
export -f correctObjectsAndFlatsMissing

divideFlatAndObjects() {
    local -n flatArray=$1
    local -n objectArray=$2

    for ((i=0; i<${#correctedArray[@]}; i+=2)); do
        file="${correctedArray[i]}"
        key="${correctedArray[i+1]}"

        if [[ "$file" == flat* ]]; then
            flatArray+=("$file" "$key")
        else
            objectArray+=("$file" "$key")
        fi
    done
}
export -f divideFlatAndObjects 

source /opt/SIE/local/glob/.bashrc_SIE
source /scratch1/sguerra/venv/bin/activate
module load gnuastro/0.22

objectImagesDir=$1
flatImagesDir=$2
flatsDir=$3
newFlatDir=$4

dateHeaderKey="DATE-OBS"

# Get relevant information for object raws
rawObjectInfo=()
getNameAndDateObs $objectImagesDir $dateHeaderKey 0 rawObjectInfo
orderFramesByDate rawObjectInfo

# Get relevant information from flat field raws
rawFlatInfo=()
getNameAndDateObs $flatImagesDir $dateHeaderKey 0 rawFlatInfo
orderFramesByDate rawFlatInfo

# Get relevant information for flats
flatInfo=()
getNameAndDateObs $flatsDir $dateHeaderKey 1 flatInfo
orderFramesByDate flatInfo

# 1.- Plot of the on off data acquisition
rawObjectInfoPythonArgument=$(printf "%s " "${rawObjectInfo[@]}")
flatFieldInfoPythonArgument=$(printf "%s " "${rawFlatInfo[@]}")
flatInfoPythonArgument=$(printf "%s " "${flatInfo[@]}")

python3 onOffAcquisitionPlot.py "$rawObjectInfoPythonArgument" "$flatFieldInfoPythonArgument" "$flatInfoPythonArgument" "./onOffAcquisitionPlot.png"


# 2.- In order to solve mapping problems (if these exists), we need to work with the flats and the object frames
# For working properly we remove the initial and last N frames, because these will always have the left and right flat
runningFlatSize=$( getRunningFlatSizeUsed flatInfo )

total_len=${#rawObjectInfo[@]}
start=$(( 2 * runningFlatSize )) # 2 * because is an array of pairs
end=$(( total_len - 4 * runningFlatSize ))
rawObjectInfoTrimmed=("${rawObjectInfo[@]:start:end}")


# Now we combine the flats and the object frames in a single data structure to later apply the operations
objectsAndFlatsCombined=()
mergeArrays rawObjectInfoTrimmed flatInfo objectsAndFlatsCombined

#  Correct flats/objects missing 
correctedArray=()
correctObjectsAndFlatsMissing objectsAndFlatsCombined correctedArray

# Recover the flat and object arrays
flatCorrectedArray=()
objectCorrectedArray=()
divideFlatAndObjects flatCorrectedArray objectCorrectedArray


# Corrected plot
rawObjectInfoPythonArgument=$(printf "%s " "${objectCorrectedArray[@]}")
flatFieldInfoPythonArgument=$(printf "%s " "${rawFlatInfo[@]}")
flatInfoPythonArgument=$(printf "%s " "${flatCorrectedArray[@]}")

python3 onOffAcquisitionPlot.py "$rawObjectInfoPythonArgument" "$flatFieldInfoPythonArgument" "$flatInfoPythonArgument" "./onOffAcquisitionPlotCorrected.png"


# Copy the flats in the format that will be correct for hardcoding in he pipeline
if ! [ -d $newFlatDir ]; then mkdir $newFlatDir; fi
counter=$(( runningFlatSize + 2 ))
lastIndex=$(( ${#flatCorrectedArray[@]} - 2 ))  # because we're stepping in pairs
for ((i=0; i<${#flatCorrectedArray[@]}; i+=2)); do
    if (( i == 0 || i == lastIndex)); then
        cp $flatsDir/${flatCorrectedArray[i]}  $newFlatDir/${flatCorrectedArray[i]}
    else
        currentFileName=${flatCorrectedArray[i]}
        newFileName=$(echo "$currentFileName" | sed -E "s/_f[0-9]+/_f$counter/")
        cp $flatsDir/$currentFileName "$newFlatDir/$newFileName"
        counter=$(( counter + 1 ))
    fi
done