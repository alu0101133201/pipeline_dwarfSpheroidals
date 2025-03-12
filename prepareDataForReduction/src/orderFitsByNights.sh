# Sorts the fits file based on the night 
# The date/time information is obtained from the keyword specified in $2

directoryWithFits=$1
keywordWithDate=$2

for i in $directoryWithFits/*.fits; do
    fileName=$( basename $i )

    date=$( gethead $i $keywordWithDate )
    date=$( echo "$date" | cut -d'.' -f1) # We trim the miliseconds
    night=$( date -d "$date" +"%Y-%m-%d")
    time=$( date -d "$date" +"%H:%M:%S")
    hour=$( date -d "$date" +"%H")


    # If the hour is before 12 then we assume it is in the night after midnight so technically
    # is associated with the previous day
    if (( 10#$hour < 12 )); then # The 10# is for decimal interpretation
        previousNight=$(date -d "$night - 1 day" +"%Y-%m-%d")
        mkdir -p $directoryWithFits/$previousNight
        mv $i $directoryWithFits/$previousNight
    else
        mkdir -p $directoryWithFits/$night
        mv $i $directoryWithFits/$night
    fi
done