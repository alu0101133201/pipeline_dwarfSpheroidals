

directoryWithFits=$1


for i in $directoryWithFits/*.fits; do
    filter=$( basename $i | awk -F'_' '{print $NF}' | awk -F'.' '{print $1}')

    if ! [ -d $directoryWithFits/$filter ]; then mkdir $directoryWithFits/$filter; fi
    mv $i $directoryWithFits/$filter
done