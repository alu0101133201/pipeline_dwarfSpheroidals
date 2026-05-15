#!/bin/bash

downloadIndex() {
    local re=$1
    local catName=$2
    local indexdir=$3

    build-astrometry-index -i $catName -e1 \
                            -P $re \
                            -S phot_g_mean_mag \
                            -E -A RA -D  DEC\
                            -o $indexdir/index_$re.fits;
}
export -f downloadIndex
