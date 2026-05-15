downloadCatalogue() {
    local surveyToUse=$1
    local ra=$2
    local dec=$3
    local radius=$4
    local catDir=$5
    local catName=$6

    case "$surveyToUse" in
        gaia)
            query_param="gaia --dataset=dr3 --center=$ra,$dec --radius=$radius --column=ra,dec,ra_error,dec_error,phot_g_mean_mag,phot_g_mean_flux,phot_g_mean_flux_error,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error"
            downloadGaiaCatalogue "$query_param" "$catDir" "$catName"
            ;;
        panstarrs)
            query_param="vizier --dataset=panstarrs1 --center=$ra,$dec --radius=$radius --column=RAJ2000,DEJ2000,gmag"
            downloadPanstarrsCatalogue "$query_param" "$catDir" "$catName"
            ;;
        *)
            echo "Unknown catalog source: $surveyToUse"
            exit 222
            ;;
    esac

}
export -f downloadCatalogue


downloadGaiaCatalogue() {
    local query=$1
    local catdir=$2
    local catName=$3

    astquery $query -o $catdir/"$objectName"_Gaia_DR3_tmp.fits
    asttable $catdir/"$objectName"_Gaia_DR3_tmp.fits -c1,2,3 -c'arith $4 abs' -c'arith $5 3 x' -c'arith $6 abs' -c'arith $7 3 x' -c'arith $8 abs' -c'arith $9 3 x' --noblank=4 -o$catdir/tmp.txt
    # I have explored 3 different ways of selecting good stars.
    # From the most restrictive to the less restrictive:

    # # Here I demand that the gaia object fulfills simultaneously that:
    # # 1.- Parallax > 3 times its error
    # # 2.- Proper motion (ra) > 3 times its error
    # # 3.- Proper motion (dec) > 3 times its error
    # asttable $catdir/tmp.txt -c1,2,3 -c'arith $4 $4 $5 gt 1000 where' -c'arith $6 $6 $7 gt 1000 where' -c'arith $8 $8 $9 gt 1000 where' -o$catdir/test_.txt
    # asttable $catdir/test_.txt -c1,2,3 -c'arith $4 $5 + $6 +' -o$catdir/test1.txt
    # asttable $catdir/test1.txt -c1,2,3 --range=ARITH_2,2999,3001 -o $catName

    # # Here I only demand that the parallax is > 3 times its error
    # asttable $catdir/tmp.txt -c1,2,3 -c'arith $4 $4 $5 gt 1000 where' -o$catdir/test_.txt
    # asttable $catdir/test_.txt -c1,2,3 --range=ARITH_2,999,1001 -o $catName

    # Here I  demand that the parallax OR a proper motion is > 3 times its error
    asttable $catdir/tmp.txt -c1,2,3 -c'arith $4 $4 $5 gt 1000 where' -c'arith $6 $6 $7 gt 1000 where' -c'arith $8 $8 $9 gt 1000 where' -o$catdir/test_.txt
    asttable $catdir/test_.txt -c1,2,3 -c'arith $4 $5 + $6 +' -o$catdir/test1.txt
    asttable $catdir/test1.txt -c1,2,3 --range=ARITH_2,999,3001 -o $catName

    # # Here we don't demand any condition
    # asttable $catdir/tmp.txt -o $catName

    rm $catdir/test1.txt $catdir/tmp.txt $catdir/test_.txt

    # This extra code is for using the catalogue for scamp (needed for LaPalma supercomputer)
    mv $catdir/"$objectName"_Gaia_DR3_tmp.fits fullCatalogueForScamp.fits
}
export -f downloadGaiaCatalogue

downloadPanstarrsCatalogue() {
    local query=$1
    local catdir=$2
    local catName=$3

    astquery $query -o $catdir/"$objectName"_Panstarrs_S1_tmp.fits
    asttable $catdir/"$objectName"_Panstarrs_S1_tmp.fits -c1,2,3  --colmetadata=1,RA,deg --colmetadata=2,DEC,deg --colmetadata=3,phot_g_mean_mag,mag -o$catName
    rm $catdir/"$objectName"_Panstarrs_S1_tmp.fits 
}
export -f downloadPanstarrsCatalogue