Important note:

I was extremely naive when I did this code, so the approach is not adequate. It should be reimplemented.
I expected to have a decent coverage of the night (on and off field) with some gaps, so I simply kept the association
but filling the potential gaps and doing small fixes. The reality is that some nights have large holes, or one of the fields stars/ends with 4 frames and the the other field stars (which has implications for how to assing the left/right flat)...

Thus, it moderately works (if you have you flat shifted is not dramatic either), but is not a good implementation. This is currently not
used frequently. But if you're going to use this a lot... simply reimplement it.

The robust approach would be to do the assignation based on the DATE-OBS, and the handle the running flat (left right flats).


The scripts contained in this folder are for reducing observations done with the on-off strategy

1.- Reduce de flat field. Thus you obtain the flats to use in your object field
2.- Run "prepareMappingBetweenObjectandFlat". This will fix errors that usually happen mainly because
    some frame was skipped. we need to fix these holes in order to get an accurate map and not to lose frames
3.- Run "hardCodeFlats" to place the flats to use in the reduction in your build folder
4.- Run the pipeline


CAVEAT: This script fixes basic problems that I have found. But the assignation can be tricky depending on the specific set
of object and flats images available (more objects than flats, the other way around, gaps in the night, etc...) so be careful
and check always that everything's fine
