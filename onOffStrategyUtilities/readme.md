The scripts contained in this folder are for reducing observations done with the on-off strategy

1.- Reduce de flat field. Thus you obtain the flats to use in your object field
2.- Run "prepareMappingBetweenObjectandFlat". This will fix errors that usually happen mainly because
    some frame was skipped. we need to fix these holes in order to get an accurate map and not to lose frames
3.- Run "hardCodeFlats" to place the flats to use in the reduction in your build folder
4.- Run the pipeline