import os,sys,glob
import numpy as np
import math
##This function will write a txt file with crops section taking into account how many frames and of which size we can combine.
##Arguments:

framesDir=sys.argv[1]
coaddSizePx=int(sys.argv[2])
availMemory=float(sys.argv[3]) #In GB
outputFile=sys.argv[4]
numberOfBlocksFile=sys.argv[5]

#Number of frames to combine
totalFrames=len(glob.glob(os.path.join(framesDir,"entirecamera*.fits")))
#Size of the folder in GB
sizeBytes=sum(os.path.getsize(f) for f in glob.glob(os.path.join(framesDir,"entirecamera*.fits")))
sizeGB=sizeBytes/(1024.0**3)

totalNumberOfPixels=coaddSizePx**2
#The idea is how much memory is needed to combine all frames at one pixel
#In other words, with $sizeGB we know how much memory is needed to combine a frame of coaddSizePx x coaddSizePx
#We want to know the minimum value of N such that we can combine N frames at once with availMemory
GbPerPx=sizeGB/totalNumberOfPixels #This is: GB needed to combine a block of one pixel
maxPixelsPerCombination=availMemory/GbPerPx
sizeBlock_pix=int(np.sqrt(maxPixelsPerCombination))
numberOfBlocks=math.ceil(coaddSizePx/sizeBlock_pix)

sections=[]
for i in range(numberOfBlocks):
    for j in range(numberOfBlocks):
        i1=i*sizeBlock_pix+1
        i2=min((i+1)*sizeBlock_pix,coaddSizePx)
        j1=j*sizeBlock_pix+1
        j2=min((j+1)*sizeBlock_pix,coaddSizePx)
        sections.append((i1,i2,j1,j2,i+1,j+1))

with open(outputFile,"w") as f:
        for (i1,i2,j1,j2,i,j) in sections:
                f.write(f"{i1}:{i2},{j1}:{j2} {i} {j}\n")


with open(numberOfBlocksFile,"w") as f2:
      f2.write(f"{numberOfBlocks}")

