import os,sys,glob
import numpy as np
import math
##This function will write a txt file with crops section taking into account how many frames and of which size we can combine.
##Arguments:

framesDir=sys.argv[1]
sizeX,sizeY=map(int,sys.argv[2].split(","))
availMemory=float(sys.argv[3]) #In GB
outputFile=sys.argv[4]
numberOfBlocksFile=sys.argv[5]

# Optional: minimum total number of blocks (numberOfBlocksX * numberOfBlocksY)
# to produce, e.g. the number of nodes available for the coaddition step, so
# a memory-driven grid smaller than the node count doesn't leave nodes idle.
# Defaults to 1 (no minimum) so existing callers are unaffected.
minNumberOfBlocks=int(sys.argv[6]) if len(sys.argv) > 6 else 1

#Number of frames to combine
totalFrames=len(glob.glob(os.path.join(framesDir,"*.fits")))


totalNumberOfPixels=sizeX*sizeY
sizeBytes=4.*totalNumberOfPixels*totalFrames
sizeGB=sizeBytes/1e9
#The idea is how much memory is needed to combine all frames at one pixel
#In other words, with $sizeGB we know how much memory is needed to combine a frame of sizeX x sizeY
#We want to know the minimum value of N such that we can combine N frames at once with availMemory
GbPerPx=sizeGB/totalNumberOfPixels #This is: GB needed to combine a block of one pixel
maxPixelsPerCombination=availMemory/GbPerPx
ratio=sizeX/sizeY
sizeBlockX_pix=max(1,int(np.sqrt(maxPixelsPerCombination*ratio)))
sizeBlockY_pix=max(1,int(np.sqrt(maxPixelsPerCombination/ratio)))
numberOfBlocksX=math.ceil(sizeX/sizeBlockX_pix)
numberOfBlocksY=math.ceil(sizeY/sizeBlockY_pix)

memoryDrivenBlocks=numberOfBlocksX*numberOfBlocksY
if memoryDrivenBlocks < minNumberOfBlocks:
        scale=np.sqrt(minNumberOfBlocks/memoryDrivenBlocks)
        numberOfBlocksX=math.ceil(numberOfBlocksX*scale)
        numberOfBlocksY=math.ceil(numberOfBlocksY*scale)
        sizeBlockX_pix=math.ceil(sizeX/numberOfBlocksX)
        sizeBlockY_pix=math.ceil(sizeY/numberOfBlocksY)

sections=[]
for i in range(numberOfBlocksX):
        for j in range(numberOfBlocksY):
                i1=i*sizeBlockX_pix+1
                i2=min((i+1)*sizeBlockX_pix,sizeX)
                j1=j*sizeBlockY_pix+1
                j2=min((j+1)*sizeBlockY_pix,sizeY)
                sections.append((i1,i2,j1,j2,i+1,j+1))

with open(outputFile,"w") as f:
        for (i1,i2,j1,j2,i,j) in sections:
                f.write(f"{i1}:{i2},{j1}:{j2} {i} {j}\n")


with open(numberOfBlocksFile,"w") as f2:
        f2.write(f"{numberOfBlocksX} {numberOfBlocksY}")

