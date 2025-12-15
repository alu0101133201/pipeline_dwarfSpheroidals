import sys,glob
import numpy as np
noiseFilesPath=sys.argv[1]
outputFile=sys.argv[2]
rootDir=sys.argv[3]

noiseFiles = glob.glob(noiseFilesPath + "/entirecamera_*.txt")
rms=[]
for file in noiseFiles:
    try:
        with open(file, 'r') as currentFile:
            lines=currentFile.readlines()
            for line in lines:
                data = line.split(' ')
                rms.append(float(data[2]))

    except Exception as e:
        print(f'Error processing file {file}: {e}')

rms_min = np.nanmin(rms)
rms=np.array(rms)
weights=(rms_min**2)/(rms**2)
sumOfWeights=np.nansum(weights)
file = open(rootDir + '/build/' + outputFile, "w")
file.write(str(rms_min))