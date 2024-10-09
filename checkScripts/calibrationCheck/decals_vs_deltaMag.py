import sys
import glob

import numpy as np
import matplotlib.pyplot as plt



def read_columns_from_file(file_path):
    data = np.loadtxt(file_path, comments='#')
    mag1  = np.array(data[:, 0].astype(float).tolist())
    mag2 = np.array(data[:, 1].astype(float).tolist())
    return mag1, mag2



directoryWithTheCatalogues = sys.argv[1]
imageName = sys.argv[2]

magDiff = np.array([])
mag1Total = np.array([])

for i in glob.glob(directoryWithTheCatalogues + "/*.cat"):
    mag1, mag2 = read_columns_from_file(i)
    mag1Total = np.append(mag1Total, mag1)
    magDiff = np.append(magDiff, np.array(mag1 - mag2))


fig, ax = plt.subplots(1, 1, figsize=(15, 15))
ax.set_xlabel("DECaLS mag", fontsize=20)
ax.set_ylabel("Magnitude differences", fontsize=20)
ax.scatter(mag1Total, magDiff, color="blue", s=50, linewidths=1.5, edgecolor="black")

plt.savefig(imageName + ".jpg")
