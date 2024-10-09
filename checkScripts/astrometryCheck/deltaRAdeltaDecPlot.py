import glob

import numpy as np
import matplotlib.pyplot as plt


def read_columns_from_file(file_path):
    data = np.loadtxt(file_path, comments='#')
    ra1  = np.array(data[:, 0].astype(float).tolist())
    dec1 = np.array(data[:, 1].astype(float).tolist())
    ra2  = np.array(data[:, 2].astype(float).tolist())
    dec2 = np.array(data[:, 3].astype(float).tolist())
    return ra1, dec1, ra2, dec2

raArrays = []
decArrays = []

for i in glob.glob("*.cat"):
    ra1, dec1, ra2, dec2 = read_columns_from_file(i)

    raArrays.append(ra1-ra2)
    decArrays.append(dec1-dec2)

for i in range(len(raArrays)):
    raArrays[i] = raArrays[i]*3600
    decArrays[i] = decArrays[i]*3600


fig, ax = plt.subplots(1, 1, figsize=(15, 15))

for i in range(len(raArrays)):
    ax.scatter(raArrays[i], decArrays[i], color="blue", s=50, linewidths=1.5, edgecolor="black")
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.hlines(y=0, xmin=-1, xmax=1, color="black", linestyle="--")
ax.vlines(x=0, ymin=-1, ymax=1, color="black", linestyle="--")
plt.savefig("deltaRadeltaDec.jpg")
