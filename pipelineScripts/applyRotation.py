import sys
import numpy as np

rotation_angle=int(sys.argv[1])
x_ccd=float(sys.argv[2])
y_ccd=float(sys.argv[3])
xaxis_size=float(sys.argv[4])


if rotation_angle==90:
    cos_th=0; sin_th=1
elif rotation_angle==-90:
    cos_th=0; sin_th=-1
else:
    rotation_angle_rad=np.radians(rotation_angle)
    cos_th=np.cos(rotation_angle_rad)
    sin_th=np.sin(rotation_angle_rad)

x_new=x_ccd*cos_th-y_ccd*sin_th
y_new=x_ccd*sin_th+y_ccd*cos_th +xaxis_size#We translate the refference point to the other corner of the CCD x-axis 


print(f"x_new={x_new}")
print(f"y_new={y_new}")