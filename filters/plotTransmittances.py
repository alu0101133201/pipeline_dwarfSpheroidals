import numpy as np
import matplotlib.pyplot as plt

def read_transmittance_data(file_path):
    try:
        data = np.loadtxt(file_path, skiprows=1)
        wavelengths = data[:, 0]
        transmittance = data[:, 1]
        return wavelengths, transmittance
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None, None

def plot_transmittance_data(data_dict):
    if not data_dict:
        print("The provided dictionary is empty. Nothing to plot.")
        return

    plt.figure(figsize=(10, 6))
    for file_name, (wavelengths, transmittance) in data_dict.items():
        plt.plot(wavelengths, transmittance, label=file_name)

    plt.title('Transmittance vs. Wavelength for Multiple Samples')
    plt.xlabel('Wavelength ($nm$)')
    plt.ylabel('Transmittance (T)')
    plt.legend()
    plt.grid(True)

    plt.savefig(f"./{instrument}_filters.png")

# file_names = ["DECaLS_g.dat", "DECaLS_i.dat", "DECaLS_r.dat"]; instrument = "DECaLS"
# file_names = ["HIPERCAM_g.dat", "HIPERCAM_i.dat", "HIPERCAM_r.dat", "HIPERCAM_z.dat"]; instrument = "HIPERCAM"
# file_names = ["OSIRIS+_g.dat", "OSIRIS+_i.dat", "OSIRIS+_r.dat", "OSIRIS+_z.dat", "OSIRIS+_u.dat"]; instrument = "OSIRIS+"
file_names = ["PANSTARRS_g.dat", "PANSTARRS_r.dat", "PANSTARRS_i.dat", "PANSTARRS_z.dat", "PANSTARRS_y.dat"]; instrument = "PANSTARRS"
# file_names = ["SDSS_g.dat", "SDSS_r.dat", "SDSS_u.dat"]; instrument = "SDSS"
# file_names = ["TST_g.dat", "TST_r.dat", "TST_i.dat", "TST_lum.dat", "TST_Ha.dat"]; instrument = "TST"
# file_names = ["TTT3_QHY_g.dat", "TTT3_QHY_r.dat", "TTT3_QHY_i.dat", "TTT3_QHY_lum.dat", "TTT3_QHY_Ha.dat"]; instrument = "TTT3_QHY"
# file_names = ["TTT3_iKon_u.dat", "TTT3_iKon_g.dat", "TTT3_iKon_r.dat", "TTT3_iKon_i.dat", "TTT3_iKon_lum.dat", "TTT3_iKon_Ha.dat", "TTT3_iKon_y.dat", "TTT3_iKon_zs.dat"]; instrument = "TTT3_iKon"
# file_names = ["WHT_TWFC_g.dat", "WHT_TWFC_r.dat", "WHT_TWFC_i.dat"]; instrument = "WHT"

all_data = {}

for file_name in file_names:
    wavelengths, transmittance = read_transmittance_data(file_name)
    if wavelengths is not None and transmittance is not None:
        all_data[file_name] = (wavelengths, transmittance)

plot_transmittance_data(all_data)