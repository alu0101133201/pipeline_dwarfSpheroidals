import os
from astropy.io import fits
import sys

mapping_file=   sys.argv[1]
folder=         sys.argv[2]

# === READ THE MAPPING FILE ===
mapping={}
with open(mapping_file,'r') as f:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue
        parts=line.strip().split()
        if len(parts) >= 2:
            new_name,original_file=parts[0],parts[1]
            mapping[original_file]=new_name

# === STEP 1: TEMPORARY RENAME ===
temp_mapping={}
for filename in os.listdir(folder):
    if not filename.endswith('.fits'):
        continue
    filepath=os.path.join(folder,filename)
    temp_filename=f"tmp_{filename}"
    temp_filepath=os.path.join(folder,temp_filename)
    os.rename(filepath,temp_filepath)
    temp_mapping[temp_filename]=filename

# === STEP 2: RENAME TO FINAL NAMES ===
for temp_filename in list(temp_mapping.keys()):
    temp_filepath=os.path.join(folder,temp_filename)

    try:
        with fits.open(temp_filepath) as hdul:
            header=hdul[1].header
            original_value=header.get('ORIGINAL_FILE', None)
        if original_value and original_value in mapping:
            new_filename=mapping[original_value]
            new_filepath=os.path.join(folder,new_filename)
            
            if os.path.exists(new_filepath):
                print(f"Warning: {new_filepath} already exists. Skipping rename of {temp_filepath}.")
                continue

            os.rename(temp_filepath,new_filepath)
            print(f"Renamed {temp_filepath} to {new_filepath}")
        else:
            print(f"Warning: ORIGINAL_FILE not found or not in mapping for {temp_filepath}. Skipping.")
    except Exception as e:
        print(f"Error processing {temp_filepath}: {e}")
        continue