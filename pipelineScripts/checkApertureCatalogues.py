import os 
import pandas as pd
import sys
from shutil import copyfile

"""
With this python script, we will check the aperture photometry catalogues
generated in prepareCalibrationData for 2 filters (typically g and r), and 
if one of them does not exist in one filter, we remove both the catalogue and the 
frame from both filters.
"""
catalogues_1=sys.argv[1]
catalogues_2=sys.argv[2]
images_1=sys.argv[3]
images_2=sys.argv[4]
filter1=sys.argv[5]
filter2=sys.argv[6]

def get_basename(file,suffix):
    """Remove filter and extension, e.g., 't021.1118+7.1450.g.cat' → 't021.1118+7.1450'"""
    return file.replace(f".{suffix}.cat", "").replace(f".{suffix}.fits", "")

cats_1 = {get_basename(f, filter1) for f in os.listdir(catalogues_1) if f.endswith(f".{filter1}.cat")}
cats_2 = {get_basename(f, filter2) for f in os.listdir(catalogues_1) if f.endswith(f".{filter2}.cat")}

common = cats_1 & cats_2
print(common)
sys.exit()
missing_1 = cats_1 - common
missing_2 = cats_2 - common

print(f"Catalogues removed from {filter1}: {len(missing_1)}")
print(f"Catalogues removed from {filter2}: {len(missing_2)}")

def delete_files(folder,missing,filter_band,ext):
    for name in missing:
        path = os.path.join(folder,f"{name}.{filter_band}{ext}")
        if os.path.exists(path):
            os.remove(path)
            print(f"Deleted {path}")

delete_files(catalogues_1,missing_1 | missing_2, filter1, ".cat")
delete_files(catalogues_2,missing_1 | missing_2, filter2, ".cat")
delete_files(catalogues_1,missing_1 | missing_2, filter1, ".fits")
delete_files(catalogues_2,missing_1 | missing_2, filter2, ".fits")

def update_brick_file(brick_file,filter_band,common_names):
    if not os.path.exists(brick_file):
        print(f"⚠️ No brick file found: {brick_file}")
        return

    # Backup original
    backup = brick_file + ".bak"
    copyfile(brick_file, backup)
    print(f"Backup created: {backup}")

    df = pd.read_csv(brick_file, delim_whitespace=True)
    df["BaseName"] = df["BrickName"].apply(lambda x: x.replace(f".{filter_band}", ""))
    df_filtered = df[df["BaseName"].isin(common_names)]

    df_filtered.drop(columns="BaseName").to_csv(brick_file, sep="\t", index=False)
    print(f"Updated {brick_file}: {len(df_filtered)} valid entries")

update_brick_file(os.path.join(images_1, "brickIdentification.txt"), filter1, common)
update_brick_file(os.path.join(images_2, "brickIdentification.txt"), filter2, common)