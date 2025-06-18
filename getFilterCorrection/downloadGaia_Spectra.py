import os
import glob

from astroquery.gaia import Gaia

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# fieldName   = "ic1613"
# ra          = 16.2086
# dec         = 2.1243

# fieldName   = "ugc5364"
# ra          = 149.8528
# dec         = 30.7493

fieldName   = "Coma"
ra          = 195.20333
dec         = 27.66072



sizeOfField = 1.5
halfSizeOfFild = sizeOfField / 2
raMin  = ra - halfSizeOfFild
raMax  = ra + halfSizeOfFild
decMin = dec - halfSizeOfFild
decMax = dec + halfSizeOfFild

query = f"SELECT * \
            FROM gaiadr3.gaia_source WHERE (ra BETWEEN {raMin} AND {raMax}) AND (dec BETWEEN {decMin} AND {decMax}) AND \
            has_xp_sampled = 'True'"
            
job     = Gaia.launch_job_async(query)
results = job.get_results()

chunk_size   = 500
ids          = results['source_id']
ids_chunks   = list(chunks(ids, chunk_size))

print(f'* Input list contains {len(ids)} source_IDs')
print(f'* This list is split into {len(ids_chunks)} chunks of <= {chunk_size} elements each')   


retrieval_type = 'XP_SAMPLED'
data_structure = 'INDIVIDUAL' 
data_release   = 'Gaia DR3'   
datalink_all   = []

ii = 0
for chunk in ids_chunks:
   ii = ii + 1
   print(f'Downloading Chunk #{ii}; N_files = {len(chunk)}')
   datalink  = Gaia.load_data(ids=chunk, data_release = data_release, retrieval_type=retrieval_type, format = 'fits', data_structure = data_structure, dump_to_file=True)
   datalink_all.append(datalink) 

datalink_out = datalink_all[0]
for inp_dict in datalink_all[1:]:
   datalink_out.update(inp_dict)

keys = list(datalink_out.keys())
print(f'* The merged dictionary contains {len(keys)} elements')  
if not os.path.exists(f'./gaiaSpectra_{fieldName}'):
    os.system(f"mkdir gaiaSpectra_{fieldName}")

os.system(f"mv ./datalink_output*.zip ./gaiaSpectra_{fieldName}")
for file in glob.glob(f"./gaiaSpectra_{fieldName}/*.zip"):
    os.system(f"unzip {file} -d ./gaiaSpectra_{fieldName}")
    os.system(f"rm {file}")
