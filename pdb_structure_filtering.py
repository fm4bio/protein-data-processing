import os
import gzip
from Bio import PDB
from tqdm import tqdm
import datetime
from multiprocessing import Pool, cpu_count

def parse_date(date_str):
    # Attempt to parse the date using different formats
    for date_format in ('%d-%b-%y', '%Y-%m-%d'):
        try:
            return datetime.datetime.strptime(date_str, date_format)
        except ValueError:
            continue
    # If no formats match, raise an error
    raise ValueError(f"Date {date_str} does not match any known format")

def parse_pdb_file_metadata(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    with gzip.open(pdb_file, 'rt') as f:
        structure = parser.get_structure('struct', f)
    
    header = structure.header
    deposition_date = parse_date(header['deposition_date'])
    resolution = header.get('resolution', None)
    experiment_method = header['structure_method']
    
    return deposition_date, resolution, experiment_method

def filter_pdb_file(args):
    pdb_file, filter_date, resolution_threshold, require_xray = args
    deposition_date, resolution, experiment_method = parse_pdb_file_metadata(pdb_file)

    if deposition_date < filter_date and (not require_xray or experiment_method.lower() == 'x-ray diffraction') and resolution and resolution < resolution_threshold:
        return pdb_file
    return None

def filter_pdb_directory(directory, filter_date, resolution_threshold=9.0, require_xray=True, num_workers=cpu_count()):
    pdb_files = []
    filtered_pdb_files = []

    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.ent.gz'):
                pdb_files.append(os.path.join(root, file))

    with Pool(processes=num_workers) as pool:
        args = [(pdb_file, filter_date, resolution_threshold, require_xray) for pdb_file in pdb_files]
        for result in tqdm(pool.imap_unordered(filter_pdb_file, args), total=len(pdb_files)):
            if result:
                filtered_pdb_files.append(result)

    return filtered_pdb_files

# Example usage
filter_date = datetime.datetime(2020, 5, 1)
resolution_threshold = 9.0
require_xray = True  # Set to False if X-ray is not required
pdb_directory = 'pdb/structures'

# Filter the PDB structures
filtered_pdb_files = filter_pdb_directory(pdb_directory, filter_date, resolution_threshold, require_xray, num_workers=128)

# Print the number of filtered structures
print("Total structures matching criteria:", len(filtered_pdb_files))
