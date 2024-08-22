import os
import tarfile
import json
import gzip
from statistics import mean
from Bio.PDB import MMCIFParser
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# Define the directory where the 'proteomes' folder is located
proteomes_directory = 'proteomes/'

# Threshold for pLDDT
plddt_threshold = 0.7
# Distance threshold for long-range contacts (in Å)
distance_threshold = 8.0
# Minimum sequence separation for long-range contacts
min_seq_distance = 12

# Function to calculate long-range contacts for a structure
def calculate_long_range_contacts(structure):
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    atoms.append((residue['CA'].get_coord(), residue.get_id()[1]))
                    
    long_range_contacts = 0
    num_atoms = len(atoms)
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            if abs(atoms[i][1] - atoms[j][1]) > min_seq_distance:
                distance = np.linalg.norm(atoms[i][0] - atoms[j][0])
                if distance < distance_threshold:
                    long_range_contacts += 1
    
    return long_range_contacts, num_atoms

# Function to process a single tar file for pLDDT and globularity filtering
def process_tar_file(tar_file_name):
    tar_file_path = os.path.join(proteomes_directory, tar_file_name)
    proteins_meeting_criteria = 0
    total_proteins = 0
    
    try:
        with tarfile.open(tar_file_path, 'r') as tar:
            # Create a dictionary to store the model and confidence files by protein name
            protein_files = {}

            # Iterate over the files in the tar archive
            for member in tar.getmembers():
                # Extract the base protein identifier from the filename
                if 'model' in member.name and member.name.endswith('.cif.gz'):
                    protein_id = member.name.split('-model')[0]
                    if protein_id not in protein_files:
                        protein_files[protein_id] = {}
                    protein_files[protein_id]['model'] = member.name
                elif 'confidence' in member.name and member.name.endswith('.json.gz'):
                    protein_id = member.name.split('-confidence')[0]
                    if protein_id not in protein_files:
                        protein_files[protein_id] = {}
                    protein_files[protein_id]['confidence'] = member.name

            # Process each protein identified in the tar file
            for protein_id, files in protein_files.items():
                # print(f"Processing protein: {protein_id}")
                total_proteins += 1
                # Ensure both model and confidence files exist for this protein
                if 'model' in files and 'confidence' in files:
                    # Process the confidence JSON file
                    confidence_json_file = tar.extractfile(files['confidence'])
                    if confidence_json_file:
                        with gzip.open(confidence_json_file, 'rt') as f:
                            confidence_data = json.load(f)
                            plddt_scores = confidence_data.get('confidenceScore', [])
                            if plddt_scores:
                                average_plddt = mean(plddt_scores)
                                if average_plddt < plddt_threshold:
                                    # print(f"{protein_id}: Fails pLDDT filtering (average pLDDT: {average_plddt})")
                                    continue  # Skip this protein, fails pLDDT filtering

                    # Process the CIF model file
                    model_cif_file = tar.extractfile(files['model'])
                    if model_cif_file:
                        with gzip.open(model_cif_file, 'rt') as f:
                            parser = MMCIFParser()
                            structure = parser.get_structure(protein_id, f)
                            long_range_contacts, L = calculate_long_range_contacts(structure)
                            if long_range_contacts < 0.5 * L:
                                # print(f"{protein_id}: Fails globularity filtering (long-range contacts: {long_range_contacts}, length: {L})")
                                continue  # Skip this protein, fails globularity filtering

                    # print(f"{protein_id}: Passes both filters")
                    proteins_meeting_criteria += 1  # Passes both filters

                else:
                    print(f"{protein_id}: Missing model or confidence file, skipping")

        return proteins_meeting_criteria, total_proteins
        
    except Exception as e:
        print(f"Error processing {tar_file_name}: {e}")
    return proteins_meeting_criteria, total_proteins

# Get a list of all tar files
tar_files = [file for file in os.listdir(proteomes_directory) if file.endswith('.tar')]

# Print the total number of files
print(f"Total number of tar files available: {len(tar_files)}")

# Initialize counters for proteins
total_proteins_processed = 0
total_proteins_meeting_criteria = 0

# Process files in parallel with tqdm for progress tracking
with ProcessPoolExecutor(max_workers=128) as executor:
    for proteins_meeting_criteria, total_proteins in tqdm(executor.map(process_tar_file, tar_files), total=len(tar_files)):
        total_proteins_processed += total_proteins
        total_proteins_meeting_criteria += proteins_meeting_criteria

# Display the summary results
print(f"Total proteins processed: {total_proteins_processed}")
print(f"Proteins meeting criteria: {total_proteins_meeting_criteria}")
print(f"Ratio of proteins meeting criteria: {total_proteins_meeting_criteria / total_proteins_processed:.2f}")
