import os
import tarfile
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# Define the directory where the 'proteomes' folder is located
proteomes_directory = 'proteomes/'

# Function to check if a file ends with '.tar' and exists within the 'proteomes/' directory
def is_tar_file(file_name):
    full_path = os.path.join(proteomes_directory, file_name)
    if file_name.endswith('.tar') and os.path.isfile(full_path):
        return full_path
    return None

# Use ThreadPoolExecutor to parallelize the filtering of tar files
with ThreadPoolExecutor(max_workers=128) as executor:
    tar_files = list(filter(None, tqdm(executor.map(is_tar_file, os.listdir(proteomes_directory)), total=len(os.listdir(proteomes_directory)))))

# Save the tar file names (with full paths) to a file for later use
with open('tar_files_list.txt', 'w') as file:
    for tar_file in tar_files:
        file.write(f"{tar_file}\n")  # Full paths are already stored

# # Now, to load the tar files from the saved file
# with open('tar_files_list.txt', 'r') as file:
#     tar_files = [line.strip() for line in file.readlines()]

# Function to count 'confidence' files within a tar file
def count_confidence_files(tar_file):
    if not os.path.exists(tar_file):
        print(f"File not found: {tar_file}")  # Logging missing files
        return 0
    try:
        with tarfile.open(tar_file, 'r') as tar:
            return sum(1 for member in tar.getmembers() if 'confidence' in member.name)
    except Exception as e:
        print(f"Error processing file {tar_file}: {e}")
        return 0

# Initialize the total confidence files counter
total_confidence_files = 0

# Use ThreadPoolExecutor to parallelize the counting task
with ThreadPoolExecutor(max_workers=128) as executor:  # Using threads for I/O-bound operations
    results = list(tqdm(executor.map(count_confidence_files, tar_files), total=len(tar_files)))
    total_confidence_files = sum(results)

# Print the total number of tar files and the total number of 'confidence' files
print(f"Total number of tar files available: {len(tar_files)}")
print(f"Total number of 'confidence' files across all tar files: {total_confidence_files}")
