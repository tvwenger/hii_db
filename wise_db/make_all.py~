import os
import subprocess
import glob
import tarfile

# Define the top-level directory where the scripts are located
root_directory = '/Users/loren/catalogs'

# Use glob to search for all Python files recursively
python_files = glob.glob(f'{root_directory}/**/*.py', recursive=True)

# Iterate through each Python file and execute it
for script_path in python_files:
    # Get the directory of the script
    script_dir = os.path.dirname(script_path)
    
    # Change the current working directory to the script's directory
    os.chdir(script_dir)
    
    print(f"Executing {script_path} in directory {script_dir}...")

    # Run the script using subprocess
    result = subprocess.run(['python', script_path], capture_output=True, text=True)

    # Output the result of the script (stdout and stderr)
    print("Output:")
    print(result.stdout)
    print("Errors:")
    print(result.stderr)

# make a tar file
os.chdir('/Users/loren/papers/wise/python/rrl_surveys/')
tarball_name = 'rrl_surveys.tar.gz'

# Open the tarball in write mode with gzip compression
with tarfile.open(tarball_name, 'w:gz') as tar:
    # Add the directory to the tarball
    tar.add('.pkl', arcname=os.path.basename('.'))

print(f"Tarball {tarball_name} created successfully!")

