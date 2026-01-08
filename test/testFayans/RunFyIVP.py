import os
import subprocess
import shutil
from pathlib import Path

# Define nuclei with their properties
nuclei = {
    'Sn120': {'Z': 50, 'A': 120},
    'Sn132': {'Z': 50, 'A': 132},
    'Pb208': {'Z': 82, 'A': 208},
    'Gd146': {'Z': 64, 'A': 146}
}

# Paths
current_dir = Path(__file__).resolve().parent
base_dir = current_dir.parent.parent
exe_path = base_dir / 'haforpa.new'
input_dir = base_dir / 'test' / 'testFayans'
output_base_dir = base_dir / 'test' / 'testFayans' / 'testNuclei'
# Input files with their corresponding folder names
input_configs = {
    'for005.fyivpFR': 'FR',  # Finite Range
    'for005.fyivpZR': 'ZR'   # Zero Range
}

def update_input_file(file_path, Z, N):
    """Update Z,N values on line 21 of input file"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Line 21 (index 20) contains Z,N values
    if len(lines) > 20:
        # Parse the current line to preserve format
        parts = lines[20].split(',')
        parts[0] = str(Z)  # Update Z (protons)
        parts[1] = str(N)  # Update N (neutrons)
        lines[20] = ','.join(parts)
        
        with open(file_path, 'w') as f:
            f.writelines(lines)

# Run for each nucleus
for nucleus, properties in nuclei.items():
    Z = properties['Z']
    A = properties['A']
    N = A - Z  # Calculate neutron number
    
    # Create main output folder for nucleus
    nucleus_dir = os.path.join(output_base_dir, nucleus)
    os.makedirs(nucleus_dir, exist_ok=True)
    
    # Run executable for each input file type
    print(f"Running {nucleus} (Z={Z}, N={N})...")
    for input_file, folder_type in input_configs.items():
        # Create specific folder for this calculation type (ZR or FR)
        calc_dir = os.path.join(nucleus_dir, folder_type)
        os.makedirs(calc_dir, exist_ok=True)
        
        # Copy and update input file
        src = os.path.join(input_dir, input_file)
        dst = os.path.join(calc_dir, input_file)
        if os.path.exists(src):
            shutil.copy(src, dst)
            # Update Z,N values for this nucleus
            update_input_file(dst, Z, N)
        
            print(f"Processing {folder_type} ({input_file})...")
            try:
                # Check if executable exists
                if not os.path.exists(exe_path):
                    print(f"Executable not found: {exe_path}")
                    continue
                
                # List files before execution
                files_before = set(os.listdir(calc_dir))
                
                # Run with input redirection: executable < input_file
                with open(dst, 'r') as f:
                    result = subprocess.run(
                        exe_path,
                        stdin=f,
                        cwd=calc_dir,
                        capture_output=False,  # Don't capture output to allow files to be written properly
                        text=True
                    )
                
                # Wait a moment for files to be written
                import time
                time.sleep(1)
                
                # List files after execution
                files_after = set(os.listdir(calc_dir))
                new_files = files_after - files_before
                
                print(f"Return code: {result.returncode}")
                if new_files:
                    print(f"New files created: {', '.join(sorted(new_files))}")
                else:
                    print(f"No new files created")
                    # Let's also check if there are any files in the directory
                    print(f"Current files in directory: {', '.join(sorted(os.listdir(calc_dir)))}")
                
                if result.returncode == 0:
                    print(f"{nucleus} {folder_type} ({input_file}) completed successfully")
                else:
                    print(f"{nucleus} {folder_type} ({input_file}) failed with return code {result.returncode}")
                    if result.stderr:
                        print(f"Error output: {result.stderr}")
                    if result.stdout:
                        print(f"Standard output: {result.stdout}")
                        
            except Exception as e:
                print(f"Error running {nucleus} {folder_type} ({input_file}): {e}")
        else:
            print(f"Input file not found: {src}")

print("All nuclei processed.")