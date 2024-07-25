import numpy as np
import re
from io import StringIO

# Step 1: Read the data from the file
file_path = "./collect_cat/collect_cat_SVbas_0.0d0_hfb.dat"
with open(file_path, "r") as file:
    data = file.read()


# Step 2: Preprocess the data to insert spaces where needed
def preprocess_line(line):
    # Add a space before negative signs that are not part of scientific notation
    line = re.sub(r"(\d)(-)", r"\1 -", line)
    # Ensure that six decimal places are not split by adding spaces correctly
    line = re.sub(r"(-?\d\.\d{6})(?=-?\d)", r"\1 ", line)
    return line


# Process each line
processed_lines = []
for line in data.splitlines():
    if line.startswith(" Z,N,r_C,E,epair,gap,sumuv="):
        # This is a data line, so we need to fix it
        line = preprocess_line(line).split()
    processed_lines.append(line)

print(processed_lines[205])
# Join the processed lines into a single string
# processed_data = "\n".join(processed_lines)

# print(len(processed_data))
# Step 3: Use StringIO to treat the processed_data string as a file for np.genfromtxt
# data_io = StringIO(processed_data)

# Skip the non-data lines and read only the data lines
# np_data = np.genfromtxt(data_io, skip_header=1, delimiter=None, usecols=range(2, 12))
# print(np_data)
