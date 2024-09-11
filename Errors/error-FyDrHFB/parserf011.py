import re
import pandas as pd
import matplotlib.pyplot as plt


def parse_table(lines, start_line, colNames):
    data = []
    # colNames = []
    for line in lines[start_line:]:
        if line.strip() == "" or re.match(r"^\s*[a-zA-Z]", line):
            break
        data.append([float(x) for x in line.split()])

    # Dynamically determine the number of columns
    num_columns = len(data[0]) if data else 0
    if len(colNames) < num_columns:
        diff = num_columns - len(colNames)
        for i in range(diff):
            colNames.append(f"col_{i+1}")
    return pd.DataFrame(data, columns=colNames)


def extract_tables(filename, table_keywords, skip_keywords):
    tables = {}
    with open(filename, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        for keyword in table_keywords:
            colNames = []
            if keyword in line:
                # Skip the table if it's in the skip_keywords list
                if any(skip_keyword in line for skip_keyword in skip_keywords):
                    break
                colNames = lines[i + 1].split()
                tables[keyword] = parse_table(lines, i + 2, colNames)
                break
        i += 1

    return tables


def plot_df(df, x, y):
    df.plot(x, y)
    plt.show()


# Define the keywords for the tables you want to extract
table_keywords = ["densities", "formfactor", "the potentials", "l*s potentials"]

# List of keywords to skip
skip_keywords = ["integrated (grad rho)**2"]

# Usage
filename = "for011"  # Replace with your file path
tables = extract_tables(filename, table_keywords, skip_keywords)

# for name, df in tables.items():
#     print(f"Table: {name}")
#     print(df)
#     print("\n")

for name, df in tables.items():
    if name == "densities":
        plot_df(df, "r", "rho-coul")
