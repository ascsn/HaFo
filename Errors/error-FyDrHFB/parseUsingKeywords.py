"""Code to parse for011 using keywords"""

import re


# Define a function to extract values based on keywords with proper column context
def extract_values(data, keyword):
    values = []
    lines = data.strip().split("\n")
    for line in lines:
        # Match keyword followed by values
        match = re.search(rf"{keyword}\s+((?:[-+]?\d*\.\d+\s+)+)", line)
        if match:
            # Extract the matched values
            extracted_values = match.group(1).split()
            values.append(extracted_values)
    return values


def main():
    inputFile = "for011"
    file_content = open(inputFile, "r").read()

    # Example keywords to extract
    keywords = ["mass radii:", "delta_p/n", "gapfc_p/n", "rho-coul", "rho-mass"]

    # Extract and print values for each keyword
    for keyword in keywords:
        extracted_values = extract_values(file_content, keyword)
        print(f"Values for '{keyword}': {extracted_values}")


if __name__ == "__main__":
    main()
