import os
import re
import pandas as pd
import fortranformat as ff


def is_float(string):
    # Compile a regular expression pattern to match valid float values
    pattern = r"^[-+]?[0-9]*\.?[0-9]+$"

    # Use re.match to check if the string matches the pattern
    # Returns a match object if there is a match, else None
    match = re.match(pattern, string)

    # Convert the match object to a boolean value
    # Returns True if there is a match, else False
    return bool(match)


def parse_tables(file_content, keywords):
    # Dictionary to store parsed tables with their keywords
    parsed_tables = {}

    # Create a regular expression pattern to match the keywords
    keyword_pattern = re.compile(
        rf'^\s*({"|".join(re.escape(kw) for kw in keywords)})\s*$'
    )

    lines = file_content.splitlines()
    current_keyword = None
    table_data = []

    for line in lines:
        line = line.strip()

        # Check if the line is a keyword line
        if keyword_pattern.match(line):
            # If we're already capturing a table, save the previous one
            if current_keyword and table_data:
                columns = [
                    "a",
                    "z",
                    "f(exp)",
                    "f(fit)",
                    "difr(%)",
                    "stdev(%)",
                    "chi",
                    "sum chi**2",
                ]
                df = pd.DataFrame(table_data, columns=columns)
                parsed_tables[current_keyword] = df
                table_data = []

            current_keyword = line.strip()

        # Check if the line is a valid table row
        elif current_keyword and re.match(r"^[\d\s\.\-E]+$", line):
            # Convert the values into a list of floats
            values = [
                float(x) for x in re.split(r"\s+", line) if re.match(r"[\d\.\-E]+", x)
            ]
            table_data.append(values)

        # If we encounter a non-table line, stop reading the table
        elif current_keyword and table_data:
            columns = [
                "a",
                "z",
                "f(exp)",
                "f(fit)",
                "difr(%)",
                "stdev(%)",
                "chi",
                "Net sum chi**2",
                "sum chi**2",
            ]
            df = pd.DataFrame(table_data, columns=columns)
            parsed_tables[current_keyword] = df
            current_keyword = None
            table_data = []

    # Save the last table if it exists
    if current_keyword and table_data:
        columns = [
            "a",
            "z",
            "f(exp)",
            "f(fit)",
            "difr(%)",
            "stdev(%)",
            "chi",
            "Net sum chi**2",
            "sum chi**2",
        ]
        df = pd.DataFrame(table_data, columns=columns)
        parsed_tables[current_keyword] = df

    return parsed_tables


def main():
    dirname = "."
    keywords = [
        "binding-energy",
        "diffr.-radius",
        "surface thickness",
        "rms radius",
        "gap-proton",
        "gap-neutr.",
    ]
    inputFile = f"for002.SVmin-softcut"
    file_content = open(inputFile, "r").read()
    writeData = parse_tables(file_content, keywords)

    # Output the parsed tables
    for keyword, df in writeData.items():
        print(f"Table: {keyword}")
        print(df)
        print("\n")


if __name__ == "__main__":
    main()
