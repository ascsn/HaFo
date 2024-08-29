import os
import re
import pandas as pd
import fortranformat as ff


def parse_tables(file_content, keywords):
    # Dictionary to store parsed tables with their keywords
    parsed_tables = {}
    # Dictionary to store the r.m.s. error (abs) values
    rms_errors = {}
    # Dictionary to store the average error (abs) values
    avg_errors = {}

    # Create a regular expression pattern to match the keywords
    keyword_pattern = re.compile(
        rf'^\s*({"|".join(re.escape(kw) for kw in keywords)})\s*$'
    )

    # Combine rms and avg patterns into one to match both on the same line
    error_pattern = re.compile(
        r"r\.m\.s\. error \(abs\)=\s*([\d\.\-E]+)\s+average error \(abs\)=\s*([\d\.\-E]+)"
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

        # If the line contains the average error (abs) value, capture it
        elif current_keyword and error_pattern.search(line):
            rms_value, avg_value = error_pattern.search(line).groups()
            avg_errors[current_keyword] = float(avg_value)
            rms_errors[current_keyword] = float(rms_value)
            # Reset current_keyword after extracting error values
            current_keyword = None

        # Check if the line is a valid table row
        elif current_keyword and re.match(r"^[\d\s\.\-E]+$", line):
            # Convert the values into a list of floats
            values = [
                float(x) for x in re.split(r"\s+", line) if re.match(r"[\d\.\-E]+", x)
            ]
            table_data.append(values)

        # If we encounter a non-table line, stop reading the table
        elif current_keyword and table_data and not re.match(r"^[\d\s\.\-E]+$", line):
            columns = [
                "a",
                "z",
                "f(exp)",
                "f(fit)",
                "difr(%)",
                "stdev(%)",
                "chi",
                "net_sum_chi**2",
                "sum_chi**2",
            ]
            df = pd.DataFrame(table_data, columns=columns)
            parsed_tables[current_keyword] = df
            # current_keyword = None
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
            "net_sum_chi**2",
            "sum_chi**2",
        ]
        df = pd.DataFrame(table_data, columns=columns)
        parsed_tables[current_keyword] = df

    return parsed_tables, rms_errors, avg_errors


def write_data(df, fileName):
    # Write dataframe to a file
    colNames = df.columns.tolist()
    format_string = "(4f8.3, a30)"
    header_line = ff.FortranRecordWriter(format_string)
    Formatted_df = df.apply(lambda x: header_line.write(x.values), axis=1)
    with open(fileName, "w") as fp:
        fp.write("# " + " ".join(colNames) + "\n")
    Formatted_df.to_csv(fileName, mode="a", index=False, header=False)


def main():
    dirname = "."
    keywords = [
        "binding-energy",
        "diffr.-radius",
        "surface thickness",
        "rms radius",
        # "gap-proton",
        # "gap-neutr.",
    ]

    inputFile = f"for002.SVmin-width15"
    names = inputFile.split(".")
    file_content = open(inputFile, "r").read()
    writeData, rms_error, avg_errors = parse_tables(file_content, keywords)

    print(rms_error)
    # Output the parsed tables
    # for keyword, df in writeData.items():
    #     print(f"Table: {keyword}")
    #     print(df)
    #     print("\n")

    # We can now extract error table as outputs:
    error_208Pb = []
    sum_chi2 = []
    rmsErrors = []
    avgErrors = []
    for keyword, df in writeData.items():
        dataVal = df.loc[df["a"] == 208]
        lastRow = df.iloc[-1]

        if not dataVal.empty:
            expVal = dataVal["f(exp)"].values[0]
            fitVal = dataVal["f(fit)"].values[0]
            error_208Pb.append(expVal - fitVal)

        if not lastRow.empty:
            sum_chi2.append(lastRow["sum_chi**2"])

        if keyword in rms_error:
            rmsErrors.append(rms_error[keyword])
            avgErrors.append(avg_errors[keyword])

    list_of_names = [
        "energy E[MeV]",
        "diffr. radius [fm]",
        "surface thickness [fm]",
        "rms radius [fm]",
    ]

    dict_of_obs = {
        "rms_error": rmsErrors,
        "average_error": avgErrors,
        "Sum_chi**2": sum_chi2,
        "error_208Pb": error_208Pb,
        "observable": list_of_names,
    }

    df = pd.DataFrame(dict_of_obs)

    write_data(df, f"python_errors.{names[1]}")


if __name__ == "__main__":
    main()
