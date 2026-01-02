import re

line = 'Z,N,r_C,E,epair,gap,sumuv=  50  70  4.646077  -1044.12      -0.000-132.601   0.000   5.401   0.000  24.552'

# Remove the text labels and the '=' sign
# This leaves us with a string containing only numbers and whitespace
clean_line = re.sub(r'Z,N,r_C,E,epair,gap,sumuv=', '', line)

# Use a regular expression to find all floating-point numbers
# This pattern looks for a sequence of digits, optional minus signs, and decimal points
numbers_as_strings = re.findall(r'-?\d+\.\d+|-?\d+', clean_line)

# Convert the list of strings to a list of floats
numbers_as_floats = [float(num) for num in numbers_as_strings]

print(numbers_as_floats)