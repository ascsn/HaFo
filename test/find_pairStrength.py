# Code to find the pairing strength for 120Sn using HaFo code
# Written: Sudhanva Lalit, August 2025
import numpy as np
import subprocess
import re


def run_simulation(input_parameter):
    """
    This function runs the simulation, and reads the output.
    Modify this function to match your specific file formats and commands.
    """

    process = subprocess.run(['sh', 'scancut_120Sn.bat', f"{input_parameter}"], check=True, capture_output=True, text=True)
    if process.stderr:
        print("Error:", process.stderr)
        return None
    # Read the output file to extract the required values
    try:
        with open('for011', 'r') as f:
            lines = f.readlines()
            for line in lines:
                values = line.strip()
                if "Z,N,r_C,E,epair,gap,sumuv=" in values:
                    clean_line = re.sub(r'Z,N,r_C,E,epair,gap,sumuv=', '', values)
                    numbers_as_strings = re.findall(r'-?\d+\.\d+|-?\d+', clean_line)
                    numbers_as_floats = [float(num) for num in numbers_as_strings]
                    print("Values:", values)
                    (
                        Z,
                        N,
                        r_C,
                        E,
                        epair_p,
                        epair_n,
                        gap_p,
                        gap_n,
                        sumuv_p,
                        sumuv_n,
                    ) = map(float, numbers_as_floats)
            return gap_n  # Return the neutron gap as the output value
                
    except FileNotFoundError:
        print("Error: output file not found.")
        return None
    except (ValueError, IndexError) as e:
        print(f"Parsing error: {e}")
        return None

def find_parameter_for_output_secant(target_output, initial_guess1, initial_guess2, tolerance=1e-5, max_iterations=100):
    """
    Finds the input parameter using the Secant Method.
    Requires two initial guesses.
    """
    x0 = initial_guess1
    x1 = initial_guess2

    f0 = run_simulation(x0) - target_output
    f1 = run_simulation(x1) - target_output
    print(f0, f1)

    if f0 is None or f1 is None:
        return None

    for i in range(max_iterations):
        # Secant method formula
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        
        f2 = run_simulation(x2) - target_output
        
        if f2 is None:
            return None

        print(f"Iteration {i+1}: Param = {x2}, Output = {f2 + target_output}")

        # Check for convergence
        if abs(f2) < tolerance:
            print(f"Found parameter: {x2} (Output: {f2 + target_output})")
            return x2

        # Update values for the next iteration
        x0 = x1
        x1 = x2
        f0 = f1
        f1 = f2

        if f1 == f0:
            print("Division by zero in Secant method; stopping iteration.")
            x1 *= 1.1  # Slightly adjust to avoid division by zero
            f1 = run_simulation(x1) - target_output
    
    print("Could not converge to the target output within the maximum number of iterations.")
    return None

def find_parameter_for_output(target_output, initial_guess, tolerance=1e-5, max_iterations=100):
    """
    Finds the input parameter that produces the target output value using a simple bisection-like method.
    """
    lower_bound = initial_guess * 0.9  # Initial guess for a range
    upper_bound = initial_guess * 1.1

    for i in range(max_iterations):
        # Calculate the middle of the current search range
        mid_param = (lower_bound + upper_bound) / 2
        
        # Run the simulation with the midpoint parameter
        output = run_simulation(mid_param)
        
        if output is None:
            return None

        # Check if the output is close enough to the target
        if abs(output - target_output) < tolerance:
            print(f"Found parameter: {mid_param} (Output: {output})")
            return mid_param

        # Adjust the search range based on the output
        if output < target_output:
            # The parameter is too low, so increase the lower bound
            lower_bound = mid_param
        else:
            # The parameter is too high, so decrease the upper bound
            upper_bound = mid_param
            
        print(f"Iteration {i+1}: Param = {mid_param}, Output = {output}")

    print("Could not converge to the target output within the maximum number of iterations.")
    return None

def brent_method(func, a, b, target, tol=1e-5, max_iter=100):
    """
    Finds the input parameter for a given target output using Brent's method.
    
    Args:
        func (function): The function to find the root of.
                         In our case, this will be run_simulation.
        a (float): A lower bound for the input parameter that brackets the root.
        b (float): An upper bound for the input parameter that brackets the root.
        target (float): The desired output value.
        tol (float): The tolerance for the solution.
        max_iter (int): Maximum number of iterations.
        
    Returns:
        float: The found input parameter.
    """

    # We are finding the root of the function f(x) - target = 0
    # So we define a helper function to wrap run_simulation
    def f(x):
        result = func(x)
        if result is None:
            raise ValueError("Simulation failed to produce a result.")
        return result - target

    fa = f(a)
    fb = f(b)

    if fa * fb >= 0:
        print("Error: The initial guesses do not bracket the root.")
        return None

    c = a
    fc = fa
    d = b - a
    e = d

    for i in range(max_iter):
        if abs(fc) < abs(fb):
            a, b, c = b, a, b
            fa, fb, fc = fb, fa, fb

        p = 0.0
        q = 0.0
        s = 0.0

        if abs(fc) > tol and abs(fb) > tol:
            # Inverse quadratic interpolation
            if fa != fc and fb != fc:
                s = a * fb * fc / ((fa - fb) * (fa - fc)) + \
                    b * fa * fc / ((fb - fa) * (fb - fc)) + \
                    c * fa * fb / ((fc - fa) * (fc - fb))
            # Secant method
            else:
                s = b - fb * (b - a) / (fb - fa)
        else:
            # Bisection method if a bracketed root is found
            s = (a + b) / 2.0

        # Conditions for bisection step
        if not ((3 * a + b) / 4.0 < s < b) or \
           (abs(s - b) >= abs(e) / 2.0) or \
           (abs(s - b) >= abs(d) / 2.0) or \
           (abs(d) < tol):
            s = (a + b) / 2.0
            d = b - a
            e = d
        else:
            d = c - b
            e = d
        
        fs = f(s)
        
        # Update bracket
        c = b
        fc = fb
        if fa * fs < 0:
            b = s
            fb = fs
        else:
            a = s
            fa = fs
        
        # Check for convergence
        if abs(fb) < tol:
            return b

    print("Could not converge within the maximum number of iterations.")
    return b


if __name__ == '__main__':
    # Define your target output and an initial guess for the input parameter
    TARGET_OUTPUT = 1.310
    INITIAL_GUESS = 20000  # Make a reasonable initial guess
    upper_bound = INITIAL_GUESS * 1.1

    final_parameter = find_parameter_for_output_secant(TARGET_OUTPUT, INITIAL_GUESS, INITIAL_GUESS * 1.1)

    # final_parameter = brent_method(run_simulation, INITIAL_GUESS, upper_bound, TARGET_OUTPUT)

    if final_parameter:
        print(f"\nFinal input parameter for an output of {TARGET_OUTPUT} is: {final_parameter}")