import numpy as np
from scipy.integrate import trapezoid

# 1-1-1 Checking point 1) Write a Python function named compute_work_isothermal.py to compute the work
# 1-1-1 Checking point 2) Use scipy.integrate.trapezoid to compute the work. The pressure in this case is given by the ideal gas law

def compute_work_isothermal(V_i, V_f, n, R, T, num_points=1000):
    """
    Computes the work done during an isothermal process for an ideal gas.
    
    Parameters:
    V_i (float): Initial volume in cubic meters.
    V_f (float): Final volume in cubic meters.
    n (float): Amount of substance in moles (mol).
    R (float): Ideal gas constant in J/mol-K.
    T (float): Temperature in Kelvin (constant for isothermal process).
    num_points (int): Number of points to use for numerical integration (default is 1000).
    
    Returns:
    float: The work done during the isothermal process in Joules.
    
    Explanation:
    The function computes the work done during an isothermal process based on the relation:
        W = - ∫ P dV, where P = nRT / V (ideal gas law for isothermal process).
    
    The work is computed using the trapezoidal rule for numerical integration over the specified 
    volume range (V_i to V_f). The pressure at each volume is determined using the ideal gas law.
    """
    # Define the volume array for integration
    V = np.linspace(V_i, V_f, num_points)
    
    # Define the pressure as a function of volume using the ideal gas law
    P = (n * R * T) / V
    
    # Compute the work using trapezoidal integration
    work = -trapezoid(P, V)
    
    return work

# Given parameters
n = 1      # mol
R = 8.314  # J/mol-K
T = 300    # K
V_i = 0.1  # m^3
gamma = 1.4
P_i = (n * R * T) / V_i  # From ideal gas law

# Define the range of final volumes
V_f_values = np.linspace(V_i, 3 * V_i, 100)


import csv
import os

# Initialize empty lists to store the computed work
work_isothermal = []

# Compute work for each final volume for the isothermal process
for V_f in V_f_values:
    w_iso = compute_work_isothermal(V_i, V_f, n, R, T)
    
    # Append the computed work values to the list
    work_isothermal.append(w_iso)

# Get the current working directory
current_directory = os.getcwd()

# Specify the CSV file path in the current working directory
csv_file_path = os.path.join(current_directory, "work_isothermal_vs_volume.csv")

# Write the results to a CSV file in the current directory
with open(csv_file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow(['Final Volume (m^3)', 'Work Isothermal (J)'])
    
    # Write data rows
    for V_f, w_iso in zip(V_f_values, work_isothermal):
        writer.writerow([V_f, w_iso])

print(f"CSV file successfully created at {csv_file_path}")
