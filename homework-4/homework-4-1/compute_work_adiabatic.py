import numpy as np
from scipy.integrate import trapezoid

# 1-1-2 Checking point 1) Write a Python function named compute_work_adiabatic.py to compute the work
# 1-1-2 Checking point 2) Compute the constant from the initial conditions and perform numerical integration to find the work

def compute_work_adiabatic(V_i, V_f, P_i, gamma, num_points=1000):
    """
    Computes the work done during an adiabatic process for an ideal gas.
    
    Parameters:
    V_i (float): Initial volume in cubic meters.
    V_f (float): Final volume in cubic meters.
    P_i (float): Initial pressure in Pascals.
    gamma (float): Adiabatic index (ratio of specific heats, Cp/Cv).
    num_points (int): Number of points to use for numerical integration (default is 1000).
    
    Returns:
    float: The work done during the adiabatic process in Joules.
    
    Explanation:
    The function computes the work done during an adiabatic process based on the relation:
        W = ∫ P dV, where P = C / V^gamma and C is a constant determined by the initial conditions.
    
    The work is computed using the trapezoidal rule for numerical integration over the specified 
    volume range (V_i to V_f). The pressure at each volume is determined using the adiabatic relation.
    """
    # Define the volume array for integration
    V = np.linspace(V_i, V_f, num_points)
    
    # Compute the constant for the adiabatic process
    C = P_i * V_i**gamma
    
    # Define the pressure as a function of volume
    P = C / (V**gamma)
    
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
work_adiabatic = []

# Compute work for each final volume for both processes
for V_f in V_f_values:
    w_adi = compute_work_adiabatic(V_i, V_f, P_i, gamma)
    
    # Append the computed work values to the lists
    work_adiabatic.append(w_adi)

# Get the current working directory
current_directory = os.getcwd()

# Specify the CSV file path in the current working directory
csv_file_path = os.path.join(current_directory, "work_adiabatic_vs_volume.csv")

# Write the results to a CSV file in the current directory
with open(csv_file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow(['Final Volume (m^3)', 'Work Adiabatic (J)'])
    
    # Write data rows
    for V_f, w_adi in zip(V_f_values, work_adiabatic):
        writer.writerow([V_f, w_adi])

print(f"CSV file successfully created at {csv_file_path}")
