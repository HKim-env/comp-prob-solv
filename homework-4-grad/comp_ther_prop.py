import numpy as np
from scipy.constants import k, h
from python_code import internal_energy
from python_code import partition_function
from python_code import T_values
import csv
import os

# I re-utilized previous codes in python_code.py. I also uploaded the new version of python_code.

# # 3-1-2 Checking point 1) Write Python functions to compute U and CV from the partition function
# Function to calculate U using np.gradient
def internal_energy(T):
    """
    Computes the internal energy U(T) using the partition function Z(T).
    
    Parameters:
    - T (float): Temperature in Kelvin.

    Returns:
    - U (float): Internal energy at the given temperature.
    """
    beta = 1 / (k * T)
    Z = partition_function(T)  # Assumed that this function doesn't save to a CSV
    delta_T = T * 0.01  # Small finite difference step for numerical differentiation
    Z_plus = partition_function(T + delta_T)
    Z_minus = partition_function(T - delta_T)
    delta_beta = (1 / (k * (T + delta_T))) - (1 / (k * (T - delta_T)))
    dlnZ_dBeta = (np.log(Z_plus) - np.log(Z_minus)) / delta_beta
    U = -(dlnZ_dBeta)
    return U

# Heat capacity C_V from internal energy U
def heat_capacity(T):
    """
    Calculates the heat capacity C_V(T) from the internal energy U(T).
    
    Parameters:
    - T (float): Temperature in Kelvin.

    Returns:
    - C_V_values (array-like): Heat capacity values for each temperature in T_values.
    - U_values (array-like): Internal energy values for each temperature in T_values.
    """
    U_values = np.array([internal_energy(T) for T in T_values])
    C_V_values = np.gradient(U_values, T_values)  # Derivative of U with respect to T
    return C_V_values, U_values

# Write data to CSV file in the current working directory
def write_csv(T_values, U_values, C_V_values):
    """
    Writes temperature, internal energy, and heat capacity data to a CSV file in the current directory.
    
    Parameters:
    - T_values (array-like): Array of temperature values.
    - U_values (array-like): Array of internal energy values.
    - C_V_values (array-like): Array of heat capacity values.
    """
    # Get the current working directory
    current_directory = os.getcwd()

    # Specify the CSV file path directly in the current directory
    csv_file_path = os.path.join(current_directory, "internal_energy_heat_capacity_vs_temperature.csv")

    # Write the results to a CSV file
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(['Temperature (K)', 'Internal Energy (J)', 'Heat Capacity (J/K)'])
        
        # Write data rows
        for T, U, C_V in zip(T_values, U_values, C_V_values):
            writer.writerow([T, U, C_V])

    print(f"CSV file successfully created at {csv_file_path}")

# Main execution
C_V_values, U_values = heat_capacity(T_values)  # Calculate U_values and C_V_values

# Write the internal energy and heat capacity data to a CSV file
write_csv(T_values, U_values, C_V_values)

