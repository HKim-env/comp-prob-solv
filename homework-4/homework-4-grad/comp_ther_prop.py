import numpy as np
from scipy.constants import k, h
from python_code import partition_function
from python_code import T_values
import csv
import os

# I re-utilized previous codes in python_code.py. I also uploaded the new version of python_code.

# # 3-1-2 Checking point 1) Write Python functions to compute U and CV from the partition function
# Function to calculate U using np.gradient

# Internal energy calculation using np.gradient
def internal_energy(T_values):
    """
    Computes the internal energy U(T) for a range of temperatures using the partition function Z(T).

    The internal energy is calculated as:
    U(T) = -d(ln(Z))/d(beta), where beta = 1 / (k_B * T)

    Parameters:
    - T_values (array-like): Array of temperature values in Kelvin.

    Returns:
    - U_values (array-like): Array of internal energy values corresponding to each temperature.
    """
    Z_values = np.array([partition_function(T) for T in T_values])
    lnZ_values = np.log(Z_values)
    beta_values = 1 / (k * T_values)
    dlnZ_dBeta = np.gradient(lnZ_values, beta_values)
    U_values = -dlnZ_dBeta  # Internal energy U = -d(lnZ)/d(beta)
    return U_values

# Heat capacity calculation using np.gradient
def heat_capacity(T_values):
    """
    Computes the heat capacity C_V(T) for a range of temperatures based on internal energy U(T).

    The heat capacity is calculated as:
    C_V(T) = dU/dT, where U is the internal energy and T is the temperature.

    Parameters:
    - T_values (array-like): Array of temperature values in Kelvin.

    Returns:
    - C_V_values (array-like): Array of heat capacity values corresponding to each temperature.
    - U_values (array-like): Array of internal energy values corresponding to each temperature.
    """
    U_values = internal_energy(T_values)
    C_V_values = np.gradient(U_values, T_values)
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

