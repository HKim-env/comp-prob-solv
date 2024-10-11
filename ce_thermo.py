import numpy as np
from scipy.constants import k  # Boltzmann constant
import matplotlib.pyplot as plt
import os
import csv

def numerical_derivative(func, T, delta=1e-5):
    """
    Computes the numerical derivative of a function `func` at temperature `T`
    using the central difference method.

    Parameters:
    - func: Function to compute the derivative of
    - T: Temperature at which the derivative is evaluated
    - delta: Small change in T used for finite differences

    Returns:
    - The numerical derivative at the given temperature T
    """
    return (func(T + delta) - func(T - delta)) / (2 * delta)

# 2-1-1 Checking point 1) Write a Python function to compute the partition function for an isolated Ce3+ ion with 14-fold degeneracy, where all states have zero energy.
# Isolated Ce3+ with 14-fold degeneracy, zero energy for all states
def partition_function_isolated(T):
    g_iso = 14  # Degeneracy (constant partition function)
    Z_iso = g_iso  # Partition function is constant
    return Z_iso

# 2-1-1 Checking point 2) Compute the thermodynamic properties 
# Thermodynamic properties based on the given expressions
def thermodynamic_properties_isolated(T):
    Z = partition_function_isolated(T)
    
    # Internal Energy: U = - d(ln Z) / d(beta), where beta = 1 / (k_B T)
    beta = 1 / (k * T)
    
    # Calculate ln(Z) for the partition function
    lnZ = np.log(Z)
    
    # Since Z is constant (degeneracy only), ln(Z) is constant
    # Derivative of ln(Z) with respect to beta (since ln(Z) is constant, d(ln Z)/d(beta) = 0)
    dlnZ_dbeta = 0
    U = -dlnZ_dbeta  # Internal energy, derived from the correct expression

    # Free Energy: F = - k_B T ln(Z)
    F = - k * T * lnZ
    
    # Use the numerical derivative of free energy for entropy
    S = - numerical_derivative(lambda T: -k * T * np.log(partition_function_isolated(T)), T)
    
    return U, F, S, Z

# 2-1-2 Checking point 1) Write a Python function to compute the partition function for Ce3+ with SOC
# Ce3+ with SOC: Two levels with degeneracies 6 and 8, separated by 0.28 eV
def partition_function_soc(T):
    g_soc = np.array([6, 8])  # Degeneracies
    E_soc = np.array([0, 0.28 * 1.60218e-19])  # Energy levels in Joules (0.28 eV to J)
    beta = 1 / (k * T)
    Z_soc = np.sum(g_soc * np.exp(-beta * E_soc))  # Partition function
    return Z_soc, g_soc, E_soc

# 2-1-2 Checking point 2) Use the partition function to compute the thermodynamic properties 
def thermodynamic_properties_soc(T):
    Z, g_soc, E_soc = partition_function_soc(T)  # Ensure g_soc and E_soc are retrieved
    
    # Internal Energy: U = - d(ln Z) / d(beta)
    beta = 1 / (k * T)
    # Compute internal energy from the energy levels
    U = np.sum((g_soc * E_soc * np.exp(-beta * E_soc)) / Z)  # Weighted sum over energy levels
    
    # Free Energy: F = - k_B T ln(Z)
    F = - k * T * np.log(Z)
    
    # Use the numerical derivative of free energy for entropy
    S = -numerical_derivative(lambda T: -k * T * np.log(partition_function_soc(T)[0]), T)
    
    return U, F, S, Z

# 2-1-3 Checking point 2) Write a Python function to compute the partition function for Ce3+ with SOC and CFS
# Ce3+ with SOC and CFS: Five distinct levels with different degeneracies
def partition_function_soc_cfs(T):
    g_soc_cfs = np.array([4, 2, 2, 4, 2])  # Degeneracies
    E_soc_cfs = np.array([0, 0.12, 0.25, 0.32, 0.46]) * 1.60218e-19  # 2-1-3 Checking point 1) Split into five individual energy levels
    beta = 1 / (k * T)
    Z_soc_cfs = np.sum(g_soc_cfs * np.exp(-beta * E_soc_cfs))  # Partition function
    return Z_soc_cfs, g_soc_cfs, E_soc_cfs

# 2-1-3 Checking point 3) Use the partition function to compute the thermodynamic properties 
def thermodynamic_properties_soc_cfs(T):
    Z, g_soc_cfs, E_soc_cfs = partition_function_soc_cfs(T)  # Ensure g_soc_cfs and E_soc_cfs are retrieved
    
    # Internal Energy: U = - d(ln Z) / d(beta)
    beta = 1 / (k * T)
    # Compute internal energy directly from the energy levels
    U = np.sum((g_soc_cfs * E_soc_cfs * np.exp(-beta * E_soc_cfs)) / Z)  # Weighted sum over energy levels
    
    # Free Energy: F = - k_B T ln(Z)
    F = - k * T * np.log(Z)
    
    # Use the numerical derivative of free energy for entropy
    S = -numerical_derivative(lambda T: -k * T * np.log(partition_function_soc_cfs(T)[0]), T)
    
    return U, F, S, Z

# Temperature range
T_values = np.linspace(300, 2000, 100)  # From 300 K to 2000 K

# Arrays to store the data
data_isolated = []
data_soc = []
data_soc_cfs = []

# Loop over the temperature values and compute properties for each case
for T in T_values:
    U_iso, F_iso, S_iso, Z_iso = thermodynamic_properties_isolated(T)
    U_soc, F_soc, S_soc, Z_soc = thermodynamic_properties_soc(T)
    U_soc_cfs, F_soc_cfs, S_soc_cfs, Z_soc_cfs = thermodynamic_properties_soc_cfs(T)
    
    data_isolated.append([T, U_iso, F_iso, S_iso])
    data_soc.append([T, U_soc, F_soc, S_soc])
    data_soc_cfs.append([T, U_soc_cfs, F_soc_cfs, S_soc_cfs])

# Define the temperature range
T_values = np.linspace(300, 2000, 100)

# Arrays to store the thermodynamic properties and partition functions
U_iso_values, F_iso_values, S_iso_values, Z_iso_values = [], [], [], []
U_soc_values, F_soc_values, S_soc_values, Z_soc_values = [], [], [], []
U_soc_cfs_values, F_soc_cfs_values, S_soc_cfs_values, Z_soc_cfs_values = [], [], [], []

# 2-1-4 Checking point 1) Compute and compare the partition functions and the thermodynamic properties 
# Compute thermodynamic properties and partition functions for each temperature
for T in T_values:
    U_iso, F_iso, S_iso, Z_iso = thermodynamic_properties_isolated(T)
    U_soc, F_soc, S_soc, Z_soc = thermodynamic_properties_soc(T)
    U_soc_cfs, F_soc_cfs, S_soc_cfs, Z_soc_cfs = thermodynamic_properties_soc_cfs(T)
    
    # Store the results for isolated Ce3+
    U_iso_values.append(U_iso)
    F_iso_values.append(F_iso)
    S_iso_values.append(S_iso)
    Z_iso_values.append(Z_iso)
    
    # Store the results for Ce3+ with SOC
    U_soc_values.append(U_soc)
    F_soc_values.append(F_soc)
    S_soc_values.append(S_soc)
    Z_soc_values.append(Z_soc)
    
    # Store the results for Ce3+ with SOC & CFS
    U_soc_cfs_values.append(U_soc_cfs)
    F_soc_cfs_values.append(F_soc_cfs)
    S_soc_cfs_values.append(S_soc_cfs)
    Z_soc_cfs_values.append(Z_soc_cfs)

# 2-1-4 Checking point 2) Provide plots of the thermodynamic properties as a function of temperature (from 300 K to 2000 K)
# Plot Partition Functions for the three cases
plt.figure(figsize=(10, 6))
plt.plot(T_values, Z_iso_values, label='Isolated Ce3+')
plt.plot(T_values, Z_soc_values, label='Ce3+ with SOC')
plt.plot(T_values, Z_soc_cfs_values, label='Ce3+ with SOC & CFS')
plt.xlabel('Temperature (K)')
plt.ylabel('Partition Function Z')
plt.title('Partition Function vs Temperature')
plt.legend()
plt.grid(True)
plt.show()

# Plot Internal Energy for the three cases
plt.figure(figsize=(10, 6))
plt.plot(T_values, U_iso_values, label='Isolated Ce3+')
plt.plot(T_values, U_soc_values, label='Ce3+ with SOC')
plt.plot(T_values, U_soc_cfs_values, label='Ce3+ with SOC & CFS')
plt.xlabel('Temperature (K)')
plt.ylabel('Internal Energy (J)')
plt.title('Internal Energy vs Temperature')
plt.legend()
plt.grid(True)
plt.show()

# Plot Free Energy for the three cases
plt.figure(figsize=(10, 6))
plt.plot(T_values, F_iso_values, label='Isolated Ce3+')
plt.plot(T_values, F_soc_values, label='Ce3+ with SOC')
plt.plot(T_values, F_soc_cfs_values, label='Ce3+ with SOC & CFS')
plt.xlabel('Temperature (K)')
plt.ylabel('Free Energy (J)')
plt.title('Free Energy vs Temperature')
plt.legend()
plt.grid(True)
plt.show()

# Plot Entropy for the three cases
plt.figure(figsize=(10, 6))
plt.plot(T_values, S_iso_values, label='Isolated Ce3+')
plt.plot(T_values, S_soc_values, label='Ce3+ with SOC')
plt.plot(T_values, S_soc_cfs_values, label='Ce3+ with SOC & CFS')
plt.xlabel('Temperature (K)')
plt.ylabel('Entropy (J/K)')
plt.title('Entropy vs Temperature')
plt.legend()
plt.grid(True)
plt.show()


# Get the current working directory
current_directory = os.getcwd()

# Specify the folder where the file will be saved, relative to the current directory
folder_name = "comp-prob-solv/homework-4-2"
directory = os.path.join(current_directory, folder_name)

# Create the directory if it does not exist
os.makedirs(directory, exist_ok=True)

# Write the results to CSV files for each case

# Isolated Ce3+
csv_file_path_isolated = os.path.join(directory, "thermo_isolated.csv")
with open(csv_file_path_isolated, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow(['Temperature (K)', 'Internal Energy (J)', 'Free Energy (J)', 'Entropy (J/K)'])
    # Write data rows
    writer.writerows(data_isolated)

# Ce3+ with SOC
csv_file_path_soc = os.path.join(directory, "thermo_soc.csv")
with open(csv_file_path_soc, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Temperature (K)', 'Internal Energy (J)', 'Free Energy (J)', 'Entropy (J/K)'])
    writer.writerows(data_soc)

# Ce3+ with SOC and CFS
csv_file_path_soc_cfs = os.path.join(directory, "thermo_soc_cfs.csv")
with open(csv_file_path_soc_cfs, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Temperature (K)', 'Internal Energy (J)', 'Free Energy (J)', 'Entropy (J/K)'])
    writer.writerows(data_soc_cfs)

print(f"CSV files successfully created in {directory}")