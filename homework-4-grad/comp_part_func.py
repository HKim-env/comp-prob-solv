import numpy as np
from scipy.constants import k, h
from scipy.integrate import trapezoid
from math import sqrt
import csv
import os

# Given Parameters
epsilon = 0.0103 * 1.60218e-19  # Lennard-Jones energy parameter in Joules (convert eV to Joule)
sigma = 3.4e-10  # Lennard-Jones length scale in meters (converted from Å)
V = 1000 * 1e-30  # Volume of the cubic box in cubic meters (converted from Å^3)
T_min = 10  # Minimum temperature in Kelvin
T_max = 1000  # Maximum temperature in Kelvin
N_points = 100  # Number of temperature points
m = 39.948 * 1.66053904e-27  # Mass of Argon in kg 
h = 6.62607015e-34  # Planck's constant in J·s

# Temperature range
T_values = np.linspace(T_min, T_max, N_points)

# Maximum relative distance in the cubic box (diagonal of the cubic box)
L_max = np.cbrt(V)

# Thermal wavelength λ
def thermal_wavelength(T):
    return sqrt(h ** 2 / (2 * np.pi * m * k * T)) 


# Lennard-Jones potential between two particles at relative distance r
def lj_potential(r):
    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)

# 3-1-1 Checking point 1) Write a Python function that numerically computes the classical partition function of two LJ particles in a cubic box using the trapezoidal rule for spatial integration

# I converted the equation as spherical coordinate system.

def partition_function(T):
    λ = thermal_wavelength(T)  # Use thermal wavelength formula with Planck's constant
    pre_factor = (4 * np.pi) ** 2 / (λ ** 6)  # Pre-factor accounting for spherical integration for both particles // 4 * pi square because there are two particles

    # Set the minimum and maximum relative distances based on the fixed volume
    r_min = 0.001 * sigma  # Avoid zero to prevent singularity
    r_max = np.cbrt(V)

    # Discretize the range of relative distances between two particles
    r_values = np.linspace(r_min, r_max, 1000)

    # Compute the integrand for the partition function over the relative distance r
    integrand = np.exp(-lj_potential(r_values) / (k * T)) * r_values**2

    # Perform the trapezoidal integration over the relative distance r for one particle
    Z_spherical = trapezoid(integrand, r_values)

    # Square the result to account for spherical integration over both particles
    Z_total = Z_spherical ** 2

    return pre_factor * Z_total  # Return the partition function result


# Compute partition function for all temperatures
partition_values = [partition_function(T) for T in T_values]

# Get the current working directory
current_directory = os.getcwd()

# Specify the folder where the file will be saved, relative to the current directory
folder_name = "comp-prob-solv/homework-4-grad"
directory = os.path.join(current_directory, folder_name)

# Create the directory if it does not exist
try:
    os.makedirs(directory, exist_ok=True)
    print(f"Directory created or exists at: {directory}")
except Exception as e:
    print(f"Error creating directory: {e}")

# Specify the CSV file path (without duplicating the folder structure)
csv_file_path = os.path.join(directory, "partition_function_vs_temperature.csv")

# Write the results to a CSV file
try:
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(['Temperature (K)', 'Partition Function'])
        
        # Write data rows
        for T, Z in zip(T_values, partition_values):
            writer.writerow([T, Z])
    
    print(f"CSV file successfully created at {csv_file_path}")
except PermissionError:
    print(f"Permission denied: Unable to write to {csv_file_path}. Check your permissions or if the file is in use.")
except Exception as e:
    print(f"An error occurred while writing the file: {e}")
