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

# Partition function for two LJ particles in a cubic box
def partition_function(T):
    λ = thermal_wavelength(T)  # Use thermal wavelength formula with Planck's constant
    # In the given equation, we can convert the Eq. (10) as spherical coordinate system. Since the system is isotropic, the energy or interaction potential depends only on the radial distance. The polar and azimuthal angle can be simplified as 4 * pi.
    # because dV = r^2 * sinϕ * drdϕdθ
    # integral range: 0 ~ 2pi for  θ, 0 ~ pi for ϕ
    # I applied the 4 * pi to the pre_factor
    pre_factor = 4 * np.pi / ((h ** 6) * (λ ** 6))  # Pre-factor from Eq. (10) includes h and 4*pi

    # # 3-1-1 Checking point 2) Assume the cubic box has a fixed volume, and the integration should be performed over the relative distances between the particles
    # Perform integration over relative distance r, using spherical coordinates
    # Set the maximum relative distance based on the fixed volume
    r_min = 0.001 * sigma  # Avoid zero to prevent singularity
    r_max = L_max
    
    # Discretize the range of relative distances
    r_values = np.linspace(r_min, r_max, 1000)
    
    # Compute the Boltzmann factor for each relative distance r
    integrand_values = np.array([np.exp(-lj_potential(r) / (k * T)) * r**2 for r in r_values])  # Multiply by r^2 for spherical coords
    
    # Perform the integration using the trapezoidal rule
    Z = trapezoid(integrand_values, r_values)
    
    return pre_factor * Z  # Return a single value



# Compute partition function for all temperatures
partition_values = [partition_function(T) for T in T_values]

# Get the current working directory
current_directory = os.getcwd()

# Specify the folder where the file will be saved, relative to the current directory
folder_name = "comp-prob-solv/homework-4-grad"
directory = os.path.join(current_directory, folder_name)

# Create the directory if it does not exist
os.makedirs(directory, exist_ok=True)

# Specify the CSV file path
csv_file_path = os.path.join(directory, "partition_function_vs_temperature.csv")

# Write the results to a CSV file
with open(csv_file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow(['Temperature (K)', 'Partition Function'])
    
    # Write data rows
    for T, Z in zip(T_values, partition_values):
        writer.writerow([T, Z])

print(f"CSV file successfully created at {csv_file_path}")