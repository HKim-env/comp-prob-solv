# compute_bond_length function from Homework 1
import numpy as np

def compute_bond_length(coord1, coord2):
    """
    Compute the bond length between two atoms using the Cartesian coordinates.
    
    parameters
    :coord1: List of [x, y, z] coordinates of the first atom.
    :coord2: List of [x, y, z] coordinates of the second atom.
    :return: The bond length between the two atoms.
    """
    # Calculate the bond length using the distance formula
    bond_length = np.sqrt((coord2[0] - coord1[0])**2 + 
                          (coord2[1] - coord1[1])**2)
    
    return bond_length


# defining compute_bond_angle from Homework 1
def compute_bond_angle(coord1, coord2, coord3):
    """
    Compute the bond angle between three atoms using their Cartesian coordinates.
    
    Parameters
    :coord1: List of [x, y] coordinates of the first atom.
    :coord2: List of [x, y] coordinates of the second atom.
    :coord3: List of [x, y] coordinates of the third atom.
    :return: The bond angle in degrees.
    """
    # Create vectors AB and BC
    AB = np.array(coord1) - np.array(coord2)
    BC = np.array(coord3) - np.array(coord2)
    
    # Calculate the dot product and magnitudes of the vectors
    dot_product = np.dot(AB, BC)
    magnitude_AB = np.linalg.norm(AB)
    magnitude_BC = np.linalg.norm(BC)
    
    # Compute the cosine of the angle using the dot product formula
    cos_theta = dot_product / (magnitude_AB * magnitude_BC)
    
    # Calculate the angle in radians and then convert to degrees
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))  # Clip to handle numerical errors
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg

def lennard_jones(r, epsilon, sigma):
    """
    Compute the lennard_jones_potential with given equation
    
    parameters
    :r : distance
    :sigma : the value of sigma
    :epsilon: the value of epsilon
    """
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

import numpy as np
from scipy.constants import k, h
from scipy.integrate import trapezoid
from math import sqrt

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

# Thermal wavelength λ
def thermal_wavelength(T):
    return sqrt(h ** 2 / (2 * np.pi * m * k * T)) 

# Maximum relative distance in the cubic box (diagonal of the cubic box)
L_max = np.cbrt(V)


# Lennard-Jones potential between two particles at relative distance r
def lj_potential(r):
    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)

# Partition function for two LJ particles in a cubic box
def partition_function(T):
    λ = thermal_wavelength(T)  # Use thermal wavelength formula with Planck's constant
    pre_factor = (4 * np.pi) ** 2 / (λ ** 6)  # Pre-factor accounting for spherical integration for both particles

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

# Function to calculate U using np.gradient
def internal_energy(T):
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
    U_values = np.array([internal_energy(T) for T in T_values])
    C_V_values = np.gradient(U_values, T_values)  # Derivative of U with respect to T
    return C_V_values, U_values

# Calculate the dissociation temperature
def find_dissociation_temperature(C_V_values):
    max_index = np.argmax(C_V_values)
    dissociation_temp = T_values[max_index]
    max_CV = C_V_values[max_index] 
    
    # Print the maximum heat capacity and corresponding temperature only once
    print(f"Maximum Heat Capacity: {max_CV:.2e} J/K")
    print(f"Dissociation Temperature (Maximum C_V): {dissociation_temp:.2f} K")
    
    return dissociation_temp, max_CV

# Main execution
C_V_values, U_values = heat_capacity(T_values)  # Calculate U_values and C_V_values


import numpy as np
from scipy.integrate import trapezoid

# 1-1-1 Checking point 1) Write a Python function named compute_work_isothermal.py to compute the work
# 1-1-1 Checking point 2) Use scipy.integrate.trapezoid to compute the work. The pressure in this case is given by the ideal gas law

def compute_work_isothermal(V_i, V_f, n, R, T, num_points=1000):
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

import numpy as np
from scipy.integrate import trapezoid

# 1-1-2 Checking point 1) Write a Python function named compute_work_adiabatic.py to compute the work
# 1-1-2 Checking point 2) Compute the constant from the initial conditions and perform numerical integration to find the work

def compute_work_adiabatic(V_i, V_f, P_i, gamma, num_points=1000):
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