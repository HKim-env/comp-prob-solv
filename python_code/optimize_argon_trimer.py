import numpy as np # define numpy as np

# 1-1-2-1. checkpoint 1) The total potential energy of Ar3 is the sum of the Lennard-Jones interactions between all three pairs of atoms

def total_potential_energy(r12, r13, r23, epsilon=0.01, sigma=3.4):
    """
    Compute total_potential_energy.
    
    parameters
    :r12 : distance between atom 1 and atom 2
    :r13 : distance between atom 1 and atom 3
    :r23 : distance between atom 2 and atom 3
    : epsilon : the epsilon value
    : sigma : the sigma value
    """
    # Calculate the Lennard-Jones potential for each pair of distances
    V12 = lennard_jones(r12, epsilon, sigma)
    V13 = lennard_jones(r13, epsilon, sigma)
    V23 = lennard_jones(r23, epsilon, sigma)
    
    # Total potential energy is the sum of the potentials for each pair
    
    return V12 + V13 + V23

# 1-1-2-2. checkpoint 1) Use scipy.optimize.minimize to find the optimal geometry of Ar3. Represent the positions of the three Ar atoms in two dimensions:
from scipy.optimize import minimize

# Redefine the Lennard-Jones potential
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    potential = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
    return potential

from python_code import compute_bond_length

# Define the total potential energy in terms of r12, x3, and y3
def total_potential_energy_trimer(variables, epsilon=0.01, sigma=3.4):
    """
    Compute the total potential energy for Argon trimer.
    
    parameters
    :variables : variable of defined function.
    :epsilon : the epsilon number
    : sigma: the sigma number in the system
    """
    r12 = variables[0]  # distance between atom 1 and atom 2
    x3 = variables[1]   # x-coordinate of atom 3
    y3 = variables[2]   # y-coordinate of atom 3

    # Positions of the three atoms in 2D
    atom1 = [0, 0]  # 1-1-2-2. checkpoint 2) The first Atom is fixed at the origin
    atom2 = [r12, 0]  # 1-1-2-2. checkpoint 3) The second Atom is placed on the x-axis at (r12, 0).
    atom3 = [x3, y3]  # 1-1-2-2. checkpoint 4) atom’s position is (x3, y3), where both x3 and y3 are variables to be optimized.

    # Compute bond lengths using the bond length function
    r12 = compute_bond_length(atom1, atom2)  # Bond length between atom 1 and atom 2
    r13 = compute_bond_length(atom1, atom3)  # Bond length between atom 1 and atom 3
    r23 = compute_bond_length(atom2, atom3)  # Bond length between atom 2 and atom 3

    return total_potential_energy(r12, r13, r23, epsilon, sigma)

# Initial guess for the trimer geometry [r12, x3, y3]
initial_guess_trimer = [4.0, 4.0, 4.0]

# Optimize the total potential energy of the trimer
result_trimer = minimize(total_potential_energy_trimer, initial_guess_trimer, bounds=[(3.0, 6.0), (0, 6.0), (0, 6.0)]) # the bound range was decided based on the previous problem condition (3 ~ 6 A)

# Extract optimized results
optimized_r12, optimized_x3, optimized_y3 = result_trimer.x
optimized_total_energy = total_potential_energy_trimer([optimized_r12, optimized_x3, optimized_y3])

print(optimized_r12, optimized_x3, optimized_y3, optimized_total_energy)
print(f"Optimized vector for Ar3:")
print(f"(0, 0, 0): (0, 0, 0)")
print(f"(r12, 0, 0): ({optimized_r12:.2f}, 0, 0)")
print(f"(x3, y3, 0): ({optimized_x3:.2f}, {optimized_y3:.2f}, 0)")

from python_code import compute_bond_angle

import math # impor math to make an access to the Python math functions provided by the C standard

# Define the atom coordinates based on the optimized values
def get_atom_coordinates(optimized_r12, optimized_x3, optimized_y3):
    """
    Compute the coordinate of atoms.
    
    parameters
    :optimized_r12 : optimized length between atom 1 and 2.
    :optimized_x3 : optimized x3 coordinate
    :optimized_y3 : optimized y3 coordinate
    """
    atom1 = [0, 0]  # Atom 1 is at the origin
    atom2 = [optimized_r12, 0]  # Atom 2 is on the x-axis at (r12, 0)
    atom3 = [optimized_x3, optimized_y3]  # Atom 3 at (x3, y3)
    return atom1, atom2, atom3

# define coordinate of atom1, atom2, and atom3
atom1, atom2, atom3 = get_atom_coordinates(optimized_r12, optimized_x3, optimized_y3)

# Use the optimized values to compute distances and angles
r12 = compute_bond_length(atom1, atom2)
r13 = compute_bond_length(atom1, atom3)
r23 = compute_bond_length(atom2, atom3)

angle_123 = compute_bond_angle(atom1, atom2, atom3)
angle_132 = compute_bond_angle(atom1, atom3, atom2)
angle_213 = compute_bond_angle(atom2, atom1, atom3)

# 1-1-2-3. checkpoint 1) Print the optimal distances r12, r13, r23 between the atoms in the trimer.
print(f"Optimal distances between atoms in Ar3:")
print(f"r12 (distance between atom 1 and 2): {r12:.2f} Å")
print(f"r13 (distance between atom 1 and 3): {r13:.2f} Å")
print(f"r23 (distance between atom 2 and 3): {r23:.2f} Å")

# 1-1-2-3. checkpoint 2) Print the optimal angles between the atoms in the trimer.
print(f"\nOptimal angles between atoms in Ar3:")
print(f"Angle at atom 2 (∠123): {angle_123:.2f} degrees")
print(f"Angle at atom 3 (∠132): {angle_132:.2f} degrees")
print(f"Angle at atom 1 (∠213): {angle_213:.2f} degrees")

# 1-1-2-3. checkpoint 3)Comment on the geometric arrangement
is_equilateral = (np.isclose(r12, r13, atol=0.01) and np.isclose(r13, r23, atol=0.01))
geometry_comment = "The atoms form an equilateral triangle." if is_equilateral else "The atoms do not form an equilateral triangle."
print(f"\nGeometric arrangement: {geometry_comment}")

import os

# Define the local path to the directory where the file will be saved
directory_path = os.path.join("C:/Users/khh38/Desktop/PhD/Class/Class [2024 Fall]/Computational-chemistry/comp-prob-solv", "homework-2-1")

# Define a function to write the XYZ file with the optimized geometry
def write_xyz_file(file_path, atom1, atom2, atom3):
    """
    Write the atom coordinates to an .xyz file.
    
    :file_path: Full path to the .xyz file
    :atom1, atom2, atom3: Coordinates of the three Argon atoms
    """
    with open(file_path, 'w') as f:
        f.write(f"3\n")  # Number of atoms
        f.write(f"Argon trimer geometry\n")  # Comment line
        f.write(f"Ar {atom1[0]:.4f} {atom1[1]:.4f} 0.0000\n")  # Atom 1 with z=0 for 2D
        f.write(f"Ar {atom2[0]:.4f} {atom2[1]:.4f} 0.0000\n")  # Atom 2 with z=0 for 2D
        f.write(f"Ar {atom3[0]:.4f} {atom3[1]:.4f} 0.0000\n")  # Atom 3 with z=0 for 2D

# Coordinates of the three Argon atoms after optimization (replace with your optimized values)
def get_atom_coordinates(optimized_r12, optimized_x3, optimized_y3):
    atom1 = [0, 0]  # Atom 1 is at the origin
    atom2 = [optimized_r12, 0]  # Atom 2 is on the x-axis at (r12, 0)
    atom3 = [optimized_x3, optimized_y3]  # Atom 3 at (x3, y3)
    return atom1, atom2, atom3

# Get the coordinates for the three Argon atoms
atom1, atom2, atom3 = get_atom_coordinates(optimized_r12, optimized_x3, optimized_y3)

# Define the path to the .xyz file inside the homework-2-1 directory
xyz_file_path = os.path.join(directory_path, "argon_trimer_geometry.xyz")

# Write the XYZ file
write_xyz_file(xyz_file_path, atom1, atom2, atom3)

print(f"XYZ file written to {xyz_file_path}")