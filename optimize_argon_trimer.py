# 1-1-2-1. checkpoint 1) The total potential energy of Ar3 is the sum of the Lennard-Jones interactions between all three pairs of atoms
def total_potential_energy(r12, r13, r23, epsilon=0.01, sigma=3.4):
    # Calculate the Lennard-Jones potential for each pair of distances
    V12 = lennard_jones(r12, epsilon, sigma)
    V13 = lennard_jones(r13, epsilon, sigma)
    V23 = lennard_jones(r23, epsilon, sigma)
    
    # Total potential energy is the sum of the potentials for each pair
    V_total = V12 + V13 + V23
    return V_total

# Example test with arbitrary distances between atoms (r12, r13, r23)
r12_example = 4.0  # Å
r13_example = 4.0  # Å
r23_example = 4.0  # Å

# Calculate the total potential energy for the Argon trimer
V_total_example = total_potential_energy(r12_example, r13_example, r23_example)

print(V_total_example)

# 1-1-2-2. checkpoint 1) Use scipy.optimize.minimize to find the optimal geometry of Ar3. Represent the positions of the three Ar atoms in two dimensions:
from scipy.optimize import minimize

# Redefine the Lennard-Jones potential
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    potential = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
    return potential

# Redefine the total potential energy for the Argon trimer based on pairwise distances
def total_potential_energy(r12, r13, r23, epsilon=0.01, sigma=3.4):
    V12 = lennard_jones(r12, epsilon, sigma)
    V13 = lennard_jones(r13, epsilon, sigma)
    V23 = lennard_jones(r23, epsilon, sigma)
    return V12 + V13 + V23

# Using compute_bond_length function from Homework 1
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

# Define the total potential energy in terms of r12, x3, and y3
def total_potential_energy_trimer(variables, epsilon=0.01, sigma=3.4):
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

import math

# Define the atom coordinates based on the optimized values
def get_atom_coordinates(optimized_r12, optimized_x3, optimized_y3):
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