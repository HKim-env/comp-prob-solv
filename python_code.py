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