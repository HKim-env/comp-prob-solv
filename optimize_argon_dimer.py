import numpy as np

# 1-1-1-1 Checking point 1) Write a Python function lennard_jones(r, epsilon=0.01, sigma=3.4) that takes the distance r between two Ar atoms and returns the potential energy V (r) using the Lennard-Jones formula.
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    # Calculate the Lennard-Jones potential using the formula
    potential = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
    return potential

# Test the function with an example value of r
r_test = 4.0  # example interatomic distance in Angstroms
potential_test = lennard_jones(r_test)

potential_test

# 1-1-1-2. Checking point 1) Use scipy.optimize.minimize to find the distance between two Ar atoms that minimizes the Lennard-Jones potential
# import scipy.optimize.minimize to find roots by minimizing the absolute value of a function
from scipy.optimize import minimize

# Define the function to optimize (minimize)
def lennard_jones_minimize(r, epsilon=0.01, sigma=3.4):
    # Since minimize works with scalar functions, we return only the potential energy
    return lennard_jones(r[0], epsilon, sigma)

# 1-1-1-2. Checking point 2) Start with an initial guess of r = 4  ̊A.
initial_guess = [4.0]  # Starting with r = 4 Å

# Perform the optimization
result = minimize(lennard_jones_minimize, initial_guess, bounds=[(0.1, None)])  # Setting lower bound for r (lower bound set at 0.1 (to avoid non-physical values such as zero or negative distances)

# Optimized distance
optimized_distance = result.x[0]
optimized_potential = lennard_jones(optimized_distance)

print(optimized_distance, optimized_potential)

import matplotlib.pyplot as plt

# 1-1-1-3. Checking point 1) Plot the Lennard-Jones potential V (r) as a function of the distance r between 3  ̊A ≤ r ≤ 6  ̊A
r_values = np.linspace(3, 6, 500)

# Calculate Lennard-Jones potential for each distance
potential_values = [lennard_jones(r) for r in r_values]

# 1-1-1-3. Checking point 2) Mark the equilibrium distance (i.e., the distance at the minimum potential) on the plot.
plt.figure(figsize=(10, 8))
plt.plot(r_values, potential_values, label='Lennard-Jones Potential')
plt.axvline(x=optimized_distance, color='r', linestyle='--', label=f'Equilibrium Distance = {optimized_distance:.2f} Å')

# Add labels and title
plt.xlabel('Distance (Å)')
plt.ylabel('Potential Energy (eV)')
plt.title('Lennard-Jones Potential Curve for Ar₂')
plt.legend()

# Show the plot
plt.grid(True)
plt.show()