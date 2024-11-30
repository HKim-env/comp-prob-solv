import numpy as np
import matplotlib.pyplot as plt

# Define the Bohr radius in atomic units
a0 = 1.0

# check point 1.1.1-1) Convert the hydrogen atom’s 2p orbital ψ2pz from polar coordinates (r, θ) to Cartesian coordinates
# check point 1.1.1-2) Write a Python function psi_2p_z(x, y, z) that computes ψ2pz (x, y, z) at a given point
# Function to compute the hydrogen 2p_z orbital wavefunction in Cartesian coordinates
def psi_2p_z(x, y, z):
    """
    Computes the hydrogen 2p_z orbital wavefunction ψ2pz at a given point (x, y, z) in Cartesian coordinates.

    Parameters:
    x (float): x-coordinate of the point.
    y (float): y-coordinate of the point.
    z (float): z-coordinate of the point.

    Returns:
    float: The value of the 2p_z wavefunction at the given point.
    """
    r = np.sqrt(x**2 + y**2 + z**2)  # radial distance
    if r == 0:
        return 0
    psi_value = (1 / (4 * np.sqrt(2 * np.pi) * a0**(3/2))) * (r / a0) * np.exp(-r / (2 * a0)) * (z / r)
    return psi_value

# Monte Carlo Integration Function
def monte_carlo_overlap_integral(R, N, L):
    # check point 1.1.2-2) Generate random points (x, y, z) within the cube [-L, L]
    x = np.random.uniform(-L, L, N)
    y = np.random.uniform(-L, L, N)
    z = np.random.uniform(-L, L, N)

    # check point 1.1.2-3) Compute the integrand at each random point
    integrand_values = np.array([
        psi_2p_z(x[i], y[i], z[i] + R / 2) * psi_2p_z(x[i], y[i], z[i] - R / 2)
        for i in range(N)
    ])

    # Estimate the average value of the integrand
    average_integrand = np.mean(integrand_values)

    # Compute the volume of the cubic region
    volume = (2 * L) ** 3

    # check point 1.1.2-4) Estimate the overlap integral S(R)
    S_R = volume * average_integrand
    return S_R

from scipy.stats import expon

# integration method #2 - importance sampling
# Define constants
a0 = 1.0  # Bohr radius in atomic units
R = 2.0   # check point 1.1.2-1) Separation distance between the orbitals in atomic units
L = 20.0  # Cube size for integration limits

# check point 1.1.2-5) Perform calculations for various numbers of random points N
N_values = [10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
S_R_values = []

# Set random seed for reproducibility
np.random.seed(42)

# Loop over different N values
for N in N_values:
    S_R = monte_carlo_overlap_integral(R, N, L)
    S_R_values.append(S_R)
    print(f"Estimated S(R) for N_random sampling = {N}: {S_R}")


# integration method #2 - importance sampling

# Define the hydrogen 2p_z orbital wavefunction in Cartesian coordinates
def psi_2p_z(x, y, z, a0=1.0):
    """
    Computes the hydrogen 2p_z orbital wavefunction ψ2pz at a given point (x, y, z) in Cartesian coordinates.

    Parameters:
    x (float): x-coordinate of the point.
    y (float): y-coordinate of the point.
    z (float): z-coordinate of the point.
    a0 (float): Bohr radius in atomic units.

    Returns:
    float: The value of the 2p_z wavefunction at the given point.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    return (1 / (4 * np.sqrt(2 * np.pi) * a0**(3/2))) * (r / a0) * np.exp(-r / (2 * a0)) * (z / r)

# check point 1.1.3-1) Define the importance sampling distribution g(x, y, z)
def g_distribution(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    return np.exp(-r)

# # check point 1.1.3-4) Set the number of points to sample
n_points_list = [100, 1000, 10000, 100000, 1000000, 10000000, 100000000]

# Create a list to store the results
averages = []

# Loop over the number of points to sample
for n_points in n_points_list:
    # check point 1.1.3-2)  Generate random points using exponential distribution for importance sampling
    x = expon.rvs(size=n_points, scale=1)  # Sampling from exponential distribution
    y = expon.rvs(size=n_points, scale=1)
    z = expon.rvs(size=n_points, scale=1)
    
    # Numerator: Overlap between psi_2p_z at (x, y, z + R/2) and (x, y, z - R/2)
    numer = psi_2p_z(x, y, z + 1.0) * psi_2p_z(x, y, z - 1.0)
    
    # Denominator: Importance sampling probability density function g(x, y, z)
    # denom = g_distribution(x, y, z)
    denom = expon.pdf(x) * expon.pdf(y) * expon.pdf(z)

    # # check point 1.1.3-3) Adjust the integrand in your calculation by dividing by g(x, y, z)
    integrand = numer / denom
    
    # Estimate the integral using the average of the integrand
    integral = 8 * np.mean(integrand)
    
    # Store the integral value
    averages.append(integral)

    # Print the estimated integral for the current number of points
    print(f"Estimated integral for N_importance sampling = {n_points}: {integral}")