import numpy as np
from scipy.stats import expon
from scipy.stats import norm

# functions for ψ1s(x, y, z), 
def psi_1s(x, y, z, Z=1, a0=1.0):
    """
    Computes the hydrogen 1s orbital wavefunction ψ1s at a given point (x, y, z) in Cartesian coordinates.

    Parameters:
    x (float): x-coordinate of the point.
    y (float): y-coordinate of the point.
    z (float): z-coordinate of the point.
    Z (int): Atomic number (default is 1 for hydrogen).
    a0 (float): Bohr radius in atomic units.

    Returns:
    float: The value of the 1s wavefunction at the given point.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    return (Z**(3/2) / np.sqrt(np.pi * a0**3)) * np.exp(-Z * r / a0)

# functions for ∇2ψ1s(x, y, z)
def laplacian_psi_1s(x, y, z, Z=1, a0=1.0):
    """
    Computes the Laplacian of the hydrogen 1s orbital wavefunction at a given point (x, y, z) in Cartesian coordinates.

    Parameters:
    x (float): x-coordinate of the point.
    y (float): y-coordinate of the point.
    z (float): z-coordinate of the point.
    Z (int): Atomic number (default is 1 for hydrogen).
    a0 (float): Bohr radius in atomic units.

    Returns:
    float: The value of the Laplacian of the 1s wavefunction at the given point.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    if r == 0:
        return 0
    # Compute the expression
    numerator = (-2 * a0 + r) * np.exp(-r / a0)
    denominator = np.sqrt(np.pi * a0**3) * r * a0**2
    # Final expression
    laplacian = numerator / denominator
    return laplacian

# 1. Compute the Diagonal Kinetic Energy Matrix Element Using Random Sampling
# Define the Bohr radius
a0 = 1.0

# Monte Carlo Integration for Kii (diagonal kinetic energy matrix element)
def monte_carlo_Kii(N, L):
    """
    Computes the diagonal kinetic energy matrix element Kii using Monte Carlo integration.

    Parameters:
    N (int): Number of random points to sample.
    L (float): Half-length of the cubic region for integration.

    Returns:
    float: The estimated value of the diagonal kinetic energy matrix element Kii.
    """
    # check point 2.1.2-2)Generate random points in the cubic region [-L, L]
    x = np.random.uniform(-L, L, N)
    y = np.random.uniform(-L, L, N)
    z = np.random.uniform(-L, L, N)
    
    # Compute the integrand -1/2 * psi_1s * Laplacian(psi_1s)
    integrand = np.array([
        (-0.5) * psi_1s(x[i], y[i], z[i]) * laplacian_psi_1s(x[i], y[i], z[i])
        for i in range(N)
    ])
    
    # check point 2.1.2-3) Compute the average value of the integrand
    integrand_avg = np.mean(integrand)
    
    # Compute the volume of the cubic region
    volume = (2 * L) ** 3
    
    # check point 2.1.2-4) Estimate Kii
    Kii = volume * integrand_avg
    return Kii

# Parameters
L = 7.0  # Cube size for integration limits

# check point 2.1.2-5) Number of samples for Monte Carlo integration
N_values = [10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8]

# check point 2.1.2-6) Array to store the results
Kii_values = []

# Set random seed for reproducibility
np.random.seed(42)

# Perform the Monte Carlo integration for each N
for N in N_values:
    Kii = monte_carlo_Kii(N, L)
    Kii_values.append(Kii)
    print(f"N = {N}: Kii = {Kii} for diagonal kinetic energy matrix element (random sampling)")

# 2. Compute the Diagonal Kinetic Energy Matrix Element Using importance Sampling
# Monte Carlo Integration for Kii (diagonal kinetic energy matrix element)
def monte_carlo_Kii(N, L):
    """
    Computes the diagonal kinetic energy matrix element Kii using Monte Carlo integration.

    Parameters:
    N (int): Number of random points to sample.
    L (float): Half-length of the cubic region for integration.

    Returns:
    float: The estimated value of the diagonal kinetic energy matrix element Kii.
    """
    # check point 2.1.2-2)Generate random points in the cubic region [-L, L]
    x = np.random.uniform(-L, L, N)
    y = np.random.uniform(-L, L, N)
    z = np.random.uniform(-L, L, N)
    
    # Compute the integrand -1/2 * psi_1s * Laplacian(psi_1s)
    integrand = np.array([
        (-0.5) * psi_1s(x[i], y[i], z[i]) * laplacian_psi_1s(x[i], y[i], z[i])
        for i in range(N)
    ])
    
    # check point 2.1.2-3) Compute the average value of the integrand
    integrand_avg = np.mean(integrand)
    
    # Compute the volume of the cubic region
    volume = (2 * L) ** 3
    
    # check point 2.1.2-4) Estimate Kii
    Kii = volume * integrand_avg
    return Kii

# Parameters
L = 7.0  # Cube size for integration limits

# check point 2.1.2-5) Number of samples for Monte Carlo integration
N_values = [10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8]

# check point 2.1.2-6) Array to store the results
Kii_values = []

# Set random seed for reproducibility
np.random.seed(42)

# Perform the Monte Carlo integration for each N
for N in N_values:
    Kii = monte_carlo_Kii(N, L)
    Kii_values.append(Kii)
    print(f"N = {N}: Kii = {Kii} for diagonal kinetic energy matrix element (random sampling)")

# 3. Compute the Off-Diagonal Kinetic Energy Matrix Element Using Random Sampling
# Monte Carlo Integration for Kij (off-diagonal kinetic energy matrix element)
def monte_carlo_Kij(N, L, Rz):
    """
    Computes the off-diagonal kinetic energy matrix element Kij using Monte Carlo integration.

    Parameters:
    N (int): Number of random points to sample.
    L (float): Half-length of the cubic region for integration.
    Rz (float): Separation distance between the two hydrogen atoms along the z-axis.

    Returns:
    float: The estimated value of the off-diagonal kinetic energy matrix element Kij.
    """
    # # check point 2.1.4-2) Generate random points in the cubic region [-L, L]
    x = np.random.uniform(-L, L, N)
    y = np.random.uniform(-L, L, N)
    z = np.random.uniform(-L, L, N)
    
    # Compute the integrand -1/2 * psi_1s(x, y, z + Rz/2) * Laplacian(psi_1s(x, y, z - Rz/2))
    integrand = np.array([
        (-0.5) * psi_1s(x[i], y[i], z[i] + Rz / 2) * laplacian_psi_1s(x[i], y[i], z[i] - Rz / 2)
        for i in range(N)
    ])
    
    # # check point 2.1.4-3) Compute the average value of the integrand
    integrand_avg = np.mean(integrand)
    
    # Compute the volume of the cubic region
    volume = (2 * L) ** 3
    
    # # check point 2.1.4-4) Estimate Kij
    Kij = volume * integrand_avg
    return Kij

# Parameters
L = 7.0  # Cube size for integration limits
Rz = 1.4 # Separation distance between the two hydrogen atoms

# # check point 2.1.4-5) Number of samples for Monte Carlo integration
N_values = [10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8]

# check point 2.1.4-6) Record Kij
Kij_values = []

# Set random seed for reproducibility
np.random.seed(42)

# Perform the Monte Carlo integration for each N
for N in N_values:
    Kij = monte_carlo_Kij(N, L, Rz)
    Kij_values.append(Kij)
    print(f"N = {N}: Kij = {Kij} for off-diagonal kinetic energy matrix element (importance sampling)")

#4. Compute the Off-Diagonal Kinetic Energy Matrix Element Using importance Sampling
# Monte Carlo Integration for Kij with Importance Sampling using Gaussian Sampling and norm.pdf for the denominator
def monte_carlo_Kij_importance(N, Z=1, Rz=1.4, a0=1.0):
    """
    Computes the diagonal kinetic energy matrix element Kii using Monte Carlo integration with importance sampling.

    Parameters:
    N (int): Number of random points to sample.
    Z (int): Atomic number (default is 1 for hydrogen).
    a0 (float): Bohr radius in atomic units.

    Returns:
    float: The estimated value of the diagonal kinetic energy matrix element Kii.
    """
    # Use Gaussian distribution to sample x, y, z values with mean shifted for z
    x = np.random.normal(0, 1, N)  # Gaussian sampling for x
    y = np.random.normal(0, 1, N)  # Gaussian sampling for y
    z = np.random.normal(Rz / 2, 1, N)  # check point 2.1.5-2) Since the orbitals are centered at different locations, consider a Gaussian distribution centered between them
    
    # Numerator: Overlap between psi_1s centered at ±Rz/2
    numer = np.array([
        (-0.5) * psi_1s(x[i], y[i], z[i], Z, a0) * laplacian_psi_1s(x[i], y[i], z[i], Z, a0) 
        for i in range(N)
    ])
    
    # Denominator: Using norm.pdf for the Gaussian distribution for x, y, z
    denom = norm.pdf(x, 0, 1) * norm.pdf(y, 0, 1) * norm.pdf(z, Rz / 2, 1)  # Gaussian PDF for x, y, z
    
    # check point 2.1.5-4) adjust the integrand
    integrand = numer / denom

    # Estimate Kij using the average integrand and multiply by 8 to account for all octants
    Kij = np.mean(integrand)
    
    return Kij

# check point 2.1.5-5) use the same N values as before
N_values = [10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8]  # Number of random samples
Rz = 1.4  # Center for z

# check point 2.1.5-6) Record the estimated values of Kij
Kij_values_importance = []

# Set random seed for reproducibility
np.random.seed(42)

# Perform Monte Carlo integration for each N using importance sampling
for N in N_values:
    Kij_importance = monte_carlo_Kij_importance(N, Rz=Rz)
    Kij_values_importance.append(Kij_importance)
    print(f"N = {N}: Kij (Importance Sampling with norm.pdf) = {Kij_importance}")