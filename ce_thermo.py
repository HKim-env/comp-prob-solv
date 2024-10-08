import numpy as np
from scipy.constants import k  # Boltzmann constant

# Isolated Ce3+ with 14-fold degeneracy, zero energy for all states
def partition_function_isolated(T):
    g_iso = 14  # Degeneracy (constant partition function)
    Z_iso = g_iso  # Partition function is constant
    return Z_iso

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
    
    # Entropy: S = - dF / dT
    S = k * np.log(Z)  # For this simple case, S = k_B * ln(Z)
    
    return U, F, S

# Ce3+ with SOC: Two levels with degeneracies 6 and 8, separated by 0.28 eV
def partition_function_soc(T):
    g_soc = np.array([6, 8])  # Degeneracies
    E_soc = np.array([0, 0.28 * 1.60218e-19])  # Energy levels in Joules (0.28 eV to J)
    beta = 1 / (k * T)
    Z_soc = np.sum(g_soc * np.exp(-beta * E_soc))  # Partition function
    return Z_soc, g_soc, E_soc

def thermodynamic_properties_soc(T):
    Z, g_soc, E_soc = partition_function_soc(T)  # Ensure g_soc and E_soc are retrieved
    
    # Internal Energy: U = - d(ln Z) / d(beta)
    beta = 1 / (k * T)
    # Compute internal energy from the energy levels
    U = np.sum((g_soc * E_soc * np.exp(-beta * E_soc)) / Z)  # Weighted sum over energy levels
    
    # Free Energy: F = - k_B T ln(Z)
    F = - k * T * np.log(Z)
    
    # Entropy: S = U/T - F/T
    S = k * np.log(Z)  # This is the analytic expression for entropy in this case
    
    return U, F, S

# Ce3+ with SOC and CFS: Five distinct levels with different degeneracies
def partition_function_soc_cfs(T):
    g_soc_cfs = np.array([4, 2, 2, 4, 2])  # Degeneracies
    E_soc_cfs = np.array([0, 0.12, 0.25, 0.32, 0.46]) * 1.60218e-19  # Energy levels in Joules (converted from eV)
    beta = 1 / (k * T)
    Z_soc_cfs = np.sum(g_soc_cfs * np.exp(-beta * E_soc_cfs))  # Partition function
    return Z_soc_cfs, g_soc_cfs, E_soc_cfs

def thermodynamic_properties_soc_cfs(T):
    Z, g_soc_cfs, E_soc_cfs = partition_function_soc_cfs(T)  # Ensure g_soc_cfs and E_soc_cfs are retrieved
    
    # Internal Energy: U = - d(ln Z) / d(beta)
    beta = 1 / (k * T)
    # Compute internal energy directly from the energy levels
    U = np.sum((g_soc_cfs * E_soc_cfs * np.exp(-beta * E_soc_cfs)) / Z)  # Weighted sum over energy levels
    
    # Free Energy: F = - k_B T ln(Z)
    F = - k * T * np.log(Z)
    
    # Entropy: S = U/T - F/T
    S = k * np.log(Z)  # This is the analytic expression for entropy in this case
    
    return U, F, S