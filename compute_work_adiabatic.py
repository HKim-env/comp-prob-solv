import numpy as np
from scipy.integrate import trapezoid

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

# Example usage
V_i = 0.1  # m^3
V_f = 0.3  # m^3 (example, can be changed)
P_i = 100000  # Pa (example initial pressure)
gamma = 1.4   # Adiabatic index

work_adiabatic = compute_work_adiabatic(V_i, V_f, P_i, gamma)
print(f"Work done in adiabatic process: {work_adiabatic} J")