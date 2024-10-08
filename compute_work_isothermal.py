import numpy as np
from scipy.integrate import trapezoid

def compute_work_isothermal(V_i, V_f, n, R, T, num_points=1000):
    # Define the volume array for integration
    V = np.linspace(V_i, V_f, num_points)
    
    # Define the pressure as a function of volume using the ideal gas law
    P = (n * R * T) / V
    
    # Compute the work using trapezoidal integration
    work = -trapezoid(P, V)
    
    return work