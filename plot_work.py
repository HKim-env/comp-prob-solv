import matplotlib.pyplot as plt
import numpy as np
from python_code import compute_work_adiabatic
from python_code import compute_work_isothermal

# 1-1-3 Checking point 1) Plot the work done in both processes as a function of the final volume Vf for values of Vf between Vi and 3Vi
# 1-1-3 Checking point 2) Assume the following parameters for the ideal gas

# Given parameters
n = 1      # mol
R = 8.314  # J/mol-K
T = 300    # K
V_i = 0.1  # m^3
gamma = 1.4
P_i = (n * R * T) / V_i  # From ideal gas law

# Define the range of final volumes
V_f_values = np.linspace(V_i, 3 * V_i, 100)

# Compute work for each final volume
work_isothermal = [compute_work_isothermal(V_i, V_f, n, R, T) for V_f in V_f_values]
work_adiabatic = [compute_work_adiabatic(V_i, V_f, P_i, gamma) for V_f in V_f_values]

# Plotting
plt.plot(V_f_values, work_isothermal, label='Isothermal Process')
plt.plot(V_f_values, work_adiabatic, label='Adiabatic Process')

plt.xlabel('Final Volume $V_f$ (m^3)')
plt.ylabel('Work (J)')
plt.title('Work Done in Isothermal and Adiabatic Processes')
plt.legend()
plt.grid(True)
plt.show()