import matplotlib.pyplot as plt
import numpy as np
from python_code import internal_energy
from python_code import heat_capacity
from python_code import T_values

# # 3-1-3 Checking point 1) Plot CV as a function of temperature, and determine the temperature at which CV reaches its maximum

# Plotting the heat capacity and adding the dissociation temperature
def plot_heat_capacity_with_dissociation(C_V_values, dissociation_temp, max_CV):
    plt.plot(T_values, C_V_values, label='Heat Capacity')
    
    # Mark the maximum heat capacity point (convert to lists to avoid errors in plt.scatter)
    plt.scatter([dissociation_temp], [max_CV], color='red', zorder=5)
    plt.text(dissociation_temp, max_CV, f'  Max C_V = {max_CV:.2e}\n  T = {dissociation_temp:.2f} K', 
             verticalalignment='bottom', horizontalalignment='right', color='red')
    
    plt.title('Heat Capacity vs. Temperature with Maximum C_V')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Heat Capacity C_V (J/K)')
    plt.grid(True)
    plt.legend()
    plt.show()

# Calculate the dissociation temperature
def find_dissociation_temperature(C_V_values):
    max_index = np.argmax(C_V_values)
    dissociation_temp = T_values[max_index]
    max_CV = C_V_values[max_index] 
    
    # # 3-1-3 Checking point 2) This maximum corresponds to the atomization temperature of the LJ dimer

    # Print the maximum heat capacity and corresponding temperature only once
    print(f"Maximum Heat Capacity: {max_CV:.2e} J/K")
    print(f"Dissociation Temperature (Maximum C_V): {dissociation_temp:.2f} K")
    
    return dissociation_temp, max_CV

# Main execution
C_V_values, U_values = heat_capacity(T_values)  # Calculate C_V and internal energy
dissociation_temp, max_CV = find_dissociation_temperature(C_V_values)  # Find and print dissociation temperature

# Plot the heat capacity with the maximum heat capacity annotated
plot_heat_capacity_with_dissociation(C_V_values, dissociation_temp, max_CV)
