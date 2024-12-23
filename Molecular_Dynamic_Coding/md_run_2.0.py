import F24_MD_code_Tim_Tan_ver2
import numpy as np
import matplotlib.pyplot as plt
# Modify epsilon before running, whether you want to calculate per particle or per mol basis. 

# Part 1. File Input
particle_lst, dimension = F24_MD_code_Tim_Tan_ver2.load_particles_from_csv('Modified_Molecular_Dynamics_Input_Data.csv')
num_particle = particle_lst.size

# Checkpoint 1
# Uncomment the following code to check to see if it had input successfully
'''
print(f"This input file takes in {num_particle} particles, the dimensionality of input is {dimension}")
print(particle_lst[0])
print(particle_lst[10])
'''
#This incorporate with checkpoint 2
'''
print("Before shift center of mass \n")
for particle in particle_lst:
    print(particle)
'''

# Part 2, Setting up parameters. 
box_size = 10.0
volume = box_size ** dimension
density = num_particle / volume
#print(f"Volume = {volume}, Density = {density}")

# Normalize positions to box-scaled units
particle_lst = F24_MD_code_Tim_Tan_ver2.normalization(particle_lst, box_size)
# Calculate the center of mass; Shift positions so that the center of mass is at the origin
particle_lst = F24_MD_code_Tim_Tan_ver2.shift_com(particle_lst)

# Checkpoint 2
'''
print("After shift \n")
for particle in particle_lst:
    print(particle)
'''

# Part 3, Setting up the simulation.
num_steps = 10000
time_step = 0.0032 
target_temperature = 0.5 
dump_frequency = 10 # Frequency to save positions to the output file

# Part 4, setting up umbrella sampling
k = 1.0
r0 = 1.0
# Part 5, running.

def main(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimensions):
    num_particles = particle_lst.size  # Number of particles
    # Initialize arrays to store simulation results
    kinetic_energy_avg = np.zeros(num_steps)
    potential_energy_avg = np.zeros(num_steps)
    temperatures = np.zeros(num_steps)
    virials = np.zeros(num_steps)
    pressures = np.zeros(num_steps)
    potential_energies = np.zeros(num_particles)
    reaction_coordinates = np.zeros(num_steps)
    # Open output file for saving positions (where we save the file)
    with open('output_R.xyz', 'w') as file:
        # Main simulation loop
        for step in range(num_steps):
            # Apply periodic boundary conditions to keep particles within the box
            particle_lst = F24_MD_code_Tim_Tan_ver2.apply_bc(particle_lst)
            # Update positions using Velocity Verlet (Step 1)
            particle_lst = F24_MD_code_Tim_Tan_ver2.apply_verlet_pos(particle_lst, time_step)

            # Evaluate Reaction Coordinate
            R = F24_MD_code_Tim_Tan_ver2.reaction_coordinate_lookup(particle_lst)
            reaction_coordinates[step] = R

            
            # Calculate temperature and kinetic energy (Step 2)
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan_ver2.calculate_temperature(particle_lst, box_size)
            # Rescale velocities to match the target temperature
            scaling_factor = np.sqrt(target_temperature / temperatures[step])

            particle_lst = F24_MD_code_Tim_Tan_ver2.scale_velocity(particle_lst, scaling_factor, time_step)
            
            # Compute forces, potential energy, and virial (Step 3)
            particle_lst, potential_energy_avg[step], virials[step] = F24_MD_code_Tim_Tan_ver2.compute_forces(particle_lst, potential_energies, box_size, dimension)
            
            # Complete the velocity update (Step 4)
            particle_lst = F24_MD_code_Tim_Tan_ver2.complete_force_update(particle_lst, time_step)

            # Calculate temperature again after velocity update
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan_ver2.calculate_temperature(particle_lst, box_size)
            
            # Calculate pressure
            density = num_particles / (box_size ** dimensions)
            pressures[step] = density * temperatures[step] + virials[step] / (box_size ** dimensions)
            
            
            # Save data to file every `dump_frequency` steps
            if step % dump_frequency == 0:
                # For kJ/mol
                #print(f"Energy {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature {temperatures[step]:.5f}\n")
                # For kJ
                print(f"Energy {kinetic_energy_avg[step] + potential_energy_avg[step]}, Temperature {temperatures[step]:.5f}\n")
                print(f"Current distance between solutes is {R}")
                file.write(f"Step {step}\n")
                file.write(f"Energy(kJ) {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature(K in reduced unit) {temperatures[step]:.5f}\n")
                file.write(f"Kinetic Energy(kJ) {kinetic_energy_avg[step]:.5f}, Potential Energy(kJ) {potential_energy_avg[step]:.5}\n")
                file.write(f"Reaction Coordinate {R}")
                for particle in particle_lst:
                    file.write(f"{particle.atom_type} ")
                    file.write(" ".join(f"{coord * box_size:.5f}" for coord in particle.position))
                    file.write("\n")


    return kinetic_energy_avg, potential_energy_avg, temperatures, pressures, particle_lst, reaction_coordinates


def main_US(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimensions, k, r0):
    num_particles = particle_lst.size  # Number of particles
    # Initialize arrays to store simulation results
    kinetic_energy_avg = np.zeros(num_steps)
    potential_energy_avg = np.zeros(num_steps)
    temperatures = np.zeros(num_steps)
    virials = np.zeros(num_steps)
    pressures = np.zeros(num_steps)
    potential_energies = np.zeros(num_particles)
    reaction_coordinates = np.zeros(num_steps)
    # Open output file for saving positions (where we save the file)
    with open('output_R.xyz', 'w') as file:
        # Main simulation loop
        for step in range(num_steps):
            # Apply periodic boundary conditions to keep particles within the box
            particle_lst = F24_MD_code_Tim_Tan_ver2.apply_bc(particle_lst)
            # Update positions using Velocity Verlet (Step 1)
            particle_lst = F24_MD_code_Tim_Tan_ver2.apply_verlet_pos(particle_lst, time_step)

            # Evaluate Reaction Coordinate
            R = F24_MD_code_Tim_Tan_ver2.reaction_coordinate_lookup(particle_lst)
            reaction_coordinates[step] = R

            
            # Calculate temperature and kinetic energy (Step 2)
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan_ver2.calculate_temperature(particle_lst, box_size)
            # Rescale velocities to match the target temperature
            scaling_factor = np.sqrt(target_temperature / temperatures[step])

            particle_lst = F24_MD_code_Tim_Tan_ver2.scale_velocity(particle_lst, scaling_factor, time_step)
            
            bias_force = F24_MD_code_Tim_Tan_ver2.bias(k, r0, R)
            #print(bias_force)

            # Compute forces, potential energy, and virial (Step 3)
            particle_lst, potential_energy_avg[step], virials[step] = F24_MD_code_Tim_Tan_ver2.compute_forces_US(particle_lst, potential_energies, box_size, dimension, bias_force)
            
            # Complete the velocity update (Step 4)
            particle_lst = F24_MD_code_Tim_Tan_ver2.complete_force_update(particle_lst, time_step)

            # Calculate temperature again after velocity update
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan_ver2.calculate_temperature(particle_lst, box_size)
            
            # Calculate pressure
            density = num_particles / (box_size ** dimensions)
            pressures[step] = density * temperatures[step] + virials[step] / (box_size ** dimensions)
            
            
            # Save data to file every `dump_frequency` steps
            if step % dump_frequency == 0:
                # For kJ/mol
                #print(f"Energy {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature {temperatures[step]:.5f}\n")
                # For kJ
                print(f"Energy {kinetic_energy_avg[step] + potential_energy_avg[step]}, Temperature {temperatures[step]:.5f}\n")
                print(f"Current distance between solutes is {R}")
                file.write(f"Step {step}\n")
                file.write(f"Energy(kJ) {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature(K in reduced unit) {temperatures[step]:.5f}\n")
                file.write(f"Kinetic Energy(kJ) {kinetic_energy_avg[step]:.5f}, Potential Energy(kJ) {potential_energy_avg[step]:.5}\n")
                file.write(f"Reaction Coordinate {R}")
                for particle in particle_lst:
                    file.write(f"{particle.atom_type} ")
                    file.write(" ".join(f"{coord * box_size:.5f}" for coord in particle.position))
                    file.write("\n")


    return kinetic_energy_avg, potential_energy_avg, temperatures, pressures, particle_lst, reaction_coordinates




# MD simulation
'''
kinetic_energy_avg, potential_energy_avg, temeperatures, pressures, particle_lst, reaction_coordinates= main(
    particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimension
)
'''
# Umbrella Sampling
kinetic_energy_avg, potential_energy_avg, temeperatures, pressures, particle_lst, reaction_coordinates= main_US(
    particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimension, k, r0
)

# ------ Plotting ------
energy_avg = kinetic_energy_avg + potential_energy_avg
steps = np.arange(num_steps)

# Energy plotting
def energy_plot_kJ(steps, energy_avg):
    plt.figure(figsize=(8, 6))
    plt.plot(steps, energy_avg, label='Average Energy')
    plt.xlabel('Steps', fontsize=12)
    plt.ylabel('Energy (Average)', fontsize=12)
    plt.title('Energy Average vs. Number of Steps', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.show()
    return

# Define a threshold for outliers
def energy_plot_kJmol(steps, energy_avg):
    threshold = 1e5 # For kJ/mol  
    # Filter energy_avg and steps
    filtered_steps = steps[energy_avg < threshold]
    filtered_energy_avg = energy_avg[energy_avg < threshold]
    # Plot the filtered data
    plt.figure(figsize=(8, 6))
    plt.plot(filtered_steps, filtered_energy_avg, label='Filtered Energy', linewidth=2)
    plt.xlabel('Steps', fontsize=12)
    plt.ylabel('Energy (Average)', fontsize=12)
    plt.title('Filtered Energy Average vs. Number of Steps', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.show()
    return

# Scatter plot of reaction coordinates and energy
def scatter(energy_avg, reaction_coordinates):
    plt.figure(figsize=(8, 6))
    plt.scatter(reaction_coordinates, energy_avg, color='blue', label='Energy vs. Distance')
    plt.xlabel('Bond Distance', fontsize=12)
    plt.ylabel('Energy', fontsize=12)
    plt.title('Energy vs. Bond Distance', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.show()

def bond_distance_plot(steps, reaction_coordinates):
    plt.figure(figsize=(8, 6))
    plt.plot(steps, reaction_coordinates, label='bond_distance')
    plt.xlabel('Steps', fontsize=12)
    plt.ylabel('bond_distance', fontsize=12)
    plt.title('Bond distance vs. Number of Steps under bias harmonic potential', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.show()
    return

energy_plot_kJ(steps, energy_avg)
bond_distance_plot(steps, reaction_coordinates)