import F24_MD_code_Tim_Tan_ver2
import numpy as np
import matplotlib.pyplot as plt
# Modify epsilon before running, whether you want to calculate per particle or per mol basis. 

# Part 1. File Input
particle_lst, dimension = F24_MD_code_Tim_Tan_ver2.load_particles_from_csv('fixed_test.csv')
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
k = 1.25
r0_lst = np.linspace(2, 6, 20)
r0_norm = r0_lst / box_size
# Part 5, running.

def md_simulation(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimensions):
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
    with open('equilibrium.xyz', 'w') as file:
        # Main simulation loop
        for step in range(num_steps):
            # Apply periodic boundary conditions to keep particles within the box
            particle_lst = F24_MD_code_Tim_Tan_ver2.apply_bc(particle_lst)
            # Update positions using Velocity Verlet (Step 1)
            particle_lst = F24_MD_code_Tim_Tan_ver2.apply_verlet_pos(particle_lst, time_step)

            # Evaluate Reaction Coordinate
            R, relative_position = F24_MD_code_Tim_Tan_ver2.reaction_coordinate_lookup(particle_lst)
            reaction_coordinates[step] = R * box_size

            
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
                print(f"Current distance between solutes is {R * box_size:.5f}")
                file.write(f"Step {step}\n")
                file.write(f"Energy(kJ) {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature(K in reduced unit) {temperatures[step]:.5f}\n")
                file.write(f"Kinetic Energy(kJ) {kinetic_energy_avg[step]:.5f}, Potential Energy(kJ) {potential_energy_avg[step]:.5}\n")
                file.write(f"Reaction Coordinate {R * box_size: .5f}")
                for particle in particle_lst:
                    file.write(f"{particle.atom_type} ")
                    file.write(" ".join(f"{coord * box_size:.5f}" for coord in particle.position))
                    file.write("\n")


    return kinetic_energy_avg, potential_energy_avg, temperatures, pressures, particle_lst, reaction_coordinates


def bias_md_simulation(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimensions, k, r0):
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
    with open(f'output_{k}_{r0}.xyz', 'w') as file:
        # Main simulation loop
        for step in range(num_steps):
            # Apply periodic boundary conditions to keep particles within the box
            particle_lst = F24_MD_code_Tim_Tan_ver2.apply_bc(particle_lst)
            # Update positions using Velocity Verlet (Step 1)
            particle_lst = F24_MD_code_Tim_Tan_ver2.apply_verlet_pos(particle_lst, time_step)

            # Evaluate Reaction Coordinate
            R, relative_position = F24_MD_code_Tim_Tan_ver2.reaction_coordinate_lookup(particle_lst)
            bias_force = F24_MD_code_Tim_Tan_ver2.bias(k, r0, R, relative_position)
            
            if R > 10.0:
                solute_1, solute_2 = F24_MD_code_Tim_Tan_ver2.reaction_coordinate_verify(particle_lst)
                print(f"Error occurred at iteration {step} for k {k} at a set target distance {r0 * box_size:.5f}")
                print(f"solute is is currently located at {solute_1.position} and {solute_2.position}")
                print(f"Acceleration is {solute_1.acceleration} and {solute_2.acceleration}")
                print(f"Energy {kinetic_energy_avg[step] + potential_energy_avg[step]}, Temperature {temperatures[step]:.5f}\n")
                break
            
            reaction_coordinates[step] = R * box_size
            # Calculate temperature and kinetic energy (Step 2)
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan_ver2.calculate_temperature(particle_lst, box_size)
            # Rescale velocities to match the target temperature
            scaling_factor = np.sqrt(target_temperature / temperatures[step])

            particle_lst = F24_MD_code_Tim_Tan_ver2.scale_velocity(particle_lst, scaling_factor, time_step)
    

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
                print(f"Current distance between solutes is {R * box_size:.5f}")
                s_vel, s_acc = F24_MD_code_Tim_Tan_ver2.monitor_the_surface(particle_lst)
                print(f"The surface molecule currently has velocity {s_vel} and acceleration {s_acc}")
                file.write(f"Step {step}\n")
                file.write(f"Energy(kJ) {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature(K in reduced unit) {temperatures[step]:.5f}\n")
                file.write(f"Kinetic Energy(kJ) {kinetic_energy_avg[step]:.5f}, Potential Energy(kJ) {potential_energy_avg[step]:.5}\n")
                file.write(f"Reaction Coordinate {R * box_size:.5f}")
                
                for particle in particle_lst:
                    file.write(f"{particle.atom_type} ")
                    file.write(" ".join(f"{coord * box_size:.5f}" for coord in particle.position))
                    file.write("\n")


    return kinetic_energy_avg, potential_energy_avg, temperatures, pressures, reaction_coordinates




# MD simulation
'''
kinetic_energy_avg, potential_energy_avg, temeperatures, pressures, particle_lst, reaction_coordinates= md_simulation(
    particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimension
)
'''
# Umbrella Sampling

def umbrella_sampling(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimension, k, r0_lst):
    # Initialize the figure and axis
    plt.figure()
    plt.xlabel("Reaction coordinate")
    plt.ylabel("Frequency")
    plt.title("Accumulating Histograms")
    plt.grid(True)
    # run an equilibrium MD to get a better intial pos, vel and acc
    eq_kinetic_energy_avg, eq_potential_energy_avg, eq_temeperatures, eq_pressures, eq_particle_lst, eq_reaction_coordinates= md_simulation(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimension,)
    num_bins = 40
    windows_histogram = []
    for r0 in r0_lst:
        kinetic_energy_avg, potential_energy_avg, temeperatures, pressures, reaction_coordinates= bias_md_simulation(
            eq_particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimension, k, r0
        )

        hist, bin_edges = np.histogram(reaction_coordinates, bins= num_bins, range=(2.0, 6.0))
        hist_norm = hist / np.sum(hist)
        windows_histogram.append(hist_norm)
        bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])  # to plot using bin centers
    
    for w, h in enumerate(windows_histogram):
        plt.plot(bin_centers, h, label=f"Window {w}")
    plt.xlabel("Reaction Coordinate (Å)")
    plt.ylabel("Probability")
    plt.title("Umbrella Window Histograms")
    #plt.legend()
    plt.show()

umbrella_sampling(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimension, k, r0_norm)