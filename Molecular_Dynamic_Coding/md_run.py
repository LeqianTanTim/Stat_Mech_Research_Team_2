import F24_MD_code_Tim_Tan
import numpy as np
import matplotlib.pyplot as plt
from IPython import display

# Part 1. File Input
particle_lst, dimension = F24_MD_code_Tim_Tan.load_particles_from_csv('Modified_Molecular_Dynamics_Input_Data.csv')
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
print(f"Volume = {volume}, Density = {density}")

# Normalize positions to box-scaled units
particle_lst = F24_MD_code_Tim_Tan.normalization(particle_lst, box_size)
# Calculate the center of mass; Shift positions so that the center of mass is at the origin
particle_lst = F24_MD_code_Tim_Tan.shift_com(particle_lst)

# Checkpoint 2
'''
print("After shift \n")
for particle in particle_lst:
    print(particle)
'''

# Part 3, Setting up the simulation.
num_steps = 1000 
time_step = 0.0032 
target_temperature = 0.5 
dump_frequency = 100 # Frequency to save positions to the output file

# Part 4, running.

def main(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimensions):
    """
    Main function for molecular dynamics simulation using the Velocity Verlet algorithm.

    Parameters:
        positions: np.ndarray
            Initial positions of particles (num_particles x dimensions).
        num_steps: int
            Number of simulation steps.
        time_step: float
            Time step for integration in reduced units.
        target_temperature: float
            Target temperature for the system in reduced units.
        dump_frequency: int
            Frequency to save positions to the output file.
        epsilon: float
            Lennard-Jones energy parameter.
        box_size: float
            Length of the cubic simulation box.
        dimensions: int
            Dimensionality of the system (2D or 3D).

    Returns:
        Tuple containing arrays of average kinetic energy, potential energy, temperature,
        pressure, and final positions.
    """
    num_particles = particle_lst.size  # Number of particles

    # Initialize arrays to store simulation results
    kinetic_energy_avg = np.zeros(num_steps)
    potential_energy_avg = np.zeros(num_steps)
    temperatures = np.zeros(num_steps)
    virials = np.zeros(num_steps)
    pressures = np.zeros(num_steps)
    potential_energies = np.zeros(num_particles)

    # Initialize velocities and accelerations with random values
    # Should we also initialize these two quantities with random values in real practice?
    # Prerun
    velocities = np.random.randn(num_particles, dimensions) - 0.5
    accelerations = np.random.randn(num_particles, dimensions) - 0.5

    # Open output file for saving positions (where we save the file)
    with open('trajectory.xyz', 'w') as file:
        # Main simulation loop
        for step in range(num_steps):
            # Apply periodic boundary conditions to keep particles within the box
            particle_lst = F24_MD_code_Tim_Tan.apply_bc(particle_lst)

            # Update positions using Velocity Verlet (Step 1)
            particle_lst = F24_MD_code_Tim_Tan.apply_verlet_pos(particle_lst, time_step)

            # Calculate temperature and kinetic energy (Step 2)
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan.calculate_temperature(particle_lst, box_size)

            # Rescale velocities to match the target temperature
            scaling_factor = np.sqrt(target_temperature / temperatures[step])
            
            particle_lst = F24_MD_code_Tim_Tan.scale_velocity(particle_lst, scaling_factor, time_step)

            # Compute forces, potential energy, and virial (Step 3)
            particle_lst, potential_energy_avg[step], virials[step] = F24_MD_code_Tim_Tan.compute_forces(
                particle_lst, potential_energies, box_size, dimension)

            # Complete the velocity update (Step 4)
            particle_lst = F24_MD_code_Tim_Tan.complete_force_update(particle_lst, time_step)
            
            # Calculate temperature again after velocity update
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan(particle_lst, box_size)

            # Calculate pressure
            density = num_particles / (box_size ** dimensions)
            pressures[step] = density * temperatures[step] + virials[step] / (box_size ** dimensions)

            # Save data to file every `dump_frequency` steps
            if step % dump_frequency == 0:
                file.write(f"{num_particles}\n")
                file.write(f"Energy {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature {temperatures[step]:.5f}\n")
                for particle in range(num_particles):
                    file.write("X ")
                    file.write(" ".join(f"{positions[particle, dim] * box_size:.5f}" for dim in range(dimensions)))
                    file.write("\n")

                # Visualization for 2D systems
                if dimensions == 2:
                    plt.cla()
                    plt.xlim(-0.5 * box_size, 0.5 * box_size)
                    plt.ylim(-0.5 * box_size, 0.5 * box_size)
                    for particle in range(num_particles):
                        plt.plot(positions[particle, 0] * box_size, positions[particle, 1] * box_size, 'o', markersize=4)
                    display.clear_output(wait=True)
                    display.display(plt.gcf())

    return kinetic_energy_avg, potential_energy_avg, temperatures, pressures, particle_lst