import F24_MD_code_Tim_Tan
import numpy as np
# Part 1. File Input
particle_lst, dimension = F24_MD_code_Tim_Tan.load_particles_from_csv('./Molecular_Dynamic_Coding/Modified_Molecular_Dynamics_Input_Data.csv')
num_particle = particle_lst.size

# Checkpoint 1
# Uncomment the following code to check to see if it had input successfully
#print(f"This input file takes in {num_particle} particles, the dimensionality of input is {dimension}")
#print(particle_lst[0])
#print(particle_lst[10])

# Part 2, Setting up parameters. 
box_size = 10.0
volume = box_size ** dimension
density = num_particle / volume
print(f"Volume = {volume}, Density = {density}")
positions = np.zeros([num_particle, dimension])

# Load positions from file; modify based on need.
positions = np.genfromtxt('output.dat', skip_header=1)
# Normalize positions to box-scaled units
positions = positions[:, :dimensions] / box_size

# Calculate the center of mass
center_of_mass = np.sum(positions, axis=0) / num_particles

# Shift positions so that the center of mass is at the origin
for dim in range(dimensions):
    positions[:, dim] -= center_of_mass[dim]
