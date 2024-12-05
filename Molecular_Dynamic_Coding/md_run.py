import F24_MD_code_Tim_Tan
import numpy as np
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