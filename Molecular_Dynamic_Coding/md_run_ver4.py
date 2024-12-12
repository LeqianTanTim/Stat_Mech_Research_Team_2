import F24_MD_code_Tim_Tan_ver3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Evaluate time for running

# Modify epsilon before running, whether you want to calculate per particle or per mol basis. 

# Part 1. File Input
particle_lst, dimension = F24_MD_code_Tim_Tan_ver3.load_particles_from_csv('surface_test_md.csv')
num_particle = particle_lst.size

# Checkpoint 1
# Uncomment the following code to check to see if it had input successfully
#print(f"This input file takes in {num_particle} particles, the dimensionality of input is {dimension}")
#print(particle_lst[0])
#print(particle_lst[10])

#This incorporate with checkpoint 2
'''
print("Before shift center of mass \n")
for particle in particle_lst:
    print(particle)
'''

# Part 2, Setting up parameters. 
box_size = 15.0
volume = box_size ** dimension
density = num_particle / volume
#print(f"Volume = {volume}, Density = {density}")

# Normalize positions to box-scaled units
particle_lst = F24_MD_code_Tim_Tan_ver3.normalization(particle_lst, box_size)
# Calculate the center of mass; Shift positions so that the center of mass is at the origin
particle_lst = F24_MD_code_Tim_Tan_ver3.shift_com(particle_lst)



# Checkpoint 2

print("After shift \n")
for particle in particle_lst:
    print(particle)

# velocity prescaling
# Part 3, Setting up the simulation.
num_steps = 5000 # 10 ps
umbrella_sampling_num_steps = 250000 # 500 ps
time_step = 2 # fs
target_temperature = 298 # K
dump_frequency = 500 # dump every ps
abnormal_count = []

# For testing purposes


# Part 4, running.

def main(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimensions):
    
    num_particles = particle_lst.size  # Number of particles
    # Initialize arrays to store simulation results
    kinetic_energy_avg = np.zeros(num_steps)
    potential_energy_avg = np.zeros(num_steps)
    temperatures = np.zeros(num_steps)
    potential_energies = np.zeros(num_particles)
    reaction_coordinates = np.zeros(num_steps)
    position_per_frame = []
    
    # Open output file for saving positions (where we save the file)
    with open('output.xyz', 'w') as file:
        # Main simulation loop
        for step in range(num_steps):
            # Apply periodic boundary conditions to keep particles within the box
            # Update positions using Velocity Verlet (Step 1)
            particle_lst = F24_MD_code_Tim_Tan_ver3.apply_bc(particle_lst)
            particle_lst = F24_MD_code_Tim_Tan_ver3.apply_verlet_pos(particle_lst, time_step)
            
            # Calculate temperature and kinetic energy (Step 2)
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan_ver3.calculate_temperature(particle_lst, box_size)
            
            # Rescale velocities to match the target temperature
            scaling_factor = np.sqrt(target_temperature / temperatures[step])
            
            particle_lst = F24_MD_code_Tim_Tan_ver3.scale_velocity(particle_lst, scaling_factor, time_step)
            
            # Compute forces, potential energy, and virial (Step 3)
            particle_lst, potential_energy_avg[step] = F24_MD_code_Tim_Tan_ver3.compute_forces(particle_lst, potential_energies, box_size, dimension)
           


            # Complete the velocity update (Step 4)
            particle_lst = F24_MD_code_Tim_Tan_ver3.complete_force_update(particle_lst, time_step)

            # Calculate temperature again after velocity update
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan_ver3.calculate_temperature(particle_lst, box_size)
            
            # Calculate reaction coordinate
            R, _ = F24_MD_code_Tim_Tan_ver3.reaction_coordinate_lookup(particle_lst)
            reaction_coordinates[step] = R * box_size
            position_per_frame.append(F24_MD_code_Tim_Tan_ver3.record_position(particle_lst, box_size))

            # Save data to file every `dump_frequency` steps
            if step % dump_frequency == 0:
                # For kJ/mol
                #print(f"Energy {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature {temperatures[step]:.5f}\n")
                # For kJ
                print(f"Kinetic Energy {kinetic_energy_avg[step]} Potential energy {potential_energy_avg[step]} Total Energy {kinetic_energy_avg[step] + potential_energy_avg[step]}, Temperature {temperatures[step]:.5f}\n")
                #for particle in particle_lst:
                    #print(particle)
                print(f"The reaction coordiate that the step {step} is {R * box_size :.5f}")
                particle_lst = F24_MD_code_Tim_Tan_ver3.apply_bc(particle_lst)
                file.write(f"Step {step}\n")
                file.write(f"Energy(kJ) {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature(K in reduced unit) {temperatures[step]:.5f}\n")
                file.write(f"Kinetic Energy(kJ) {kinetic_energy_avg[step]:.5f}, Potential Energy(kJ) {potential_energy_avg[step]:.5f}\n")
                for particle in particle_lst:
                    file.write(f"{particle.atom_type} ")
                    file.write(" ".join(f"{coord * box_size:.5f}" for coord in particle.position))
                    file.write("\n")


    return kinetic_energy_avg, potential_energy_avg, temperatures, particle_lst, reaction_coordinates, position_per_frame




kinetic_energy_avg, potential_energy_avg, temperatures, eq_particle_lst, reaction_coordinates, position_per_frame = main(
    particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimension
)
#print(f"The longest reaction coordinate record in the simulation is {np.max(reaction_coordinates)}")
energy_avg = kinetic_energy_avg + potential_energy_avg
steps = np.arange(num_steps)

# Animation part. Uncomment to access


surface = [0,1,2,3,4,5,6,7,8,9,10]
solute = [11]
fig, ax = plt.subplots()
ax.set_xlim(-7.5, 7.5)
ax.set_ylim(-7.5, 7.5)
scatter = ax.scatter([],[])
highlight_surface = ax.scatter([], [], color='red', s=100, label='surface')
highlight_solute = ax.scatter([], [], color='violet', s=100, label='solute')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))  # Adjust bbox_to_anchor to place legend outside

# Function to initialize the animation
def init():
    scatter.set_offsets(np.empty((0,2)))
    highlight_surface.set_offsets(np.empty((0, 2)))  # Empty 2D array    
    highlight_solute.set_offsets(np.empty((0, 2)))
    return scatter, highlight_surface, highlight_solute

# Function to update the animation for each frame
def update(frame):
    # Update all points
    scatter.set_offsets(frame)
    
    # Highlight specific particles
    marked_surface = frame[surface]  # Get positions of marked particles
    marked_solute = frame[solute]
    highlight_surface.set_offsets(marked_surface)
    highlight_solute.set_offsets(marked_solute)
    
    return scatter, highlight_surface, highlight_solute

# Create the animation
ani = FuncAnimation(fig, update, frames=position_per_frame, init_func=init, blit=True, interval=500)

# Show the animation
plt.show()

writer = FFMpegWriter(fps=30, codec='h264', extra_args=['-pix_fmt', 'yuv420p'])

# Save the animation using FFMpegWriter
ani.save('animation_surface.mp4', writer=writer)

