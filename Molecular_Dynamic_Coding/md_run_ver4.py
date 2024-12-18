import F24_MD_code_Tim_Tan_ver3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import time

# Evaluate time for running
start_time = time.time()
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
umbrella_sampling_num_steps = 50000 # 500 ps
time_step = 2 # fs
target_temperature = 298 # K
dump_frequency = 500 # dump every ps
abnormal_count = []

# Umbrella sampling parameter
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




##### Umbrella Sampling
def umbrella_sampling(particle_lst, num_steps, time_step, target_temperature, dump_frequency, box_size, dimensions, k, r0, cnt):
    num_particles = particle_lst.size  # Number of particles
    # Initialize arrays to store simulation results
    kinetic_energy_avg = np.zeros(num_steps)
    potential_energy_avg = np.zeros(num_steps)
    temperatures = np.zeros(num_steps)
    potential_energies = np.zeros(num_particles)
    reaction_coordinates = np.zeros(num_steps)
    position_per_frame = []
    
    # Open output file for saving positions (where we save the file)
    with open(f'umbrella_sampling_{cnt}.xyz', 'w') as file:
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
            # Determine the current R and relative position here for bias force evaluation
            R, relative_position = F24_MD_code_Tim_Tan_ver3.reaction_coordinate_lookup(particle_lst)
            R = R * box_size # Scale it correctly (in angstrom)
            reaction_coordinates[step] = R
            # k is in unit of N/m 
            relative_position = relative_position * 1e-10 * box_size
            bias_force = -k * ((R - r0)/R) * relative_position 
            particle_lst, potential_energy_avg[step] = F24_MD_code_Tim_Tan_ver3.compute_forces_umbrella_sampling(particle_lst, potential_energies, box_size, dimensions, bias_force)
           


            # Complete the velocity update (Step 4)
            particle_lst = F24_MD_code_Tim_Tan_ver3.complete_force_update(particle_lst, time_step)

            # Calculate temperature again after velocity update
            kinetic_energy_avg[step], temperatures[step] = F24_MD_code_Tim_Tan_ver3.calculate_temperature(particle_lst, box_size)
            
            # Calculate reaction coordinate
            
            position_per_frame.append(F24_MD_code_Tim_Tan_ver3.record_position(particle_lst, box_size))

            # Save data to file every `dump_frequency` steps
            if step % dump_frequency == 0:
                # For kJ/mol
                #print(f"Energy {kinetic_energy_avg[step] + potential_energy_avg[step]:.5f}, Temperature {temperatures[step]:.5f}\n")
                # For kJ
                print(f"Umbrella sampling {step} window {cnt}")
                print(f"Kinetic Energy {kinetic_energy_avg[step]} Potential energy {potential_energy_avg[step]} Total Energy {kinetic_energy_avg[step] + potential_energy_avg[step]}, Temperature {temperatures[step]:.5f}\n")
                #for particle in particle_lst:
                    #print(particle)
                print(f"The bias reaction coordiate that the step {step} is {R :.5f}")
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

'''
filter = mask = (temperatures > 400) | (temperatures < 100)
filter_reaction_coordinates = reaction_coordinates[filter]
hist, bin_edges = np.histogram(filter_reaction_coordinates, bins=50, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.figure(figsize=(8,6))
plt.plot(bin_centers, hist, drawstyle='steps-mid', color='blue')
plt.xlabel('Reaction Coordinate')
plt.ylabel('Probability Density')
plt.title('reaction_coordinates')
plt.grid(True)
plt.tight_layout()
plt.savefig('reaction_coordinate_probability.png')
plt.show()

'''
k_lst = [500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500]
bins = 20
windows_histogram = []
bin_centers_list = []  # To store bin centers for each window
cnt = 0

# custom cutoff
r0_lst = [1.50, 2.39, 3.28, 4.17, 5.06, 5.94, 6.4, 6.83, 7.2, 7.72, 8.61, 9.5]

for idx, (r0, k) in enumerate(zip(r0_lst, k_lst)):
    # Debug: Print original r_range to inspect its structure
    print(f"\nProcessing Window {idx+1}: r0 = {r0}, k = {k}")
    r_range = (0.0, 11.0)
    # Ensure r_range is a tuple with two elements (min, max)
    if isinstance(r_range, (list, np.ndarray)):
        if len(r_range) >= 2:
            r_min, r_max = r_range[:2]
            try:
                r_min = float(r_min)
                r_max = float(r_max)
            except (ValueError, TypeError) as e:
                print(f"Error converting r_min and r_max to float for window {idx+1}: {e}")
                continue  # Skip this window
            r_range_fixed = (r_min, r_max)
            print(f"Fixed r_range: {r_range_fixed}, type = {type(r_range_fixed)}")
        else:
            print(f"r_range for window {idx+1} does not have enough elements: {r_range}")
            continue  # Skip this window
    elif isinstance(r_range, tuple) and len(r_range) == 2:
        r_min, r_max = r_range
        try:
            r_min = float(r_min)
            r_max = float(r_max)
        except (ValueError, TypeError) as e:
            print(f"Error converting r_min and r_max to float for window {idx+1}: {e}")
            continue  # Skip this window
        r_range_fixed = (r_min, r_max)
        print(f"Fixed r_range: {r_range_fixed}, type = {type(r_range_fixed)}")
    else:
        print(f"r_range for window {idx+1} is not in a recognized format: {r_range}, type = {type(r_range)}")
        continue  # Skip this window

    # Perform umbrella sampling (ensure that the function and variables are defined)
    cnt += 1
    try:
        us_kinetic_energy_avg, us_potential_energy_avg, us_temperatures, us_particle_lst, us_reaction_coordinates, us_position_per_frame = umbrella_sampling(
            eq_particle_lst, umbrella_sampling_num_steps, time_step, target_temperature,
            dump_frequency, box_size, dimension, k, r0, cnt
        )
    except Exception as e:
        print(f"Umbrella sampling failed for window {idx+1} (r0={r0}): {e}")
        continue  # Skip to the next window

    # Filter based on temperature (uncomment if filtering is desired)
    temp_filter = (us_temperatures > 100) & (us_temperatures < 400)  # Corrected logical condition
    filter_reaction_coordinates = us_reaction_coordinates[temp_filter]

    # Ensure there are data points after filtering
    if filter_reaction_coordinates.size == 0:
        print(f"No data points in window {idx+1} after temperature filtering. Skipping histogram.")
        continue

    # Compute histogram within the specified range
    try:
        hist, bin_edges = np.histogram(filter_reaction_coordinates, bins=bins, range=r_range_fixed)
        print(f"Histogram computed successfully for window {idx+1}")
    except Exception as e:
        print(f"Histogram computation failed for window {idx+1} (r0={r0}): {e}")
        continue  # Skip to the next window

    # Normalize the histogram
    hist_norm = hist / np.sum(hist)
    windows_histogram.append(hist_norm)

    # Calculate bin centers and store them
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_centers_list.append(bin_centers)

# Check if any histograms were generated
if not windows_histogram:
    raise RuntimeError("No histograms were generated. Check your data and filtering criteria.")

# Plotting all histograms with their corresponding bin centers
plt.figure(figsize=(10, 6))
for idx, (hist_norm, bin_centers) in enumerate(zip(windows_histogram, bin_centers_list)):
    plt.plot(bin_centers, hist_norm, label=f"Window {idx+1} (r0={r0_lst[idx]:.2f})")

plt.xlabel("Reaction Coordinate (Å)")
plt.ylabel("Probability")
plt.title("Umbrella Window Histograms")
plt.legend()
plt.tight_layout()
plt.savefig(f"umbrella_window_histograms_{umbrella_sampling_num_steps}_100.png")


# Animation part. Uncomment to access

'''
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
'''

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Operation time {elapsed_time}")
plt.show()
plt.close()

# Compute two trials of umbrella sampling with 5000 steps to show the inconsistency
# Compute 10000, 50000 and 250000 (representing 20 ps, 100 ps, and 500 ps) 
