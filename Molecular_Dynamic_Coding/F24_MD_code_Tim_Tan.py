import csv
import numpy as np
class Particle:
    def __init__(self, atom_type, position, fixed_or_not, charge):
        self.atom_type = atom_type
        self.position = np.array(position, dtype=float)
        self.fixed_or_not = fixed_or_not
        self.charge = charge
        self.dimension = self.position.size
        # Initialize epsilon and sigma based on the atom type
        self.epsilon, self.sigma = self._assign_parameters()
        # Initialize velocities and accelerations with random values
        if self.fixed_or_not:
            self.velocity = np.zeros(self.dimension)
            self.acceleration = np.zeros(self.dimension)
        else:
            np.random.seed(1)
            self.velocity = np.random.randn(self.dimension) - 0.5
            self.acceleration = np.random.randn(self.dimension) - 0.5

    def _assign_parameters(self):
        # epsilon (KJ)
        # sigma (angstrom)
        if self.atom_type == "CH4":
            epsilon = 0.15 / 6.022e23
            sigma = 3.7
        elif self.atom_type == "H2O":
            epsilon = 0.65 / 6.022e23
            sigma = 3.15
        elif self.atom_type == "surface":
            epsilon = 0.15
            sigma = 1.3
        else:
            # Default values for unknown atom types
            print('Invalid particle type was input, please update the current force field before running.')
            epsilon = 0.05
            sigma = 1.2

        return epsilon, sigma


    def __repr__(self):
        return (f"Particle(atom_type={self.atom_type}, position={self.position}, dimension= {self.dimension}, "
                f"fixed_or_not={self.fixed_or_not}, charge={self.charge}, "
                f"epsilon={self.epsilon}, sigma={self.sigma}, "
                f"velocity={self.velocity}, acceleration={self.acceleration},  )")
    
    def __sub__(self, other):
        if not isinstance(other, Particle):
            raise TypeError("Subtraction only supported between Particle objects.")
        else:
            pos_1 = self.position
            pos_2 = other.position
            relative_pos = pos_1 - pos_2
            sigma = (self.sigma + other.sigma)/2
            epsilon = np.sqrt(self.epsilon * other.epsilon)
        return relative_pos, sigma, epsilon
    
    def update_pos(self, new_value):
        old_pos = self.position
        self.position = new_value

        # Test statement
        #print(f"Position has been updated from {old_pos} to new pos {self.position}")

    def update_vel(self, new_value):
        old_vel = self.velocity
        self.velocity = new_value
        # Test statement
        #print(f"Velocity has been updated from {old_vel} to new velocity {self.velocity}")

    def update_acc(self, new_value):
        old_acc = self.acceleration
        self.acceleration = new_value
        # Test statement
        #print(f"Acceleration has been updated from {old_acc} to new acceleration {self.acceleration}")

    def zero_acc(self):
        self.acceleration = np.zeros(self.dimension)

    def get_vel(self):
        return self.velocity

# load particles from csv: returned a list of particles where each particle is loaded info from a csv file
def load_particles_from_csv(file_path):
    particle_lst = [] 

    with open(file_path, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            # Parse values from the CSV
            index = int(row['index'])
            atom_type = row['atom_type']
            position = tuple(map(float, row['position'].strip("()").split(',')))  # Convert position string to tuple
            dimension = len(position)
            fixed_or_not = row['fixed'].strip().lower() == 'true'  # Convert string to boolean
            charge = float(row['charge'])

            # Create Particle instance and append it to the list
            particle = Particle(atom_type, position, fixed_or_not, charge)
            particle_lst.append(particle)

    # Convert list to numpy array
    return np.array(particle_lst), dimension

# Normalization: Normalize the particles with the box size, return the particle list with normalized
def normalization(particle_lst, box_size):
    for particle in particle_lst:
        old_pos = particle.position
        particle.update_pos(particle.position/box_size)
    return particle_lst

# Shift position to center of mass: return the particle list that is shift to center of mass as origin
def shift_com(particle_lst):
    num_particle = particle_lst.size
    sum_pos = np.zeros(particle_lst[0].dimension)
    for particle in particle_lst:
        sum_pos += particle.position
    center_of_mass = sum_pos / num_particle
    for particle in particle_lst:
        old_pos = particle.position
        particle.update_pos(old_pos - center_of_mass)
    return particle_lst

# Calculate Temperature: return average kinetic energy, temperature
def calculate_temperature(particle_lst, box_size):
    total_kinetic_energy = 0.0
    num_particles = particle_lst.size
    for particle in particle_lst:
        scaled_velocity = box_size * particle.get_vel()
        total_kinetic_energy += 0.5 * np.dot(scaled_velocity, scaled_velocity)
    average_kinetic_energy = total_kinetic_energy / num_particles
    temperature = (2.0 * average_kinetic_energy) / particle_lst[0].dimension
    return average_kinetic_energy, temperature

# relative particle: Help determine the relationship between a particle and the rest of the particles in the list.
def relative_particle(particle, particle_lst):
    pos_lst = []
    sigma_lst = []
    epsilon_lst = []
    for i in range(len(particle_lst)):
        # I overload my subtraction operator in the class. 
        # You can check back the class for what happened when two particle class substract each other.
        relative_position, sigma_new, epsilon_new = particle - particle_lst[i]
        pos_lst.append(relative_position)
        sigma_lst.append(sigma_new)
        epsilon_lst.append(epsilon_new)
    pos_array = np.array(pos_lst)
    sigma_array = np.array(sigma_lst)
    epsilon_array = np.array(epsilon_lst)
    adjusted_pos = []
    for i in range(sigma_array.size):
        adjusted_pos.append(pos_array[i]/sigma_array[i])
    adjusted_pos_array = np.array(adjusted_pos)
    return adjusted_pos_array, epsilon_array

# Reset acceleration for all particles in the list: return a list where acceleration is reset. 
def reset_acc(particle_lst):
    for particle in particle_lst:
        particle.zero_acc()
    return particle_lst

# Calculate forces: Updated accelerations, average potential energy, and virial coefficient.
def compute_forces(particle_lst, potential_energies, box_size, dimensions):
    #Initialization
    potential_energies.fill(0.0)
    particle_lst = reset_acc(particle_lst)
    virial_coefficient = 0.0
    num_particles = len(particle_lst)
    cutoff_radius = 2.0
    # Precompute the square of the cutoff radius
    
    cutoff_radius_squared = cutoff_radius ** 2
    '''
    shifted_potential_cutoff = epsilon * (
        4.0 * ((1.0 / cutoff_radius_squared) ** 6 - (
            1.0 / cutoff_radius_squared) ** 3))
    '''
    #Cutoff Radius is pre-normalized by sigma
    #Loop over particles
    for i in range(num_particles - 1):
      # Compute relative positions for all other particles
      relative_positions, epsilon_lst = relative_particle(particle_lst[i], particle_lst[i+1:])
      #Apply periodic boundary conditions
      #adjusts relative positions for particles that cross the boundary
      relative_positions -= np.rint(relative_positions)
      #Convert to real units (Scale from dimensionless to box dimension)
      relative_positions *= box_size
      
      #Compute squared distance for all pairs
      squared_distances = np.sum(relative_positions ** 2, axis=1)

      #Apply cutoff Radius
      # Boolean mask for pairs within the cutoff radius
      within_cutoff = squared_distances < cutoff_radius_squared
      # Filter positions within cutoff
      relative_positions = relative_positions[within_cutoff]
      # Filter distances within cutoff
      squared_distances = squared_distances[within_cutoff]
      # Filter the epsilon value within cutoff
      epsilon_lst = epsilon_lst[within_cutoff]
      # Skip if no pairs within cutoff
      if len(squared_distances) == 0:
        continue

      #Lennard-Jones Potential and Force calculation
      squared_distances[squared_distances == 0] = np.inf  # prevent divisionbyzero error
      # Compute 1 / r^2 for all pairs
      inverse_squared_distances = 1.0 / squared_distances
      # Compute (1 / r^2)^3 = 1 / r^6
      inverse_sixth_distances = inverse_squared_distances ** 3
      # Compute (1 / r^6)^2 = 1 / r^12
      inverse_twelfth_distances = inverse_sixth_distances ** 2
      # Lennard Jones potential
      # As it was shown earlier, shifted potential improve cutoff smoothness

      shifted_potential_cutoff = epsilon_lst * (
        4.0 * ((1.0 / cutoff_radius_squared) ** 6 - (
            1.0 / cutoff_radius_squared) ** 3))

      potential = epsilon_lst * (
          4.0 * (
              inverse_twelfth_distances
            - inverse_sixth_distances)
            - shifted_potential_cutoff)
      #Compute the magnitude of the force from the LJ formula
      #F = -dV/dr
      force_magnitude = epsilon_lst  * 24.0 * inverse_squared_distances * (
          2.0 * inverse_twelfth_distances
          - inverse_sixth_distances)
      # Equally splitting between two interacting particles
      # Add half of the potential energy to particle i
      potential_energies[i] += np.sum(0.5 * potential)
      # Add half to the interacting particles
      potential_energies[i + 1:][within_cutoff] += 0.5 * potential
      # Virial coefficient: Used for pressure calculations
      virial_coefficient += np.sum(
          force_magnitude * np.sqrt(squared_distances))
      # Force vectors are computed and added to accelerations [i]
      forces = (force_magnitude[:, np.newaxis] *
                relative_positions) / np.sqrt(
                    squared_distances)[:, np.newaxis]
      # Add forces to particle i
      if not particle_lst[i].fixed_or_not:
        current_acc = particle_lst[i].acceleration
        particle_lst[i].update_acc(current_acc + np.sum(forces, axis=0))
        for k in range(particle_lst[i+1:][within_cutoff].size):
            if not particle_lst[i+1:][within_cutoff][k].fixed_or_not: # Not fixed
                current_acc = particle_lst[i+1:][within_cutoff][k].acceleration
                new_acc = current_acc - forces[k]
                particle_lst[i+1:][within_cutoff][k].update_acc(new_acc)
    # Final Calculations and returns
    # Average potential energy per particle
    avg_potential_energy = np.sum(potential_energies) / num_particles
    # Normalize virial coefficient by dimensions
    virial_coefficient /= -dimensions

    return particle_lst, avg_potential_energy, virial_coefficient

# Apply Boundary conditions: return a particle list with position fixed
def apply_bc(particle_lst):
    for particle in particle_lst:
        current_position = particle.position
        for i in range(current_position.size): # upper bound
            pos = current_position[i]
            if pos > 0.5:
                pos = pos - 1.0
                current_position[i] = pos
        for j in range(current_position.size): # lower bound
            pos = current_position[i]
            if pos < -0.5:
                pos = pos + 1.0
                current_position[i] = pos
    return particle_lst

# update pos with verlet velocity: return a particle list with position updated
def apply_verlet_pos(particle_lst, time_step):
    for particle in particle_lst:
        if particle.fixed_or_not:
            continue
        else:
            current_position = particle.position
            new_position = current_position + time_step * particle.velocity + 0.5 * (time_step ** 2) * particle.acceleration
            particle.update_pos(new_position)
    return particle_lst

# scale velocity with temperature: return a particle list with velocity scaled
def scale_velocity(particle_lst, scale, time_step):
    for particle in particle_lst:
        old_vel = particle.velocity
        new_vel = scale * old_vel + 0.5 * time_step * particle.acceleration
        particle.update_vel(new_vel)
    return particle_lst

# update velocity based on acc computed by force: return a particle list with velocity updated
def complete_force_update(particle_lst, time_step):
    for particle in particle_lst:
        if not particle.fixed_or_not:
            current_vel = particle.velocity
            particle.update_vel(current_vel + 0.5 * time_step * particle.acceleration)
    return particle_lst

# look up the particle list to see the current position and distance between the two solute particle
# returned the reaction coordinate at that step. 
# We might be able to map a scatter plot for R and energy (Let's see how we can implement this)
def reaction_coordinate_lookup(particle_lst):
    solute_lst = []
    solvent_lst = [] # In case anyone wants to modify this.
    R = 0.0
    for particle in particle_lst:
        if particle.atom_type == "H2O":
            solvent_lst.append(particle)
        else:
            solute_lst.append(particle)
    solute_array = np.array(solute_lst)
    if solute_array.size > 2:
        print("More than one solute available for pairwise analysis")
        return
    else:
        solute_1 = solute_array[0].position
        solute_2 = solute_array[1].position
        R = np.linalg.norm(solute_1 - solute_2)
    return R