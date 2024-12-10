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
        self.epsilon, self.sigma, self.mass = self._assign_parameters()
        # Initialize velocities and accelerations with random values
        if self.fixed_or_not:
            self.velocity = np.zeros(self.dimension)
            self.acceleration = np.zeros(self.dimension)
        else:
            self.velocity = np.random.randn(self.dimension) - 0.5
            self.acceleration = np.random.randn(self.dimension) - 0.5

    def _assign_parameters(self):
        # epsilon (J per particle)
        # sigma (angstrom)
        if self.atom_type == "CH4":
            epsilon = 1.23 * 1e3 / 6.022e23 
            sigma = 3.73
            mass = 2.66404e-26 # 16 amu, in kg
        elif self.atom_type == "H2O":
            epsilon = 0.636 * 1e3 / 6.022e23 
            sigma = 3.15
            mass = 2.99151e-26 # 18.015 amu, in kg
        elif self.atom_type == "surface": # assume graphene lattice model, approximated for 10 angstrom, 57.8 unit cells, each unit cell contain 2 carbon atoms
            epsilon = 0.23 * 1e3 / 6.022e23 
            sigma = 3.4
            mass = 2.3e-24
        else:
            # Default values for unknown atom types
            print('Invalid particle type was input, please update the current force field before running.')
            epsilon = 0.0
            sigma = 0.0
            mass = 0.0
        return epsilon, sigma, mass


    def __repr__(self):
        return (f"Particle(atom_type={self.atom_type}, position={self.position}, dimension= {self.dimension}, "
                f"fixed_or_not={self.fixed_or_not}, charge={self.charge}, "
                f"epsilon={self.epsilon}, sigma={self.sigma}, "
                f"velocity={self.velocity}, acceleration={self.acceleration}, mass={self.mass} )")
    
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
    def get_acc(self):
        return self.acceleration
    def get_pos(self):
        return self.position

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
        old_pos = particle.get_pos()
        particle.update_pos(old_pos/box_size)
        old_vel = particle.get_vel()
        old_acc = particle.get_acc()
        particle.update_vel(old_vel/box_size)
        particle.update_acc(old_acc/box_size)
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

# Calculate Temperature: return average kinetic energy in J, temperature in K
def calculate_temperature(particle_lst, box_size):
    kb = 1.3806e-23 # J/K
    total_kinetic_energy = 0.0
    num_particles = particle_lst.size
    for particle in particle_lst:
        scaled_velocity = box_size * particle.get_vel() * 5e4 # convert from angstrom/2fs to m/s
        total_kinetic_energy += 0.5 * particle.mass * np.dot(scaled_velocity, scaled_velocity) # 1/2 mv^2
    average_kinetic_energy = total_kinetic_energy / num_particles # in J
    temperature = (average_kinetic_energy) / (particle_lst.size * kb) # in K
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
    return pos_array, epsilon_array, sigma_array

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

    #Loop over particles
    for i in range(num_particles - 1):
      # Compute relative positions for all other particles
      relative_positions, epsilon_lst, sigma_lst = relative_particle(particle_lst[i], particle_lst[i+1:])
      cutoff_radius = 4e-10 # 4 angstrom
      sigma_lst = sigma_lst * 1e-10
      cutoff_radius_squared = cutoff_radius ** 2
      #print(f"the cutoff_radius_squared was {cutoff_radius_squared} in m")
      #Apply periodic boundary conditions
      #adjusts relative positions for particles that cross the boundary
      relative_positions -= np.rint(relative_positions)
      #Convert to real units (Scale from dimensionless to box dimension)
      relative_positions *= (box_size * 1e-10) # convert angstrom to m
      
      #Compute squared distance for all pairs
      squared_distances = np.sum(relative_positions ** 2, axis=1)
      #print(squared_distances)
      #Apply cutoff Radius
      # Boolean mask for pairs within the cutoff radius
      within_cutoff = squared_distances < cutoff_radius_squared
      # Filter positions within cutoff
      #print(f"Before filter {relative_positions.size}")
      relative_positions = relative_positions[within_cutoff]
      #print(f"After cutoff {relative_positions.size}")
      # Filter distances within cutoff
      squared_distances = squared_distances[within_cutoff]
      #print(f"Squared distance \n {squared_distances}")
      # Filter the epsilon value within cutoff
      epsilon_lst = epsilon_lst[within_cutoff]
      sigma_lst = sigma_lst[within_cutoff]
      #print(f"Sigma list \n {sigma_lst}")
      # Skip if no pairs within cutoff
      if len(squared_distances) == 0:
        continue
      sigma_squared = sigma_lst **2
      #Lennard-Jones Potential and Force calculation
      squared_distances[squared_distances == 0] = np.inf  # prevent divisionbyzero error
      # compute 1/r
      inverse_distance = 1 / (np.sqrt(squared_distances))
      # Compute sigma^2 / r^2 for all pairs
      inverse_squared_distances = sigma_squared / squared_distances
      #print(f"(sigma / r)^2 = {inverse_squared_distances}")
      # Compute (sigma^2 / r^2)^3 = sigma^6 / r^6
      inverse_sixth_distances = inverse_squared_distances ** 3
      # Compute (sigma^6 / r^6)^2 = sigma^12 / r^12
      inverse_twelfth_distances = inverse_sixth_distances ** 2
      # Lennard Jones potential
      # As it was shown earlier, shifted potential improve cutoff smoothness

      shifted_potential_cutoff = 4.0 * ((sigma_squared / cutoff_radius_squared) ** 6 - (sigma_squared / cutoff_radius_squared) ** 3)
      potential = epsilon_lst * (
          4.0 * (
              inverse_twelfth_distances
            - inverse_sixth_distances)
            - shifted_potential_cutoff)
      #print(f"Potential is \n {potential}")
      #Compute the magnitude of the force from the LJ formula
      #F = -dV/dr
      force_magnitude = epsilon_lst  * 24.0 * inverse_distance * (
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
      # relative position is in m
      forces = (force_magnitude[:, np.newaxis] *
                relative_positions) / np.sqrt(
                    squared_distances)[:, np.newaxis]
      # Add forces to particle i
      if not particle_lst[i].fixed_or_not:
        current_acc = particle_lst[i].acceleration
        #print(f"previously acceleration for {particle_lst[i].atom_type} is {current_acc} angstrom/4fs^2")
        new_acc = current_acc * 0.25 * 1e20 + np.sum(forces, axis=0)/particle_lst[i].mass # convert the acceleration to m/s^2 before adding F/m = a
        #print(f"new acceleration for {particle_lst[i].atom_type} is {new_acc} m/s^2")
        particle_lst[i].update_acc(new_acc * 0.25 * 1e-20) # convert back to angstrom/4fs^2
            
      for k in range(particle_lst[i+1:][within_cutoff].size):
        if not particle_lst[i+1:][within_cutoff][k].fixed_or_not: # Not fixed
            current_acc = particle_lst[i+1:][within_cutoff][k].acceleration * 0.25 * 1e20
            new_acc = current_acc - forces[k]/particle_lst[i+1:][within_cutoff][k].mass
            particle_lst[i+1:][within_cutoff][k].update_acc(new_acc * 0.25 * 1e-20)
    # Final Calculations and returns
    # Average potential energy per particle
    avg_potential_energy = np.sum(potential_energies) / num_particles
    # Normalize virial coefficient by dimensions
    virial_coefficient /= -dimensions

    return particle_lst, avg_potential_energy, virial_coefficient

# Apply Boundary conditions: return a particle list with position fixed
def apply_bc(particle_lst):
    for particle in particle_lst:
        current_position = particle.get_pos()
        new_position = current_position - np.rint(current_position)
        particle.update_pos(new_position)
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
        if not particle.fixed_or_not:
            old_vel = particle.get_vel()
            new_vel = scale * old_vel
            particle.update_vel(new_vel)
    return particle_lst

# update velocity based on acc computed by force: return a particle list with velocity updated
def complete_force_update(particle_lst, time_step):
    for particle in particle_lst:
        if not particle.fixed_or_not:
            current_vel = particle.get_vel()
            particle.update_vel(current_vel + time_step * particle.acceleration)
    return particle_lst