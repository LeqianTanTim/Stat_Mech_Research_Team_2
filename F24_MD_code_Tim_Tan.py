import csv
import numpy as np
class Particle:
    def __init__(self, atom_type, position, fixed_or_not, charge):
        self.atom_type = atom_type
        self.position = position
        self.fixed_or_not = fixed_or_not
        self.charge = charge

        # Initialize epsilon and sigma based on the atom type
        self.epsilon, self.sigma = self._assign_parameters()
        self.velocity = (0,0)
        self.acceleration = (0,0)

    def _assign_parameters(self):
        # epsilon (KJ/mol)
        # sigma (angstrom)
        if self.atom_type == "CH4":
            epsilon = 0.15
            sigma = 3.7
        elif self.atom_type == "H2O":
            epsilon = 0.65
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
        return (f"Particle(atom_type={self.atom_type}, position={self.position}, "
                f"fixed_or_not={self.fixed_or_not}, charge={self.charge}, "
                f"epsilon={self.epsilon}, sigma={self.sigma})")
    
    def __sub__(self, other):
        if not isinstance(other, Particle):
            raise TypeError("Subtraction only supported between Particle objects.")
        else:
            pos_1 = np.array(self.position)
            pos_2 = np.array(other.position)
            relative_pos = pos_1 - pos_2
            sigma = (self.sigma + other.sigma)/2
            epsilon = np.sqrt(self.epsilon * other.epsilon)
        return relative_pos, sigma, epsilon

# load particles from csv: returned a list of particles where each particle is loaded info from a csv file
def load_particles_from_csv(file_path):
    particle_list = []

    with open(file_path, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            # Parse values from the CSV
            index = int(row['index'])
            atom_type = row['atom_type']
            position = tuple(map(float, row['pos'].strip("()").split(',')))  # Convert position string to tuple
            fixed_or_not = row['fixed_or_not'].strip().lower() == 'true'  # Convert string to boolean
            charge = float(row['charge'])

            # Create Particle instance and append it to the list
            particle = Particle(atom_type, position, fixed_or_not, charge)
            particle_list.append(particle)

    # Convert list to numpy array
    return np.array(particle_list)
# Calculate Temperature: return average kinetic energy, temperature
def calculate_temperature(velocities, box_size, dimensions, num_particles):
    total_kinetic_energy = 0.0
    for i in range(num_particles):
        scaled_velocity = box_size * velocities[i, :]
        total_kinetic_energy += 0.5 * np.dot(scaled_velocity, scaled_velocity)
    average_kinetic_energy = total_kinetic_energy / num_particles
    temperature = (2.0 * average_kinetic_energy) / dimensions
    return average_kinetic_energy, temperature
# relative particle: Help determine the relationship between a particle and the rest of the particles in the list.
def relative_particle(particle, particle_lst):
    pos_lst = []
    sigma_lst = []
    epsilon_lst = []
    for i in range(len(particle_lst)):
        relative_position, sigma_new, epsilon_new = particle - particle_lst[i]
        pos_lst.append(relative_position)
        sigma_lst.append(sigma_new)
        epsilon_lst.append(epsilon_new)
    adjusted_pos = pos_lst/sigma_lst
    return adjusted_pos, epsilon_lst
# Calculate forces: Updated accelerations, average potential energy, and virial coefficient.
def compute_forces(positions, accelerations, potential_energies, box_size, dimensions, particle_list):
    #Initialization
    potential_energies.fill(0.0)
    accelerations.fill(0.0)
    virial_coefficient = 0.0
    num_particles = len(particle_list)
    cutoff_radius = 2.5
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
      relative_positions, epsilon_lst = relative_particle(particle_list[i], particle_list[i+1:])
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
      # Compute 1 / r^2 for all pairs
      inverse_squared_distances = 1.0 / squared_distances
      # Compute (1 / r^2)^3 = 1 / r^6
      inverse_sixth_distances = inverse_squared_distances ** 3
      # Compute (1 / r^6)^2 = 1 / r^12
      inverse_twelfth_distances = inverse_sixth_distances ** 2
      # Lennard Jones potential
      # As it was shown earlier, shifted potential improve cutoff smoothness
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
      accelerations[i] += np.sum(forces, axis=0)
      # Apply equal and opposite forces
      accelerations[i + 1:][within_cutoff] -= forces

    # Final Calculations and returns
    # Average potential energy per particle
    avg_potential_energy = np.sum(potential_energies) / num_particles
    # Normalize virial coefficient by dimensions
    virial_coefficient /= -dimensions

    return accelerations, avg_potential_energy, virial_coefficient