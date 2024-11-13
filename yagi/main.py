import multiprocessing as mp
from pathlib import Path
import random
import numpy as np
from tqdm import tqdm
from copy import deepcopy
from simulation import run_simulation, check_parameters, FIELDS, MAX_DIRECTORS, MIN_DIRECTORS, random_configuration

def objective_function(sim_data, design_params):
    # Extract metrics
    f0 = design_params['f0']
    resonant_freq = sim_data.get('f_resonance', 0)
    s11_at_f0 = sim_data.get('s11_mag', 0)
    directivity = sim_data.get('directivity', 0)

    # Calculate cost components
    freq_error = abs(resonant_freq - f0) / f0
    s11_cost = max(s11_at_f0 + 10, 0)
    directivity_cost = -directivity

    return 10.0 * freq_error + 5.0 * s11_cost + 1.0 * directivity_cost

def new_population(population_size, best_params, temperature, max_length):
    population = []
    
    # Calculate how many random vs mutated solutions based on temperature
    # Higher temperature = more random solutions
    random_count = int(population_size * temperature)
    mutated_count = population_size - random_count
    
    # Add some completely random configurations
    if random_count > 0:
        population.extend(random_configuration(max_length, random_count))
    
    # Generate mutated versions of the best solution
    while len(population) < population_size:
        new_params = mutate_parameters(best_params, temperature, max_length)
        if new_params is not None:
            population.append(new_params)
    
    return population

def mutate_parameters(params, temperature, max_length):
    max_attempts = 10
    for _ in range(max_attempts):
        new_params = deepcopy(params)
        
        # Mutation strength based on temperature
        mutation_strength = 0.3 * temperature  # 30% maximum variation at highest temperature
        
        # Mutate basic parameters
        for f, (fmin, fmax) in FIELDS.items():
            if not 'director' in f:
                current_val = new_params[param]
                    
                # Calculate allowed variation range
                variation = (fmax - fmin) * mutation_strength
                new_val = current_val + random.gauss(0, variation)
                new_params[param] = np.clip(new_val, fmin, fmax)
        
        # Handle director mutations
        if random.random() < temperature:
            # Potentially add or remove directors
            current_directors = len(new_params['director_lengths'])
            if random.random() < 0.5 and current_directors > MIN_DIRECTORS:
                # Remove a random director
                idx = random.randint(0, current_directors - 1)
                new_params['director_lengths'].pop(idx)
                new_params['director_spacings'].pop(idx)
            elif current_directors < MAX_DIRECTORS:
                # Add a new director
                min_len, max_len = FIELDS['director_lengths']
                min_space, max_space = FIELDS['director_spacings']
                new_params['director_lengths'].append(
                    random.uniform(min_len, max_len)
                )
                new_params['director_spacings'].append(
                    random.uniform(min_space, max_space)
                )
        
        # Mutate existing directors
        for i in range(len(new_params['director_lengths'])):
            # Mutate lengths
            min_len, max_len = FIELDS['director_lengths']
            variation = (max_len - min_len) * mutation_strength
            new_len = new_params['director_lengths'][i] + random.gauss(0, variation)
            new_params['director_lengths'][i] = np.clip(new_len, min_len, max_len)
            
            # Mutate spacings
            min_space, max_space = FIELDS['director_spacings']
            variation = (max_space - min_space) * mutation_strength
            new_space = new_params['director_spacings'][i] + random.gauss(0, variation)
            new_params['director_spacings'][i] = np.clip(new_space, min_space, max_space)
        
        # Verify the new configuration is valid
        if check_parameters(new_params, max_length):
            return new_params
    
    # If we couldn't generate a valid configuration after max_attempts, return None
    return None

def acceptance_probability(old_cost, new_cost, temperature):
    if new_cost < old_cost:
        return 1.0
    return np.exp((old_cost - new_cost) / temperature)

# Helper function to prepare the multiprocessing inputs
def prepare_mp_arguments(base_dir, design_params, parameter_sets, scripts="./matlab"):
    return [{'base_dir': base_dir, 'design_params': design_params, 'params': params, 'scripts': scripts} for params in parameter_sets]

def print_result(params, result, score):
    print("\nParameters:")
    # Print wire_rad if present
    if 'wire_rad' in params:
        print(f"Wire Radius: {params['wire_rad']:.2f} mm")
    
    # Print mandatory parameters
    print(f"Reflector Length: {params['reflector_length']:.2f} mm")
    print(f"Feed Point: {params['feed_point']:.2f} mm")
    print(f"Feed Gap: {params['feed_gap']:.2f} mm")
    print(f"Feed Length: {params['feed_length']:.2f} mm")
    
    # Print director parameters
    print("\nDirector Spacings:")
    for i, spacing in enumerate(params['director_spacings'], 1):
        print(f"  D{i}: {spacing:.2f} mm")
    
    print("\nDirector Lengths:")
    for i, length in enumerate(params['director_lengths'], 1):
        print(f"  D{i}: {length:.2f} mm")

    print("\nPerformance:")
    if 'f_resonance' in result:
        print(f"Resonant Frequency: {result['f_resonance']/1e6:.2f} MHz")
    if 's11_mag' in result:
        s11_db = 20 * np.log10(result['s11_mag'])
        print(f"S11: {s11_db:.1f} dB")
    if 'directivity' in result:
        print(f"Directivity: {10*np.log10(result['directivity']):.1f} dBi")
    if 'efficiency' in result:
        print(f"Efficiency: {result['efficiency']*100:.1f}%")
    if 'HPBW' in result:
        print(f"HPBW: {result['HPBW']:.1f}°")
   
    print(f"\nObjective Function Score: {score:.3f}")

if __name__ == "__main__":
    OCTAVE_SCRIPTS = Path("./matlab")
    SIM_DIR = Path("./working")

    # Design parameters
    DESIGN_PARAMS = {'f0': 437.25e6, 'fc': 50e6, 'R': 50}
    MAX_LENGTH = 2000 # in millimeters

    # Rough starting point
    INITIAL_PARAMS = {
        #'wire_rad': 2.0, 
        'reflector_length': 334.49, 
        'feed_point': 164.55, 
        'feed_gap': 10.0,
        'feed_length': 330.47,
        'director_spacings': [51.42, 123.41, 147.41, 171.41, 191.98, 205.69, 215.97, 226.26, 236.54, 246.83],
        'director_lengths': [315.72, 312.75, 310.03, 307.54, 305.29, 303.24, 301.40, 299.74, 298.26, 296.94],
    }

    if not check_parameters(INITIAL_PARAMS, MAX_LENGTH):
        print("Invalid starting config!!")
    else:
        print("Evaluating starting point")
        params, res = run_simulation(prepare_mp_arguments(SIM_DIR, DESIGN_PARAMS, [INITIAL_PARAMS], scripts=OCTAVE_SCRIPTS)[0])
        if res['success']:
            print_result(params, res['result'], objective_function(res['result'], DESIGN_PARAMS))
        else:
            print(f"Simulation failed: {res.get('error', 'Unknown error')}")
    
    # Optimizer params
    INITIAL_TEMP = 1.0
    FINAL_TEMP = 0.01
    COOLING_RATE = 0.95
    SIM_PARALLEL = 40
    SIM_POPULATION = 32 # power of 2
    SIM_ITERATIONS = 20

    # Create the working directory if it does not yet exist
    SIM_DIR.mkdir(exist_ok=True)

    # Simulated annealing temperature
    temperature = INITIAL_TEMP
    iteration = 0

    population = new_population(SIM_POPULATION, INITIAL_PARAMS, temperature, MAX_LENGTH)

    pool = mp.Pool(processes=SIM_PARALLEL)

    # Store the absolute best costs
    overall_best_cost = float('inf')
    overall_best_params = INITIAL_PARAMS.copy()
    overall_best_results = None

    # Accpeted solution for simulated annealing
    anneal_best_cost = float('inf')
    anneal_best_params = INITIAL_PARAMS.copy()

    while temperature > FINAL_TEMP and iteration < SIM_ITERATIONS:
        print(f"\nIteration {iteration + 1}/{SIM_ITERATIONS}, Temperature: {temperature:.4f}")

        # Evaluate the current population
        simulation_params = prepare_mp_arguments(SIM_DIR, DESIGN_PARAMS, population, scripts=OCTAVE_SCRIPTS)
        results = []
        with tqdm(total=len(simulation_params), desc="Running simulations", leave=False) as sim_pbar:
            for result in pool.imap_unordered(run_simulation, simulation_params):
                results.append(result)
                sim_pbar.update(1)

        # Discard the members of the population where the simulation
        # failed to successfully run
        processed_results = []
        processed_population = []
        for params, res in results:
            if res['success']:
                processed_population.append(params)
                processed_results.append(res['result'])
            else:
                print(f"Simulation failed: {res.get('error', 'Unknown error')}")

        costs = [objective_function(res, DESIGN_PARAMS) for res in processed_results]
        # Update absolute best solution
        min_cost_idx = np.argmin(costs)
        print(f"Best cost this iteration: {costs[min_cost_idx]}")
        if costs[min_cost_idx] < overall_best_cost:
            overall_best_cost = costs[min_cost_idx]
            overall_best_params = processed_population[min_cost_idx]
            overall_best_results = processed_results[min_cost_idx]
            print(f"New best cost: {overall_best_cost}")

        # Choose best solution with simulated annealing acceptance probability
        for param, cost in zip(processed_population, costs):
            if acceptance_probability(anneal_best_cost, cost, temperature) > np.random.random():
                anneal_best_params = params
                anneal_best_cost = cost
        
        population = new_population(SIM_POPULATION, anneal_best_params, temperature, MAX_LENGTH)

        # Cool down
        temperature *= COOLING_RATE
        iteration += 1

        if not overall_best_results is None:
            print("Best result this iteration:")
            print_result(overall_best_params, overall_best_results, overall_best_cost)
        print()

    pool.close()

    print("Best result:")
    print_result(overall_best_params, overall_best_results, overall_best_cost)
