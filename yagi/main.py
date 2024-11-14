import multiprocessing as mp
from pathlib import Path
import random
import numpy as np
from tqdm import tqdm
from copy import deepcopy
import json
import math
from datetime import datetime
from simulation import run_simulation, check_parameters, FIELDS, MAX_DIRECTORS, MIN_DIRECTORS, random_configuration, adjust_to_min_difference

def save_result_to_json(params, result, score, filename):
    """Save a single result to the JSON file"""
    data = {
        'timestamp': datetime.now().isoformat(),
        'parameters': params,
        'simulation_results': result,
        'objective_score': score
    }
    
    try:
        # Load existing data if file exists
        if Path(filename).exists():
            with open(filename, 'r') as f:
                existing_data = json.load(f)
        else:
            existing_data = {'configurations': []}
        
        # Append new data
        existing_data['configurations'].append(data)
        
        # Save updated data
        with open(filename, 'w') as f:
            json.dump(existing_data, f, indent=2)
    except Exception as e:
        print(f"Error saving to JSON: {e}")

def objective_function(sim_data, design_params):
    # Extract metrics
    f0 = design_params['f0']
    if math.isnan(sim_data.get('f_resonance', 0)):
        return float('inf')
    resonant_freq = sim_data.get('f_resonance', 0)
    s11_at_f0 = sim_data.get('s11_mag', 0)
    directivity = sim_data.get('directivity', 0)

    # Calculate cost components
    freq_error = abs(resonant_freq - f0) / f0
    s11_cost = max(s11_at_f0 + 10, 0)
    directivity_cost = -directivity

    return 10.0 * freq_error + 5.0 * s11_cost + 1.0 * directivity_cost


def new_population(population_size, top_params, temperature, max_length, min_res):
    # Generate a new population from the current bests, top_params is a ordered
    # array containing the best solutions [(params, score), ...], from best to worst
    population = []
    
    # Calculate how many random vs mutated solutions based on temperature
    # Higher temperature = more random solutions
    random_count = int(population_size * temperature)
    
    # Add some completely random configurations
    if random_count > 0:
        population.extend(random_configuration(max_length, random_count, min_res))
    
    # Generate mutated versions of the top solutions
    while len(population) < population_size:
        # Select one of the elite solutions weighted by their rank
        weights = [1/(i+1) for i in range(len(top_params))]  # Higher weight for better solutions
        total_weight = sum(weights)
        weights = [w/total_weight for w in weights]  # Normalize weights
        chosen = random.choices(top_params, weights=weights, k=1)[0][0]  # Select params only

        new_params = mutate_parameters(chosen, temperature, max_length, min_res)
        if new_params is not None:
            population.append(new_params)
    
    return population

def mutate_parameters(params, temperature, max_length, min_res):
    max_attempts = 10
    min_diff = min_res*2
    for _ in range(max_attempts):
        new_params = deepcopy(params)
        
        # Mutation strength based on temperature
        mutation_strength = 0.3 * temperature  # 30% maximum variation at highest temperature
        
        # Keep track of all lengths for minimum difference checking
        all_lengths = []
        
        # Mutate basic parameters
        for param, (fmin, fmax) in FIELDS.items():
            if not 'director' in param:
                current_val = new_params[param]
                    
                # Calculate allowed variation range
                variation = (fmax - fmin) * mutation_strength
                new_val = current_val + random.gauss(0, variation)
                new_val = np.clip(new_val, fmin, fmax)
                
                # Apply minimum difference check for length parameters
                if 'length' in param:
                    new_val = adjust_to_min_difference(new_val, all_lengths, min_diff)
                    all_lengths.append(new_val)
                    
                new_params[param] = new_val
        
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
                new_length = random.uniform(min_len, max_len)
                new_length = adjust_to_min_difference(new_length, all_lengths, min_diff)
                all_lengths.append(new_length)
                new_params['director_lengths'].append(new_length)
                new_params['director_spacings'].append(
                    random.uniform(min_space, max_space)
                )
        
        # Mutate existing directors
        for i in range(len(new_params['director_lengths'])):
            # Mutate lengths
            min_len, max_len = FIELDS['director_lengths']
            variation = (max_len - min_len) * mutation_strength
            new_len = new_params['director_lengths'][i] + random.gauss(0, variation)
            new_len = np.clip(new_len, min_len, max_len)
            new_len = adjust_to_min_difference(new_len, all_lengths, min_diff)
            all_lengths.append(new_len)
            new_params['director_lengths'][i] = new_len
        

        # Mutate spacings 
        # Shift the elements around
        for i in range(len(new_params['director_spacings'])-1):
            min_space, _ = FIELDS['director_spacings']
            max_space = new_params['director_spacings'][i+1]/2
            variation = (max_space - min_space) * mutation_strength
            new_space = new_params['director_spacings'][i] + random.gauss(0, variation)
            new_space = np.clip(new_space, min_space, max_space)
            difference = new_space - new_params['director_spacings'][i]
            # Keep overall length of antenna the same if it will make it too long
            if (sum(new_params['director_spacings']) + new_params['feed_point'] > max_length):
                new_params['director_spacings'][i+1] -= difference
            new_params['director_spacings'][i] = new_space

        # Mutate last element
        spare_length = max_length - sum(new_params['director_spacings']) + new_params['feed_point'];
        if (spare_length > 0):
            min_space, _ = FIELDS['director_spacings']
            max_space = spare_length
            variation = (max_space - min_space) * mutation_strength
            new_space = new_params['director_spacings'][-1] + random.gauss(0, variation)
            new_params['director_spacings'][-1] = np.clip(new_space, min_space, max_space)
        
        # Verify the new configuration is valid
        if check_parameters(new_params, max_length, min_res=min_res):
            return new_params
    
    # If we couldn't generate a valid configuration after max_attempts, return None
    return None

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
    
    # Print director parameters in aligned columns
    num_directors = len(params['director_spacings'])
    print("\nDirectors:")
    print("    " + "".join([f"{f'D{i}':>8} " for i in range(1, num_directors + 1)]))
    print("Spa " + "".join([f"{spacing:8.2f} " for spacing in params['director_spacings']]))
    print("Len " + "".join([f"{length:8.2f} " for length in params['director_lengths']]))

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
    
    # Create results directory if it doesn't exist
    RESULTS_DIR = Path("./results")
    RESULTS_DIR.mkdir(exist_ok=True)
    
    # Create unique filename for this optimization run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_file = RESULTS_DIR / f"optimization_results_{timestamp}.json"

    # Design parameters
    DESIGN_PARAMS = {'f0': 437.25e6, 'fc': 50e6, 'R': 50}
    MAX_LENGTH = 2000 # in millimeters
    MIN_RES=5.0

    # Save design parameters to JSON
    with open(results_file, 'w') as f:
        json.dump({
            'meta': {
                'timestamp': datetime.now().isoformat(),
                'design_parameters': DESIGN_PARAMS,
                'max_length': MAX_LENGTH,
                'min_resolution': MIN_RES
            },
            'configurations': []
        }, f, indent=2)

    # [Initial parameters remain unchanged]
    INITIAL_PARAMS = {
        'reflector_length': 340, 
        'feed_point': 164.55, 
        'feed_gap': 10.0,
        'feed_length': 330,
        'director_spacings': [51.42, 123.41, 147.41, 171.41, 191.98, 205.69, 215.97, 226.26, 236.54, 246.83],
        'director_lengths': [320, 310, 310, 310, 310, 300, 300, 290, 290, 290]
    }

    if not check_parameters(INITIAL_PARAMS, MAX_LENGTH, min_res=MIN_RES):
        print("Invalid starting config!!")
        exit()
    else:
        print("Evaluating starting point")
        params, res = run_simulation(prepare_mp_arguments(SIM_DIR, DESIGN_PARAMS, [INITIAL_PARAMS], scripts=OCTAVE_SCRIPTS)[0])
        if res['success']:
            score = objective_function(res['result'], DESIGN_PARAMS)
            print_result(params, res['result'], score)
            save_result_to_json(params, res['result'], score, results_file)
        else:
            print(f"Simulation failed: {res.get('error', 'Unknown error')}")
            exit()

    # [Optimizer parameters remain unchanged]
    INITIAL_TEMP = 0.75
    FINAL_TEMP = 0.01
    COOLING_RATE = 0.95
    TOP_N_COUNT = 15
    SIM_PARALLEL = 32
    SIM_POPULATION = 64
    SIM_ITERATIONS = 50

    # Create the working directory if it does not yet exist
    SIM_DIR.mkdir(exist_ok=True)

    temperature = INITIAL_TEMP
    iteration = 0

    pool = mp.Pool(processes=SIM_PARALLEL)

    overall_best_cost = score
    overall_best_params = INITIAL_PARAMS.copy()
    overall_best_results = res['result']

    best_solutions = [(INITIAL_PARAMS, score)] 

    # Starting population
    population = new_population(SIM_POPULATION, best_solutions, temperature, MAX_LENGTH, MIN_RES)

    while temperature > FINAL_TEMP and iteration < SIM_ITERATIONS:
        print(f"\nIteration {iteration + 1}/{SIM_ITERATIONS}, Temperature: {temperature:.4f}")

        processed_results = []
        processed_population = []
        costs = []

        # Track metrics
        error_count = 0
        curr_best = float('inf')
        curr_best_direc = 1
        sim_times = []

        # Run the simulation for the current population multithreaded
        simulation_params = prepare_mp_arguments(SIM_DIR, DESIGN_PARAMS, population, scripts=OCTAVE_SCRIPTS)
        with tqdm(total=len(simulation_params), desc="Running simulations", 
                 bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, '
                           'avg: {rate_fmt}{postfix}]') as sim_pbar:
            for params, res in pool.imap_unordered(run_simulation, simulation_params):
                if res['success']:
                    # Calculate score if successfull
                    score = objective_function(res['result'], DESIGN_PARAMS)
                    # Save to processed as progenitor for next generation
                    processed_population.append(params)
                    processed_results.append(res['result'])
                    costs.append(score)
                    # Log each successful simulation
                    save_result_to_json(params, res['result'], score, results_file)

                    # Update metrics
                    sim_times.append(res['time'])
                    if score < curr_best:
                        curr_best = score
                        curr_best_direc = res['result']['directivity']
                else:
                    # Log error count
                    error_count += 1

                # Update progress bar with metrics
                avg_time = np.mean(sim_times)
                sim_pbar.set_postfix({
                    'sim_time': f'{avg_time:.1f}s',
                    'best_cost': f'{curr_best:.3f}',
                    'best_gain': f'{10*math.log10(curr_best_direc):.1f}dB',
                    'errors': error_count
                })
                sim_pbar.update(1)

        # Update absolute best solution record
        min_cost_idx = np.argmin(costs)
        print(f"Best cost this iteration: {costs[min_cost_idx]}")
        if costs[min_cost_idx] < overall_best_cost:
            overall_best_cost = costs[min_cost_idx]
            overall_best_params = processed_population[min_cost_idx]
            overall_best_results = processed_results[min_cost_idx]
            print(f"New best cost: {overall_best_cost}")

        # Generate list of top solutions for current iteration
        current_solutions = list(zip(processed_population, costs))
        current_solutions.sort(key=lambda x: x[1])  # Sort by cost

        # Add TOP_N_COUNT to the list of best solutions
        best_solutions = best_solutions + current_solutions
        best_solutions.sort(key=lambda x: x[1])
        best_solutions = best_solutions[:min(len(best_solutions), TOP_N_COUNT)]

        # Generate a new population based on the best solutions
        population = new_population(SIM_POPULATION, best_solutions, temperature, MAX_LENGTH, MIN_RES)

        temperature *= COOLING_RATE
        iteration += 1

        if not overall_best_results is None:
            print("Best result this iteration:")
            print_result(overall_best_params, overall_best_results, overall_best_cost)
        print()

    pool.close()

    print("Best result:")
    print_result(overall_best_params, overall_best_results, overall_best_cost)