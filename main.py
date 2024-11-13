import multiprocessing as mp
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from tqdm import tqdm
from simulation import run_simulation, check_parameters, FIELDS

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

def generate_points(n, mesh_spacing, axial_length, feed_height_thres=4, center_point=None):
    # Initialize Sobol sampler for the parameter space
    sampler = qmc.Sobol(
        d=len(FIELDS),
        scramble=True
    )

    l_bounds = [param[0] for param in FIELDS.values()]
    u_bounds = [param[1] for param in FIELDS.values()]
    param_sets = []

    while len(param_sets) < n:
        sample = sampler.random(n=n)

        if center_point is not None:
            # Convert center point to [0, 1] scale
            center_scaled = np.array([(center_point[k] - v[0]) / (v[1] - v[0]) for k, v in FIELDS.items()])
            
            # Mix with center point (weighted average)
            mix_factor = 0.3  # How much to cluster around center_point
            sample = sample * (1 - mix_factor) + center_scaled * mix_factor
    
        # Scale to parameter bounds
        scaled_samples = qmc.scale(sample, l_bounds, u_bounds)

        # Convert to list of parameter dictionaries
        for point in scaled_samples:
            params = {}
            for i, name in enumerate(FIELDS.keys()):
                value = point[i]
                if name.startswith('taper_turns'):
                    value = round(value)
                params[name] = value
            
            # Only add if constraints are met
            if check_parameters(params, mesh_spacing, axial_length, feed_height_thres):
                param_sets.append(params)
            
            if len(param_sets) == n:
                break
    
    return param_sets

def new_population(n, params, temperature, mesh_spacing, axial_length, feed_height_thres=4):
    population = []

    while len(population) < n:
        new_params = params.copy()
        for f, (fmin, fmax) in FIELDS.items():
            # Temperature-scaled step size
            step_size = temperature * (fmax - fmin) * 0.1
            delta = np.random.normal(0, step_size)
            new_params[f] = np.clip(params[f] + delta, fmin, fmax)
            
            if f.startswith('taper_turns'):
                new_params[f] = round(new_params[f])
        
        if check_parameters(new_params, mesh_spacing, axial_length, feed_height_thres):
            population.append(new_params)
    
    return population

def acceptance_probability(old_cost, new_cost, temperature):
    if new_cost < old_cost:
        return 1.0
    return np.exp((old_cost - new_cost) / temperature)

# Helper function to prepare the multiprocessing inputs
def prepare_mp_arguments(base_dir, design_params, parameter_sets, scripts="./matlab"):
    return [{'base_dir': base_dir, 'design_params': design_params, 'params': params, 'scripts': scripts} for params in parameter_sets]

def print_result(params, result, score):
    print("\nParameters:")
    print(f"Pitch: {params['pitch']:.1f} mm")
    print(f"Feed Height: {params['feed_height']:.1f} mm")
    print(f"Reflector Radius: {params['reflector_radius']:.1f} mm")
    print(f"Taper Radii: {params['taper_radius1']:.1f}, {params['taper_radius2']:.1f}, {params['taper_radius3']:.1f} mm")
    print(f"Taper Turns: {params['taper_turns1']}, {params['taper_turns2']}, {params['taper_turns3']}")
    
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
    RADIUS_SPACING_THRES = 3 # in millimeters
    FEED_HEIGHT_THRES = 5 # in millimeters
    MAX_AXIAL_LENGTH = 2000 # in millimeters

    # Rough starting point
    INITIAL_PARAMS = {
        'pitch': 133,
        'feed_height': -8,
        'reflector_radius': 109,
        'taper_radius1': 120,
        'taper_radius2': 109,
        'taper_radius3': 98,
        'taper_turns1': 3,
        'taper_turns2': 4,
        'taper_turns3': 10
    }
    
    # Optimizer params
    INITIAL_TEMP = 1.0
    FINAL_TEMP = 0.01
    COOLING_RATE = 0.95
    SIM_PARALLEL = 40
    SIM_POPULATION = 32 # power of 2
    SIM_ITERATIONS = 20

    # Create the working directory if it does not yet exist
    SIM_DIR.mkdir(exist_ok=True)

    population = generate_points(SIM_POPULATION, RADIUS_SPACING_THRES, MAX_AXIAL_LENGTH, FEED_HEIGHT_THRES, center_point=INITIAL_PARAMS)

    pool = mp.Pool(processes=SIM_PARALLEL)

    # Store the absolute best costs
    overall_best_cost = float('inf')
    overall_best_params = INITIAL_PARAMS.copy()
    overall_best_results = None

    # Accpeted solution for simulated annealing
    anneal_best_cost = float('inf')
    anneal_best_params = INITIAL_PARAMS.copy()

    # Simulated annealing temperature
    temperature = INITIAL_TEMP
    iteration = 0

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
        
        population = new_population(SIM_POPULATION, anneal_best_params, temperature, RADIUS_SPACING_THRES, MAX_AXIAL_LENGTH, FEED_HEIGHT_THRES)

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
