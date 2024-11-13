import os
import tempfile
from oct2py import Oct2Py
import atexit
import multiprocessing as mp
from datetime import datetime
import shutil
from pathlib import Path
import uuid
import logging


# Helper functions to change INPUTRC deal with this issue
# https://github.com/blink1073/oct2py/issues/308
# Create a temporary inputrc and reset it after the execution is done
def setup_octave_environment():
    fd, inputrc_path = tempfile.mkstemp(prefix='octave_', suffix='.inputrc')
    with os.fdopen(fd, 'w') as f:
        f.write("set enable-bracketed-paste off\n")
    
    original_env = {
        'INPUTRC': os.environ.get('INPUTRC'),
        'TERM': os.environ.get('TERM')
    }
    
    os.environ['INPUTRC'] = inputrc_path
    os.environ['TERM'] = 'dumb'
    
    return inputrc_path, original_env

# Function to reset environment to original state after execution
def cleanup_environment(inputrc_path, original_env):
    try:
        os.remove(inputrc_path)
    except:
        pass
    
    for key, value in original_env.items():
        if value is not None:
            os.environ[key] = value
        elif key in os.environ:
            del os.environ[key]

# Function to run a single instance of the simulation
def run_simulation(sim_params):
    base_dir = sim_params['base_dir']
    design_parameters = sim_params['design_params']
    params = sim_params['params']
    sim_uuid = str(uuid.uuid4())
    process_id = mp.current_process().name

    folder = base_dir / sim_uuid

    # Create simulation directory if it doesn't exist
    os.makedirs(folder, exist_ok=True)

    # Set up environment for this process
    inputrc_path, original_env = setup_octave_environment()
    atexit.register(lambda: cleanup_environment(inputrc_path, original_env))

    result = {'success': False, 'result': {}, 'error': '', 'folder': folder}
    
    try:
         # Create Oct2Py instance with minimal logging
        oct_logger = logging.getLogger(f'oct2py_{sim_uuid}')
        oct_logger.setLevel(logging.WARNING)
        octave = Oct2Py(logger=oct_logger, temp_dir=str(folder))
        octave.addpath('./matlab')

        print(f"{process_id}: Starting simulation in folder: {folder}")

        # Run simulation
        result = octave.simulate_helix(
            str(folder),
            params,
            design_parameters['f0'],
            design_parameters['fc'],
            design_parameters['R'],
            nout=1
        )
        
        print(f"{process_id}: Simulation completed successfully in folder: {folder}")
        result['success'] = True
        result['result'] = result
        
    except Exception as e:
        print(f"{process_id}: An error occurred in folder {folder}: {str(e)}")
        result['error'] = str(e)
    finally:
        octave.exit()
        # Clean up the simulation directory
        try:
            shutil.rmtree(folder)
            print(f"{process_id}: Cleaned up folder: {folder}")
        except Exception as e:
            print(f"{process_id}: Error cleaning up folder {folder}: {e}")
    
    return result

# Helper function to prepare the multiprocessing inputs
def prepare_mp_arguments(base_dir, design_params, parameter_sets):
    return [{'base_dir': base_dir, 'design_params': design_params, 'params': params} for params in parameter_sets]

if __name__ == "__main__":
    SIM_DIR = Path("./working")
    SIM_PARALLEL = 20
    DESIGN_PARAMS = {'f0': 437.25e6, 'fc': 50e6, 'R': 50}

    # Create the working directory if it does not yet exist
    SIM_DIR.mkdir(exist_ok=True)

    pool = mp.Pool(processes=SIM_PARALLEL)

    parameter_sets = [
        {
            'pitch': 133 + i*5,  # Vary pitch
            'feed_height': -8,
            'reflector_radius': 109,
            'taper_radius1': 120,
            'taper_radius2': 109,
            'taper_radius3': 98,
            'taper_turns1': 3,
            'taper_turns2': 4,
            'taper_turns3': 10
        }
        for i in range(20)  # Create 4 different parameter sets
    ]

    simulation_params = prepare_mp_arguments(SIM_DIR, DESIGN_PARAMS, parameter_sets)
    results = pool.map(run_simulation, simulation_params)

    processed_results = []
    for res in results:
        if res['success']:
            processed_results.append(res['result'])
        else:
            processed_results.append(None)
            print(f"Simulation failed: {res.get('error', 'Unknown error')}")

    pool.close()

    # Process results
    print("\nSimulation Results:")
    for i, result in enumerate(results):
        if result is not None:
            print(f"\nSimulation {i} (pitch={parameter_sets[i]['pitch']}):")
