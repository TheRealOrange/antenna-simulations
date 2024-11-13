import os
import multiprocessing as mp
from simulation import util
from oct2py import Oct2Py
import atexit
import shutil
import uuid
import logging

class MeshCheckHandler(logging.Handler):
    def __init__(self):
        super().__init__()
        self.logs = []
        self.mesh_check_found = False

    def emit(self, record):
        msg = self.format(record)
        self.logs.append(msg)
        
        # If we find "check your mesh" for the first time
        if not self.mesh_check_found and "check your mesh" in msg.lower():
            self.mesh_check_found = True
            # Print all accumulated logs immediately
            print("\n=== Previous logs ===")
            for log in self.logs:
                print(log)
            print("=" * 50)
            # Clear the logs since we've printed them
            self.logs = []
        # After finding "check your mesh", print all subsequent logs immediately
        elif self.mesh_check_found:
            print(msg)

# Function to run a single instance of the simulation
def run_simulation(sim_params):
    base_dir = sim_params['base_dir']
    scripts_dir = sim_params['scripts']
    design_parameters = sim_params['design_params']
    params = sim_params['params']
    sim_uuid = str(uuid.uuid4())
    process_id = mp.current_process().name
    folder = base_dir / sim_uuid
    
    # Create simulation directory if it doesn't exist
    os.makedirs(folder, exist_ok=True)
    
    # Set up environment for this process
    inputrc_path, original_env = util.setup_octave_environment()
    atexit.register(lambda: util.cleanup_environment(inputrc_path, original_env))
    
    result = {'success': False, 'result': {}, 'error': '', 'folder': folder}
    
    try:
        # Create Oct2Py instance with immediate logging
        oct_logger = logging.getLogger(f'oct2py_{sim_uuid}')
        oct_logger.setLevel(logging.INFO)
        
        # Create and configure handler
        handler = MeshCheckHandler()
        formatter = logging.Formatter('%(message)s')
        handler.setFormatter(formatter)
        oct_logger.addHandler(handler)
        
        octave = Oct2Py(logger=oct_logger, temp_dir=str(folder))
        octave.addpath(str(scripts_dir))
        
        # Run simulation
        out = octave.simulate_yagi(
            str(folder),
            params,
            design_parameters['f0'],
            design_parameters['fc'],
            design_parameters['R'],
            nout=1
        )
        
        result['success'] = True
        result['result'] = out
        
    except Exception as e:
        print(f"{process_id}: An error occurred in folder {folder}: {str(e)}")
        result['error'] = str(e)
        
    finally:
        octave.exit()
        # Clean up the simulation directory
        try:
            shutil.rmtree(folder)
        except Exception as e:
            print(f"{process_id}: Error cleaning up folder {folder}: {e}")
    
    return params, result