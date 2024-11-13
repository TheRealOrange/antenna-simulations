import os
import tempfile

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