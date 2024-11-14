import random

MIN_DIRECTORS = 4
MAX_DIRECTORS = 15

# Parameter bounds (as per the octave script)
FIELDS = {
    #'wire_rad': (2, 4), 
    'reflector_length': (100, 400), 
    'feed_point': (100, 300), 
    'feed_gap': (5, 25),
    'feed_length': (100, 400),
    'director_spacings': (40, 300),
    'director_lengths': (100, 400),
}

# Parameter sanity check
def check_parameters(params, max_length, element_spacing=25.0, min_res=1):
    for f, (fmin, fmax) in FIELDS.items():
        if not f in params:
            return False
        elif (not 'director' in f) and (params[f] > fmax or params[f] < fmin):
            return False
        
    # Different length of director spacings and lengths, invalid
    if len(params['director_spacings']) != len(params['director_lengths']):
        return False
    
    # Too many or too little directors
    if len(params['director_spacings']) > MAX_DIRECTORS or len(params['director_spacings']) < MIN_DIRECTORS:
        return False
    
    # Check that none of the feed, reflector or directors positions
    # are too close together (will cause meshing issues with OpenEMS)
    all_pos = [params['feed_point']] + params['director_spacings']
    for pos in all_pos:
        if pos < element_spacing:
            return False
        
    # Check that none of the feed, reflector or directors edge lengths
    # are too close together (will cause meshing issues with OpenEMS)
    all_len = set([params['reflector_length']/2] + [params['feed_length']/2] + [length/2 for length in params['director_lengths']])
    all_len = sorted(all_len)
    for r in range(len(all_len) - 1):
        if abs(all_len[r] - all_len[r+1]) < min_res:
            return False
        
    # Check that the generated antenna length is not too long
    total_length = sum(params['director_spacings']) + params['feed_point']
    if total_length > max_length:
        return False
    
    # Check that the feed height is not too small
    if params['feed_gap'] > params['feed_length']/3:
        return False
    
    return True

def adjust_to_min_difference(value, existing_values, min_diff=2.0):
    for existing in existing_values:
        diff = abs(value - existing)
        if 0 < diff < min_diff:
            # If closer to existing value, snap to it, otherwise move away by min_diff
            if diff < min_diff/2:
                return existing
            else:
                return existing + min_diff if value > existing else existing - min_diff
    return value

def random_configuration(max_length, n=1, min_res=1):
    configs = []
    min_diff = min_res*2
    while len(configs) < n:
        params = dict()
        all_lengths = []  # Keep track of all lengths for minimum difference checking
        
        # Generate non-director parameters
        for f, (fmin, fmax) in FIELDS.items():
            if not 'director' in f:
                raw_value = random.uniform(fmin, fmax)
                if 'length' in f:  # Only adjust length parameters
                    raw_value = adjust_to_min_difference(raw_value, all_lengths, min_diff)
                    all_lengths.append(raw_value)
                params[f] = raw_value
        
        # Generate director elements
        director_count = random.randint(MIN_DIRECTORS, MAX_DIRECTORS)
        minlen, maxlen = FIELDS['director_lengths']
        minspace, maxspace = FIELDS['director_spacings']
        
        params['director_lengths'] = []
        params['director_spacings'] = []
        
        for _ in range(director_count):
            # Generate length with minimum difference check
            raw_length = random.uniform(minlen, maxlen)
            adjusted_length = adjust_to_min_difference(raw_length, all_lengths, min_diff)
            params['director_lengths'].append(adjusted_length)
            all_lengths.append(adjusted_length)
            
            # Generate spacing (no adjustment needed)
            params['director_spacings'].append(random.uniform(minspace, maxspace))
        
        if check_parameters(params, max_length, min_res=min_res):
            configs.append(params)
    
    return configs