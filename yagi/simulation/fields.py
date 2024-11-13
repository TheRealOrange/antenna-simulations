import random

MIN_DIRECTORS = 4
MAX_DIRECTORS = 15

# Parameter bounds (as per the octave script)
FIELDS = {
    'wire_rad': (2, 4), 
    'reflector_length': (100, 400), 
    'feed_point': (100, 300), 
    'feed_gap': (5, 25),
    'feed_length': (100, 400),
    'director_spacings': (40, 300),
    'director_lengths': (100, 400),
}

# Parameter sanity check
def check_parameters(params, max_length):
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
        if pos < params['wire_rad']*7:
            return False
        
    # Check that none of the feed, reflector or directors edge lengths
    # are too close together (will cause meshing issues with OpenEMS)
    all_len = set([params['reflector_length']/2] + [params['feed_length']/2] + [length/2 for length in params['director_lengths']])
    all_len = sorted(all_len)
    for r in range(len(all_len) - 1):
        if abs(all_len[r] - all_len[r+1]) < 0.5:
            return False
        
    # Check that the generated antenna length is not too long
    total_length = sum(params['director_spacings']) + params['feed_point']
    if total_length > max_length:
        return False
    
    # Check that the feed height is not too small
    if params['feed_gap'] > params['feed_length']/3:
        return False
    
    return True

def random_configuration(max_length, n=1):
    configs = []
    while len(configs) < n:
        params = dict()
        for f, (fmin, fmax) in FIELDS.items():
            if not 'director' in f:
                # For director elements, generate differently
                params[f] = random.random() * (fmax - fmin) + fmin
        
        # Decide on how many director elements to place
        director_count = random.randint(MIN_DIRECTORS, MAX_DIRECTORS)
        minlen, maxlen = FIELDS['director_lengths']
        minspace, maxspace = FIELDS['director_spacings']
        params['director_lengths'] = []
        params['director_spacings'] = []
        # For each director element, generate a spacing and a length
        for _ in range(director_count):
            params['director_lengths'].append(random.random() * (maxlen - minlen) + minlen)
            params['director_spacings'].append(random.random() * (maxspace - minspace) + minspace)

        if check_parameters(params, max_length):
            configs.append(params)
    return configs