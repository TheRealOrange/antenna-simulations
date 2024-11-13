# Parameter bounds (as per the octave script)
FIELDS = {
    'pitch': (50, 200), 
    'feed_height': (-50, 50), 
    'reflector_radius': (50, 160), 
    'taper_radius1': (50, 160),
    'taper_radius2': (50, 160),
    'taper_radius3': (50, 160),
    'taper_turns1': (1, 20),
    'taper_turns2': (1, 20),
    'taper_turns3': (1, 20)
}

# Parameter sanity check
def check_parameters(params, spacing, max_axial_length, min_feed_height=4):
    for f, (fmin, fmax) in FIELDS.items():
        if not f in params:
            return False
        elif params[f] > fmax or params[f] < fmin:
            return False
    
    # Check that the reflector radius is not the same
    # as the starting taper radius (the reflector loop will collide with the helix)
    if params['reflector_radius'] == params['taper_radius1']:
        return False
    
    # Check that none of the taper radii and reflector
    # radii are too close together (will cause meshing issues with OpenEMS)
    all_radii = set([params['reflector_radius'], params['taper_radius1'], params['taper_radius2'], params['taper_radius3']])
    all_radii = sorted(all_radii)
    for r in range(len(all_radii) - 1):
        if abs(all_radii[r] - all_radii[r+1]) < spacing:
            return False
        
    # Check that the generated antenna length is not too long
    total_turns = params['taper_turns1'] + params['taper_turns2'] + params['taper_turns3']
    axial_length = total_turns * params['pitch']
    if axial_length > max_axial_length:
        return False
    
    # Check that the feed height is not too small
    if abs(params['feed_height']) < min_feed_height:
        return False
    
    return True