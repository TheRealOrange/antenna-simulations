folder = 'hi';
params.pitch = 133;
params.feed_height = -8;
params.reflector_radius = 109;
params.taper_radius1 = 120;
params.taper_radius2 = 109;
params.taper_radius3 = 98;
params.taper_turns1 = 3;
params.taper_turns2 = 4;
params.taper_turns3 = 10;

f0 = 437.25e6;  % Target frequency
fc = 50e6;      % 20 dB corner frequency
R = 50;

simulate_helix(folder, params, f0, fc, R);