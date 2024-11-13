folder = 'hi';
params.wire_rad = 3;
params.reflector_length = 350;
params.feed_point = 254;
params.feed_gap = 22;
params.feed_length = 259;
params.director_lengths = [134.1856492532784, 230.46467311914367];
params.director_spacings = [204.415302940417, 296.41089357929064];

f0 = 437.25e6;  % Target frequency
fc = 50e6;      % 20 dB corner frequency
R = 50;

simulate_yagi(folder, params, f0, fc, R);