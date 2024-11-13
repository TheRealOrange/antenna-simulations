folder = 'hi';
params.wire_rad = 2;
params.reflector_length = 334.49;
params.feed_point = 164.55;
params.feed_gap = 10.0;
params.feed_length = 330.47;
params.director_lengths = [315.72, 312.75, 310.03, 307.54, 305.29, 303.24, 301.40, 299.74, 298.26, 296.94];
params.director_spacings = [51.42, 123.41, 147.41, 171.41, 191.98, 205.69, 215.97, 226.26, 236.54, 246.83];

f0 = 437.25e6;  % Target frequency
fc = 50e6;      % 20 dB corner frequency
R = 50;

simulate_yagi(folder, params, f0, fc, R);