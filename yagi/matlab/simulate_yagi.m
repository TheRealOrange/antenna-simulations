function result = simulate_yagi(folder, params, f0, fc, R)
    % Optimization function for helical antenna with constrained taper turns

    disp([mfilename ': SIMULATING...']);
    % Load required signal processing package for findpeaks
    pkg('load', 'signal');

    % Run function with default parameters if none supplied
    if nargin == 0
        folder = 'tmp';
        params.wire_rad = 2;
        params.reflector_length = 300;
        params.feed_point = 150;
        params.feed_gap = 5;
        params.feed_length = 140;
        params.director_spacings = [80 80 80];
        params.director_lengths = [130 130 130];

        % Design parameters
        f0 = 437.25e6;  % Target frequency
        fc = 50e6;      % 20 dB corner frequency
        R = 50;
    end

    % Setup simulation
    physical_constants;
    unit = 1e-3;    % all length in mm
    wire_rad = params.wire_rad;

    % Calculate wavelength
    lambda0 = round(c0/f0/unit);
    disp(num2str(lambda0));

    % Create reflector section
    yagi.reflector.length = params.reflector_length;

    % Create feed section
    yagi.feed.position = params.feed_point;
    yagi.feed.gap = params.feed_gap;
    yagi.feed.length = params.feed_length;

    % Create director sections
    yagi.director.positions = [];
    for space=params.director_spacings
      if length(yagi.director.positions) == 0
        yagi.director.positions = [space];
      else
        yagi.director.positions = [yagi.director.positions space+yagi.director.positions(end)];
      end
    end
    yagi.director.positions += yagi.feed.position;
    yagi.director.lengths = params.director_lengths;
    yagi.total_length = max(yagi.director.positions);

    % Grid bounds and resolution
    simBox.xy_max = max([yagi.reflector.length yagi.feed.length yagi.director.lengths]) + lambda0/4;
    simBox.z_min = -lambda0/2;
    simBox.z_max = yagi.total_length+lambda0;

    disp([mfilename ': MESH SETUP']);

    % Calculate base cell size
    mesh.maxres   = floor(c0 / (f0+fc) / unit / 40);
    mesh.fineres  = params.wire_rad / 2;
    mesh.transres = mesh.fineres * 2;

    %% Y axis mesh
    % Y mesh for wires
    mesh.y = -wire_rad*2:mesh.fineres:wire_rad*2;
    mesh.y = [mesh.y -wire_rad*2-10:2:-wire_rad*2];
    % Y axis simulation bounds
    mesh.y = SmoothMeshLines([mesh.y -simBox.xy_max simBox.xy_max], mesh.maxres, 1.4);

    %% Z axis mesh
    % Z mesh for reflector element
    reflector_mesh = -wire_rad*2:mesh.fineres:wire_rad*2;
    % Z mesh for feed element
    feed_mesh = [-wire_rad*2+yagi.feed.position:mesh.fineres:wire_rad*2+yagi.feed.position yagi.feed.position];
    % Add remaining Z mesh lines for director elements
    director_meshes = [];
    for wires=yagi.director.positions
        curr_director_mesh = -wire_rad*2+wires:mesh.fineres:wire_rad*2+wires;
        director_meshes = [director_meshes curr_director_mesh];
    end
    % Z axis simulation bounds
    mesh.z = SmoothMeshLines([reflector_mesh feed_mesh director_meshes simBox.z_min simBox.z_max], mesh.maxres, 1.3);

    %% X axis mesh
    % X axis mesh for feed
    mesh.x = SmoothMeshLines([-yagi.feed.gap/2 0 yagi.feed.gap/2], mesh.fineres);
    mesh.x = SmoothMeshLines([mesh.x  yagi.feed.gap*2  yagi.director.lengths/2  yagi.reflector.length/2  yagi.feed.length/2], mesh.transres, 1.4);
    mesh.x = SmoothMeshLines([mesh.x -yagi.feed.gap*2 -yagi.director.lengths/2 -yagi.reflector.length/2 -yagi.feed.length/2], mesh.transres, 1.4);
    % X axis simulation bounds
    mesh.x = SmoothMeshLines([mesh.x -simBox.xy_max simBox.xy_max yagi.director.lengths/2 -yagi.director.lengths/2], mesh.maxres, 1.4);


    mesh = AddPML(mesh, [8 8 8 8 8 16]);  % add PML lines in both z-directions

    %% Print some information about the antenna and the simulation
    disp( ['Yagi-Uda Parameters']);
    disp( ['Total Length       : '  num2str(yagi.total_length)   ]);
    disp( ['Wire Radius        : '  num2str(wire_rad)            ]);
    disp( ['Feed Position      : '  num2str(yagi.feed.position)  ]);
    disp( ['Feed Gap           : '  num2str(yagi.feed.gap)       ]);
    disp( ['Feed Length        : '  num2str(yagi.feed.length)    ]);
    disp( ['Director Positions : (' num2str(yagi.director.positions) ')']);
    disp( ['Director Lengths   : (' num2str(yagi.director.lengths)   ')']);
    disp( ['Reflector Length   : '  num2str(yagi.reflector.length)]);
    disp('');
    disp( ['Simulation Parameters']);
    disp( ['Mesh Resolution (XYZ) : '  num2str(numel(mesh.x)) ' x ' num2str(numel(mesh.y)) ' x ' num2str(numel(mesh.z))]);
    disp( ['Number of Cells       : '  num2str(numel(mesh.x)*numel(mesh.y)*numel(mesh.z))]);
    disp( ['Bounds X              : (' num2str(min(mesh.x)) ', ' num2str(max(mesh.x)) ')']);
    disp( ['Bounds Y              : (' num2str(min(mesh.y)) ', ' num2str(max(mesh.y)) ')']);
    disp( ['Bounds Z              : (' num2str(min(mesh.z)) ', ' num2str(max(mesh.z)) ')']);
    disp( ['Cell Sizes:']);
    disp( ['  Maximum Cell Size   : ' num2str(mesh.maxres)   ]);
    disp( ['  Fine Cell Size      : ' num2str(mesh.fineres)  ]);
    disp( ['  Transition Cell Size: ' num2str(mesh.transres) ]);

    disp([mfilename ': SIMULATION SETUP']);

    %% setup FDTD parameter & excitation function
    FDTD = InitFDTD( );
    FDTD = SetGaussExcite( FDTD, f0, fc );
    BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_16'}; % boundary conditions
    FDTD = SetBoundaryCond( FDTD, BC );

    CSX = InitCSX();
    CSX = DefineRectGrid( CSX, unit, mesh );

    disp([mfilename ': ANTENNA GEOMETRY SETUP']);

    % Add materials
    CSX = AddMaterial( CSX, 'copper' );
    CSX = SetMaterialProperty(CSX,'copper','Kappa',56e6);
    %% create the central carbon fiber support
    CSX = AddMaterial( CSX, 'support' ); % create carbon fiber material
    CSX = SetMaterialProperty(CSX,'support','Kappa',3.3e2, 'Epsilon', 5000);

    %% Create the reflector wire
    start = [-yagi.reflector.length/2 0 0];
    stop  = [ yagi.reflector.length/2 0 0];
    CSX = AddWire(CSX, 'copper', 10, [start' stop'], wire_rad);

    %% Create the feed wires (two wires)
    % Wire -X
    start = [-yagi.feed.length/2 0 yagi.feed.position];
    stop  = [-yagi.feed.gap/2    0 yagi.feed.position];
    CSX = AddWire(CSX, 'copper', 0, [start' stop'], wire_rad);
    % Wire +X
    start = [ yagi.feed.length/2 0 yagi.feed.position];
    stop  = [ yagi.feed.gap/2    0 yagi.feed.position];
    CSX = AddWire(CSX, 'copper', 10, [start' stop'], wire_rad);

    %% Create the director wires
    for n=1:length(yagi.director.positions)
      start = [-yagi.director.lengths(n)/2 0 yagi.director.positions(n)];
      stop  = [ yagi.director.lengths(n)/2 0 yagi.director.positions(n)];
      CSX = AddWire(CSX, 'copper', 10, [start' stop'], wire_rad);
    end

    %% Create the central support
    start = [0 -wire_rad*2-5 0                 ];
    stop  = [0 -wire_rad*2-5 yagi.total_length ];
    CSX=AddCylindricalShell(CSX, 'support', 5, start, stop, 4.5, 2);

    disp([mfilename ': EXCITATION SOURCE SETUP']);

    %% apply the excitation & resist as a current source
    start = [-yagi.feed.gap/2 0 yagi.feed.position];
    stop  = [ yagi.feed.gap/2 0 yagi.feed.position];
    [CSX port] = AddLumpedPort(CSX, 50 , 1, R, start, stop, [1 0 0], true);

    disp([mfilename ': NF2FF BOX SETUP']);

    start = [mesh.x(11)      mesh.y(11)     mesh.z(11)];
    stop  = [mesh.x(end-10) mesh.y(end-10) mesh.z(end-10)];
    [CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'OptResolution', lambda0/15);

    disp([mfilename ': RUNNING SIMULATION']);

    % Write and run simulation
    WriteOpenEMS([folder '/helix.xml'], FDTD, CSX);
    RunOpenEMS(folder, 'helix.xml', '--engine=multithreaded --numThreads=11');

    %% Process frequency domain results
    freq = linspace(f0-fc, f0+fc, 501);  % Increased number of points for better resolution
    port = calcPort(port, folder, freq);

    % get accepted antenna power at design frequency
    result.p_in = interp1(freq, port.P_acc, f0);

    % Calculate S-parameters and impedance
    Zin = port.uf.tot ./ port.if.tot;
    s11 = port.uf.ref ./ port.uf.inc;
    s11_dB = 20 * log10(abs(s11));

    % Calculate impedance and S11 at design frequency
    Zin_f0 = interp1(freq, Zin, f0);
    s11_f0 = interp1(freq, s11, f0);
    % Return the real and imaginary components of impedance at design frequency
    result.Zin_real  = real(Zin_f0);
    result.Zin_imag  = imag(Zin_f0);
    % Return the S11 magnitude and phase at the design frequency
    result.s11_mag   = abs(s11_f0);
    result.s11_phase = angle(s11_f0);

    %% Resonance Analysis
    [minima_values, minima_indices] = findpeaks(-s11_dB);
    minima_values = -minima_values;

    if isempty(minima_values)
        % No resonant frequencies found within frequency range
        % Skip calculating the s11 at resonance
        warning('No resonance points found in frequency range');
        result.f_resonance = NaN;
        result.s11_resonance = NaN;
    else
        minima_freqs = freq(minima_indices);
        [~, closest_idx] = min(abs(minima_freqs - f0));

        % Return the resonance frequency closest to the design frequency
        result.f_resonance = minima_freqs(closest_idx);
        result.s11_resonance = minima_values(closest_idx);

        % Calculate the power and directivity at the resonance frequency
        nf2ff_res = calculate_nf2ff(folder, result.f_resonance, nf2ff, port, freq);

        % Save the power and directivity
        result.rad_power_res   = nf2ff_res.Prad;
        result.directivity_res = nf2ff_res.Dmax;
        result.efficiency_res  = nf2ff_res.efficiency;
        result.HPBW_res        = nf2ff_res.theta_HPBW;
    end

    % Calculate the power and directivity at the design frequency f0
    nf2ff_f0 = calculate_nf2ff(folder, f0, nf2ff, port, freq);

    % Save the power and directivity
    result.rad_power   = nf2ff_f0.Prad;
    result.directivity = nf2ff_f0.Dmax;
    result.efficiency  = nf2ff_f0.efficiency;
    result.HPBW        = nf2ff_f0.theta_HPBW;

    disp( ['radiated power: Prad = ' num2str(result.rad_power) ' Watt']);
    disp( ['directivity: Dmax    = ' num2str(result.directivity) ' (' num2str(10*log10(result.directivity)) ' dBi)'] );
    disp( ['efficiency: nu_rad   = ' num2str(100*result.efficiency) ' %']);
    disp( ['theta_HPBW           = ' num2str(result.HPBW) ' °']);
end

function nf2ff_res = calculate_nf2ff(folder, f_calc, nf2ff_inp, port_inp, freq_inp)
    % calculate the far field at phi=0 degrees and at phi=90 degrees
    thetaRange = unique([0:0.5:90 90:180]);
    phiRange = (0:2:360) - 180;

    nf2ff_res = CalcNF2FF(nf2ff_inp, folder, f_calc, thetaRange*pi/180, phiRange*pi/180,'Mode',1,'Outfile','3D_Pattern.h5','Verbose',1);

    calc_p_in = interp1(freq_inp, port_inp.P_acc, f_calc);

    nf2ff_res.theta_HPBW = interp1(nf2ff_res.E_norm{1}(:,1)/max(nf2ff_res.E_norm{1}(:,1)),thetaRange,1/sqrt(2))*2;
    nf2ff_res.efficiency = nf2ff_res.Prad./calc_p_in;
end