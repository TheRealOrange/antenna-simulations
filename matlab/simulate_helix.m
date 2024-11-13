function result = simulate_helix(folder, params, f0, fc, R)
    % Optimization function for helical antenna with constrained taper turns
    
    disp([mfilename ': SIMULATING...']);
    % Load required signal processing package for findpeaks
    pkg('load', 'signal');

    % Run function with default parameters if none supplied
    if nargin == 0
        folder = 'tmp';
        params.pitch = 133;
        params.feed_height = -8;
        params.reflector_radius = 109;
        params.taper_radius1 = 120;
        params.taper_radius2 = 109;
        params.taper_radius3 = 98;
        params.taper_turns1 = 3;
        params.taper_turns2 = 4;
        params.taper_turns3 = 10;

        % Design parameters
        f0 = 437.25e6;  % Target frequency
        fc = 50e6;      % 20 dB corner frequency
        R = 50;
    end

    % Setup simulation
    physical_constants;
    unit = 1e-3;    % all length in mm

    % Calculate wavelength
    lambda0 = round(c0/f0/unit);

    % Create tapered sections
    helix.taper_turns1 = params.taper_turns1;
    helix.taper_turns2 = helix.taper_turns1+params.taper_turns2;
    helix.taper_turns3 = helix.taper_turns2+params.taper_turns3;
    helix.taper_radius = [params.taper_radius1  params.taper_radius2  params.taper_radius3];
    helix.taper_turns  = [helix.taper_turns1    helix.taper_turns2    helix.taper_turns3  ];
    helix.radius = max(helix.taper_radius);

    % Calculate helix dimensions
    helix.turns = max(helix.taper_turns);
    axial_length = helix.turns * params.pitch;
    helix.segments = 21;

    % Set the feed location to be from z=0 to z=feed_height
    feed.start = [helix.taper_radius(1) 0 0                 ];
    feed.stop  = [helix.taper_radius(1) 0 params.feed_height];

    % Grid bounds and resolution
    simBox.xy_max = params.reflector_radius + lambda0/4;
    simBox.z_min = -lambda0;
    simBox.z_max = axial_length+lambda0;

    disp([mfilename ': MESH SETUP']);

    %% Mesh setup
    % Calculate base cell size
    mesh.maxres   = floor(c0 / (f0+fc) / unit / 20);
    mesh.fineres  = mesh.maxres / 4;
    mesh.transres = mesh.maxres / 2;

    mesh.maxres_z   = min(params.pitch / 6 ,  mesh.maxres             );
    mesh.fineres_z  = min(mesh.maxres_z / 4,  abs(params.feed_height) );
    mesh.transres_z = min(mesh.maxres_z / 2,  mesh.fineres_z*2        );

    % Setup XY direction mesh to be finer around helix region
    bounding_radius = max(helix.radius, params.reflector_radius);
    mesh.x = SmoothMeshLines([-bounding_radius*1.3 -params.reflector_radius -helix.taper_radius 0 helix.taper_radius params.reflector_radius bounding_radius*1.3], mesh.fineres, 1.4);
    mesh.x = SmoothMeshLines([mesh.x  bounding_radius*1.3  bounding_radius*1.5], mesh.transres, 1.4);
    mesh.x = SmoothMeshLines([mesh.x -bounding_radius*1.3 -bounding_radius*1.5], mesh.transres, 1.4);
    mesh.x = SmoothMeshLines([mesh.x -simBox.xy_max simBox.xy_max], mesh.maxres, 1.4);
    % Copy x mesh to y axis
    mesh.y = mesh.x;

    % Z-direction mesh
    z_fine_end_pos = max(0, params.feed_height) + 1.6*mesh.fineres_z;
    z_fine_end_neg = min(0, params.feed_height) - 1.6*mesh.fineres_z;
    z_trans_end_pos = z_fine_end_pos + lambda0/4;
    z_trans_end_neg = z_fine_end_neg - lambda0/4;
    mesh.z = SmoothMeshLines([z_fine_end_neg 0 params.feed_height z_fine_end_pos], mesh.fineres_z);
    mesh.z = SmoothMeshLines([mesh.z  z_trans_end_neg z_fine_end_pos+mesh.fineres_z z_trans_end_pos], mesh.transres_z);
    mesh.z = SmoothMeshLines([mesh.z  z_trans_end_pos+mesh.transres_z axial_length], mesh.maxres_z);
    mesh.z = [mesh.z simBox.z_min simBox.z_max];
    mesh.z = SmoothMeshLines(mesh.z, mesh.maxres_z, 1.4);

    mesh = AddPML(mesh, [8 8 8 8 8 16]);  % add 8 lines in both z-directions

    %% Print some information about the antenna and the simulation
    disp( ['Helix Parameters']);
    disp( ['Axial Length     : '  num2str(axial_length) ]);
    disp( ['Turns            : '  num2str(helix.turns)  ]);
    disp( ['Pitch            : '  num2str(params.pitch) ]);
    disp( ['Taper Sections   : (' num2str(helix.taper_radius) ')']);
    disp( ['Taper Turns      : (' num2str(helix.taper_turns)  ')']);
    disp( ['Reflector Radius : '  num2str(params.reflector_radius)]);
    disp( ['Feed Height      : '  num2str(params.feed_height)]);
    disp('');
    disp( ['Simulation Parameters']);
    disp( ['Mesh Resolution (XYZ) : '  num2str(numel(mesh.x)) ' x ' num2str(numel(mesh.y)) ' x ' num2str(numel(mesh.z))]);
    disp( ['Number of Cells       : '  num2str(numel(mesh.x)*numel(mesh.y)*numel(mesh.z))]);
    disp( ['Bounds X              : (' num2str(min(mesh.x)) ', ' num2str(max(mesh.x)) ')']);
    disp( ['Bounds Y              : (' num2str(min(mesh.y)) ', ' num2str(max(mesh.y)) ')']);
    disp( ['Bounds Z              : (' num2str(min(mesh.z)) ', ' num2str(max(mesh.z)) ')']);
    disp( ['Cell Sizes:']);
    disp( ['  Maximum Cell Size   : ' num2str(mesh.maxres)   ' x ' num2str(mesh.maxres)   ' x ' num2str(mesh.maxres_z)  ]);
    disp( ['  Fine Cell Size      : ' num2str(mesh.fineres)  ' x ' num2str(mesh.fineres)  ' x ' num2str(mesh.fineres_z) ]);
    disp( ['  Transition Cell Size: ' num2str(mesh.transres) ' x ' num2str(mesh.transres) ' x ' num2str(mesh.transres_z)]);

    disp([mfilename ': SIMULATION SETUP']);

    %% setup FDTD parameter & excitation function
    FDTD = InitFDTD( );
    FDTD = SetGaussExcite( FDTD, f0, fc );
    BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_16'}; % boundary conditions
    FDTD = SetBoundaryCond( FDTD, BC );

    CSX = InitCSX();
    CSX = DefineRectGrid( CSX, unit, mesh );

    disp([mfilename ': ANTENNA GEOMETRY SETUP']);

    %% create helix using the wire primitive
    CSX = AddMetal( CSX, 'helix' ); % create a perfect electric conductor (PEC)

    ang = linspace(0, 2*pi, helix.segments);
    coil_z = ang/2/pi * params.pitch;

    helix.x=[];
    helix.y=[];
    helix.z=[];
    zpos = params.feed_height;
    section = 1;
    for n=1:helix.turns
        if (n==helix.taper_turns(section) && section != length(helix.taper_turns))
            % Helix taper change boundary, vary the radius smoothly between
            % the two different radii
            radius = linspace(helix.taper_radius(section),helix.taper_radius(section+1),helix.segments);
            section++;
        else
            % Within a section, constant radius
            radius = ones(1, helix.segments) * helix.taper_radius(section);
        end
        helix.x = [helix.x radius.*cos(ang)];
        helix.y = [helix.y radius.*sin(ang)];
        helix.z = [helix.z coil_z+zpos];
        zpos = zpos + params.pitch;
    end
    clear p
    p(1,:) = helix.x;
    p(2,:) = helix.y;
    p(3,:) = helix.z;
    CSX = AddCurve(CSX, 'helix', 0, p);

    %% create a circular loop reflector
    CSX = AddMetal( CSX, 'reflector' ); % create a perfect electric conductor (PEC)
    ang = linspace(0,2*pi,helix.segments);
    loop.x = [feed.start(1) params.reflector_radius*cos(ang)];
    loop.y = [feed.start(2) params.reflector_radius*sin(ang)];
    loop.z = [feed.start(3) zeros(1, helix.segments) ];
    clear l
    l(1,:) = loop.x;
    l(2,:) = loop.y;
    l(3,:) = loop.z;
    CSX = AddCurve(CSX, 'reflector', 0, l);

    disp([mfilename ': EXCITATION SOURCE SETUP']);

    %% apply the excitation & resist as a current source
    [CSX port] = AddLumpedPort(CSX, 50 , 1, R, feed.start, feed.stop, [0 0 sign(params.feed_height)], true);

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