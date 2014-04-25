%% Params for Electrochemical Model
%   Created June 3, 2011 by Scott Moura

%% Geometric Params
% Thickness of each layer
p.L_n = 2.885e-3;     % Thickness of negative electrode [cm]
p.L_s = 1.697e-3;     % Thickness of separator [cm]
p.L_p = 6.521e-3;     % Thickness of positive electrode [cm]

% Particle Radii
p.R_s_n = 3.596e-4;   % Radius of solid particles in negative electrode [cm]
p.R_s_p = 1.637e-5;   % Radius of solid particles in positive electrode [cm]

% Volume fractions
p.epsilon_s_n = 0.3810;   % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.4800;   % Volume fraction in solid for pos. electrode

p.epsilon_e_n = 0.6190;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 0.3041;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.5200;   % Volume fraction in electrolyte for pos. electrode

% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [cm^2/cm^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [cm^2/cm^3]

%% Transport Params
% Diffusion coefficient in solid
p.D_s_n = 1.736e-10;  % Diffusion coeff for solid in neg. electrode, [cm^2/s]
p.D_s_p = 8.256e-10;  % Diffusion coeff for solid in pos. electrode, [cm^2/s]

% Diffusion coefficient in electrolyte
p.D_e = 6.911e-6;    % Diffusion coeff for electrolyte, [cm^2/s]

p.brug = 1.452;       % Bruggeman porosity
p.D_e_n = p.D_e * p.epsilon_e_n^p.brug; % Effective diffusion coef. in neg. electrode, [cm^2/s]
p.D_e_s = p.D_e * p.epsilon_e_s^p.brug; % Effective diffusion coef. in neg. electrode, [cm^2/s]
p.D_e_p = p.D_e * p.epsilon_e_p^p.brug; % Effective diffusion coef. in neg. electrode, [cm^2/s]

% Conductivity of solid
p.sig_n = 1;    % Conductivity of solid in neg. electrode, [1/Ohms*cm]
p.sig_p = 1;    % Conductivity of solid in pos. electrode, [1/Ohms*cm]

p.sig_eff_n = p.sig_n * p.epsilon_s_n;    % Eff. conductivity in neg. electrode, [1/Ohms*cm]
p.sig_eff_p = p.sig_p * p.epsilon_s_p;    % Eff. conductivity in pos. electrode, [1/Ohms*cm]

% Conductivity of electrolyte

% Miscellaneous
p.t_plus = 0.2495;    % Transference number
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
p.Area = 0.3108*1e4;  % Electrode current collector area [cm^2]

%% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]
p.alph = 0.5;         % Charge transfer coefficients
p.R_SEI = 3.391e1;   % Resistivity of SEI layer, [Ohms*cm^2]

% Reaction rates
p.k_n = 8.696e-11; % Reaction rate in neg. electrode, [(A/cm^2)*(mol^3/mol)^(1+alpha)]
p.k_p = 1.127e-11; % Reaction rate in pos. electrode, [(A/cm^2)*(mol^3/mol)^(1+alpha)]

%% Thermodynamic Params
% Equilibrium potentials (ws/Uref_LiFEPO4.mat)
load('ws/Upn_LiFePO4');
p.theta_vec = theta_vec;
p.Upn_vec = Upn_vec;

% Thermal dynamics
p.C_p = 75;   % Heat capacity, [J/K]
p.h = 12.4;   % Heat transfer coefficient, [W/K]

% Ambient Temperature
p.T_amp = 298.15; % [K]

% Entropy coefficients
p.dUref_dT = -0.4e-3; % [V/K] approx. from Al Hallaj et al 2000, JPS

%% Concentrations
% Maxima
p.c_s_n_max = 2.948e-2;    % Max concentration in anode, [mol/cm^3]
p.c_s_p_max = 1.035e-2;    % Max concentration in cathode, [mol/cm^3]

p.c_e = 1e-3;              % Fixed electrolyte concentration for SPM, [mol/cm^3]

%% Discretization parameters
% Discrete time step
p.delta_t = 0.1;

% Pade Order
p.PadeOrder = 3;

% Finite difference points along r-coordinate
p.Nr = 10;
p.delta_r_n = p.R_s_n / p.Nr;
p.delta_r_p = p.R_s_p / p.Nr;

% Finite difference points along x-coordinate
p.Nxn = 20;
p.Nxs = 20;
p.Nxp = 20;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

%% Bessel Function Regression Params
p.poly_I1 = [0.06994, -0.00641, 0.5016, -8.554e-5];
p.poly_I2 = [0.02194, 0.1105, 0.003257, -0.0001637];

