%% Electrolyte Properties from DUALFOIL
%   Created July 5, 2012 by Scott Moura

%   nprop = 11; LiPF6 in EC:DMC (liquid)

c_e = linspace(500,1500,101);

% Diffusion coefficient in electrolyte
D_e = 5.34e-10*exp(-0.65*c_e/1000.0);

% Electrolyte conductivity
kappa = 0.0911+1.9101*c_e/1e3 - 1.052*(c_e/1e3).^2 + 0.1554*(c_e/1e3).^3;

dkappa = 1.9101/1e3 - 2*1.052*c_e/1e3/1e3 + 0.1554*3*(c_e/1e3)^2/1e3;