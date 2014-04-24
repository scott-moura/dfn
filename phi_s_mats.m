%% Matrices for Electric Potential in Solid Phase, phi_s(x,t)
%   Created May 8, 2012 by Scott Moura

function [F1n,F1p,F2n,F2p,Gn,Gp,varargout] = phi_s_mats(p)

% Conductivity and FD Length Coefficients
alpha_n = p.sig_n / (2 * p.L_n * p.delta_x_n);
alpha_p = p.sig_p / (2 * p.L_p * p.delta_x_p);

%% Block Matrices
% M1 : phi_s x
M1n = alpha_n * (diag(ones(p.Nxn-2,1),1) + diag(-1*ones(p.Nxn-2,1),-1));
M1p = alpha_p * (diag(ones(p.Nxp-2,1),1) + diag(-1*ones(p.Nxp-2,1),-1));

% M2 : phi_s z
M2n = zeros(p.Nxn-1,2);
M2n(1,1) = -alpha_n;
M2n(end,end) = alpha_n;

M2p = zeros(p.Nxp-1,2);
M2p(1,1) = -alpha_p;
M2p(end,end) = alpha_p;

% M3 : i_e
M3n = [zeros(p.Nxn-1,1), diag(-1*ones(p.Nxn-1,1)), zeros(p.Nxn-1,1)];
M3p = [zeros(p.Nxp-1,1), diag(-1*ones(p.Nxp-1,1)), zeros(p.Nxp-1,1)];

% N1 : phi_s x
N1n = zeros(2,p.Nxn-1);
N1n(1,1) = 2*alpha_n;
N1n(end,end) = -2*alpha_n;

N1p = zeros(2,p.Nxp-1);
N1p(1,1) = 2*alpha_p;
N1p(end,end) = -2*alpha_p;

% N2 : phi_s z
N2n = [-2*alpha_n, 0; 0, 2*alpha_n];
N2p = [-2*alpha_p, 0; 0, 2*alpha_p];

% N4 : I
N4n = [1; 0];
N4p = [0; 1];

%% Form F, G matrices
F1n = M1n - M2n*(N2n\N1n);
F2n = M3n;
Gn = 1 - M2n*(N2n\N4n);

F1p = M1p - M2p*(N2p\N1p);
F2p = M3p;
Gp = 1 - M2p*(N2p\N4p);

%% Compute C,D matrices for boundary values
Cn = -(N2n\N1n);
Dn = -(N2n\N4n);

Cp = -(N2p\N1p);
Dp = -(N2p\N4p);

varargout{1} = Cn;
varargout{2} = Cp;
varargout{3} = Dn;
varargout{4} = Dp;
