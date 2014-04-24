%% DAEs for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura

function [f, g, varargout] = dae_dfn(x,z,Cur,p)


%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;

% Solid Concentration
c_s_n = x(1:Ncsn);
c_s_p = x(Ncsn+1:Ncsn+Ncsp);

% Reformat into matrices
c_s_n_mat = reshape(c_s_n,p.PadeOrder,p.Nxn-1);
c_s_p_mat = reshape(c_s_p,p.PadeOrder,p.Nxp-1);

% Electrolyte concentration
c_e = x((Ncsn+Ncsp + 1):(Nc));

% Temperature
T = x(end);

% Solid Potential
phi_s_n = z(1:Nn);
phi_s_p = z(Nn+1:Nnp);

% Electrolyte Current
i_en = z(Nnp+1 : Nnp+Nn);
i_ep = z(Nnp+Nn+1 : 2*Nnp);

% Electrolyte Potential
phi_e = z(2*Nnp+1:2*Nnp+Nx);

% Molar ionic flux
jn = z(2*Nnp+Nx+1 : 2*Nnp+Nx+Nn);
jp = z(2*Nnp+Nx+Nn+1 : end);

% Temperature
T = x(end);

%% Li Diffusion in Solid Phase: c_s(x,r,t)
% Preallocate matrices for derivatives
c_s_n_dot_mat = zeros(size(c_s_n_mat));
c_s_p_dot_mat = zeros(size(c_s_p_mat));
c_ss_n = zeros(Nn,1);
c_ss_p = zeros(Np,1);
c_avg_n = zeros(Nn,1);
c_avg_p = zeros(Np,1);

% Loop through each "comb tooth" in anode
for idx = 1:Nn
    c_s_n_dot_mat(:,idx) = p.A_csn * c_s_n_mat(:,idx) + p.B_csn * jn(idx);
    y_csn = p.C_csn * c_s_n_mat(:,idx);
    c_ss_n(idx) = y_csn(1);
    c_avg_n(idx) = y_csn(2);
end
c_s_n_dot = reshape(c_s_n_dot_mat,numel(c_s_n_dot_mat),1);

% Loop through each "comb tooth" in cathode
for idx = 1:Np
    c_s_p_dot_mat(:,idx) = p.A_csp * c_s_p_mat(:,idx) + p.B_csp * jp(idx);
    y_csp = p.C_csp * c_s_p_mat(:,idx);
    c_ss_p(idx) = y_csp(1);
    c_avg_p(idx) = y_csp(2);
end
c_s_p_dot = reshape(c_s_p_dot_mat,numel(c_s_p_dot_mat),1);

%% Li Diffusion in Electrolyte Phase: c_e(x,t)
% Electrolyte current across all three regions
i_ex = [0; i_en; Cur*ones(p.Nxs+1,1); i_ep; 0];

c_e_bcs = p.C_ce*c_e;
c_ex = [c_e_bcs(1); c_e(1:Nn); c_e_bcs(2); c_e(Nn+1:Nn+p.Nxs-1); ...
        c_e_bcs(3); c_e(Nn+p.Nxs : end); c_e_bcs(4)];

% System Matrices
[A_ce, B_ce, ~] = c_e_mats(p,c_ex);
    
% Compute derivative
c_e_dot = A_ce*c_e + B_ce*i_ex;

%% Potential in Solid Phase: phi_s(x,t)
% Algebraic eqns (semi-explicit form)
i_enn = [0; i_en; Cur];
i_epp = [Cur; i_ep; 0];

phi_sn_dot = p.F1_psn*phi_s_n + p.F2_psn*i_enn + p.G_psn*Cur;
phi_sp_dot = p.F1_psp*phi_s_p + p.F2_psp*i_epp + p.G_psp*Cur;

% Terminal Voltage
phi_s_n_bcs = p.C_psn * phi_s_n + p.D_psn * Cur;
phi_s_p_bcs = p.C_psp * phi_s_p + p.D_psp * Cur;

Volt = phi_s_p_bcs(2) - phi_s_n_bcs(1);

%% Electrolyte Current: i_e(x,t)
i_en_dot = p.F1_ien*i_en + p.F2_ien*jn + p.F3_ien*Cur;
i_ep_dot = p.F1_iep*i_ep + p.F2_iep*jp + p.F3_iep*Cur;

%% Potential in Electrolyte Phase: phi_e(x,t)
% System matrices
[F1_pe,F2_pe,F3_pe,C1_pe,C2_pe,C3_pe] = phi_e_mats(p,c_ex);

% Algebraic eqns (semi-explicit form)
phi_e_dot = F1_pe*phi_e + F2_pe*i_ex + F3_pe*log(c_ex);

%% Butler-Volmer Equation
aFRT = (p.alph*p.Faraday)/(p.R*p.T_amp);

% Exchange Current Density, i_0^{\pm}
[i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e);

% Equilibrium Potential, U^{\pm}(c_ss)
Unref = refPotentialAnode(p, c_ss_n / p.c_s_n_max);
Upref = refPotentialCathode(p, c_ss_p / p.c_s_p_max);

% Overpotential, \eta
eta_n = phi_s_n - phi_e(1:Nn) - Unref - p.Faraday*p.R_f_n*jn;
eta_p = phi_s_p - phi_e(end-Np+1:end) - Upref - p.Faraday*p.R_f_p*jp;

% Algebraic eqns (semi-explicit form)
jn_dot = 2/p.Faraday * i_0n .* sinh(aFRT * eta_n) - jn;
jp_dot = 2/p.Faraday * i_0p .* sinh(aFRT * eta_p) - jp;

%% Temperature
% Equilibrium Potential and Gradient wrt bulk concentration
[Unb,~,dUnbdT] = refPotentialAnode(p, c_avg_n / p.c_s_n_max);
[Upb,~,dUpbdT] = refPotentialCathode(p, c_avg_p / p.c_s_p_max);

% Heat generated from intercalation (w/o boundaries for NOW)
Q_nx = p.a_s_n*p.Faraday * jn .* (Unb - T*dUnbdT);
Q_n = sum(Q_nx) * p.delta_x_n * p.L_n;

Q_px = p.a_s_p*p.Faraday * jp .* (Upb - T*dUpbdT);
Q_p = sum(Q_px) * p.delta_x_p * p.L_p;

Q_inter = Q_n + Q_p;

% Temperature ODE
T_dot = (p.h*(p.T_amp - T) - Cur*Volt - Q_inter) / (p.rho_avg*p.C_p);

%% Concatenate Time Derivatives
f = [c_s_n_dot; c_s_p_dot; c_e_dot; T_dot];
g = [phi_sn_dot; phi_sp_dot; i_en_dot; i_ep_dot; ...
     phi_e_dot; jn_dot; jp_dot];
     
%% Optional Output Vars

% Conservation of Li-ion matters
nLi = sum(c_avg_n) * p.delta_x_n * Nn ...
     + sum(c_avg_p) * p.delta_x_p * Np ...
     + sum(c_e(1:Nn)) * p.delta_x_n * Nn ...
     + sum(c_e(Nn+1:end-Np)) * p.delta_x_s * p.Nxs-1 ...
     + sum(c_e(end-Np+1:end)) * p.delta_x_p * Np;
 
nLidot = sum(jn) * p.delta_x_n * Nn + sum(jp) * p.delta_x_p * Np;

c_e0n = c_ex(1);
c_e0p = c_ex(end);

phi_e_bcs = C1_pe*phi_e + C2_pe*i_ex + C3_pe*log(c_ex);
eta_s_Ln = phi_s_n_bcs(2) - phi_e_bcs(2);

% Aggregate Outputs
varargout{1} = [c_ss_n; c_ss_p; c_avg_n; c_avg_p; c_ex; eta_n; eta_p;...
    c_e0n; c_e0p; eta_s_Ln;...
    Volt; nLi; nLidot];
