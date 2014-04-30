%% Jacobian w.r.t. Parameters for Doyle-Fuller-Newman Model
%   Created Feb 18, 2014 by Scott Moura
%
%   INPUTS
%   x   : States      c_s_n, c_s_p, c_e, T
%   z   : Alg. vars   phi_s_n, phi_s_p, i_en, i_ep, phi_e, jn, jp
%   Cur : Applied Current
%   p   : Model parameter structure
%
%   OUTPUTS
%   JF = \partial F / \partial \theta    Jacobian of x ODEs w.r.t. theta
%   JFb = \partial Fb / \partial \theta    Jacobian of x BCs w.r.t. theta
%   JG = \partial G / \partial \theta    Jacobian of z eqns w.r.t. theta
%   JGb = \partial Gb / \partial \theta    Jacobian of z BCs w.r.t. theta
%
%   UNCERTAIN PARAMETERS, theta
%   1  : D_s_n
%   2  : D_s_p
%   3  : D_e_n
%   4  : D_e_s
%   5  : D_e_p
%   6  : (1-t_plus)
%   7  : 1/sig_n
%   8  : 1/sig_p
%   9  : 1/kappa
%   10 : (1 + d ln f_ca / d ln c_e)
%   11 : k_n
%   12 : k_p
%   13 : R_f_n
%   14 : R_f_p
%   15 : epsilon_e_n
%   16 : epsilon_e_s
%   17 : epsilon_e_p
%   18 : c_s_n_max
%   19 : c_s_p_max
%   20 : h
%   21 : 1/(rho_avg * C_p)

function [JF, JFb, JG, JGb] = jac_p_dfn_federico(x,z,Cur,p)


%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Ns = p.Nxs - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx;
Nt = 21;    % Number of uncertain params

ind_csn = 1:Ncsn;
ind_csp = Ncsn+1:Ncsn+Ncsp;

ind_cs = 1:Ncsn+Ncsp;
ind_ce = Ncsn+Ncsp+1:Nc;

ind_T = Nc+1;

ind_phi_s_n = 1:Nn;
ind_phi_s_p = Nn+1:Nnp;

ind_ien = Nnp+1:Nnp+Nn;
ind_iep = Nnp+Nn+1:2*Nnp;

ind_phi_e = 2*Nnp+1 : 2*Nnp+Nx;

ind_jn = 2*Nnp+Nx+1 : 2*Nnp+Nx+Nn;
ind_jp = 2*Nnp+Nx+Nn+1 : Nz;

% Solid Concentration
c_s_n = x(1:Ncsn);
c_s_p = x(Ncsn+1:Ncsn+Ncsp);

% Reformat into matrices
c_s_n_mat = reshape(c_s_n,p.PadeOrder,p.Nxn-1);
c_s_p_mat = reshape(c_s_p,p.PadeOrder,p.Nxp-1);

c_ss_n = zeros(Nn,1);
c_ss_p = zeros(Np,1);
c_avg_n = zeros(Nn,1);
c_avg_p = zeros(Np,1);

% Loop through each "comb tooth" in anode
for idx = 1:Nn
    y_csn = p.C_csn * c_s_n_mat(:,idx);
    c_ss_n(idx) = y_csn(1);
    c_avg_n(idx) = y_csn(2);
end

% Loop through each "comb tooth" in cathode
for idx = 1:Np
    y_csp = p.C_csp * c_s_p_mat(:,idx);
    c_ss_p(idx) = y_csp(1);
    c_avg_p(idx) = y_csp(2);
end

% Electrolyte concentration
c_e = x(ind_ce);
c_en = c_e(1:Nn);
c_es = c_e(Nn+1:Nn+p.Nxs-1);
c_ep = c_e(end-Np+1:end);

% Temperature
T = x(end);

% Solid Potential
phi_s_n = z(ind_phi_s_n);
phi_s_p = z(ind_phi_s_p);

% Terminal Voltage
phi_s_n_bcs = p.C_psn * phi_s_n + p.D_psn * Cur;
phi_s_p_bcs = p.C_psp * phi_s_p + p.D_psp * Cur;
Volt = phi_s_p_bcs(2) - phi_s_n_bcs(1);

% Electrolyte current
i_en = z(ind_ien);
i_ep = z(ind_iep);

% Electrolyte potential
phi_e = z(ind_phi_e);

% Molar ionic flux
jn = z(ind_jn);
jp = z(ind_jp);

%% Preallocate Jacobian
JF = zeros(Nc+1,Nt);
JFb = zeros(10,Nt);
JG = zeros(Nz,Nt);
JGb = zeros(8,Nt);

% Leading Coefficients for Laplacian operator
alpha_n = 1 / (p.L_n^2 * p.delta_x_n^2);
alpha_s = 1 / (p.L_s^2 * p.delta_x_s^2);
alpha_p = 1 / (p.L_p^2 * p.delta_x_p^2);

% Leading Coefficients for Gradient operator
beta_n = 1 / (2 * p.L_n * p.delta_x_n);
beta_s = 1 / (2 * p.L_s * p.delta_x_s);
beta_p = 1 / (2 * p.L_p * p.delta_x_p);

%% [DONE] Terms of JF
%%% [DONE] Li Diffusion in Solid Phase: c_s(x,r,t)

%%%%%%%%%%%%%%%%%% Start commented by Federico %%%%%%%%%%%%%%%%
%An = zeros(3);
%An(3,2) = -6930*p.D_s_n/p.R_s_n^4;
%An(3,3) = -189/p.R_s_n^2;
%%%%%%%%%%%%%%%%%%% End commented by Federico %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start added by Federico %%%%%%%%%%%%%%%%%%
An=p.A_csn_normalized;
%%%%%%%%%%%%%%%%%%%%% End added by Federico %%%%%%%%%%%%%%%%%%%

% Loop through each "comb tooth" in anode
JF_csn = zeros(3,Nn);
for idx = 1:Nn
    JF_csn(:,idx) = An*c_s_n_mat(:,idx);
end
JF(ind_csn,1) = reshape(JF_csn,[numel(JF_csn),1]);

%%%%%%%%%%%%%%%%%% Start commented by Federico %%%%%%%%%%%%%%%%
%Ap = zeros(3);
%Ap(3,2) = -6930*p.D_s_p/p.R_s_p^4;
%Ap(3,3) = -189/p.R_s_p^2;
%%%%%%%%%%%%%%%%%%% End commented by Federico %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start added by Federico %%%%%%%%%%%%%%%%%%
Ap=p.A_csp_normalized;
%%%%%%%%%%%%%%%%%%%%% End added by Federico %%%%%%%%%%%%%%%%%%%

% Loop through each "comb tooth" in cathode
JF_csp = zeros(3,Np);
for idx = 1:Np
    JF_csp(:,idx) = Ap*c_s_p_mat(:,idx);
end
JF(ind_csp,2) = reshape(JF_csp,[numel(JF_csp),1]);


%%% [DONE] Li Diffusion in Electrolyte Phase: c_e(x,t)
% % Laplacians
% Lap_n = -2*eye(Nn) + (diag(ones(Nn-1,1),1) + diag(ones(Nn-1,1),-1));
% Lap_s = -2*eye(Ns) + (diag(ones(Ns-1,1),1) + diag(ones(Ns-1,1),-1));
% Lap_p = -2*eye(Np) + (diag(ones(Np-1,1),1) + diag(ones(Np-1,1),-1));
% 
% % Gradients
 Grad_n = diag(ones(Nn-1,1),1) - diag(ones(Nn-1,1),-1);
 Grad_s = diag(ones(Ns-1,1),1) - diag(ones(Ns-1,1),-1);
 Grad_p = diag(ones(Np-1,1),1) - diag(ones(Np-1,1),-1);

c_e_bcs = p.C_ce*c_e;
c_ex = [c_e_bcs(1); c_e(1:Nn); c_e_bcs(2); c_e(Nn+1:Nn+p.Nxs-1); ...
        c_e_bcs(3); c_e(Nn+p.Nxs : end); c_e_bcs(4)];
[A_ce, B_ce, trash_var] = c_e_mats_federico(p,c_ex);
i_ex = [0; i_en; Cur*ones(p.Nxs+1,1); i_ep; 0];

JF_ce_De = A_ce*c_e;
JF_ce_ie = B_ce*i_ex;

%%% ADD BC
% wrt D_e_n
% JF(ind_ce(1:Nn),3) = alpha_n * Lap_n * c_en;
JF(ind_ce(1:Nn),3) = JF_ce_De(1:Nn);

%%% ADD BC
% wrt D_e_s
% JF(ind_ce(Nn+1:Nn+p.Nxs-1),4) = alpha_s * Lap_s * c_es;
JF(ind_ce(Nn+1:Nn+p.Nxs-1),4) = JF_ce_De(Nn+1:Nn+p.Nxs-1);

%%% ADD BC
% wrt D_e_p
% JF(ind_ce(end-Np+1:end),5) = alpha_p * Lap_p * c_ep;
JF(ind_ce(end-Np+1:end),5) = JF_ce_De(end-Np+1:end);

% wrt (1-t_plus)
 JF(ind_ce(1:Nn),6) = beta_n/(p.epsilon_e_n * p.Faraday) * Grad_n * i_en;
 JF(ind_ce(end-Np+1:end),6) = beta_p/(p.epsilon_e_p * p.Faraday) * Grad_p * i_ep;
%JF(ind_ce,6) = 0; %B_ce*i_ex / (1 - p.t_plus);

%%% ADD BC
% wrt epsilon_e_n
% JF(ind_ce(1:Nn),15) = -beta_n*(1-p.t_plus)/(p.epsilon_e_n^2*p.Faraday) * Grad_n * i_en;
JF(ind_ce(1:Nn),15) = -(1/p.epsilon_e_n) * JF_ce_ie(1:Nn);

%%% ADD BC
% wrt epsilon_e_s

%%% ADD BC
% wrt epsilon_e_p
% JF(ind_ce(end-Np+1:end),17) = -beta_p*(1-p.t_plus)/(p.epsilon_e_p^2*p.Faraday) * Grad_p * i_ep;
JF(ind_ce(end-Np+1:end),17) = -(1/p.epsilon_e_p) * JF_ce_ie(end-Np+1:end);

%%% [DONE] Temperature: T(t)
% Equilibrium Potential and Gradient wrt bulk concentration
[Unb,trash_var,dUnbdT] = refPotentialAnode(p, c_avg_n / p.c_s_n_max);
[Upb,trash_var,dUpbdT] = refPotentialCathode(p, c_avg_p / p.c_s_p_max);

% Heat generated from intercalation (w/o boundaries for NOW)
Q_nx = p.a_s_n*p.Faraday * jn .* (Unb - T*dUnbdT);
Q_n = sum(Q_nx) * p.delta_x_n * p.L_n;

Q_px = p.a_s_p*p.Faraday * jp .* (Upb - T*dUpbdT);
Q_p = sum(Q_px) * p.delta_x_p * p.L_p;

Q_inter = Q_n + Q_p;

% Temperature Jacobian w.r.t. params
JF(ind_T,20) = 1/(p.rho_avg*p.C_p) * (p.T_amp - T);         % wrt h
JF(ind_T,21) = (p.h*(p.T_amp - T) - Cur*Volt - Q_inter);    % wrt 1/(rho_avg*C_p)

%% Terms of JFb
% For now, disregard parameter variations in the BCs

%% [DONE] Terms of JG

%%% [DONE] Potential in Solid Phase: phi_s(x,t)
% wrt 1/sig_n
JG(ind_phi_s_n,7) = -i_en + Cur;
% wrt 1/sig_n
JG(ind_phi_s_p,8) = -i_ep + Cur;


%%% [DONE] Electrolyte Current: i_e(x,t)


%%% [DONE] Potential in Electrolyte Phase: phi_e(x,t)
[trash_var,F2_pe,F3_pe,trash_var,trash_var,trash_var] = phi_e_mats(p,c_ex);

% wrt (1-t_plus)
% JG(ind_phi_e(1:Nn),6) = -(2*p.R*p.T_amp)/(p.alph*p.Faraday) * (1 + 0) ...
%                     * Grad_n * log(c_en);
% JG(ind_phi_e(Nn+1:Nn+p.Nxs-1),6) = -(2*p.R*p.T_amp)/(p.alph*p.Faraday) * (1 + 0) ...
%                     * Grad_s * log(c_es);
% JG(ind_phi_e(end-Np+1:end),6) = -(2*p.R*p.T_amp)/(p.alph*p.Faraday) * (1 + 0) ...
%                     * Grad_p * log(c_ep);
JG(ind_phi_e,6) = 1/(1-p.t_plus) * F3_pe*log(c_ex);


% wrt (1/kappa)
% JG(ind_phi_e(1:Nn),9) = i_en;
% JG(ind_phi_e(Nn+1:Nn+p.Nxs-1),9) = Cur;
% JG(ind_phi_e(end-Np+1:end),9) = i_ep;
JG(ind_phi_e,9) = F2_pe*i_ex;

% wrt (1 + d ln f_ca / d ln c_e)
% JG(ind_phi_e(1:Nn),10) = -(2*p.R*p.T_amp)/(p.alph*p.Faraday) * (1 - p.t_plus) ...
%                     * Grad_n * log(c_en);
% JG(ind_phi_e(Nn+1:Nn+p.Nxs-1),10) = -(2*p.R*p.T_amp)/(p.alph*p.Faraday) * (1 - p.t_plus) ...
%                     * Grad_s * log(c_es);
% JG(ind_phi_e(end-Np+1:end),10) = -(2*p.R*p.T_amp)/(p.alph*p.Faraday) * (1 - p.t_plus) ...
%                     * Grad_p * log(c_ep);
JG(ind_phi_e,10) = 1/(1+0) * F3_pe*log(c_ex);


%%% [DONE] Butler-Volmer Equation
aFRT = (p.alph*p.Faraday)/(p.R*p.T_amp);

% Exchange Current Density, i_0^{\pm}
[i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e);

% Equilibrium Potential, U^{\pm}(c_ss)
Unref = refPotentialAnode(p, c_ss_n / p.c_s_n_max);
Upref = refPotentialCathode(p, c_ss_p / p.c_s_p_max);

% Overpotential, \eta
eta_n = phi_s_n - phi_e(1:Nn) - Unref - p.Faraday*p.R_f_n*jn;
eta_p = phi_s_p - phi_e(end-Np+1:end) - Upref - p.Faraday*p.R_f_p*jp;

% Jacobians w.r.t. params
% wrt k-
JG(ind_jn,11) = 2/p.Faraday .* sinh(aFRT * eta_n) .* ((p.c_s_n_max - c_ss_n) .* c_ss_n .* c_en).^p.alph;
% wrt k+
JG(ind_jp,12) = 2/p.Faraday .* sinh(aFRT * eta_p) .* ((p.c_s_p_max - c_ss_p) .* c_ss_p .* c_ep).^p.alph;
% wrt R_f_n
JG(ind_jn,13) = 2/p.Faraday * i_0n .* cosh(aFRT * eta_n) .* aFRT .* (-p.Faraday*jn);
% wrt R_f_p
JG(ind_jp,14) = 2/p.Faraday * i_0p .* cosh(aFRT * eta_p) .* aFRT .* (-p.Faraday*jp);
% wrt c_s_n_max
JG(ind_jn,18) = 2/p.Faraday .* sinh(aFRT * eta_n) ...
    .* p.k_n .* p.alph .* (p.c_s_n_max - c_ss_n).^(p.alph-1) .* (c_ss_n .* c_en).^p.alph;
% wrt c_s_p_max
JG(ind_jp,19) = 2/p.Faraday .* sinh(aFRT * eta_p) ...
    .* p.k_p .* p.alph .* (p.c_s_p_max - c_ss_p).^(p.alph-1) .* (c_ss_p .* c_ep).^p.alph;

%% [DONE] Terms of JGb
%  All terms are zero.

%% Sparsify Jacobian
JF = sparse(JF * diag(p.theta0));
JFb = sparse(JFb);
JG = sparse(JG * diag(p.theta0));
JGb = sparse(JGb);

