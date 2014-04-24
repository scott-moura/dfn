%% DAEs for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura

function [A, B, varargout] = dae_dfn_lin(x,z,Cur,p)


%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;

ind_csn = 1:Ncsn;
ind_csp = Ncsn+1:Ncsn+Ncsp;

ind_cs = 1:Ncsn+Ncsp;
ind_ce = Ncsn+Ncsp+1:Nc;

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

% Electrolyte concentration
c_e = x((Ncsn+Ncsp + 1):(Nc));

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

%% Preallocate Jacobian
A11 = zeros(Nc+1);
A12 = zeros(Nc+1,Nz);
A21 = zeros(Nz,Nc+1);
A22 = zeros(Nz);

%% Li Diffusion in Solid Phase: c_s(x,r,t)
% Preallocate matrices for derivatives
Cell_Acsn = cell(Nn,1);
Cell_Bcsn = cell(Nn,1);
Cell_Acsp = cell(Np,1);
Cell_Bcsp = cell(Np,1);

c_ss_n = zeros(Nn,1);
c_ss_p = zeros(Np,1);
c_avg_n = zeros(Nn,1);
c_avg_p = zeros(Np,1);

% Loop through each "comb tooth" in anode
for idx = 1:Nn
    Cell_Acsn{idx} = p.A_csn;
    Cell_Bcsn{idx} = p.B_csn;
    
    y_csn = p.C_csn * c_s_n_mat(:,idx);
    c_ss_n(idx) = y_csn(1);
    c_avg_n(idx) = y_csn(2);
end

% Loop through each "comb tooth" in cathode
for idx = 1:Np
    Cell_Acsp{idx} = p.A_csp;
    Cell_Bcsp{idx} = p.B_csp;
    
    y_csp = p.C_csp * c_s_p_mat(:,idx);
    c_ss_p(idx) = y_csp(1);
    c_avg_p(idx) = y_csp(2);
end


A11(ind_csn,ind_csn) = blkdiag(Cell_Acsn{:});
A12(ind_csn,ind_jn) = blkdiag(Cell_Bcsn{:});

A11(ind_csp,ind_csp) = blkdiag(Cell_Acsp{:});
A12(ind_csp,ind_jp) = blkdiag(Cell_Bcsp{:});

%% Li Diffusion in Electrolyte Phase: c_e(x,t)
% Electrolyte current across all three regions
i_ex = [0; i_en; Cur*ones(p.Nxs+1,1); i_ep; 0];

c_e_bcs = p.C_ce*c_e;
c_ex = [c_e_bcs(1); c_e(1:Nn); c_e_bcs(2); c_e(Nn+1:Nn+p.Nxs-1); ...
        c_e_bcs(3); c_e(Nn+p.Nxs : end); c_e_bcs(4)];

% System Matrices
[A_ce, B_ce, ~] = c_e_mats(p,c_ex);

A11(ind_ce,ind_ce) = A_ce;
A12(ind_ce(1:Nn),ind_ien) = B_ce(1:Nn,2:Nn+1);
A12(ind_ce(end-Np+1:end),ind_iep) = B_ce(end-Np+1:end, end-Np:end-1);

%% Potential in Solid Phase: phi_s(x,t)
A22(ind_phi_s_n,ind_phi_s_n) = p.F1_psn;
A22(ind_phi_s_p,ind_phi_s_p) = p.F1_psp;

A22(ind_phi_s_n,ind_ien) = p.F2_psn(:,2:end-1);
A22(ind_phi_s_p,ind_iep) = p.F2_psp(:,2:end-1);

%% Electrolyte Current: i_e(x,t)
A22(ind_ien,ind_ien) = p.F1_ien;
A22(ind_iep,ind_iep) = p.F1_iep;

A22(ind_ien,ind_jn) = p.F2_ien;
A22(ind_iep,ind_jp) = p.F2_iep;

%% Potential in Electrolyte Phase: phi_e(x,t)
% System matrices
[F1_pe,F2_pe,F3_pe] = phi_e_mats(p,c_ex);

A22(ind_phi_e,ind_phi_e) = F1_pe;

dpedie = F2_pe(:,[2:Nn+1, end-Np:end-1]);
A22(ind_phi_e, [ind_ien, ind_iep]) = dpedie;

% Derivatives w.r.t. c_e here
pe_ce_on_n = [diag(ones(Nn,1)), zeros(Nn, p.Nxs-1+Np)];
pe_ce_on_s = [zeros(p.Nxs-1, Nn), diag(ones(p.Nxs-1,1)), zeros(p.Nxs-1, Np)];
pe_ce_on_p = [zeros(Np, p.Nxs-1+Nn), diag(ones(Np,1))];

F3_pe_bcs = [p.C_ce(1,:); pe_ce_on_n ;p.C_ce(2,:); pe_ce_on_s; p.C_ce(3,:); pe_ce_on_p; p.C_ce(4,:)];
A21(ind_phi_e,ind_ce) = F3_pe * F3_pe_bcs * diag(1./c_e);

%% Butler-Volmer Equation
aFRT = (p.alph*p.Faraday)/(p.R*p.T_amp);

% Exchange Current Density
[i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,p.c_e);

% Equilibrium Potential
[Unref,dUnref] = refPotentialAnode(p, c_ss_n / p.c_s_n_max);
[Upref,dUpref] = refPotentialCathode(p, c_ss_p / p.c_s_p_max);

% Overpotential
eta_n = phi_s_n - phi_e(1:Nn) - Unref - p.Faraday*p.R_SEI*jn;
eta_p = phi_s_p - phi_e(end-Np+1:end) - Upref;

% Algebraic eqns (semi-explicit form)
jn_dot = 2/p.Faraday * i_0n .* sinh(aFRT * eta_n) - jn;
jp_dot = 2/p.Faraday * i_0p .* sinh(aFRT * eta_p) - jp;

A

%% Concatenate Time Derivatives
f = [c_s_n_dot; c_s_p_dot; c_e_dot];
g = [phi_sn_dot; phi_sp_dot; i_en_dot; i_ep_dot; ...
     phi_e_dot; jn_dot; jp_dot];
     
%% Optional Output Vars
% Terminal Voltage
phi_s_n_bcs = p.C_psn * phi_s_n + p.D_psn * Cur;
phi_s_p_bcs = p.C_psp * phi_s_p + p.D_psp * Cur;

Volt = phi_s_p_bcs(2) - phi_s_n_bcs(1);

% Conservation of Li-ion matters
nLi = sum(c_avg_n) * p.delta_x_n * Nn ...
     + sum(c_avg_p) * p.delta_x_p * Np ...
     + sum(c_e(1:Nn)) * p.delta_x_n * Nn ...
     + sum(c_e(Nn+1:end-Np)) * p.delta_x_s * p.Nxs-1 ...
     + sum(c_e(end-Np+1:end)) * p.delta_x_p * Np;
 
nLidot = sum(jn) * p.delta_x_n * Nn + sum(jp) * p.delta_x_p * Np;

c_e0p = c_ex(end);

phi_e_bcs = C1_pe*phi_e + C2_pe*i_ex + C3_pe*log(c_ex);
eta_s_Ln = phi_s_n_bcs(2) - phi_e_bcs(2);

% Aggregate Outputs
varargout{1} = [c_ss_n; c_ss_p; c_avg_n; c_avg_p; c_ex; eta_n; eta_p;...
    c_e0p; eta_s_Ln;...
    Volt; nLi; nLidot];
