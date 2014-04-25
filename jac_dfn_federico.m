%% Jacobian for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura
%   State-dependent elements

function [f_x_full, f_z_full, g_x_full, g_z_full, varargout] = jac_dfn_federico(x,z,Cur,f_x,f_z,g_x,g_z,p)


%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx;

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
c_e = x(ind_ce);
c_en = c_e(1:Nn);
c_es = c_e(Nn+1:Nn+p.Nxs-1);
c_ep = c_e(end-Np+1:end);

% Solid Potential
phi_s_n = z(ind_phi_s_n);
phi_s_p = z(ind_phi_s_p);

% Electrolyte current
i_en = z(ind_ien);
i_ep = z(ind_iep);

% Electrolyte potential
phi_e = z(ind_phi_e);

% Molar ionic flux
jn = z(ind_jn);
jp = z(ind_jp);

% Temperature
T = x(end);

% Preallocate B matrices
B1 = zeros(Nc+1,1);
B2 = zeros(Nz,1);

%% Li Diffusion in Electrolyte Phase: c_e(x,t)
% c_e across entire sandwich
c_e_bcs = p.C_ce*c_e;
c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];

% System matrices
[A_ce, B_ce, trash_var] = c_e_mats_federico(p,c_ex);

f_x(ind_ce,ind_ce) = A_ce;
f_z(ind_ce(1:Nn),ind_ien) = B_ce(1:Nn,2:Nn+1);
f_z(ind_ce(end-Np+1:end),ind_iep) = B_ce(end-Np+1:end, end-Np:end-1);

if(nargout > 4)
    diexdI = [0; zeros(Nn,1); ones(p.Nxs+1,1); zeros(Np,1); 0];
    B1(ind_ce) = B_ce * diexdI;
end

%% Temperature: T(t)
% Bulk concentration
c_avg_n = zeros(Nn,1);
c_avg_p = zeros(Np,1);

% Loop through each "comb tooth" in anode
for idx = 1:Nn
    y_csn = p.C_csn * c_s_n_mat(:,idx);
    c_avg_n(idx) = y_csn(2);
end

% Loop through each "comb tooth" in cathode
for idx = 1:Np
    y_csp = p.C_csp * c_s_p_mat(:,idx);
    c_avg_p(idx) = y_csp(2);
end

% Equilibrium Potential and Gradient wrt bulk concentration
[Unb,dUnb] = refPotentialAnode(p, c_avg_n / p.c_s_n_max);
[Upb,dUpb] = refPotentialCathode(p, c_avg_p / p.c_s_p_max);

% Derivatives wrt c_s
dfTdcsn_int = -(p.a_s_n*p.Faraday * jn .* dUnb) / (p.rho_avg*p.C_p);
dfTdcsn = dfTdcsn_int * p.C_csn(2,:) * p.delta_x_n * p.L_n;

dfTdcsp_int = -(p.a_s_p*p.Faraday * jp .* dUpb) / (p.rho_avg*p.C_p);
dfTdcsp = dfTdcsp_int * p.C_csp(2,:) * p.delta_x_p * p.L_p;

% Derivatives wrt T
dfTdT = -p.h/(p.rho_avg*p.C_p);

% Derivatives wrt jn
dfTdjn_int = -(p.a_s_n*p.Faraday * Unb)/(p.rho_avg*p.C_p);
dfTdjn = dfTdjn_int * p.delta_x_n * p.L_n;

dfTdjp_int = -(p.a_s_p*p.Faraday * Upb)/(p.rho_avg*p.C_p);
dfTdjp = dfTdjp_int * p.delta_x_p * p.L_p;

f_x(end,ind_csn) = reshape(dfTdcsn',1,numel(dfTdcsn));
f_x(end,ind_csp) = reshape(dfTdcsp',1,numel(dfTdcsp));
f_x(end,end) = dfTdT;
f_z(end,ind_jn) = dfTdjn';
f_z(end,ind_jp) = dfTdjp';

if(nargout > 4)
    
    % Terminal Voltage
    phi_s_n_bcs = p.C_psn * phi_s_n + p.D_psn * Cur;
    phi_s_p_bcs = p.C_psp * phi_s_p + p.D_psp * Cur;
    Volt = phi_s_p_bcs(2) - phi_s_n_bcs(1);
    
    B1(end) = (Volt + Cur*(p.D_psp(end)-p.D_psn(1))) / (p.rho_avg*p.C_p);
end

%% Potential in Solid Phase: phi_s(x,t)
if(nargout > 4)
    B2(ind_phi_s_n) = p.F2_psn(:,end) + p.G_psn;
    B2(ind_phi_s_p) = p.F2_psp(:,1) + p.G_psp;
end

%% Electrolyte Current: i_e(x,t)
if(nargout > 4)
    B2(ind_ien) = p.F3_ien;
    B2(ind_iep) = p.F3_iep;
end

%% Potential in Electrolyte Phase: phi_e(x,t)
% System matrices
[F1_pe,F2_pe,F3_pe] = phi_e_mats(p,c_ex);

g_z(ind_phi_e,ind_phi_e) = F1_pe;

dpedie = F2_pe(:,[2:Nn+1, end-Np:end-1]);
g_z(ind_phi_e, [ind_ien, ind_iep]) = dpedie;

% Derivatives w.r.t. c_e here
pe_ce_on_n = [diag(ones(Nn,1)), zeros(Nn, p.Nxs-1+Np)];
pe_ce_on_s = [zeros(p.Nxs-1, Nn), diag(ones(p.Nxs-1,1)), zeros(p.Nxs-1, Np)];
pe_ce_on_p = [zeros(Np, p.Nxs-1+Nn), diag(ones(Np,1))];

F3_pe_bcs = [p.C_ce(1,:); pe_ce_on_n ;p.C_ce(2,:); pe_ce_on_s; p.C_ce(3,:); pe_ce_on_p; p.C_ce(4,:)];
g_x(ind_phi_e,ind_ce) = F3_pe * F3_pe_bcs * diag(1./c_e);

% g_x(ind_phi_e,ind_ce) = F3_pe(:,[2:Nn+1,Nn+3:Nn+p.Nxs+1,Nn+p.Nxs+3:end-1]) * diag(1./c_e);

if(nargout > 4)
    B2(ind_phi_e) = F2_pe * diexdI;
end

%% Butler-Volmer Equation
% Surface concentration
c_ss_n = (p.C_csn(1,:) * c_s_n_mat)';
c_ss_p = (p.C_csp(1,:) * c_s_p_mat)';

% Param
aFRT = (p.alph*p.Faraday)/(p.R*p.T_amp);

% Exchange Current Density
[i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e);

di0dcssn = p.k_n*c_en.*(p.c_s_n_max - 2*c_ss_n) ./ ...
            (2 * sqrt(c_en.*c_ss_n.*(p.c_s_n_max - c_ss_n)));
di0dcssp = p.k_p*c_ep.*(p.c_s_p_max - 2*c_ss_p) ./ ...
            (2 * sqrt(c_ep.*c_ss_p.*(p.c_s_p_max - c_ss_p)));

di0dcen = p.k_n*sqrt(c_en.*c_ss_n.*(p.c_s_n_max - c_ss_n)) ./ (2*c_en);
di0dcep = p.k_p*sqrt(c_ep.*c_ss_p.*(p.c_s_p_max - c_ss_p)) ./ (2*c_ep);

% Equilibrium Potential
[Unref,dUnref] = refPotentialAnode(p, c_ss_n / p.c_s_n_max);
[Upref,dUpref] = refPotentialCathode(p, c_ss_p / p.c_s_p_max);

% Overpotential
eta_n = phi_s_n - phi_e(1:Nn) - Unref - p.Faraday*p.R_f_n*jn;
eta_p = phi_s_p - phi_e(end-Np+1:end) - Upref - p.Faraday*p.R_f_p*jp;

% Components of Jacobian
df6dcssn = 2/p.Faraday * di0dcssn .* sinh(aFRT*eta_n) ...
          - 2/p.Faraday * i_0n .* cosh(aFRT*eta_n) * aFRT .* dUnref;
df6dcssp = 2/p.Faraday * di0dcssp .* sinh(aFRT*eta_p) ...
          - 2/p.Faraday * i_0p .* cosh(aFRT*eta_p) * aFRT .* dUpref;
      
df6dcen = 2/p.Faraday * di0dcen .* sinh(aFRT*eta_n);
df6dcep = 2/p.Faraday * di0dcep .* sinh(aFRT*eta_p);

df6psn = 2/p.Faraday*i_0n.*cosh(aFRT*eta_n) * aFRT;
df6psp = 2/p.Faraday*i_0p.*cosh(aFRT*eta_p) * aFRT;

df6pen = -2/p.Faraday*i_0n.*cosh(aFRT*eta_n) * aFRT;
df6pep = -2/p.Faraday*i_0p.*cosh(aFRT*eta_p) * aFRT;

df6Jn = -2/p.Faraday*i_0n.*cosh(aFRT*eta_n) * aFRT*p.R_f_n*p.Faraday - 1;
df6Jp = -2/p.Faraday*i_0p.*cosh(aFRT*eta_p) * aFRT*p.R_f_p*p.Faraday - 1;

% Input into Jacobian

% Loop through each "comb tooth" in anode
Cell_df6n = cell(Nn,1);
Cell_df6p = cell(Np,1);
for idx = 1:Nn
    Cell_df6n{idx} = df6dcssn(idx) * p.C_csn(1,:);
end
for idx = 1:Np
    Cell_df6p{idx} = df6dcssp(idx) * p.C_csp(1,:);
end

rsn = ones(Nn,1);
rsp = ones(Np,1);
csn = p.PadeOrder * ones(1,Nn);
csp = p.PadeOrder * ones(1,Np);
g_x(ind_jn,ind_csn) = blkdiagFast(rsn, csn, Cell_df6n{:});
g_x(ind_jp,ind_csp) = blkdiagFast(rsp, csp, Cell_df6p{:});

g_x(ind_jn,ind_ce(1:Nn)) = diag(df6dcen);
g_x(ind_jp,ind_ce(end-Np+1:end)) = diag(df6dcep);

g_z(ind_jn,ind_phi_s_n) = diag(df6psn);
g_z(ind_jp,ind_phi_s_p) = diag(df6psp);

g_z(ind_jn,ind_phi_e(1:Nn)) = diag(df6pen);
g_z(ind_jp,ind_phi_e(end-Np+1:end)) = diag(df6pep);

g_z(ind_jn,ind_jn) = diag(df6Jn);
g_z(ind_jp,ind_jp) = diag(df6Jp);

%%
f_x_full = f_x;
f_z_full = f_z;
g_x_full = g_x;
g_z_full = g_z;

% spy(Jac)
% pause;

if(nargout > 4)
    varargout{1} = B1;
    varargout{2} = B2;
end
