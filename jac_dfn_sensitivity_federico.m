%% Jacobian for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura
%   State-dependent elements

function [A_corr_p, A_corr_n, A_corr_s] = jac_dfn_sensitivity_federico(x,z,Cur,f_x,f_z,g_x,g_z,p)


%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;

ind_ce = Ncsn+Ncsp+1:Nc;


% Electrolyte concentration
c_e = x(ind_ce);
c_en = c_e(1:Nn);
c_es = c_e(Nn+1:Nn+p.Nxs-1);
c_ep = c_e(end-Np+1:end);

A_corr_p=zeros(size(f_x,1),1);
A_corr_n=zeros(size(f_x,1),1);
A_corr_s=zeros(size(f_x,1),1);



%% Li Diffusion in Electrolyte Phase: c_e(x,t)
% c_e across entire sandwich
c_e_bcs = p.C_ce*c_e;
c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];

% System matrices
[trash_var, trash_var, trash_var, A_corr_sens_p, A_corr_sens_n, A_corr_sens_s] = c_e_mats_sensitivity_federico(p,c_ex);

%Correction Terms for boundary conditions in the sensitivity equations
A_corr_p(ind_ce)=A_corr_sens_p*c_e;
A_corr_n(ind_ce)=A_corr_sens_n*c_e;
A_corr_s(ind_ce)=A_corr_sens_s*c_e;


