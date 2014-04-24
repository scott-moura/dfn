%% Crank-Nicolson Eqns for Doyle-Fuller-Newman Model
%   Created June 5, 2012 by Scott Moura
%   Block Jacobi Method

function [x_nxtf, z_nxtf, varargout] = cn_dfn_jacobi(x,z,Cur_vec,p)

Cur_prv = Cur_vec(1);
Cur = Cur_vec(2);
Cur_nxt = Cur_vec(3);

%% Parameters

% Newton params
maxIters = 50;
tol = 1e-5;

% Sparse Linear System Solver Params
linEqnTol = 1e-5; % Default is 1e-6
linEqnMaxIter = 100;  % Default is 20

relres(1) = 1;

%% Index Parameters
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

%% Preallocate
c_sn_nxt = zeros(Ncsn,maxIters);
c_sp_nxt = zeros(Ncsp,maxIters);
c_e_nxt = zeros(Nce,maxIters);
phi_sn_nxt = zeros(Nn,maxIters);
phi_sp_nxt = zeros(Np,maxIters);
i_en_nxt = zeros(Nn,maxIters);
i_ep_nxt = zeros(Nn,maxIters);
phi_e_nxt = zeros(Nx,maxIters);
jn_nxt = zeros(Nn,maxIters);
jp_nxt = zeros(Np,maxIters);

Del_c_sn_nxt = zeros(Ncsn,maxIters);
Del_c_sp_nxt = zeros(Ncsp,maxIters);
Del_c_e_nxt = zeros(Nce,maxIters);
Del_phi_sn_nxt = zeros(Nn,maxIters);
Del_phi_sp_nxt = zeros(Np,maxIters);
Del_i_en_nxt = zeros(Nn,maxIters);
Del_i_ep_nxt = zeros(Nn,maxIters);
Del_phi_e_nxt = zeros(Nx,maxIters);
Del_jn_nxt = zeros(Nn,maxIters);
Del_jp_nxt = zeros(Np,maxIters);

%% Solve for consistent ICs at k
if(Cur ~= Cur_prv)
    
    disp('Solving for consistent ICs');
    
    % Preallocate
    z_cons = zeros(Nz,maxIters);
    z_cons(:,1) = z;
    
    for idx = 1:(maxIters-1)
        
        % DAE eqns for current time-step
        [~,g] = dae_dfn(x,z_cons(:,idx),Cur,p);
        
        % Jacobian of DAE eqns
        [~,g_z] = jac_dfn(x,z_cons(:,idx),p.g_x,p.g_z,p);
        
        % Newton Iteration
        [LL,UU] = ilu(g_z,struct('type','ilutp','droptol',1e-6,'udiag',1));
        [Delta_z,flag] = bicg(g_z,-g,linEqnTol,linEqnMaxIter,LL,UU);
        z_cons(:,idx+1) = z_cons(:,idx) + Delta_z;
        
        % Check stopping criterion
        relres_z = norm(Delta_z,inf) / norm(z,inf);
        if(relres_z < tol)
            break;
        elseif(idx == (maxIters-1))
            fprintf(1,'Warning: Max Newton Iters Reached | RelChange = %3.2f%%\n',relres_z*100);
        end
        
    end
    
    z = z_cons(:,idx+1);
    
end

%% Solve Nonlinear System using Newton's Method
% DAE eqns for current time-step
[f,~] = dae_dfn(x,z,Cur,p);

% Initialize next x,z
c_sn_nxt(:,1) = x(ind_csn);
c_sp_nxt(:,1) = x(ind_csp);
c_e_nxt(:,1) = x(ind_ce);
phi_sn_nxt(:,1) = z(ind_phi_s_n);
phi_sp_nxt(:,1) = z(ind_phi_s_p);
i_en_nxt(:,1) = z(ind_ien);
i_ep_nxt(:,1) = z(ind_iep);
phi_e_nxt(:,1) = z(ind_phi_e);
jn_nxt(:,1) = z(ind_jn);
jp_nxt(:,1) = z(ind_jp);

x_nxt = x;
z_nxt = z;

% Iterate Newton's Method
for idx = 1:(maxIters-1)
    
    % DAE eqns for next time-step
    [f_nxt, g_nxt] = dae_dfn(x_nxt, z_nxt, Cur_nxt, p);

    if(~isreal(g_nxt))
        Cur_nxt
        disp('z_nxt')
        z_nxt
        disp('g_nxt')
        g_nxt(ind_jp)
    end
    
    % Nonlinear system of equations
    F1 = x - x_nxt + p.delta_t/2 * (f + f_nxt);
    F2 = g_nxt;
    
    % Jacobian of DAE eqns
    [g_x, g_z] = jac_dfn(x_nxt,z_nxt,p.g_x,p.g_z,p);

    %%% Jacobi method for solving Newton's method iteration
    % Solid concentration
    D_csn = -eye(Ncsn) + p.delta_t/2 * p.f_x(ind_csn,ind_csn);
    F_csn = -F1(ind_csn) - p.delta_t/2 * p.f_z(ind_csn,ind_jn) * Del_jn_nxt(:,idx);
    
    D_csp = -eye(Ncsp) + p.delta_t/2 * p.f_x(ind_csp,ind_csp);
    F_csp = -F1(ind_csp) - p.delta_t/2 * p.f_z(ind_csp,ind_jp) * Del_jp_nxt(:,idx);
    
    Del_c_sn_nxt(:,idx+1) = bicg(D_csn,F_csn,[],maxIters);
    Del_c_sp_nxt(:,idx+1) = bicg(D_csp,F_csp,[],maxIters);

    c_sn_nxt(:,idx+1) = c_sn_nxt(:,idx) + Del_c_sn_nxt(:,idx+1);
    c_sp_nxt(:,idx+1) = c_sp_nxt(:,idx) + Del_c_sp_nxt(:,idx+1);
    
    % Electrolyte concentration
    D_ce = -eye(Nx) + p.delta_t/2*p.A_ce;
    Del_ie_nxt = [0; Del_i_en_nxt(:,idx); zeros(p.Nxs+1,1); Del_i_ep_nxt(:,idx); 0];
    F_ce = -F1(ind_ce) - p.delta_t/2 * p.B_ce * Del_ie_nxt;
    
    Del_c_e_nxt(:,idx+1) = bicg(D_ce,F_ce,[],maxIters);
    c_e_nxt(:,idx+1) = c_e_nxt(:,idx) + Del_c_e_nxt(:,idx+1);
    
    % Solid Potential
    F_psn = -F2(ind_phi_s_n) - p.F2_psn * [0; Del_i_en_nxt(:,idx); 0];
    F_psp = -F2(ind_phi_s_p) - p.F2_psp * [0; Del_i_ep_nxt(:,idx); 0];
    
    Del_phi_sn_nxt(:,idx+1) = bicg(p.F1_psn,F_psn,[],maxIters);
    Del_phi_sp_nxt(:,idx+1) = bicg(p.F1_psp,F_psp,[],maxIters);
    
    phi_sn_nxt(:,idx+1) = phi_sn_nxt(:,idx) + Del_phi_sn_nxt(:,idx+1);
    phi_sp_nxt(:,idx+1) = phi_sp_nxt(:,idx) + Del_phi_sp_nxt(:,idx+1);
    
    % Electrolyte Current
    F_ien = -F2(ind_ien) - p.F2_ien * Del_jn_nxt(:,idx);
    F_iep = -F2(ind_iep) - p.F2_iep * Del_jp_nxt(:,idx);
    
    Del_i_en_nxt(:,idx+1) = bicg(p.F1_ien,F_ien,[],maxIters);
    Del_i_ep_nxt(:,idx+1) = bicg(p.F1_iep,F_iep,[],maxIters); 

    i_en_nxt(:,idx+1) = i_en_nxt(:,idx) + Del_i_en_nxt(:,idx+1);
    i_ep_nxt(:,idx+1) = i_ep_nxt(:,idx) + Del_i_ep_nxt(:,idx+1);
    
    if(~isreal(i_en_nxt(:,idx+1)))
        disp(i_en_nxt(:,idx+1));
        pause;
    end

    if(~isreal(i_ep_nxt(:,idx+1)))
        disp('i_ep_nxt')
        F2(ind_iep)
        Del_jp_nxt(:,idx)
        pause;
    end
    
    % Electrolyte Potential
    D_pe = g_x(ind_phi_e,ind_phi_e);
    dpedie = [zeros(Nx,1), g_z(ind_phi_e,ind_ien), zeros(Nx,p.Nxs+1), g_z(ind_phi_e,ind_iep), zeros(Nx,1)];
    Del_ie_nxtkp1 = [0; Del_i_en_nxt(:,idx+1); zeros(p.Nxs+1,1); Del_i_ep_nxt(:,idx+1); 0];
    
    F_pe = -F2(ind_phi_e) - dpedie * Del_ie_nxtkp1 ...
                          - g_x(ind_phi_e,ind_ce) * Del_c_e_nxt(:,idx+1);
                      
    Del_phi_e_nxt(:,idx+1) = bicg(D_pe,F_pe,[],maxIters);
    phi_e_nxt(:,idx+1) = phi_e_nxt(:,idx) + Del_phi_e_nxt(:,idx+1);
    
    % Molar Ionic Flux
    df6dcssn = g_x(ind_jn,ind_csn);
    df6dcen = g_x(ind_jn,ind_ce(1:Nn));
    df6dpsn = g_z(ind_jn,ind_phi_s_n);
    df6dpen = g_z(ind_jn,ind_phi_e(1:Nn));

    df6dcssp = g_x(ind_jp,ind_csp);
    df6dcep = g_x(ind_jp,ind_ce(end-Np+1:end));
    df6dpsp = g_z(ind_jp,ind_phi_s_p);
    df6dpep = g_z(ind_jp,ind_phi_e(end-Np+1:end));
    
    D_jn = g_z(ind_jn,ind_jn);
    D_jp = g_z(ind_jp,ind_jp);
    
    F_jn = -F2(ind_jn) - df6dcssn * Del_c_sn_nxt(:,idx+1) - df6dcen * Del_c_e_nxt(1:Nn,idx+1)...
           -df6dpsn * Del_phi_sn_nxt(:,idx+1) - df6dpen * Del_phi_e_nxt(1:Nn,idx+1);
    F_jp = -F2(ind_jp) - df6dcssp * Del_c_sp_nxt(:,idx+1) - df6dcep * Del_c_e_nxt(end-Np+1:end,idx+1)...
           -df6dpsp * Del_phi_sp_nxt(:,idx+1) - df6dpep * Del_phi_e_nxt(end-Np+1:end,idx+1);
       
    if(Cur_nxt > 0.01)
        F_jp
    end
       
    Del_jn_nxt(:,idx+1) = D_jn \ F_jn;
    Del_jp_nxt(:,idx+1) = D_jp \ F_jp;
    
    jn_nxt(:,idx+1) = jn_nxt(:,idx) + Del_jn_nxt(:,idx+1);
    jp_nxt(:,idx+1) = jp_nxt(:,idx) + Del_jp_nxt(:,idx+1);
    
    if(~isreal(jp_nxt(:,idx+1)))
        disp('jp_nxt')
        F_jp
        F2(ind_jp)
    end
    
    % Compile into vectors
    x_nxt = [c_sn_nxt(:,idx+1); c_sp_nxt(:,idx+1); c_e_nxt(:,idx+1)];
    z_nxt = [phi_sn_nxt(:,idx+1); phi_sp_nxt(:,idx+1); ...
             i_en_nxt(:,idx+1); i_ep_nxt(:,idx+1); ...
             phi_e_nxt(:,idx+1); jn_nxt(:,idx+1); jp_nxt(:,idx+1)];
       
    % Check stopping criterion
%     y_nxt = [x_nxt(:,idx+1); z_nxt(:,idx+1)];
%     relres(idx+1) = norm(Delta_y,inf) / norm(y_nxt,inf);
%     if(relres(idx+1) < tol)
%         break;
%     elseif(idx == (maxIters-1))
%         fprintf(1,'Warning: Max Newton Iters Reached | RelChange = %1.4e%%\n',relres(end)*100);
%     end
    
end

%% Output next x,z
x_nxtf = x_nxt;
z_nxtf = z_nxt;

newtonStats.iters = idx;
newtonStats.relres = relres;
varargout{1} = newtonStats;

