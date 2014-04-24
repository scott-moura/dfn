%% Matrices for Li Diffusion in Electrolyte Phase, c_e(x,t)
%   Created July 15, 2011 by Scott Moura

function [F1,F2,A1,A2] = c_e_mats_implicit(p)

% Effective Diffusivities
D_eff_n = p.D_e_n;
D_eff_s = p.D_e_s;
D_eff_p = p.D_e_p;

% Leading Coefficients
alpha_n = (D_eff_n * p.delta_t) / (2 * p.L_n^2 * p.delta_x_n^2);
alpha_s = (D_eff_s * p.delta_t) / (2 * p.L_s^2 * p.delta_x_s^2);
alpha_p = (D_eff_p * p.delta_t) / (2 * p.L_p^2 * p.delta_x_p^2);

beta_n = (1 - p.t_plus) / (p.epsilon_e_n * p.Faraday * 4 * p.L_n * p.delta_x_n);
beta_s = (1 - p.t_plus) / (p.epsilon_e_s * p.Faraday * 4 * p.L_s * p.delta_x_s);
beta_p = (1 - p.t_plus) / (p.epsilon_e_p * p.Faraday * 4 * p.L_p * p.delta_x_p);

%% Block Matrices
% M1 : c_e x k+1
M1n = diag((1+2*alpha_n)*ones(p.Nxn-1,1))...
         + diag(-alpha_n*ones(p.Nxn-2,1),-1)...
         + diag(-alpha_n*ones(p.Nxn-2,1),1);
     
M1s = diag((1+2*alpha_s)*ones(p.Nxs-1,1))...
         + diag(-alpha_s*ones(p.Nxs-2,1),-1)...
         + diag(-alpha_s*ones(p.Nxs-2,1),1);
     
M1p = diag((1+2*alpha_p)*ones(p.Nxp-1,1))...
         + diag(-alpha_p*ones(p.Nxp-2,1),-1)...
         + diag(-alpha_p*ones(p.Nxp-2,1),1);

M1 = sparse(blkdiag(M1n,M1s,M1p));

% M2 : c_e z k+1
M2n = zeros(p.Nxn-1,2);
M2n(1,1) = -alpha_n;
M2n(end,end) = -alpha_n;

M2s = zeros(p.Nxs-1,2);
M2s(1,1) = -alpha_s;
M2s(end,end) = -alpha_s;

M2p = zeros(p.Nxp-1,2);
M2p(1,1) = -alpha_p;
M2p(end,end) = -alpha_p;

M2 = [M2n, zeros(p.Nxn-1,2); ...
      zeros(p.Nxs-1,1), M2s, zeros(p.Nxs-1,1);...
      zeros(p.Nxp-1,2), M2p];
M2 = sparse(M2);
  
% M3 : i_e k+1
M3 = zeros(p.Nx-3,p.Nx+1);
for idx = 1:p.Nx-3
    
    if(idx <= p.Nxn-1)
        M3(idx,idx) = beta_n;
        M3(idx,idx+2) = -beta_n;
        
    elseif(idx <= p.Nxn+p.Nxs-2)
        M3(idx,idx+1) = beta_s;
        M3(idx,idx+3) = -beta_s;
        
    else
        M3(idx,idx+2) = beta_p;
        M3(idx,idx+4) = -beta_p;
        
    end
    
end

% M4 : c_e x k
M4n = diag((1-2*alpha_n)*ones(p.Nxn-1,1))...
         + diag(alpha_n*ones(p.Nxn-2,1),-1)...
         + diag(alpha_n*ones(p.Nxn-2,1),1);
     
M4s = diag((1-2*alpha_s)*ones(p.Nxs-1,1))...
         + diag(alpha_s*ones(p.Nxs-2,1),-1)...
         + diag(alpha_s*ones(p.Nxs-2,1),1);
     
M4p = diag((1-2*alpha_p)*ones(p.Nxp-1,1))...
         + diag(alpha_p*ones(p.Nxp-2,1),-1)...
         + diag(alpha_p*ones(p.Nxp-2,1),1);

M4 = sparse(blkdiag(M4n,M4s,M4p));

% M5 : c_e z k
M5 = -M2;

% M6 : i_e k
M6 = -M3;

%% Boundary Conditions
N1 = zeros(4,p.Nx-3); % c_e x
N2 = zeros(4);        % c_e z

% BC1
N1(1,1) = 1/(p.L_n * p.delta_x_n);  % c_{n,1}
N2(1,1) = -1/(p.L_n * p.delta_x_n); % c_{n,0}

% BC2
N1(2,p.Nxn-1) = -(p.epsilon_e_n*D_eff_n)/(p.L_n * p.delta_x_n);  % c_{n,N-1}
N1(2,p.Nxn) = -(p.epsilon_e_s*D_eff_s)/(p.L_s * p.delta_x_s);    % c_{s,1}
N2(2,2) = (p.epsilon_e_n*D_eff_n)/(p.L_n * p.delta_x_n) + (p.epsilon_e_s*D_eff_s)/(p.L_s * p.delta_x_s); % c_{sp}

% BC3
N1(3,p.Nxn+p.Nxs-2) = -(p.epsilon_e_s*D_eff_s)/(p.L_s * p.delta_x_s); % c_{s,N-1}
N1(3,p.Nxn+p.Nxs-1) = -(p.epsilon_e_p*D_eff_p)/(p.L_p * p.delta_x_p); % c_{p,1}
N2(3,3) = (p.epsilon_e_s*D_eff_s)/(p.L_s * p.delta_x_s) + (p.epsilon_e_p*D_eff_p)/(p.L_p * p.delta_x_p); % c_{sp}

% BC4
N1(4,end) = -1/(p.L_p * p.delta_x_p); % c_{p,N-1}
N2(4,4) = 1/(p.L_p * p.delta_x_p); % c_{p,N}

%% F,A,B Matrics
F1 = sparse(M1 - M2*(N2\N1));
F2 = M3;
A1 = sparse(M4 - M5*(N2\N1));
A2 = M6;


