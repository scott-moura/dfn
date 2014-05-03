%% Linear Dependence of Parameter Sensitivity for DFN Model
%   Created Apr 30, 2014 by Scott Moura
%
%   This code analyzes linear dependence in the DFN model parameter
%   sensitivity equations. The goal is to find which subset of parameters
%   are uniquely identifiable (i.e. linearly dependent)
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
%
%   OUTPUTS OF INTEREST, y
%   1  : Volt
%   2  : SOC
%   3  : T

clc;
clear;

%% Load Sensitivities
fn = 'data/sensitivity/sensitivity_new_c_s_UDDS_500s.mat';
load(fn);
disp(['Loaded Sensitivity data file:  ' fn]); 

% Parse output sensitivity data (S3)
dfn_fn = out.fn;
S3 = out.S3;
clear out;

% Vector Sizes
Nt = 21;
NT = size(S3,3);

% Parse out sensitivity for each output
S_volt = squeeze(S3(1,:,:))';
S_soc = squeeze(S3(2,:,:))';
S_temp = squeeze(S3(3,:,:))';
clear S3;

%% Parameter Ranking by successive orthogonalization

% Use QR decomposition in economy mode
[Q,R,E] = qr(S_volt,0);

% Extract Diagonal
D = diag(diag(R));
% D(D < 0) = -D(D<0);
Rbar = D\R;

%% Output Ranking
params = {'D_s^-','D_s^+','D_e^-','D_e^{sep}','D_e^+','1-t_c^0',...
    '1/sigma^-','1/sigma^+','1/kappa','1+(d ln f_c/a)(d ln c_e)','k^-','k^+',...
    'R_f^-','R_f^+','epsilon_e^-','epsilon_e^sep','epsilon_e^+',...
    'c_s_max^-','c_s_max^+','h','(rho_avg*C_p)^-1'};

fprintf(1,'-----------------------------\n');
fprintf(1,'PARAMETER SENSITIVITY RANKING\n');
fprintf(1,'-----------------------------\n');
for idx = 1:Nt
    fprintf(1,['Rank :%2.0f | Param : ' params{E(idx)} '\n'],idx); 
end

%% Plot Linear Dependence Matrix
fs = 16;

cutoff = 0.1;

figure(1); clf;
set(gcf,'Position',[330    39   710   659]);
for ii = 1:Nt
    for jj = (ii+1):Nt
        
        % If NaN, then don't include
        if(isnan(Rbar(ii,jj)))
            col = [0 0 0];
        
        % (LINEARLY DEPENDENT) If dot product of sensitivity is sufficiently close to one, then
        % only plot the border and make interior white
        elseif(abs(Rbar(ii,jj)) > cutoff)
            col = [abs(Rbar(ii,jj)), 0, 0];
            plot(ii,jj,'s','MarkerSize',24,'MarkerEdgeColor',col,'LineWidth',2);
        % (LINEARLY INDEPENDENT) If dot product of sensitivity is sufficiently close to one, then
        % only plot the border and make interior white
        else
            col = [abs(Rbar(ii,jj)), 0, 0];
            plot(ii,jj,'s','MarkerSize',24,'MarkerEdgeColor',col,...
            'MarkerFaceColor',col);
        end
        hold on;
        
    end
end

% colorbar;

xlim([0 22])
ylim([0 22])

set(gca,'FontSize',fs);
set(gca,'XTick',get(gca,'YTick'));
grid on;
xlabel('$$\textrm{Param} \ \ \ \theta_i \ \ \textrm{index}$$','FontSize',fs+2,'interpreter','latex')
ylabel('$$\textrm{Param} \ \ \ \theta_j \ \ \textrm{index}$$','FontSize',fs+2,'interpreter','latex')
tit = ['$$<S_i,S_j>/ \|S_i\| \|S_j\| \qquad', sprintf(' \\textrm{Cutoff = %1.1f}$$',cutoff)];
title(tit,'interpreter','latex');