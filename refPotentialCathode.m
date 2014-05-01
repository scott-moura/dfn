%% Reference Potential for Pos. Electrode: Unref(theta_n)
%   Created July 12, 2011 by Scott Moura

function [Uref,varargout] = refPotentialCathode(p,theta)

% % Polynomail Fit
% Uref = ppvalFast(p.Uppp,theta);

%%% DUALFOIL: CoO2 (Cobalt Dioxide) 0.5 < y < 0.99
% Uref = 2.16216+0.07645*tanh(30.834-54.4806*theta) ...
%  + 2.1581*tanh(52.294-50.294*theta) ...
%  - 0.14169*tanh(11.0923-19.8543*theta) ...
%  + 0.2051*tanh(1.4684-5.4888*theta) ...
%  + 0.2531*tanh((-theta+0.56478)/0.1316) ...
%  - 0.02167*tanh((theta-0.525)/0.006);

%% Bosch Manganese Chemsitry (0.2 <= SOC <= 1.0)
Uref = 4.06279 + 0.0677504*tanh(-21.8502*theta + 12.8268) ...
      -0.105734 * ((1.00167 - theta).^(-0.379571) - 1.575994) ...
      -0.045 * exp(-71.69*theta.^8) + 0.01 * exp(-200.0*(theta - 0.19));

% Gradient of OCP wrt theta
if(nargout == 2)
    
%     % Polynomail Fit
%     dUref = ppvalFast(p.dUppp,theta);
%     varargout{1} = dUref / p.c_s_p_max;

%%% DUALFOIL: CoO2 (Cobalt Dioxide) 0.5 < y < 0.99
% dUref = 0.07645*(-54.4806/p.c_s_p_max)* ...
% ((1.0./cosh(30.834-54.4806*theta)).^2) ...
% +2.1581*(-50.294/p.c_s_p_max)*((cosh(52.294-50.294*theta)).^(-2)) ...
% +0.14169*(19.854/p.c_s_p_max)*((cosh(11.0923-19.8543*theta)).^(-2)) ...
% -0.2051*(5.4888/p.c_s_p_max)*((cosh(1.4684-5.4888*theta)).^(-2)) ...
% -0.2531/0.1316/p.c_s_p_max*((cosh((-theta+0.56478)/0.1316)).^(-2)) ...
%  - 0.02167/0.006/p.c_s_p_max*((cosh((theta-0.525)/0.006)).^(-2));

%% Bosch Manganese Chemsitry (0.2 <= SOC <= 1.0)
dUref = -0.0677504*21.8502*(cosh(-21.8502*theta + 12.8268)).^(-2) ...
    -0.105734*0.379571*(1.00167-theta).^(-0.379571-1) ...
    +0.045*71.69*8*theta.^7 .* exp(-71.69*theta.^8) ...
    -0.01*200.0*exp(-200.0*(theta-0.19));

varargout{1} = dUref;
    
end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end