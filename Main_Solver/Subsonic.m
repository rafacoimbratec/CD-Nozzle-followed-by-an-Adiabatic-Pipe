function [M1, Me, flag_area] = Subsonic(T0,p0,gamma,AR,TubeLenght,epsilon,pe)
% SUBSONICSOLVER Computes M1 and Me for a subsonic flow through a converging-diverging nozzle 
% followed by a pipe with friction, solving isentropic and Fanno relations.
%
% Inputs:
%   p0  - Stagnation pressure (Pa)
%   gamma - Specific heat ratio (dimensionless)
%   Ae  - Exit area (m²)
%   At  - Throat area (m²)
%   L   - Pipe length (m)
%   D   - Pipe diameter (m)
%   f   - Darcy friction factor (dimensionless)
%   pe  - Exit pressure (Pa)
%
% Outputs:
%   M1  - Mach number at the inlet of the pipe
%   Me  - Mach number at the exit of the pipe
%   flag_area - 1 if At is too small (flow chokes), 0 otherwise

% ---- Define Initial Guess for Solver ----
y_guess = [p0, 0.2, 0.1]; % Initial guesses for [p0e, Me, M1]

% ---- Solve Nonlinear Equations Using fsolve ----
[y, exitflag] = fsolve(@(x) SubsonicEquations(x, T0,p0, gamma, pe, epsilon, TubeLenght, AR), y_guess);

% ---- Extract Solutions ----
if exitflag <= 0
    error('fsolve did not converge. Try different inputs.');
end

p0e = y(1); % Stagnation pressure at pipe exit (not used in outputs)
Me = y(2);  % Mach number at pipe exit
M1 = y(3); % Mach number at pipe inlet

% ---- Check Critical Throat Area Condition ----
A_critical =(1/M1 * ((2 + (gamma - 1) * M1^2) / (gamma + 1))^((gamma + 1) / (2 * (gamma - 1))));

if AR > A_critical
    flag_area = 1; % Throat area too small → choking
else
    flag_area = 0; % Flow remains subsonic
end

end

% ---- Function for Subsonic Equations ---- %
function F = SubsonicEquations(x, T0, p0, gamma, pe, epsilon, TubeLenght, AR)
% Defines the system of nonlinear equations to solve for M1 and Me.

p0e = x(1); % Stagnation pressure at the exit
Me = x(2);  % Mach number at the exit
M1 = x(3);  % Mach number at the inlet

% Equation (1): Isentropic relation for stagnation pressures
F(1) = (p0e / p0) - (M1 / Me) * ((2 + (gamma - 1) * Me^2) / (2 + (gamma - 1) * M1^2))^((gamma + 1) / (2 * (gamma - 1)));

% Equation (2): Isentropic relation for exit static pressure
F(2) = (p0e / pe) - (1 + (gamma - 1) / 2 * Me^2)^(gamma / (gamma - 1));

T1 = T0 * (1+((gamma-1)/2)*M1^2)^-1;
p1 = p0 * (T0/T1)^(-gamma/(gamma-1));
TE = T0 * (1+((gamma-1)/2)*Me^2)^-1;
[ReM1,fM1]=Churchill(gamma,T1,p1,M1,epsilon,AR);
[ReMe,fMe]=Churchill(gamma,TE,pe,Me,epsilon,AR);

% Equation (3): Fanno flow relation (accounting for friction)
F(3) = (fMe * TubeLenght / (2*AR)) - ...
       ((1 - M1^2) / (gamma * M1^2) + (gamma + 1) / (2 * gamma) * log(((gamma + 1) * M1^2) / (2 + (gamma - 1) * M1^2))) + ...
       ((1 - Me^2) / (gamma * Me^2) + (gamma + 1) / (2 * gamma) * log(((gamma + 1) * Me^2) / (2 + (gamma - 1) * Me^2)));

end
