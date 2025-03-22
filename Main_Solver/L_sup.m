function [x_sh, M_B] = L_sup(p0, T0, gamma, pe, AR, TubeLength, epsilon, yT, GT, radius)
%LONG TUBE Calculates the position of the shock wave for a long tube
%given that the inlet flow is supersonic.

R = 287; % Specific gas constant for air

% ---- Compute **Critical Mass Flow Rate** ----
m_crit = p0 / sqrt(T0) * sqrt(gamma / R) * (1 + ((gamma - 1) / 2))^(-0.5 * (gamma + 1) / (gamma - 1)) * 2 * yT;
m_crit = double(m_crit);

% ---- Solve for Supersonic Mach Number at Tube Inlet (M1) ----
M1 = fzero(@(M) -AR + 1/M * (2 / (gamma + 1) * (1 + ((gamma - 1) / 2) * M^2))^(0.5 * (gamma + 1) / (gamma - 1)), [1, 100]);
M1 = double(M1);

%STAGNATION PRESSURE AT THE OUTLET OF THE TUBE (p0e_crit) for Me=1:
p0e_crit_lmax = m_crit*sqrt(T0)/((2*AR)*sqrt(gamma/R)*1*(1+(gamma-1)/2*1^2)^(-0.5*(gamma+1)/(gamma-1)));
%PRESSURE AT THE OUTLET OF THE TUBE (pe_crit) for Me=1:
pe_crit_lmax=p0e_crit_lmax/( (1+(gamma-1)/2*1^2)^(gamma/(gamma-1)));
fprintf('critical_5 = %.15f\n', pe_crit_lmax);

%Lets compare with the user input
if (pe_crit_lmax<=pe)
    %Subsonic outlet in the tube, in the case pe_crit_lmax=pe ME=1
    % ---- Solve for Mach Number at Tube Outlet (Me) ----
    Me = fzero(@(M) -m_crit + pe * 2 * AR * M * sqrt(gamma / (R * (T0 / (1 + ((gamma - 1) / 2) * M^2)))), 0.5);
    Me = double(Me);

    % ---- Solve for Shock Location and Mach Number After the Shock ----
    shock_guess = [0.5, 1.5, 0.5]; % Initial guess: [x_shock, M_A, M_B]

    
    shock_lmax = fsolve(@(g) M_shock(g, T0,p0,pe,gamma, M1, TubeLength, epsilon, AR, Me), shock_guess);

    % ---- Assign Shock Location and Post-Shock Mach Number ----
    x_sh = shock_lmax(1) + 2 * yT;
    M_B = shock_lmax(3);

else
    %No solution, because it means the flow would be chocked at the throath
    %and chocked on another place after the shock doubly chocked
    x_sh=-1;
    M_B=0;
end
end

%% ---- Function to Solve for Shock Relations ----
function H = M_shock(g,T0,p0,pe, gamma, M1, L, epsilon, AR, Me)
% Solves for shock relations inside the tube using Fanno flow and shock conditions
% Inputs:
%   - g       : Array containing [x_shock, M_A, M_B]
%   - gamma   : Specific heat ratio
%   - M1      : Mach number at inlet
%   - L       : Tube length
%   - f       : Darcy friction factor
%   - AR      : Area ratio
%   - Me      : Exit Mach number

T1 = T0 * (1+((gamma-1)/2)*M1^2)^-1;
p1 = p0 * (T0/T1)^(-gamma/(gamma-1));
TE = T0 * (1+((gamma-1)/2)*Me^2)^-1;
[ReM1,fM1]=Churchill(gamma,T1,p1,M1,epsilon,AR);
[ReMe,fMe]=Churchill(gamma,TE,pe,Me,epsilon,AR);

% Fanno flow equation up to the shock location
H(1) = fM1 * g(1) / (2*AR) - (1 - M1^2) / (gamma * M1^2) - ((gamma + 1) / (2 * gamma)) * log(((gamma + 1) * M1^2) / (2 * (1 + (gamma - 1) / 2 * M1^2))) + ...
       (1 - g(2)^2) / (gamma * g(2)^2) + ((gamma + 1) / (2 * gamma)) * log(((gamma + 1) * g(2)^2) / (2 * (1 + (gamma - 1) / 2 * g(2)^2)));

% Fanno flow equation after the shock
H(2) = fMe * (L - g(1)) / (2*AR) - (1 - g(3)^2) / (gamma * g(3)^2) - ((gamma + 1) / (2 * gamma)) * log(((gamma + 1) * g(3)^2) / (2 * (1 + (gamma - 1) / 2 * g(3)^2))) + ...
       (1 - Me^2) / (gamma * Me^2) + ((gamma + 1) / (2 * gamma)) * log(((gamma + 1) * Me^2) / (2 * (1 + (gamma - 1) / 2 * Me^2)));

% Normal shock relation (Mach number transition)
H(3) = g(3) - sqrt(((gamma - 1) * g(2)^2 + 2) / (2 * gamma * g(2)^2 - (gamma - 1)));

end