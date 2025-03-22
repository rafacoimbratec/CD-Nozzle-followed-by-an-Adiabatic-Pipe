function [x_sh, M_B, outshock, flag] = L_shock(p0, T0, gamma, pe, AR, TubeLength, epsilon, yT, GT, radius)
% L_INF Determines the location of the shock inside the pipe or at the exit.
%
% Inputs:
%   - p0         : Stagnation pressure (Pa)
%   - T0         : Stagnation temperature (K)
%   - gamma      : Specific heat ratio
%   - pe         : Exit static pressure (Pa)
%   - AR         : Area ratio (Ae/At)
%   - TubeLength : Length of the pipe after nozzle (m)
%   - f          : Darcy friction factor
%   - yT         : Throat radius (m)
%   - GT         : Throat grid index
%   - radius     : Nozzle radius array
%
% Outputs:
%   - x_sh     : Shock position (m)
%   - M_B      : Mach number after the shock
%   - outshock : Boolean (true if the shock occurs outside the pipe)
%   - flag     : Indicates different flow conditions:
%       1 = Normal shockwave at the exit
%       2 = Overexpanded flow with oblique shock outside
%       3 = No shock, perfectly expanded flow
%       4 = Underexpanded flow with expansion waves
%       0 = Shock occurs inside the tube

R = 287; % Specific gas constant for air

% ---- Compute Critical Mass Flow Rate (m_crit) ----
m_crit = p0 / sqrt(T0) * sqrt(gamma / R) * (1 + ((gamma - 1) / 2))^(-0.5 * (gamma + 1) / (gamma - 1)) * 2 * yT;
m_crit = double(m_crit);

% ---- Solve for Supersonic Mach Number at Tube Inlet (M1) ----
M1 = fzero(@(M) -AR + 1/M * (2 / (gamma + 1) * (1 + ((gamma - 1) / 2) * M^2))^(0.5 * (gamma + 1) / (gamma - 1)), [1, 100]);
M1 = double(M1);

T1 = T0 * (1+((gamma-1)/2)*M1^2)^-1;
p1 = p0 * (T0/T1)^(-gamma/(gamma-1));
[Rem1,fM1]=Churchill(gamma,T1,p1,M1,epsilon,AR);

% ---- Solve for Mach Number Before the Shock (M_A) ----
M_A = fzero(@(M) -fM1 * TubeLength / (2 * AR) + ...
    ((1 - M1^2) / (gamma * M1^2) + (gamma + 1) / (2 * gamma) * log(((gamma + 1) * M1^2) / (2 + (gamma - 1) * M1^2))) - ...
    ((1 - M^2) / (gamma * M^2) + (gamma + 1) / (2 * gamma) * log(((gamma + 1) * M^2) / (2 + (gamma - 1) * M^2))), [1, 100]);
M_A = double(M_A);

% ---- Compute Static Pressure Before the Shock (p_A) ----
p_A = fzero(@(p) -m_crit + p * 2 * AR * M_A * sqrt(gamma / (R * (T0 / (1 + ((gamma - 1) / 2) * M_A^2)))), [0, 1e9]);
p_A = double(p_A);

% ---- Compute Critical Exit Pressure for Shock at the Outlet (pe_sk) ----
pe_sk = p_A * (1 + (2 * gamma / (gamma + 1)) * (M_A^2 - 1));% Eq. 3.57 Anderson
fprintf('critical_3 = %.10f\n', pe_sk);
fprintf('critical_4 = %.15f\n', p_A);
% ---- Check if Shock is at the Exit or Inside the Tube ----
if pe <= pe_sk
    % The shock occurs at the pipe exit
    x_sh = -1;  % No internal shock
    M_B = 0;
    outshock = true;

    % Determine type of shock interaction at the exit
    if pe == pe_sk
        flag = 1; % Normal shockwave at the tube outlet
    elseif (pe < pe_sk) && (pe > p_A)
        flag = 2; % Overexpanded flow with oblique shock outside
    elseif pe==p_A
        flag = 3; % No shock, perfectly expanded flow
    elseif pe < p_A
        flag = 4; % Underexpanded flow with expansion waves outside
    end

else
    % The shock occurs inside the tube

    % ---- Solve for Mach Number at Tube Outlet (Me) ----
    Me = fzero(@(M) -m_crit + pe * 2 * AR * M * sqrt(gamma / (R * (T0 / (1 + ((gamma - 1) / 2) * M^2)))), 0.5);
    Me = double(Me);

    % ---- Solve for Shock Location and Mach Number After the Shock ----
    shock_guess = [0.5, 1.5, 0.5]; % Initial guess: [x_shock, M_A, M_B]

    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-6);
    shock_lmax = fsolve(@(g) M_shock(g, T0,p0,pe, gamma, M1, TubeLength, epsilon, AR, Me), shock_guess, options);
    % ---- Assign Shock Location and Post-Shock Mach Number ----
    x_sh = shock_lmax(1) + 2 * yT;
    M_B = shock_lmax(3);
    outshock = false;
    flag = 0; % Shock inside the tube
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