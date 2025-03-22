function [Re,f] = Churchill(gamma,T,P,M,epsilon,AR)

R = 287;    % J/(kg*K) Specific gas constant
mu_ref = 1.716e-5; % Reference viscosity at Tref (Ns/m²)
T_ref = 273.15;    % Reference temperature (K)
S = 110.4;         % Sutherland's constant (K)

% Compute dynamic viscosity using the Sutherland equation
mu = mu_ref * (T / T_ref)^(3/2) * (T_ref + S) / (T + S);
rho = P / (R * T); % Inlet density [kg/m³] %Assuming ideal gas
a = sqrt(gamma * R * T); % Speed of sound at inlet [m/s]
V= M*a;

Re = (rho * 2 * AR * V) / (mu);
A = (-2.457*log((7/Re)^0.9+0.27*epsilon/(2*AR)))^16;
B = (37530/Re)^16;
f = 8*((8/Re)^12+(A+B)^(-1.5))^(1/12); %f=4Cf
end