function GUI_Fanno_Flow_pipe(Dpipe, epsilon, L, P1, T1, V1, gamma, hpNozzle, hpPressure, hpMach, hpTemperature)
% Group introduction, all functions used
% The flow is adiabatic but friction along the tube walls is not negligible
% The Darcy friction factor f is approximated as constant based on
% conditions at the inlet and the Churchill equation is used to calculate f
%https://www.me.psu.edu/cimbala/me320/Lesson_Notes/Fluid_Mechanics_Lesson_15H.pdf
% Define constants
R = 287;    % J/(kg*K) Specific gas constant
mu_ref = 1.716e-5; % Reference viscosity at Tref (Ns/m²)
T_ref = 273.15;    % Reference temperature (K)
S = 110.4;         % Sutherland's constant (K)

% Compute dynamic viscosity using the Sutherland equation
mu = mu_ref * (T1 / T_ref)^(3/2) * (T_ref + S) / (T1 + S);

% Compute inlet properties
A = pi * Dpipe^2 / 4; % Pipe cross-sectional area
rho1 = P1 / (R * T1); % Inlet density [kg/m³]
a1 = sqrt(gamma * R * T1); % Speed of sound at inlet [m/s]
M1 = V1 / a1; % Inlet Mach number

% Compute Reynolds number
Re = (rho1 * Dpipe * V1) / (mu);

% Compute Darcy friction factor using Churchill equation
A = (-2.457*log((7/Re)^0.9+0.27*epsilon/Dpipe))^16;
B = (37530/Re)^16;
f = 8*((8/Re)^12+(A+B)^(-1.5))^(1/12); %f=4Cf

% Discretize pipe into grid points
N = 30000; % Number of grid points
x = linspace(0, L, N); % Pipe length discretized
Dx = L / (N - 1); % Grid spacing

% Initialize arrays for Mach number, pressure, temperature, density
M = zeros(1, N);
P = zeros(1, N);
T = zeros(1, N);
rho = zeros(1, N);
V = zeros(1,N);

% Set initial conditions at inlet
M(1) = M1;
P(1) = P1;
T(1) = T1;
rho(1) = rho1;
V(1) = V1;

% Fanno flow relations as anonymous functions of Mach number
M_fcn = @(M) 1 + (gamma-1)/2 * M.^2;
p_ratio = @(M) sqrt((gamma+1)/2 ./ M_fcn(M)) ./ M;
T_ratio = @(M) (gamma+1)/2 ./ M_fcn(M);
delta_s_R = @(M) log(M .* ((gamma+1)/2 ./ M_fcn(M)).^((gamma+1)/2/(gamma-1)));
fanno = @(M) (1 - M.^2)./(gamma*M.^2) + (gamma+1)/(2*gamma) * log((gamma+1)/2 .* M.^2 ./ M_fcn(M));
% Sonic outlet pressure [Pa]
p_star = P1/p_ratio(M1);
% Sonic outlet temperature [K]
T_star = T1/T_ratio(M1);
% Sonic outlet specific enthalpy [J/kg]
h_star = (gamma*R/(gamma-1))*T_star;
% Sonic outlet density [kg/m^3]
rho_star = p_star/(R*T_star);
% Pipe length to reach sonic outlet condition from inlet condition [m]
L_star = fanno(M(1)) * Dpipe / f;
tol=0.5;
% Solve Fanno flow equations iteratively along the pipe
for i = 1:N-1
    A = (-2.457*log((7/Re)^0.9+0.27*epsilon/Dpipe))^16;
    B = (37530/Re)^16;
    f = 8*((8/Re)^12+(A+B)^(-1.5))^(1/12);
   % Compute the term in the differential equation
    dM2_dx = (gamma * M(i)^4 * (1 + ((gamma - 1) / 2) * M(i)^2)) / (1 - M(i)^2);
    
    % Update Mach number squared using discretized form
    M_squared = M(i)^2 + (f * Dx / Dpipe) * dM2_dx;
    
    % Ensure M(i+1) is real and positive
    if M_squared < 0
        hThroat = uicontrol('Style','text','String','Unphysical Mach number calculated. Check input values.','Position',[10,100,400,20]);
        set(hThroat,'FontSize',10);
        error("Unphysical Mach number calculated. Check input values.");
    end
    M(i+1) = sqrt(M_squared);  % Solve for M(i+1)    
    
    % Compute temperature using ratio equation
    T(i+1)= ((1 + ((gamma - 1) / 2) * M(i)^2) / (1 + ((gamma - 1) / 2) * M(i+1)^2)) * T(i);    
    a2=sqrt(gamma * R * T(i+1));
    V(i+1) = a2*M(i+1);
    % Compute density using ideal gas law
    rho(i+1) = rho(i)*(V(i)/V(i+1));
    Re = (rho(i+1) * Dpipe * V(i+1)) / (mu);
    % Compute pressure using isentropic relation
    P(i+1) = rho(i+1) * R * T(i+1);
end
    M_range = linspace(0.15, 3, 100);
    % Compute entropy and enthalpy for the Fanno line
    entropy = delta_s_R(M_range) * R;
    enthalpy = T_ratio(M_range) * T_star * (gamma * R / (gamma - 1));
    
    % Compute the entry and exit points
    entry_point = [delta_s_R(M1) * R, T_ratio(M1) * T_star * (gamma * R / (gamma - 1))];
    exit_point = [delta_s_R(M(end)) * R, T_ratio(M(end)) * T_star * (gamma * R / (gamma - 1))];

    % Plot the Fanno Line
    figure;
    plot(entropy, enthalpy, 'b', 'LineWidth', 1.5);
    hold on;
    scatter(entry_point(1), entry_point(2), 80, 'ro', 'filled'); % Entry point (red)
    scatter(exit_point(1), exit_point(2), 80, 'go', 'filled');   % Exit point (green)
    grid on;
    xlabel("Specific Entropy [J/(kg*K)]");
    ylabel("Specific Enthalpy [J/kg]");
    title("Fanno Line");
    legend("Fanno Line", "Entry Point", "Exit Point");
    hold off;

 % Plot results
    axes(hpNozzle);
    hold on;
    plot(x, Dpipe/2 * ones(size(x)), 'b', 'LineWidth', 2); % Linha superior
    plot(x, -(Dpipe/2) * ones(size(x)), 'b', 'LineWidth', 2); % Linha inferior
    xlabel('Pipe Length (m)');
    ylabel('Diameter (m)');
    title('Pipe Representation');
    ylim([-Dpipe, Dpipe]); % Ajusta os eixos para melhor visualização
    grid on;
    hold off;

    axes(hpPressure);
    plot(x, P, 'r', 'LineWidth', 2);
    xlabel('Pipe Length (m)');
    ylabel('Pressure (Pa)');
    ylim([0, P(1)]); % Ajusta os eixos para melhor visualização
    grid on;

    axes(hpMach);
    plot(x, M, 'g', 'LineWidth', 2);
    xlabel('Pipe Length (m)');
    ylabel('Mach Number');
    ylim([0, 1.2*M(end)]); % Ajusta os eixos para melhor visualização
    grid on;

    axes(hpTemperature);
    plot(x, T, 'm', 'LineWidth', 2);
    xlabel('Pipe Length (m)');
    ylabel('Temperature (K)');
    ylim([0.99*T(end) , 1.01*T(1)]); % Ajusta os eixos para melhor visualização
    grid on;

% Create a text box in the GUI to display L_star
hLstar = uicontrol('Style', 'text', ...
                   'String', sprintf('Sonic Length (L*): %.3f meters', L_star), ...
                   'Position', [10, 80, 250, 20], ...
                   'FontSize', 10, ...
                   'BackgroundColor', [0.95 0.95 0.95]);
set(hLstar,'FontSize',10);
end