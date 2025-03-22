function [x_skl, M1, p01] = Shock_Div_Nozzle(yT, GT,GE, radius, G, p0, T0, gamma, AR, TubeLenght, epsilon, pe)
    % NOZZLE_SHOCK Determines the position of the shock in the nozzle,
    % when pe < pc_1 && pe > pc_2.
    %
    % Inputs:
    %   - yT         : Throat radius
    %   - GT         : Throat grid index
    %   - radius     : Nozzle radius array
    %   - G          : Grid positions
    %   - p0, T0     : Stagnation pressure and temperature
    %   - gamma      : Specific heat ratio
    %   - AR         : Area Ratio (Ae/At)
    %   - TubeLength : Pipe length after the nozzle
    %   - f          : Fanno friction factor
    %   - pe         : Exit pressure
    %
    % Outputs:
    %   - x_sh  : Shock position
    %   - M1    : Mach number before the shock
    %   - p01   : Stagnation pressure before the shock
    %   - MA    : Mach number before the shock wave
    %   - MB    : Mach number after the shock wave
    
    R = 287; % Specific gas constant for air

    % ---- Compute **Critical Mass Flow Rate** ----
    m_crit = p0 / sqrt(T0) * sqrt(gamma / R) * (1 + ((gamma - 1) / 2))^(-0.5 * (gamma + 1) / (gamma - 1)) * 2 * yT;

    % ---- Solve for Exit Mach Number (Me) ----
    Me = fzero(@(M) -m_crit + pe * 2 * AR * M * sqrt(gamma / (R * (T0 / (1 + ((gamma - 1) / 2) * M^2)))), 0.5);
    Me = double(Me);

    % ---- Compute **Stagnation Pressure at Exit (p0e)** ----
    p0e = pe * (1 + ((gamma - 1) / 2) * Me^2)^(gamma / (gamma - 1));
    
    TE = T0 * (1+((gamma-1)/2)*Me^2)^-1;
    [ReM1,fMe]=Churchill(gamma,TE,pe,Me,epsilon,AR);

    % ---- Solve for Mach Number at Inlet of the Pipe (M1) ----
    M1 = fzero(@(M) -fMe * TubeLenght / (2 * AR) + ...
        ((1 - M^2) / (gamma * M^2) + (gamma + 1) / (2 * gamma) * log(((gamma + 1) * M^2) / (2 + (gamma - 1) * M^2))) - ...
        ((1 - Me^2) / (gamma * Me^2) + (gamma + 1) / (2 * gamma) * log(((gamma + 1) * Me^2) / (2 + (gamma - 1) * Me^2))), 0.5);
    M1 = double(M1);

    % ---- Compute **Stagnation Pressure Before (after?) the Shock (p01)** ----
    p01 = fzero(@(p) -p0e / p + M1 / Me * ((2 + (gamma - 1) * Me^2) / (2 + (gamma - 1) * M1^2))^((gamma + 1) / (2 * (gamma - 1))), p0);
    p01 = double(p01);

    % ---- Compute **Mach Number Before the Shock (M_A)** ----
    MA = fzero(@(M) (((gamma + 1) * M^2) / ((gamma - 1) * M^2 + 2))^(gamma / (gamma - 1)) * ...
                (((gamma + 1) / (2 * gamma * M^2 - (gamma - 1)))^(1 / (gamma - 1))) - (p01 / p0), 2);
    MA = double(MA);

    % ---- Compute **Mach Number After the Shock (M_B)** using Normal Shock Relations ----
    MB = sqrt(((gamma - 1) * MA^2 + 2) / (2 * gamma * MA^2 - (gamma - 1)));
    MB = double(MB);

    % ---- Determine the Shock Position (x_sh) ----
    % The shock occurs where the area ratio supports the upstream Mach number (M_A)
    A_skl = fzero(@(A) -A/yT + 1/MA * (2/(gamma+1)*(1+(gamma-1)/2*MA^2)) ^ (0.5*(gamma+1)/(gamma-1)), [0 10000]);
    Glenght=length(G);
    tol=10e-4;
    for i=GT:GE
        if abs(A_skl-radius(i))<tol
            x_skl=i;
        end
    end
end
