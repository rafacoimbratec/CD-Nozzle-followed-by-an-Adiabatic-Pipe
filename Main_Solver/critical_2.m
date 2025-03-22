function [pc_2, flag_2] = critical_2(T0, p0, gamma, AR, TubeLenght, epsilon)
    % **Function Description:**
    % Computes the **critical outside pressure (pc_2)** that causes a shock 
    % wave at the inlet of the tube.
    % 
    % Inputs:
    %   - p0        : Total stagnation pressure (Pa)
    %   - gamma     : Specific heat ratio
    %   - AR        : Area Ratio (Ae/At)
    %   - TubeLenght: Length of the pipe (m)
    %   - f         : Fanno friction factor
    % Outputs:
    %   - pc_2      : Critical pressure (Pa) to initiate shock wave
    %   - flag_2    : 0 (if shock is at the inlet), 1 (if no valid solution)

    % ---- Compute the **Mach Number Before Shock (M_A)** ----
    f_inlet = @(M) -AR + (1/M) * (2/(gamma+1) * (1 + ((gamma-1)/2)*M^2))^(0.5*(gamma+1)/(gamma-1));

    % Solve for M_S using Bisection method (supersonic range)
    M_S = Bissection(f_inlet, 1, 100, 100);
    M_S = double(M_S);

    % ---- Compute the **Mach Number After the Shock (M1)** ----
    M1 = sqrt(((gamma-1)*M_S^2 + 2) / (2*gamma*M_S^2 - (gamma-1)));
    T1 = T0 * (1+((gamma-1)/2)*M1^2)^-1;
    P1 = p0 * (T0/T1)^(-gamma/(gamma-1));
    [ReM1,f] = Churchill(gamma,T1,P1,M1,epsilon,AR);
    % ---- Compute the **Mach Number at the Outlet (Me)** using Fanno Flow Equation ----
    f_exit = @(M) -f*TubeLenght/(2*AR) + ...
                  ((1-M1^2)/(gamma*M1^2) + (gamma+1)/(2*gamma) * log(((gamma+1)*M1^2)/(2+(gamma-1)*M1^2))) - ...
                  ((1-M^2)/(gamma*M^2) + (gamma+1)/(2*gamma) * log(((gamma+1)*M^2)/(2+(gamma-1)*M^2)));

    % Solve for Me using Bisection method (valid range: subsonic)
    c = f_exit(0.01); %Set c equal to function at a
    d = f_exit(1); %Set d equal to function at b
    if c*d > 0.0 %If c times d is positive, bisection method will not work.
        Me=10; %Prevents interrupts
    else
    % Use Bisection Method to solve for Me (Adjust range for subsonic/supersonic)
    Me = Bissection(f_exit, 0.01, 1, 100);
    Me = double(Me);
    end

    % ---- Check for **Valid Mach Solution** ----
    if (Me > 1) || isnan(Me)
        cond = 1;  % **Invalid conditions (choked flow)**
    else
        cond = 0;  % **Valid solution exists**
    end
    % ---- Compute **Critical Pressure (pc_2)** ----
    if cond == 0
        flag_2 = 0;  % **Shock occurs at the inlet**
        
        % Compute stagnation pressure after the shock (p02)
        p02 = p0 * (((gamma+1)/2*M_S^2) / (1+(gamma-1)/2*M_S^2))^(gamma/(gamma-1)) * ...
                    (1 / (2*gamma/(gamma+1) * M_S^2 - (gamma-1)/(gamma+1)))^(1/(gamma-1));

        % Temperature ratio Te/T1 (Using Fanno Relations)
        Te_T1 = (1 + ((gamma-1)/2)*M1^2) / (1 + ((gamma-1)/2)*Me^2);

        % Compute stagnation pressure at the outlet (p0e)
        p0e = p02 * M1 ./ Me * (Te_T1).^((gamma+1)/(2*(1-gamma))); %Eq 3.102

        % Compute **Critical Exit Pressure (pc_2)**
        pc_2 = p0e ./ (1 + ((gamma-1)/2)*Me^2).^(gamma/(gamma-1));

    else
        flag_2 = 1;  % **No valid solution (doubly choked or unrealistic inputs)**

        % Compute stagnation pressure after the shock (p01)
        p01 = p0 * (((gamma+1)/2*M1^2) / (1+(gamma-1)/2*M1^2))^(gamma/(gamma-1)) * ...
                    (1 / (2*gamma/(gamma+1) * M1^2 - (gamma-1)/(gamma+1)))^(1/(gamma-1));

        % Compute Temperature ratio assuming exit Mach = 1 (choked condition)
        Te_T1 = (1 + ((gamma-1)/2)*M1^2) / (1 + ((gamma-1)/2)*1^2);

        % Compute stagnation pressure at the outlet assuming Me = 1
        p0e = p01 * M1 / 1 * (Te_T1).^((gamma+1)/(2*(1-gamma)));

        % Compute **Minimum Allowable Static Pressure (pc_2) assuming Me = 1**
        pc_2 = p0e ./ (1 + ((gamma-1)/2)*1^2).^(gamma/(gamma-1));
    end

end
