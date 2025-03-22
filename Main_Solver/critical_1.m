function [pc_1, flag_1] = critical_1(T0, p0, gamma, AR, TubeLenght, epsilon)
    % **Function Description:**
    % Computes the **critical pressure (pc_1)** at the **throat** and determines if 
    % the flow is **choked or subsonic** at the exit of the pipe.
    % Inputs:
    %   - p0        : Total stagnation pressure (Pa)
    %   - gamma     : Specific heat ratio
    %   - AR        : Area Ratio (Ae/At)
    %   - TubeLenght: Length of the pipe (m)
    %   - f         : Fanno friction factor
    % Outputs:
    %   - pc_1      : Critical pressure (Pa) at the throat
    %   - flag_1    : 0 (if flow is subsonic), 1 (if the flow is doubly choked)

    % ---- Compute the **Subsonic Inlet Mach Number (M1)** ----
    f_inlet = @(M) ((1/M^2) * ((2/(gamma+1)) * (1+((gamma-1)/2)*M^2))^((gamma+1)/(gamma-1)) - (AR)^2);
    % Use Bisection Method to solve for M1 (Ensure correct bounds)
    M1 = Bissection(f_inlet, 0.01, 1, 100);  % M1 should be **subsonic**
    M1 = double(M1);
    T1 = T0 * (1+((gamma-1)/2)*M1^2)^-1;       % eq 3.28
    P1 = p0 * (T0/T1)^(-gamma/(gamma-1));   % eq 3.29
    [ReM1,f] = Churchill(gamma,T1,P1,M1,epsilon,AR);
    % ---- Compute the **Exit Mach Number (Me)** using Fanno Flow Equation ----
    f_exit = @(M) -f*TubeLenght/(2*AR) + ...
                  ((1-M1^2)/(gamma*M1^2) + (gamma+1)/(2*gamma) * log(((gamma+1)*M1^2)/(2+(gamma-1)*M1^2))) - ...
                  ((1-M^2)/(gamma*M^2) + (gamma+1)/(2*gamma) * log(((gamma+1)*M^2)/(2+(gamma-1)*M^2)));  %eq 3.107
    
    c = f_exit(0.01); %Set c equal to function at a
    d = f_exit(1); %Set d equal to function at b
    if c*d > 0.0 %If c times d is positive, bisection method will not work.
        Me=10; %Prevents interrupts
    else
    % Use Bisection Method to solve for Me (Adjust range for subsonic/supersonic)
    Me = Bissection(f_exit, 0.01, 1, 100);
    Me = double(Me);
    end

    % ---- Check if Flow is **Subsonic or Choked** at Pipe Exit ----
    if (Me > 1) || isnan(Me)
        cond = 1; % **Doubly choked flow**
    else
        cond = 0; % **Subsonic exit**
    end

    % ---- Compute **Critical Pressure (pc_1)** ----
    if cond == 0
        flag_1 = 0;  % **Subsonic exit**
        
        % Temperature ratio: Te/T1 (Using Fanno Relations) (eq 3.98)
        Te_T1 = (2+(gamma-1)*M1.^2) / (2+(gamma-1)*Me.^2);

        % Compute exit stagnation pressure p0e (Anderson Eq. 3.102)
        p0e = p0 .* M1 ./ Me .* (1/Te_T1).^((gamma+1)/(2*(gamma-1)));

        % Compute **Minimum Allowable Static Pressure (pc_1)**  (eq 3.30)
        pc_1 = p0e ./ (1+((gamma-1)./2).*Me.^2).^(gamma/(gamma-1));

    else
        flag_1 = 1;  % **Doubly choked flow warning**

        % In this case, the exit Mach is **assumed to be choked (Me = 1)**
        Te_T1 = (1+((gamma-1)/2)*M1.^2) / (1+((gamma-1)/2)*1.^2);
        p0e = p0 .* M1 ./ 1 .* (1/Te_T1).^((gamma+1)/(2*(gamma-1)));

        % Compute **Minimum Allowable Static Pressure (pc_1) with Me = 1**
        pc_1 = p0e ./ (1+((gamma-1)./2).*1.^2).^(gamma/(gamma-1)); %The lowest pressure the user can input
    end

end
