function GUI_Fanno_Flow_Nozzle2(AR, p0, T0, pe, TubeLenght, epsilon, gamma, hpNozzle, hpPressure, hpMach, hpTemperature)
NLenght=2;
R=287;
G = (0:1e-3:NLenght)';
GLenght = length(G);
[yT, GT, radius,G, GE] = NozzleDraw(AR,G,GLenght,NLenght,1,TubeLenght);
GLenght = length(G);
Dx = (2+TubeLenght)/length(G);
%% First Boundaries for outside pressure
    %CRITICAL PRESSURE FOR M=1 AT THE THROAT
    [pc_1, flag_1] = critical_1(T0, p0, gamma, AR, TubeLenght, epsilon);
    fprintf('critical_1 = %d\n', pc_1); % Using fprintf for better formatting

    % CRITICAL PRESSURE FOR SHOCK WAVE AT THE INLET
    [pc_2, flag_2] = critical_2(T0, p0, gamma, AR, TubeLenght, epsilon);
    fprintf('critical_2 = %d\n', pc_2); % Displaying pc_2 properly
%%
    hShock = uicontrol('Style','text','String',' ','Position',[10,100,200,20]);
    set(hShock,'FontSize',10);   
    if(pe>pc_1) || ((pe==pc_1) && (flag_1==0))
       [M(1), M(2), flag_area]=Subsonic(T0,p0,gamma,AR,TubeLenght,epsilon,pe);
       if flag_area==1
             msgbox('The specified Aspect Ratio is too small for accurate computation. Please increase the Aspect Ratio and try again.')
            return
       end
         %MASS FLOW RATE, CONSTANT
         m= 2*AR * p0/sqrt(T0) * sqrt(gamma/R) * M(1) * (1+(gamma-1)/2 * M(1)^2) ^ (-0.5 * (gamma+1)/(gamma-1));
         M = zeros(1, GLenght); % Mach number along the entire domain
         T = zeros(1, GLenght); % Temperature along the entire domain
         p = zeros(1, GLenght); % Pressure along the entire domain
         f = zeros(1,GLenght);
         Re = zeros(1,GLenght);
            %Plots
         for n = 1:GE
            % Compute Mach using fzero() to solve mass flow equation
            M(n) = fzero(@(Ma) -m + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), 0.6);
            if isnan(M(n))
             M(n) = fzero(@(Ma) -m + 2* radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), 0.9);
            end
            % Compute Temperature using isentropic relations
            T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;

            % Compute Pressure using isentropic relations
            p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
         end
          for n = GE:GLenght-1
            [Re(n),f(n)]=Churchill(gamma,T(n),p(n),M(n),epsilon,AR);
            dM2_dx = (gamma * M(n)^4 * (1 + ((gamma - 1) / 2) * M(n)^2)) / (1 - M(n)^2);
             % Update Mach number squared using discretized form
            M_squared = M(n)^2 + (f(n) * Dx / (2*AR)) * dM2_dx;
             % Ensure M(i+1) is real and positive
            if M_squared < 0
            error("Unphysical Mach number calculated. Check input values.");
            end
            M(n+1) = sqrt(M_squared);  % Solve for M(i+1)
            % Compute temperature using ratio equation
            T(n+1)= ((1 + ((gamma - 1) / 2) * M(n)^2) / (1 + ((gamma - 1) / 2) * M(n+1)^2)) * T(n);    
            % Compute pressure using Fanno Flow properties
            p(n+1) = p(n) * (M(n)/M(n+1)) * (T(n+1)/T(n))^(0.5);
          end
            %Nozzle Plot
            axes(hpNozzle)
            hold on;
            cla(hpNozzle);
            % Plot the upper half of the nozzle
            plot(G, radius, 'b', 'LineWidth', 2);
            % Plot the lower half for symmetry
            plot(G, -radius, 'b', 'LineWidth', 2);
            % Add labels and title
            plot(G(GE)*ones(1,10), linspace(-AR,AR,10), 'k', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Radius (m)');
            if radius(end) > 3
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-1.5*radius(end), 1.5*radius(end)]); % Y-axis from -10 to 10
            else
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-4, 4]); % Y-axis from -10 to 10
            end
            hold off;

            hold on;
            % Mach Number Plot
            axes(hpMach);
            cla(hpMach);
            plot(G, M, 'r', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Mach Number');
            hold off;

            % Temperature Plot
            hold on;
            axes(hpTemperature);
            cla(hpTemperature);
            plot(G, T, 'b', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Temperature (K)');
            hold off;

            % Pressure Plot
            hold on;
            axes(hpPressure);
            cla(hpPressure);
            plot(G, p, 'g', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Pressure (Pa)');
            hold off;

            figure;
            hold on;
            plot(G(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Friction Factor (f)');
            title('Friction Factor Variation');
            hold off;


            figure;
            hold on;
            plot(Re(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Reynolds Number (Re)');
            ylabel('Friction factor(f)');
            title('Friction factor(f) x Reynolds Number (Re)');
            hold off;

            %End plots
    elseif (pe<=pc_1) && (flag_1==1)
        msgbox('No valid solution for the given inputs. The pipe is experiencing double choking. Please adjust the parameters.')
        return
    elseif (pe<pc_1) && (pe>=pc_2)
        [x_skl, M1, p01] = Shock_Div_Nozzle(yT, GT,GE, radius, G, p0, T0, gamma, AR, TubeLenght, epsilon, pe);
        M = zeros(1, GLenght); % Mach number along the entire domain
        T = zeros(1, GLenght); % Temperature along the entire domain
        p = zeros(1, GLenght); % Pressure along the entire domain
        f = zeros(1,GLenght);
        Re = zeros(1,GLenght);
        m_crit = p0 / sqrt(T0) * sqrt(gamma / R) * (1 + ((gamma - 1) / 2))^(-0.5 * (gamma + 1) / (gamma - 1)) * 2 * yT;
        m_crit = double(m_crit);
        for n = 1:GT
            % Compute Mach using fzero() to solve mass flow equation
            M(n) = fzero(@(Ma) -m_crit + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), [0.0001,1]);
            % Compute Temperature using isentropic relations
            T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;

            % Compute Pressure using isentropic relations
            p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
        end
        for n = GT:x_skl
            % Compute Mach using fzero() to solve mass flow equation
            M(n) = fzero(@(Ma) -m_crit + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), [1,1000]);
            % Compute Temperature using isentropic relations
            T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;

            % Compute Pressure using isentropic relations
            p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
        end
        for n = x_skl:GE
            AR_new = sqrt((1/(M1^2))*((2/(gamma+1))*(1+((gamma-1)/2)*M1^2))^((gamma+1)/(gamma-1)));
            A_2star= yT*AR/AR_new;
            % Compute Mach using fzero() to solve mass flow equation
            fun = @(M) (1/M^2)*((2/(gamma+1))*(1+((gamma-1)/2)*M^2))^((gamma+1)/(gamma-1))-(radius(n)/A_2star)^2;
            M(n) = Bissection(fun,0.001,1,100);
            % Compute Temperature using isentropic relations
            T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;
            % Compute Pressure using isentropic relations
            p(n) = p01 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
        end
        for n = GE:GLenght-1
            [Re(n),f(n)] = Churchill(gamma,T(n),p(n),M(n),epsilon,AR);
            dM2_dx = (gamma * M(n)^4 * (1 + ((gamma - 1) / 2) * M(n)^2)) / (1 - M(n)^2);
             % Update Mach number squared using discretized form
            M_squared = M(n)^2 + (f(n) * Dx / (2*AR)) * dM2_dx;
             % Ensure M(i+1) is real and positive
            if M_squared < 0
            error("Unphysical Mach number calculated. Check input values.");
            end
            M(n+1) = sqrt(M_squared);  % Solve for M(i+1)
            % Compute temperature using ratio equation
            T(n+1)= ((1 + ((gamma - 1) / 2) * M(n)^2) / (1 + ((gamma - 1) / 2) * M(n+1)^2)) * T(n);    
            % Compute pressure using Fanno Flow properties
            p(n+1) = p(n) * (M(n)/M(n+1)) * (T(n+1)/T(n))^(0.5);
        end
             %Nozzle Plot
            axes(hpNozzle)
            hold on;
            cla(hpNozzle);
            % Plot the upper half of the nozzle
            plot(G, radius, 'b', 'LineWidth', 2);
            % Plot the lower half for symmetry
            plot(G, -radius, 'b', 'LineWidth', 2);
            % Add labels and title
            plot(G(x_skl)*ones(1,10),linspace(-radius(x_skl),radius(x_skl),10),'r', 'LineWidth', 2);
            plot(G(GE)*ones(1,10), linspace(-AR,AR,10), 'k', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Radius (m)');
            if radius(end) > 3
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-1.5*radius(end), 1.5*radius(end)]); % Y-axis from -10 to 10
            else
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-4, 4]); % Y-axis from -10 to 10
            end
            hold off;
            hold on;
            % Mach Number Plot
            axes(hpMach);
            cla(hpMach);
            plot(G, M, 'r', 'LineWidth', 2);
            plot(G(GE)*ones(1,10), linspace(-AR,AR,10), 'k', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Mach Number');
            hold off;
            hold on;
            % Temperature Plot
            axes(hpTemperature);
            cla(hpTemperature);
            plot(G, T, 'b', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Temperature (K)');
            hold off;
            hold on;
            % Pressure Plot
            axes(hpPressure);
            cla(hpPressure);
            plot(G, p, 'g', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Pressure (Pa)');
            hold off;

            figure;
            hold on;
            plot(G(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Friction Factor (f)');
            title('Friction Factor Variation');
            hold off;

            figure;
            hold on;
            plot(Re(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Reynolds Number (Re)');
            ylabel('Friction factor(f)');
            title('Friction factor(f) x Reynolds Number (Re)');
            hold off;
    elseif (pe < pc_2) && (flag_2==1)
        msgbox('No valid solution for the given inputs. The pipe is experiencing double choking. Please adjust the parameters.')
        return
    else
        %Supersonic flow at the tube inlet, lets first calculate the Mach1
        %inlet
        %and then deal with tube Lenght
        Minlet = fzero(@(M) -AR + 1/M * (2/(gamma+1)*(1+((gamma-1)/2)*M^2))^(0.5*(gamma+1)/(gamma-1)), [1 100]);
        Minlet=double(Minlet);
        %Now lets calculate LMAX so that with Minlet we get at the end of
        %the pipe M=1
        var=(1-Minlet^2)/(gamma*Minlet^2) + (gamma+1)/(2*gamma) * log(((gamma+1)*Minlet^2)/(2+(gamma-1)*Minlet^2));
        T1 = T0 * (1+((gamma-1)/2)*Minlet^2);
        p1 = p0 * (T0/T1)^(-gamma/(gamma-1));
        [Re1,fM1]=Churchill(gamma,T1,p1,Minlet,epsilon,AR);
        L_max = (var*2*AR)/(fM1);
        L_max
        if (TubeLenght<L_max) %Short tube
        [x_sh_tube, M_B, outshock,flag]=L_shock(p0,T0,gamma,pe,AR,TubeLenght,epsilon,yT,GT,radius);
            if flag==1
                 %Normal Shock
                hShock = uicontrol('Style','text','String','NSW tube outlet','Position',[10,100,200,20]);
                set(hShock,'FontSize',10);   
            elseif flag==2
                hShock = uicontrol('Style','text','String','Oblique shock outside, overexpanded flow','Position',[10,100,200,20]);
                set(hShock,'FontSize',10);
            elseif flag==3
                 hShock = uicontrol('Style','text','String','No shock, perfectly expanded flow','Position',[10,100,200,20]);
                set(hShock,'FontSize',10);
            elseif flag==4
                hShock = uicontrol('Style','text','String','Expansion fans outside, underexpanded flow','Position',[10,100,200,20]);
                set(hShock,'FontSize',10);
            end
            if (x_sh_tube>0) && (x_sh_tube<=TubeLenght+2*yT)
                    %PLOT %SHOCK INSIDE THE SHORT TUBE
                    M = zeros(1, GLenght); % Mach number along the entire domain
                    T = zeros(1, GLenght); % Temperature along the entire domain
                    p = zeros(1, GLenght); % Pressure along the entire domain
                    f = zeros(1, GLenght);
                    Re = zeros(1,GLenght);
                    m_crit = p0 / sqrt(T0) * sqrt(gamma / R) * (1 + ((gamma - 1) / 2))^(-0.5 * (gamma + 1) / (gamma - 1)) * 2 * yT;
                    m_crit = double(m_crit);

                    for n = 1:GT
                    % Compute Mach using fzero() to solve mass flow equation
                    M(n) = fzero(@(Ma) -m_crit + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), [0.0001,1]);
                    % Compute Temperature using isentropic relations
                    T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;

                    % Compute Pressure using isentropic relations
                    p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
                    end

                    for n = GT:GE
                    % Compute Mach using fzero() to solve mass flow equation
                    M(n) = fzero(@(Ma) -m_crit + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), [1,1000]);
                    % Compute Temperature using isentropic relations
                    T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;

                    % Compute Pressure using isentropic relations
                    p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
                    end
                    l=int16((x_sh_tube*GLenght)/(NLenght+TubeLenght));
                    for n=GE:l
                    [Re(n),f(n)] = Churchill(gamma,T(n),p(n),M(n),epsilon,AR);
                    dM2_dx = (gamma * M(n)^4 * (1 + ((gamma - 1) / 2) * M(n)^2)) / (1 - M(n)^2);
                     % Update Mach number squared using discretized form
                    M_squared = M(n)^2 + (f(n) * Dx / (2*AR)) * dM2_dx;
                     % Ensure M(i+1) is real and positive
                    if M_squared < 0
                    error("Unphysical Mach number calculated. Check input values.");
                    end
                    M(n+1) = sqrt(M_squared);  % Solve for M(i+1)
                    % Compute temperature using ratio equation
                    T(n+1)= ((1 + ((gamma - 1) / 2) * M(n)^2) / (1 + ((gamma - 1) / 2) * M(n+1)^2)) * T(n);    
                    % Compute pressure using Fanno Flow properties
                    p(n+1) = p(n) * (M(n)/M(n+1)) * (T(n+1)/T(n))^(0.5);
                    end
                    M(l+1)=M_B;
                    T(l+1)=T(l)*(2 + (gamma - 1) * M(l)^2) * ((2 * gamma * M(l)^2 - (gamma - 1)) / ((gamma + 1)^2 * M(l)^2));
                    p(l+1)= p(l)*(1+(2*gamma)/(gamma+1)*(M(l)^2-1));

                    for n=l+1:GLenght-1
                    [Re(n),f(n)] = Churchill(gamma,T(n),p(n),M(n),epsilon,AR);
                    dM2_dx = (gamma * M(n)^4 * (1 + ((gamma - 1) / 2) * M(n)^2)) / (1 - M(n)^2);
                     % Update Mach number squared using discretized form
                    M_squared = M(n)^2 + (f(n) * Dx / (2*AR)) * dM2_dx;
                     % Ensure M(i+1) is real and positive
                    if M_squared < 0
                    error("Unphysical Mach number calculated. Check input values.");
                    end
                    M(n+1) = sqrt(M_squared);  % Solve for M(i+1)
                    % Compute temperature using ratio equation
                    T(n+1)= ((1 + ((gamma - 1) / 2) * M(n)^2) / (1 + ((gamma - 1) / 2) * M(n+1)^2)) * T(n);    
                    % Compute pressure using Fanno Flow properties
                    p(n+1) = p(n) * (M(n)/M(n+1)) * (T(n+1)/T(n))^(0.5);
                    end
            %Nozzle Plot 
            axes(hpNozzle)
            hold on;
            cla(hpNozzle);
            % Plot the upper half of the nozzle
            plot(G, radius, 'b', 'LineWidth', 2);
            % Plot the lower half for symmetry
            plot(G, -radius, 'b', 'LineWidth', 2);
            % Add labels and title
            plot(x_sh_tube*ones(1,10),linspace(-AR,AR,10),'r', 'LineWidth', 2);
            plot(G(GE)*ones(1,10), linspace(-AR,AR,10), 'k', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Radius (m)');
            if radius(end) > 3
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-1.5*radius(end), 1.5*radius(end)]); % Y-axis from -10 to 10
            else
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-4, 4]); % Y-axis from -10 to 10
            end
            hold off;
            hold on;
             % Mach Number Plot
            axes(hpMach);
            cla(hpMach);
            plot(G, M, 'r', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Mach Number');
            hold off;
            hold on;
            % Temperature Plot
            axes(hpTemperature);
            cla(hpTemperature);
            plot(G, T, 'b', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Temperature (K)');
            hold off;
            hold on;
            % Pressure Plot
            axes(hpPressure);
            cla(hpPressure);
            plot(G, p, 'g', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Pressure (Pa)');
            hold off;

            figure;
            hold on;
            plot(G(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Friction Factor (f)');
            title('Friction Factor Variation');
            hold off;

            figure;
            hold on;
            plot(Re(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Reynolds Number (Re)');
            ylabel('Friction factor(f)');
            title('Friction factor(f) x Reynolds Number (Re)');
            hold off;

            else
            %No shock inside the tube
            M = zeros(1, GLenght); % Mach number along the entire domain
            T = zeros(1, GLenght); % Temperature along the entire domain
            p = zeros(1, GLenght); % Pressure along the entire domain
            m_crit = p0 / sqrt(T0) * sqrt(gamma / R) * (1 + ((gamma - 1) / 2))^(-0.5 * (gamma + 1) / (gamma - 1)) * 2 * yT;
            m_crit = double(m_crit);
            for n = 1:GT
            % Compute Mach using fzero() to solve mass flow equation
            M(n) = fzero(@(Ma) -m_crit + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), [0.0001,1]);
            % Compute Temperature using isentropic relations
            T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;

            % Compute Pressure using isentropic relations
            p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
            end

            for n = GT:GE
            % Compute Mach using fzero() to solve mass flow equation
            M(n) = fzero(@(Ma) -m_crit + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), [1,1000]);
            % Compute Temperature using isentropic relations
            T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;

            % Compute Pressure using isentropic relations
            p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
            end

            for n = GE:GLenght-1
            [Re(n),f(n)] = Churchill(gamma,T(n),p(n),M(n),epsilon,AR);
            dM2_dx = (gamma * M(n)^4 * (1 + ((gamma - 1) / 2) * M(n)^2)) / (1 - M(n)^2);
             % Update Mach number squared using discretized form
            M_squared = M(n)^2 + (f(n) * Dx / (2*AR)) * dM2_dx;
             % Ensure M(i+1) is real and positive
            if M_squared < 0
            error("Unphysical Mach number calculated. Check input values.");
            end
            M(n+1) = sqrt(M_squared);  % Solve for M(i+1)
            % Compute temperature using ratio equation
            T(n+1)= ((1 + ((gamma - 1) / 2) * M(n)^2) / (1 + ((gamma - 1) / 2) * M(n+1)^2)) * T(n);    
            % Compute pressure using Fanno Flow properties
            p(n+1) = p(n) * (M(n)/M(n+1)) * (T(n+1)/T(n))^(0.5);
            end

            cla(hShock);
            axes(hpNozzle)
            hold on;
            cla(hpNozzle);
            % Plot the upper half of the nozzle
            plot(G, radius, 'b', 'LineWidth', 2);
            % Plot the lower half for symmetry
            plot(G, -radius, 'b', 'LineWidth', 2);
            plot(G(GE)*ones(1,10), linspace(-AR,AR,10), 'k', 'LineWidth', 2);
            % Add labels and title
            plot(G(end)*ones(1,10),linspace(-radius(end),radius(end),10),'r', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Radius (m)');
            if radius(end) > 3
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-1.5*radius(end), 1.5*radius(end)]); % Y-axis from -10 to 10
            else
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-4, 4]); % Y-axis from -10 to 10
            end
            hold off;

            hold on;
            % Mach Number Plot
            axes(hpMach);
            cla(hpMach);
            plot(G, M, 'r', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Mach Number');
            hold off;

            hold on;
            % Temperature Plot
            axes(hpTemperature);
            cla(hpTemperature);
            plot(G, T, 'b', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Temperature (K)');
            hold off;

            hold on;
            % Pressure Plot
            axes(hpPressure);
            cla(hpPressure);
            plot(G, p, 'g', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Pressure (Pa)');
            hold off;
            
            figure;
            hold on;
            plot(G(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Friction Factor (f)');
            title('Friction Factor Variation');
            hold off;
            figure;
            hold on;
            plot(Re(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Reynolds Number (Re)');
            ylabel('Friction factor(f)');
            title('Friction factor(f) x Reynolds Number (Re)');
            hold off;
            end
        else %Long Tube
        [x_sh_l_tube, M_B]=L_sup(p0,T0,gamma,pe,AR,TubeLenght,epsilon,yT,GT,radius);
        if (x_sh_l_tube>2*yT) && (x_sh_l_tube<=TubeLenght+2*yT)
            %Normal shock inside the tube 
            %PLOT %SHOCK INSIDE THE SHORT TUBE
                    M = zeros(1, GLenght); % Mach number along the entire domain
                    T = zeros(1, GLenght); % Temperature along the entire domain
                    p = zeros(1, GLenght); % Pressure along the entire domain
                    f = zeros(1, GLenght);
                    Re = zeros(1, GLenght);
                    m_crit = p0 / sqrt(T0) * sqrt(gamma / R) * (1 + ((gamma - 1) / 2))^(-0.5 * (gamma + 1) / (gamma - 1)) * 2 * yT;
                    m_crit = double(m_crit);

                    for n = 1:GT
                    % Compute Mach using fzero() to solve mass flow equation
                    M(n) = fzero(@(Ma) -m_crit + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), [0.0001,1]);
                    % Compute Temperature using isentropic relations
                    T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;

                    % Compute Pressure using isentropic relations
                    p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
                    end
                    for n = GT:GE
                    % Compute Mach using fzero() to solve mass flow equation
                    M(n) = fzero(@(Ma) -m_crit + 2*radius(n) * p0/sqrt(T0) * sqrt(gamma/R) * Ma * (1+(gamma-1)/2 * Ma^2) ^ (-0.5 * (gamma+1)/(gamma-1)), [1,1000]);
                    % Compute Temperature using isentropic relations
                    T(n) = T0 * (1 + (gamma - 1)/2 * M(n)^2)^-1;
                    % Compute Pressure using isentropic relations
                    p(n) = p0 * (1 + (gamma - 1)/2 * M(n)^2)^(-gamma/(gamma-1));
                    end

                    l=int16((x_sh_l_tube*GLenght)/(NLenght+TubeLenght));

                    for n=GE:l
                    [Re(n),f(n)] = Churchill(gamma,T(n),p(n),M(n),epsilon,AR);
                    dM2_dx = (gamma * M(n)^4 * (1 + ((gamma - 1) / 2) * M(n)^2)) / (1 - M(n)^2);
                     % Update Mach number squared using discretized form
                    M_squared = M(n)^2 + (f(n)* Dx / (2*AR)) * dM2_dx;
                     % Ensure M(i+1) is real and positive
                    if M_squared < 0
                    error("Unphysical Mach number calculated. Check input values.");
                    end
                    M(n+1) = sqrt(M_squared);  % Solve for M(i+1)
                    % Compute temperature using ratio equation
                    T(n+1)= ((1 + ((gamma - 1) / 2) * M(n)^2) / (1 + ((gamma - 1) / 2) * M(n+1)^2)) * T(n);    
                    % Compute pressure using Fanno Flow properties
                    p(n+1) = p(n) * (M(n)/M(n+1)) * (T(n+1)/T(n))^(0.5);
                    end

                    M(l+1)=M_B;
                    T(l+1)=T(l)*(2 + (gamma - 1) * M(l)^2) * ((2 * gamma * M(l)^2 - (gamma - 1)) / ((gamma + 1)^2 * M(l)^2));
                    p(l+1)= p(l)*(1+(2*gamma)/(gamma+1)*(M(l)^2-1));

                    for n=l+1:GLenght-1
                    [Re(n),f(n)] = Churchill(gamma,T(n),p(n),M(n),epsilon,AR);
                    dM2_dx = (gamma * M(n)^4 * (1 + ((gamma - 1) / 2) * M(n)^2)) / (1 - M(n)^2);
                     % Update Mach number squared using discretized form
                    M_squared = M(n)^2 + (f(n) * Dx / (2*AR)) * dM2_dx;
                     % Ensure M(i+1) is real and positive
                    if M_squared < 0
                    error("Unphysical Mach number calculated. Check input values.");
                    end
                    M(n+1) = sqrt(M_squared);  % Solve for M(i+1)
                    % Compute temperature using ratio equation
                    T(n+1)= ((1 + ((gamma - 1) / 2) * M(n)^2) / (1 + ((gamma - 1) / 2) * M(n+1)^2)) * T(n);    
                    % Compute pressure using Fanno Flow properties
                    p(n+1) = p(n) * (M(n)/M(n+1)) * (T(n+1)/T(n))^(0.5);
                    end

            %Nozzle Plot 
            axes(hpNozzle)
            hold on;
            cla(hpNozzle);
            % Plot the upper half of the nozzle
            plot(G, radius, 'b', 'LineWidth', 2);
            % Plot the lower half for symmetry
            plot(G, -radius, 'b', 'LineWidth', 2);
            % Add labels and title
            plot(x_sh_l_tube*ones(1,10),linspace(-AR,AR,10),'r', 'LineWidth', 2);
            plot(G(GE)*ones(1,10), linspace(-AR,AR,10), 'k', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Radius (m)');
            if radius(end) > 3
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-1.5*radius(end), 1.5*radius(end)]); % Y-axis from -10 to 10
            else
            xlim([0, TubeLenght+2]); % X-axis from 0 to TubeLenght
            ylim([-4, 4]); % Y-axis from -10 to 10
            end
            hold off;
            hold on;
             % Mach Number Plot
            axes(hpMach);
            cla(hpMach);
            plot(G, M, 'r', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Mach Number');
            hold off;
            hold on;
            % Temperature Plot
            axes(hpTemperature);
            cla(hpTemperature);
            plot(G, T, 'b', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Temperature (K)');
            hold off;
            hold on;
            % Pressure Plot
            axes(hpPressure);
            cla(hpPressure);
            plot(G, p, 'g', 'LineWidth', 2);
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Pressure (Pa)');
            
            figure;
            hold on;
            plot(G(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Axial Position (m)');
            ylabel('Friction Factor (f)');
            title('Friction Factor Variation');
            hold off;
            figure;
            hold on;
            plot(Re(GE:GLenght-1), f(GE:GLenght-1), 'm', 'LineWidth', 2); % Plotting f(n) along the pipe
            grid on;
            xlabel('Reynolds Number (Re)');
            ylabel('Friction factor(f)');
            title('Friction factor(f) x Reynolds Number (Re)');
            hold off;

             
        else
            %Doubly chocked
            msgbox('Outside pressure too low, required a new fanno line. Model is not valid. Please introduce different inputs.');
                 return
        end
        end
               
    end
end
