function GUI_Analytical(PR,gamma,AR,hpNozzle,hpPressure,hpMach,hpTemperature)
% Group introduction, all functions used
q1dAR(PR,gamma,AR)
function q1dAR(PR, gamma, AR)
    PR=1/PR;
    AreaRatio=AR;
    NLenght=2;
    G = (0:1e-3:NLenght)';
    GLenght = length(G);
    [yT, GT, radius] = NozzleDraw(AR,G,GLenght,NLenght,0,0);
    Mach = nan(GLenght,1);
    M = Mach;
    for i = 1:GLenght
    f = @(M) (1/M^2)*((2/(gamma+1))*(1+((gamma-1)/2)*M^2))^((gamma+1)/(gamma-1))-(radius(i)/yT)^2; %Area/Astar formula
    if i < GT
       Mach(i) = Bissection(f,0.01,1,100); %Subsonic
       M(i) = Mach(i); %Subsonic
    elseif i == GT
       Mach(i) = 1; %Throat
       M(i) = 1; %Throat
    else
       Mach(i) = Bissection(f,1,10,100); %Supersonic solution limit is 10
       M(i) = Bissection(f,0.01,1,100); %Subsonic solution limit is 0.1
    end       
    end
  
    PresRatiosub = nan(GLenght,1);
    % Pressure Ratio for subsonic choked solution
    PresRatiosuper = nan(GLenght,1);
    % Pressure Ratio for supersonic solution (perfectly expanded)
    for i = 1:GLenght
        PresRatiosub(i) = (1+((gamma-1)/2)*M(i)^2)^(-gamma/(gamma-1));
        PresRatiosuper(i) = (1+((gamma-1)/2)*Mach(i)^2)^(-gamma/(gamma-1)); 
    end
        Subsonic = PresRatiosub(end) %Pressure ratio at the exit
        Supersonic = PresRatiosuper(end) %Pressure ratio at the exit
        1/Subsonic
        1/Supersonic
%Here will be the plot functions
    function Plot_Nozzle(varargin)
    type = varargin{1,1}; %type of shocks present    
        if type == 0
        hNozzle = plot(hpNozzle,G,radius, ... 
            G,-(radius), ...  
            G(GT)*ones(1,10),linspace(-yT,yT,10));
        for NozzleIndex = 1:length(hNozzle)-1
            set(hNozzle(NozzleIndex),'LineWidth',3,'Color','b');
        end
            set(hNozzle(length(hNozzle)),'LineWidth',1.2,'Color','k');

        elseif type == 1 %Normal Shock
            GShock = varargin{1,2};
            hNozzle = plot(hpNozzle,G,radius, ... 
            G,-(radius), ... 
            G(GT)*ones(1,10),linspace(-yT,yT,10),GShock*ones(1,10),linspace(-(radius(ShockL)),radius(ShockL),10));
        for NozzleIndex = 1:length(hNozzle)-2
                set(hNozzle(NozzleIndex),'LineWidth',3,'Color','b');
        end
            set(hNozzle(length(hNozzle)-1),'LineWidth',1.2,'Color','k');
            set(hNozzle(length(hNozzle)),'LineWidth',5,'Color','r');

        elseif type == 2 %Oblique Shock
            %slope = varargin{1,2};
            Beta = varargin{1,2};
            NewGrid = (G(end):1e-3:3)';
            NewGridRadius = nan(length(NewGrid),1);
            NewGridRadius(1) = radius(end);
            %Inclination = tan(slope);
            Inclination = tan(Beta);
            NGLenght = length(NewGrid);

            for p = 2:NGLenght
                if NewGrid(p)<=2.5
                NewGridRadius(p) = NewGridRadius(p-1) - (NewGrid(p)-NewGrid(p-1))*Inclination;
                else 
                NewGridRadius(p) = nan;
                end

            end
            hNozzle = plot(hpNozzle, G, radius, G, -radius, G(GT)*ones(1,100), linspace(-yT, yT, 100), NewGrid, NewGridRadius, NewGrid, -(NewGridRadius));
            
            for j = 1:length(hNozzle)
                if j <= 2
                    set(hNozzle(j),'LineWidth',3,'Color','b');
                elseif j == 3
                    set(hNozzle(j),'LineWidth',1.2,'Color','k');
                else
                    set(hNozzle(j),'LineWidth',2,'Color','r');
                end
            end

        elseif type == 3 %Expansion Fans
            NewGrid = (G(end):1e-3:3)';
            NewGridRadius = nan(length(NewGrid),1);
            NewGridRadius(1) = radius(end);
            NewGridRadius2 = NewGridRadius;
            NewGridRadius3 = NewGridRadius;

            Angles = varargin{1,2};
            Inclination = tan(Angles(2));
            Inclination2 = tan(Angles(1));
            Inclination3 = tan(Angles(3));

            NGLenght = length(NewGrid);

            for p = 2:NGLenght
                if NewGrid(p)<=3
                NewGridRadius(p) = NewGridRadius(p-1) - (NewGrid(p)-NewGrid(p-1))*Inclination;  
                NewGridRadius2(p) = NewGridRadius2(p-1) + (NewGrid(p)-NewGrid(p-1))*Inclination2;
                NewGridRadius3(p) = NewGridRadius3(p-1) - (NewGrid(p)-NewGrid(p-1))*Inclination3;
                else 
                NewGridRadius(p) = nan;
                NewGridRadius2(p) = nan;
                NewGridRadius3(p) = nan;
                end
            end
            inicio=NewGridRadius(2);
            fim=NewGridRadius(NGLenght);
            slope1=atan(inicio-fim)*180/pi
            inicio2=NewGridRadius3(2);
            fim2=NewGridRadius3(NGLenght);
            slope2=atan(inicio2-fim2)*180/pi
            
            hNozzle = plot(hpNozzle,G,radius, ... 
            G, -(radius), ...
            G(GT)*ones(1,10),linspace(-yT,yT,10), ...
            NewGrid,NewGridRadius, ...
            NewGrid,-NewGridRadius, ...
            NewGrid,NewGridRadius2, ...
            NewGrid,-NewGridRadius2, ...
            NewGrid,NewGridRadius3, ...
            NewGrid,-NewGridRadius3);          
        for j = 1:length(hNozzle)
            if j <= 2
                set(hNozzle(j),'LineWidth',3,'Color','b');
            elseif j == 3
                set(hNozzle(j),'LineWidth',1.2,'Color','k');
            elseif j == 6 || j == 7
                set(hNozzle(j),'LineWidth',1.2,'Color','k');
            else
                set(hNozzle(j),'LineWidth',0.5,'Color','g')
            end

        end
        else % No shocks or expansion fans 
        end
        T = sprintf('CONVERGING-DIVERGING NOZZLE for PR = %4.3f',1/PR);
        title(hpNozzle,T);
        if radius(end) > 3
            set(hpNozzle,'XLim',[G(1) 1.5*G(end)],'YLim',[-1.5*radius(end) 1.5*radius(end)]);
        else
            set(hpNozzle,'XLim',[G(1) 1.5*G(end)],'YLim',[-4 4]);
        end
        hThroat = uicontrol('Style','text','String','Thin black line is the throat','Position',[10,100,200,20]);
        set(hThroat,'FontSize',10);
        hShock  = uicontrol('Style','text','String','Red lines are shocks','Position',[10,75,200,20]);
        set(hShock,'FontSize',10);
        hEF = uicontrol('Style','text','String','Green lines are expansion fans','Position',[10,50,200,20]);
        set(hEF,'FontSize',10);
    end

    function Plot_Pressure(Entries,PR,Opt)
    %The parameter can take validity values, pressure values or Mach values
    % Parameter contains the Mach number to be plotted 
    % (if the solution contains shocks or expansion fans, or 
    % is subsonic everywhere), OR it contains the scalar value
    % 0 or 1, indicating the subsonic or supersonic solution.
    if length(Entries) ~= 1 && Opt == 1
        %Normal Shock, %Entries=Pressure Ratio
        Gnew = [G;(G(end)+0.001:0.001:3)']; %This is to draw after the exit
        Glengthnew=length(Gnew);
        Entries = [Entries;Entries(end)*ones(Glengthnew-GLenght,1)];
        hPressure = plot(hpPressure,Gnew,Entries,G(end)*ones(10,1),(linspace(0,1.1,10))');
        set(hPressure(1),'LineWidth',3,'Color','b');
        set(hPressure(2),'LineWidth',2,'Color','k');
    elseif length(Entries) ~= 1 && Opt == 2
        %Oblique Shock %Entries=Pressure Ratio supersonic
        Gnew = [G;(G(end)+0.001:0.001:3)']; %This is to draw after the exit
        Glengthnew=length(Gnew);
        Entries = [Entries;Entries(end)*ones(100,1);PR*ones(Glengthnew-GLenght-100,1)];
        hPressure = plot(hpPressure,Gnew,Entries,G(end)*ones(10,1),(linspace(0,1.1,10))');
        set(hPressure(1),'LineWidth',3,'Color','b');
        set(hPressure(2),'LineWidth',2,'Color','k');
    elseif isscalar(Entries) && Opt == 2 %Parameter lenght == 1
            %Expansion Fan
            Gnew = [G;(G(end)+0.001:0.001:3)'];
            Gnewlength = length(Gnew);
            LSpline = 400;
            S1 = [G(end-10:end);Gnew(end-LSpline:end)];
            S2 = [PresRatiosuper(end-10:end);PR*ones(LSpline+1,1)];
            A = nan(Gnewlength-GLenght,1);
            for z = 1:length(A)
                A(z) = spline(S1,S2,Gnew(GLenght+z));
            end
            PresRatio2new = [PresRatiosuper;A];
            hPressure = plot(hpPressure,Gnew,PresRatio2new,G(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');
    elseif Opt == -1
            % No flow
            Gnew = [G;(G(end)+0.001:0.001:3)']; %This is to draw after the exit
            Glengthnew=length(Gnew);
            PressureZero = ones(Glengthnew,1);
            hPressure = plot(hpPressure,Gnew,PressureZero,G(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');
     elseif length(Entries) ~= 1 && Opt == 0
            % Subsonic but not choked
            %Entries= MachNumberDistribution
            PresRatioSS = nan(GLenght,1);
            for s = 1:GLenght
                PresRatioSS(s) = (1+((gamma-1)/2)*Entries(s)^2)^(-gamma/(gamma-1));
            end
            Gnew = [G;(G(end)+0.001:0.001:3)']; %This is to draw after the exit
            Glengthnew=length(Gnew);
            PresRatioSS = [PresRatioSS;PresRatioSS(end)*ones(Glengthnew-GLenght,1)];
            hPressure = plot(hpPressure,Gnew,PresRatioSS,G(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');
      elseif Entries == true && Opt == 0
            % Plot supersonic solution
            Gnew = [G;(G(end)+0.001:0.001:3)'];
            Gnewlength = length(Gnew);
            PresRatio2new = [PresRatiosuper;PresRatiosuper(end)*ones(Gnewlength-GLenght,1)];
            hPressure = plot(hpPressure,Gnew,PresRatio2new,G(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k');
       elseif Entries == false && Opt == 0
            % Plot subsonic solution 
            Gnew = [G;(G(end)+0.001:0.001:3)'];
            Gnewlenght = length(Gnew);
            PresRationew = [PresRatiosub;PresRatiosub(end)*ones(Gnewlenght-GLenght,1)];
            hPressure = plot(hpPressure,Gnew,PresRationew,G(end)*ones(10,1),(linspace(0,1.1,10))');
            set(hPressure(1),'LineWidth',3,'Color','b');
            set(hPressure(2),'LineWidth',2,'Color','k'); 
       end
        title(hpPressure,'Pressure ratio plot');
        set(hpPressure,'XLim',[G(1) 1.5*G(end)],'YLim',[0 1.1]);
        grid(hpPressure);
    end
    function varargout = Plot_Mach(Entries,PR,Opt)
        % Entries contains the Mach number to be plotted 
        % (if the solution contains shocks or expansion fans, or 
        % is subsonic everywhere), OR it contains the scalar value
        % 0 or 1, indicating the subsonic or supersonic solution.
        % Opt = 1: Shock
        % Opt = 2: Expansion fan
        % Opt = 0: Neither
        if length(Entries) ~= 1 && Opt == 1
            %Normal Shock, feed in the shock location and M1    
            Gnew = [G;(G(end)+0.001:0.001:3)'];
            Gnewlenght = length(Gnew);
            Entries = [Entries;Entries(end)*ones(Gnewlenght-GLenght,1)];
            hMach = plot(hpMach,Gnew,Entries,G(end)*ones(10,1),(linspace(0,1.5*max(Entries),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif length(Entries) ~= 1 && Opt == 2
            % Oblique shock
            % Feed in supersonic solution,shock away from the nozzle exit
            Gnew = [G;(G(end)+0.001:0.001:3)'];
            Gnewlenght = length(Gnew);
            Entries = [Entries;Entries(end)*ones(100,1);PR*ones(Gnewlenght-GLenght-100,1)];
            hMach = plot(hpMach,Gnew,Entries,G(end)*ones(10,1),(linspace(0,1.5*max(Entries),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
         elseif length(Entries) == 1 && Opt == 2
             % Expansion fan
             % Expansion fans are gradual, so the solution should be 
             % a smooth curve.
             % Solve for downstream Mach number (M2)
             F = @(M2) ((1+((gamma-1)/2)*Mach(end)^2)/(1+((gamma-1)/2)*M2^2))^(gamma/(gamma-1)) - (PR/PresRatiosuper(end)); %This final is p01/p02
             M2 = abs(fzero(F,1.1*Mach(end)));
            varargout{1,1} = M2;
            xnew = [G;(G(end)+0.001:0.001:3)'];
            Z = length(xnew);
            LSpline = 400;
            S1 = [G(end-10:end);xnew(end-LSpline:end)];
            S2 = [Mach(end-10:end);M2*ones(LSpline+1,1)];
            A = nan(Z-GLenght,1);
            for z = 1:length(A)
                A(z) = spline(S1,S2,xnew(GLenght+z));
            end
            Machnew = [Mach;A];
            hMach = plot(hpMach,xnew,Machnew,G(end)*ones(10,1),(linspace(0,1.5*max(Machnew),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif  length(Entries) ~= 1 && Opt == 0
            %Subsonic not Choked
            NewG = [G;(G(end)+0.001:0.001:3)']; %Extends the Grid
            NewGLength = length(NewG);
            Entries = [Entries;Entries(end)*ones(NewGLength-GLenght,1)] ;%Mach on the extendend grid
            hMach = plot(hpMach,NewG,Entries,G(end)*ones(10,1),(linspace(0,1.5*max(Entries),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif Entries == true
            %Supersonic Solution
            NewG = [G;(G(end)+0.001:0.001:3)']; %Extends the Grid
            NewGLength = length(NewG);
            Machnew = [Mach;Mach(end)*ones(NewGLength-GLenght,1)];
            hMach = plot(hpMach,NewG,Machnew,G(end)*ones(10,1),(linspace(0,1.1*max(Mach),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        elseif Opt == -1
            % No flow
            xnew = [G;(G(end)+0.001:0.001:3)'];
            Z = length(xnew);
            MachZero = zeros(Z,1);
            hMach = plot(hpMach,xnew,MachZero,G(end)*ones(10,1),(linspace(0,0.01,10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');
        else
            % Plot subsonic choked solution 
            xnew = [G;(G(end)+0.001:0.001:3)'];
            Z = length(xnew);
            Mnew = [M;M(end)*ones(Z-GLenght,1)];
            hMach = plot(hpMach,xnew,Mnew,G(end)*ones(10,1),(linspace(min(M),1.1*max(M),10))');
            set(hMach(1),'LineWidth',3,'Color','b');
            set(hMach(2),'LineWidth',2,'Color','k');   
        end
        title(hpMach,'Mach number plot');
        grid(hpMach);
    end
    
    function Plot_Temperature(Entries, PR, Opt)
    % Entries contains the Mach number or a scalar value indicating the flow type
    % Opt = 1: Shock
    % Opt = 2: Expansion fan
    % Opt = 0: Neither
        if length(Entries) ~= 1 && Opt == 1
        % Normal Shock
        %Feed the Mach
            Gnew = [G; (G(end)+0.001:0.001:3)'];
            GnewLength = length(Gnew);
            Entries = [Entries; Entries(end)*ones(GnewLength-GLenght,1)];
            TemperatureRatio = (1 + ((gamma-1)/2) * Entries.^2).^-1; % Temperature Ratio T/T0
            hTemp = plot(hpTemperature, Gnew, TemperatureRatio, G(end)*ones(10,1), (linspace(0,1.5,10))');
            set(hTemp(1),'LineWidth',3,'Color','b');
            set(hTemp(2),'LineWidth',2,'Color','k');
        elseif length(Entries) == 1 && Opt == 3
            Gnew = [G;(G(end)+0.001:0.001:3)'];
            Gnewlenght = length(Gnew);
            Entries = [Entries;Entries(end)*ones(100,1);PR*ones(Gnewlenght-GLenght-100,1)];
            TemperatureRatio = (1 + ((gamma-1)/2) * Entries.^2).^-1; % Temperature Ratio T/T0
            hTemp = plot(hpTemperature, Gnew, TemperatureRatio, G(end)*ones(10,1), (linspace(0,1.5,10))');
            set(hTemp(1),'LineWidth',3,'Color','b');
            set(hTemp(2),'LineWidth',2,'Color','k');
        elseif length(Entries) == 1 && Opt == 2
            % Expansion fan
            F = @(M2) ((1+((gamma-1)/2)*Mach(end)^2)/(1+((gamma-1)/2)*M2^2))^(gamma/(gamma-1)) - (PR/PresRatiosuper(end)); %This final is p01/p02
            M2 = abs(fzero(F,1.1*Mach(end)));
            xnew = [G;(G(end)+0.001:0.001:3)'];
            Z = length(xnew);
            LSpline = 400;
            S1 = [G(end-10:end);xnew(end-LSpline:end)];
            S2 = [Mach(end-10:end);M2*ones(LSpline+1,1)];
            A = nan(Z-GLenght,1);
            for z = 1:length(A)
                A(z) = spline(S1,S2,xnew(GLenght+z));
            end
            Machnew = [Mach;A];
            TemperatureRatioNew = (1+((gamma-1)/2)*Machnew.^2).^-1;
            hTemp = plot(hpTemperature, xnew, TemperatureRatioNew, G(end)*ones(10,1), (linspace(0,1.5,10))');
            set(hTemp(1),'LineWidth',3,'Color','b');
            set(hTemp(2),'LineWidth',2,'Color','k');
        elseif isscalar(Entries) && Opt == -1
            % No flow condition
            Gnew = [G; (G(end)+0.001:0.001:3)'];
            GnewLength = length(Gnew);
            TempZero = ones(GnewLength,1);
            hTemp = plot(hpTemperature, Gnew, TempZero, G(end)*ones(10,1), (linspace(0,1.5,10))');
            set(hTemp(1),'LineWidth',3,'Color','b');
            set(hTemp(2),'LineWidth',2,'Color','k');
        else
        % Regular temperature distribution
            Gnew = [G; (G(end)+0.001:0.001:3)'];
            GnewLength = length(Gnew);
            TemperatureRatioNew = [(1+((gamma-1)/2)*Entries.^2).^-1; (1+((gamma-1)/2)*Entries(end)^2).^-1*ones(GnewLength-GLenght,1)];
            hTemp = plot(hpTemperature, Gnew, TemperatureRatioNew, G(end)*ones(10,1), (linspace(0,1.5,10))');
            set(hTemp(1),'LineWidth',3,'Color','b');
            set(hTemp(2),'LineWidth',2,'Color','k');
        end
    title(hpTemperature, 'Temperature Ratio (T/T_0) Across Nozzle');
    grid(hpTemperature);
    end
   

%Start here the Analysis
p0 = 200; %Some value (Analyse)
zero = 1e-3;

if PR>Subsonic && PR~=1 
    message = uicontrol('Style','text','String','Subsonic-Subsonic Flow','Position',[250,605,120,40]);
    set(message,'FontSize',10);
    % In this case A* doesn´t correspond to the throat so we need to find
    % the area  at which the flow would reach Mach 1
    Me = sqrt((2/(gamma-1))*((PR)^(-(gamma-1)/gamma)-1));
    AreaRatio2 = sqrt((1/(Me^2))*((2/(gamma+1))*(1+((gamma-1)/2)*Me^2))^((gamma+1)/(gamma-1)));
    % AreaRatio2 is A/A**
    A_star2 = yT*AreaRatio/AreaRatio2;  %A** (hipotetic area at which the flow would reach Mach 1)
    for i = 1:GLenght
        f = @(M) (1/M^2)*((2/(gamma+1))*(1+((gamma-1)/2)*M^2))^((gamma+1)/(gamma-1))-(radius(i)/A_star2)^2;
        MachSubSonic(i) = Bissection(f,0.001,1,100);
    end
    Plot_Pressure(MachSubSonic,PR,0);
    Plot_Nozzle(0);
    Plot_Mach(MachSubSonic',PR,0);
    Plot_Temperature(MachSubSonic',PR,0);
    %Plot;
elseif PR == 1 
    % No flow
    message = uicontrol('Style','text','String','No flow.','Position',[250,605,120,40]);
    set(message,'FontSize',15);
    Plot_Pressure(2,PR,-1);
    Plot_Nozzle(0);
    Plot_Mach(2,PR,-1);
    Plot_Temperature(2,PR,-1);
    %Plot;
elseif abs(PR-Subsonic) < zero 
    % Subsonic solution 
    message = uicontrol('Style','text','String','Subsonic and choked','Position',[250,605,120,40]);
    set(message,'FontSize',11);
    Plot_Pressure(0,PR,0);
    Plot_Nozzle(0);
    Plot_Mach(0,PR,0);
    Plot_Temperature(M,PR,0);
    %Plot;
elseif PR < Subsonic && PR > Supersonic && AR > 1
    %This solution is where we will find shocks, either normal
    %shocks or oblique shocks
    %All equations from the book "Modern Compressible flow"
    %Lets first find the shock location 
    Mexit = sqrt((-1/(gamma-1))+sqrt((1/((gamma-1)^2))+(2/(gamma-1))*(2/(gamma+1))^((gamma+1)/(gamma-1))*(1/PR)^2*(1/AreaRatio)^2)); %Eq 5.28
    pexit = PR*p0;
    pstagexit = pexit*(1+((gamma-1)/2)*Mexit^2)^(gamma/(gamma-1)); %Second step of the "Direct Method" po2
    Ratio = (pstagexit/p0); %p02/p01
    %Formulario Aero II
     S = @(M1) (( ( (gamma+1) * M1.^2 ) ./ ( (gamma-1) * M1.^2 + 2 ) ).^(gamma/(gamma-1)).* ( (gamma+1) ./ ( 2*gamma*M1.^2 - (gamma-1) ) ).^(1/(gamma-1))) - Ratio;
     MachCritical = Bissection(S,0.9,5,100);
     ShockL = nan;
       for i = 1:GLenght
           if abs(Mach(i)-MachCritical) < zero
                 ShockL = i;
           break
           end
        end
            if isnan(ShockL) == 1 && MachCritical > Mach(end) || Mexit > 1 
               %Oblique case shock
               message = uicontrol('Style','text','String','Oblique shock outside.','Position',[250,605,120,40]);
               set(message,'FontSize',11);
               NewPR = @(B) 1+(2*gamma/(gamma+1))*(Mach(end)^2*sin(B)^2-1) - PR*(1/PresRatiosuper(end)); %Eq 4.9 %This last part is p2/p1;
               beta = Bissection(NewPR,0,pi/2,100) %Normal shocks are oblique shocks when beta=pi/2
               theta = atan(2*cot(beta)*((Mach(end)^2*sin(beta)^2-1)/(Mach(end)^2*(gamma+cos(2*beta))+2))) %Eq 4.17 Turned Flow
               %slope = beta + theta;
               %if slope >= pi/2
               %    slope = pi/2;
               %end
               MeOS = (1/sin(beta-theta))*sqrt((1+((gamma-1)/2)*Mach(end)^2*sin(beta)^2)/(gamma*Mach(end)^2*sin(beta)^2-((gamma-1)/2))); %Eq 4.12 + Eq 4.10
               Plot_Pressure(PresRatiosuper,PR,2);
               Plot_Mach(Mach,MeOS,2);
               %Plot_Nozzle(2, slope);
               Plot_Nozzle(2, beta);
               Plot_Temperature(Mach,MeOS,3);
               %Plot
            else %Normal shock case
               message = uicontrol('Style','text','String','Normal shock inside.','Position',[250,605,120,40]);
               set(message,'FontSize',11);
               SLocation=G(ShockL);
               Plot_Nozzle(1, SLocation, ShockL);
               Mach_NS = nan(GLenght,1);
               Mach_NS(1:ShockL) = Mach(1:ShockL); %Isentropic Solution until ShockL
               %Find out now the solution after the shock based on the Me
               %(A**)
               AR_new = sqrt((1/(Mexit^2))*((2/(gamma+1))*(1+((gamma-1)/2)*Mexit^2))^((gamma+1)/(gamma-1)));
               A_2star= yT*AreaRatio/AR_new;
               %Now for the pressure
               for u = ShockL+1:GLenght
                    f = @(M) (1/M^2)*((2/(gamma+1))*(1+((gamma-1)/2)*M^2))^((gamma+1)/(gamma-1))-(radius(u)/A_2star)^2;
                    Mach_NS(u) = Bissection(f,0.001,1,100);
               end
               PresRatioNS = nan(GLenght,1);
               for w = 1:GLenght
                    if w > ShockL
                        PresRatioNS(w) = (1+((gamma-1)/2)*Mach_NS(w)^2)^(-gamma/(gamma-1))*Ratio; %p/p02
                    else
                        PresRatioNS(w) = (1+((gamma-1)/2)*Mach_NS(w)^2)^(-gamma/(gamma-1)); %p/p01
                    end
                end
                Plot_Mach(Mach_NS,PR,1);
                Plot_Pressure(PresRatioNS,PR,1);    
                Plot_Temperature(Mach_NS,PR,1);
            end
elseif abs(PR-Supersonic) < zero && AR > 1
    % Perfectly expanded
    message = uicontrol('Style','text','String','Perfectly expanded jet.','Position',[250,605,120,40]);
    set(message,'FontSize',11);
    Plot_Pressure(1,PR,0);
    Plot_Nozzle(0);
    Plot_Mach(1,PR,0);
    Plot_Temperature(Mach,PR,0);
    %Plot
else
    % Underexpanded jet
    % Prandtl-Meyer expansion fan
    message = uicontrol('Style','text','String','Underexpanded jet.','Position',[250,605,120,40]);
    set(message,'FontSize',11);
    Mexit = Plot_Mach(true,PR,2);
    Mach(end)
    PMteta = (Prandtl_Meyer(Mexit)-Prandtl_Meyer(Mach(end)))
    miu1=asin(1/Mach(end))*180/pi
    miu2=asin(1/Mexit)*180/pi
    Angles=[PMteta, asin(1/Mach(end)), asin(1/Mexit)];
    Plot_Pressure(true,PR,2);
    Plot_Nozzle(3,Angles);
    Plot_Temperature(true,PR,2);
    %END OF IF
end    
  %END OF Q1DAR
end
%END OF ALL
end