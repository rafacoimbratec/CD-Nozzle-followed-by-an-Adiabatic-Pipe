function [yT, GT, radius,G,GE] = NozzleDraw(AR, G, GLenght, NLenght,Parameter,TubeLength)
    %Nozzle Making
if Parameter==0
        zero = 1e-3; %Aproximate Zero
    for i = 1:GLenght %Iterates through every grid point
        if abs(G(i)-(NLenght/8))<zero
            GR = i; %Point in the Grid at which the Reservoir ends
        elseif abs(G(i)-(NLenght/2))<zero
            GT = i; %Point in the Grid in which the Throat is
            break
        end
    end
    yR = 10; %Reservoir Radius
    yT = yR/3; %Throat Radius
    % Construct converging cosine wave
    ampC = (yR-yT)/2; %Amplitude wave
    wC = (pi)/(G(GT)-G(GR)); %Angular frequency (half-period)
    psC = -wC*G(GR); %Phase shift (GR max value)
    CW = @(G) ampC*cos(wC*G+psC)+yT+ampC;
    % Construct diverging cosine wave
    EA = AR*yT;
    ampD = (EA-yT)/2;
    wD = pi/(G(end)-G(GT));
    psD = -wD*G(end);
    DW = @(G) ampD*cos(wD*G+psD)+yT+ampD;  
    radius = yR*ones(GR,1);
    for i = 1:GLenght
        if i > GR && i < GT
            radius(i) = CW(G(i));
        elseif i == GT
            radius(i) = yT;
        elseif i > GT
            radius(i) = DW(G(i));
        end
    end
    elseif Parameter == 1
    zero = 1e-3; % Approximate Zero
    % Find Grid Points for Reservoir, Throat, and Exit
    for i = 1:GLenght  % Iterates through every grid point
        if abs(G(i) - (NLenght/8)) < zero
            GR = i; % Grid point where Reservoir ends
        elseif abs(G(i) - (NLenght/2)) < zero
            GT = i; % Grid point where the Throat is
        end
    end
    
    GE = GLenght;  % Define Grid Exit (Assume last point of original grid)
    
    % Define Radius at Key Points
    yR = 3; % Reservoir Radius
    yT = yR / 3; % Throat Radius
    EA = AR * yT; % Exit Area (AR is Area Ratio)

    % Construct Converging Section using a Cosine Function
    ampC = (yR - yT) / 2;
    wC = pi / (G(GT) - G(GR));
    psC = -wC * G(GR);
    CW = @(x) ampC * cos(wC * x + psC) + yT + ampC;

    % Construct Diverging Section using a Cosine Function
    ampD = (EA - yT) / 2;
    wD = pi / (G(GE) - G(GT));
    psD = -wD * G(GE);
    DW = @(x) ampD * cos(wD * x + psD) + yT + ampD;

    % Create Extended Pipe Section (Straight Tube)
    Gtube = (G(GE) + 1e-3 : 1e-3 : TubeLength + NLenght)';  % New tube grid points
    G = [G; Gtube];  % Append tube grid to the existing grid
    GLenght = length(G);  % Update total grid length

    % Initialize Radius Array
    radius = yR * ones(GLenght, 1); % Set reservoir radius by default

    % Assign Radii Based on Region
    for i = 1:GLenght
        if i > GR && i < GT
            radius(i) = CW(G(i)); % Apply Converging Section
        elseif i == GT
            radius(i) = yT; % Throat radius
        elseif i > GT && i <= GE
            radius(i) = DW(G(i)); % Apply Diverging Section
        elseif i > GE
            radius(i) = EA; % Pipe section (constant radius)
        end
    end
end
end    
