   function GUI_Main
    close all force;
    global PR gamma AR ModelChoice Dpipe epsilon L P1 T1 V1;
    PR = nan; gamma = nan; AR = nan; ModelChoice = 1;
    Dpipe = nan; epsilon = nan; L = nan; P1 = nan; T1 = nan; V1 = nan;
    
    GUI = figure('Visible', 'off', 'Position', [200, 100, 1200, 600], 'Color', [0.95 0.95 0.95]);
    
    % Panel for Model Selection
    modelPanel = uipanel('Title', 'Select Model', 'FontSize', 12, 'Position', [0.05 0.86 0.3 0.1]);

    uicontrol(modelPanel, 'Style', 'popupmenu', 'String', {'Analytical', 'Friction in a Pipe'}, ...
        'Position', [10, 10, 150, 30], 'Callback', @SelectModel);
    
    % Enlarged Input Panel
    inputPanel = uipanel('Title', 'Input Parameters', 'FontSize', 10, 'Position', [0.05 0.4 0.3 0.455]);

    function SelectModel(source, ~)
        ModelChoice = get(source, 'Value');
        delete(allchild(inputPanel)); % Clear previous inputs
        
        switch ModelChoice
            case 1 % Analytical Model
                uicontrol(inputPanel, 'Style', 'text', 'String', 'Pressure Ratio (p0/pe):', 'Position', [10, 110, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 110, 60, 25], 'Callback', @Activate_Pressure_Ratio);
                
                uicontrol(inputPanel, 'Style', 'text', 'String', 'Area Ratio:', 'Position', [10, 80, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 80, 60, 25], 'Callback', @Activate_AR);
                
                uicontrol(inputPanel, 'Style', 'text', 'String', 'Gamma:', 'Position', [10, 50, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 50, 60, 25], 'Callback', @Activate_Gamma);
                
            case 2 % Fanno Flow Model

                uicontrol(inputPanel, 'Style', 'text', 'String', 'Pipe Diameter (m):', 'Position', [10, 190, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 200, 60, 25], 'Callback', @Activate_Dpipe);
                
                uicontrol(inputPanel, 'Style', 'text', 'String', 'Roughness (m):', 'Position', [10, 160, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 170, 60, 25], 'Callback', @Activate_Epsilon);
                
                uicontrol(inputPanel, 'Style', 'text', 'String', 'Pipe Length (m):', 'Position', [10, 130, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 140, 60, 25], 'Callback', @Activate_L);
               
                uicontrol(inputPanel, 'Style', 'text', 'String', 'Inlet Pressure (Pa):', 'Position', [10, 100, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 110, 60, 25], 'Callback', @Activate_P1);
                
                uicontrol(inputPanel, 'Style', 'text', 'String', 'Inlet Temperature (K):', 'Position', [10, 70, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 80, 60, 25], 'Callback', @Activate_T1);
                
                uicontrol(inputPanel, 'Style', 'text', 'String', 'Inlet Velocity (m/s):', 'Position', [10, 40, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 50, 60, 25], 'Callback', @Activate_V1);

                uicontrol(inputPanel, 'Style', 'text', 'String', 'Gamma:', 'Position', [10, 20, 150, 25]);
                uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 20, 60, 25], 'Callback', @Activate_Gamma);
        end
    end

    % Buttons Panel
    controlPanel = uipanel('Title', 'Controls', 'FontSize', 12, 'Position', [0.05 0.2 0.3 0.2]);
    uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Plot', 'Position', [20, 5, 180, 95], 'Callback', @Plot);
    uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Restart', 'Position', [220, 55, 130, 45], 'Callback', @Restart);
    uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Exit', 'Position', [220, 5, 130, 45], 'Callback', @Button);
    % Criar um botão na GUI para abrir o arquivo de ajuda
    uicontrol('Style', 'pushbutton', ...
                            'String', 'Help', ...
                            'Position', [10, 10, 80, 30], ...
                            'FontSize', 10, ...
                            'Callback', @openHelpFile);

    % Função para abrir o arquivo de ajuda
    function openHelpFile(~, ~)
        helpFile = 'help.txt'; % Nome do arquivo de ajuda

        % Verifica se o arquivo existe antes de tentar abrir
        if exist(helpFile, 'file')
            winopen(helpFile); % Abre o arquivo no Windows
        else
            errordlg('O arquivo de ajuda não foi encontrado.', 'Erro');
        end
    end
    
    % Plot Panels
    plotPanel = uipanel('Title', 'Plots', 'FontSize', 12, 'Position', [0.4 0.1 0.55 0.9]);
    hpNozzle = axes('Parent', plotPanel, 'Units', 'normalized', 'Position', [0.1, 0.7, 0.8, 0.25]);
    hpPressure = axes('Parent', plotPanel, 'Units', 'normalized', 'Position', [0.1, 0.48, 0.8, 0.15]);
    hpMach = axes('Parent', plotPanel, 'Units', 'normalized', 'Position', [0.1, 0.25, 0.8, 0.15]);
    hpTemperature = axes('Parent', plotPanel, 'Units', 'normalized', 'Position', [0.1, 0.03, 0.8, 0.15]);

    % Make GUI visible
    set(GUI, 'Visible', 'on');

    function Activate_Pressure_Ratio(source, ~)
        PR = str2double(get(source, 'String'));
    end

    function Activate_Gamma(source, ~)
        gamma = str2double(get(source, 'String'));
    end

    function Activate_AR(source, ~)
        AR = str2double(get(source, 'String'));
    end

    function Activate_Dpipe(source, ~)
        Dpipe = str2double(get(source, 'String'));
    end

    function Activate_Epsilon(source, ~)
        epsilon = str2double(get(source, 'String'));
    end

    function Activate_L(source, ~)
        L = str2double(get(source, 'String'));
    end

    function Activate_P1(source, ~)
        P1 = str2double(get(source, 'String'));
    end

    function Activate_T1(source, ~)
        T1 = str2double(get(source, 'String'));
    end

    function Activate_V1(source, ~)
        V1 = str2double(get(source, 'String'));
    end

    function Button(~, ~)
        fprintf('Closing.\n');
        close all force;
    end

    function Restart(~, ~)
        clc;
        GUI_Main;
    end

    function Plot(~, ~)
        switch ModelChoice
            case 1
                fprintf('Running Analytical Model...\n');
                GUI_Analytical(PR, gamma, AR, hpNozzle, hpPressure, hpMach, hpTemperature);
            case 2
                fprintf('Running Fanno Flow Model in a Pipe...\n');
                GUI_Fanno_Flow_pipe(Dpipe, epsilon, L, P1, T1, V1, gamma, hpNozzle, hpPressure, hpMach, hpTemperature);
        end
    end
end
