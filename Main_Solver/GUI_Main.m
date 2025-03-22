function GUI_Main
    close all force;
    global AR p0 T0 pe TubeLenght epsilon gamma;
    AR = nan; p0 = nan; T0 = nan; pe = nan; TubeLenght = nan; epsilon = nan; gamma = nan;

    % Create GUI Figure
    GUI = figure('Visible', 'off', 'Position', [200, 100, 1200, 600], 'Color', [0.95 0.95 0.95]);

     % Panel for Model Selection
    modelPanel = uipanel('Title', '', 'FontSize', 12, 'Position', [0.05 0.86 0.3 0.1]);
    uicontrol(modelPanel, 'Style', 'text','FontSize', 13, 'String', 'Solver: Converging-Diverging Nozzle Followed by Adiabatic Pipe with Friction', 'Position', [10, 10, 350, 50]);
   
    % Enlarged Input Panel
    inputPanel = uipanel('Title', 'Input Parameters', 'FontSize', 10, 'Position', [0.05 0.4 0.3 0.455]);

    % Input Fields
    uicontrol(inputPanel, 'Style', 'text', 'String', 'Area Ratio (Ae/At):', 'Position', [10, 220, 150, 25]);
    uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 220, 60, 25], 'Callback', @Activate_AR);
    
    uicontrol(inputPanel, 'Style', 'text', 'String', 'Total Pressure p0 (Pa):', 'Position', [10, 190, 150, 25]);
    uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 190, 60, 25], 'Callback', @Activate_p0);
    
    uicontrol(inputPanel, 'Style', 'text', 'String', 'Total Temperature T0 (K):', 'Position', [10, 160, 150, 25]);
    uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 160, 60, 25], 'Callback', @Activate_T0);

    uicontrol(inputPanel, 'Style', 'text', 'String', 'Exit Pressure pe (Pa):', 'Position', [10, 130, 150, 25]);
    uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 130, 60, 25], 'Callback', @Activate_pe);

    uicontrol(inputPanel, 'Style', 'text', 'String', 'Tube Length (m):', 'Position', [10, 100, 150, 25]);
    uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 100, 60, 25], 'Callback', @Activate_L);

    uicontrol(inputPanel, 'Style', 'text', 'String', 'Roughness:', 'Position', [10, 70, 150, 25]);
    uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 70, 60, 25], 'Callback', @Activate_Epsilon);

    uicontrol(inputPanel, 'Style', 'text', 'String', 'Gamma:', 'Position', [10, 40, 150, 25]);
    uicontrol(inputPanel, 'Style', 'edit', 'Position', [170, 40, 60, 25], 'Callback', @Activate_Gamma);

    % Buttons Panel
    controlPanel = uipanel('Title', 'Controls', 'FontSize', 12, 'Position', [0.05 0.2 0.3 0.2]);
    uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Plot', 'Position', [20, 5, 180, 95], 'Callback', @Plot);
    uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Restart', 'Position', [220, 55, 130, 45], 'Callback', @Restart);
    uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Exit', 'Position', [220, 5, 130, 45], 'Callback', @Button);
    % Criar um botão na GUI para abrir o arquivo de ajuda
    uicontrol('Style', 'pushbutton', ...
                            'String', 'Help', ...
                            'Position', [1110, 10, 80, 30], ...
                            'FontSize', 10, ...
                            'Callback', @openHelpFile);

    % Função para abrir o arquivo de ajuda
    function openHelpFile(~, ~)
        helpFile = 'help.txt'; % Nome do arquivo de ajuda

        % Verifica se o arquivo existe antes de tentar abrir
        if exist(helpFile, 'file')
            winopen(helpFile); % Abre o arquivo no Windows
        else
            errordlg('File help.txt not found.', 'Error');
        end
    end

    % Plot Panels
    plotPanel = uipanel('Title', 'Plots', 'FontSize', 10, 'Position', [0.4 0.1 0.55 0.9]);
    hpNozzle = axes('Parent', plotPanel, 'Units', 'normalized', 'Position', [0.1, 0.7, 0.8, 0.25]);
    hpPressure = axes('Parent', plotPanel, 'Units', 'normalized', 'Position', [0.1, 0.48, 0.8, 0.15]);
    hpMach = axes('Parent', plotPanel, 'Units', 'normalized', 'Position', [0.1, 0.25, 0.78, 0.15]);
    hpTemperature = axes('Parent', plotPanel, 'Units', 'normalized', 'Position', [0.1, 0.03, 0.8, 0.15]);

    % Make GUI visible
    set(GUI, 'Visible', 'on');

    % Activation Functions
    function Activate_AR(source, ~), AR = str2double(get(source, 'String')); end
    function Activate_p0(source, ~), p0 = str2double(get(source, 'String')); end
    function Activate_T0(source, ~), T0 = str2double(get(source, 'String')); end
    function Activate_pe(source, ~), pe = str2double(get(source, 'String')); end
    function Activate_L(source, ~), TubeLenght = str2double(get(source, 'String')); end
    function Activate_Epsilon(source, ~), epsilon = str2double(get(source, 'String')); end
    function Activate_Gamma(source, ~), gamma = str2double(get(source, 'String')); end


    % Exit Function
    function Button(~, ~)
        fprintf('Closing.\n');
        close all force;
    end

    % Restart Function
    function Restart(~, ~)
        clc;
        GUI_Main;
    end

    % **Main Solver Function**
    function Plot(~, ~)
        % **Validation & Warnings**
        if T0 <= 0 || p0 <= 0 || pe <= 0 || gamma <= 1 || gamma > 2
            msgbox('Please input different parameters.', 'Input Error', 'error');
            return;
        elseif AR <= 1
            msgbox('Throat area cannot be greater than the exit area. Please input different parameters.', 'Input Error', 'error');
            return;
        elseif TubeLenght > 10000
            msgbox('Model is not valid for the given parameters. Please introduce different inputs.', 'Input Error', 'error');
            return;
        elseif pe>=p0
            msgbox('No flow or inverse flow pe>=po', 'Input Error', 'error');
            return;
        end
        
        fprintf('Running Converging-Diverging Nozzle Solver with Adiabatic Pipe...\n');
        GUI_Fanno_Flow_Nozzle2(AR, p0, T0, pe, TubeLenght, epsilon,gamma, hpNozzle, hpPressure, hpMach, hpTemperature);
    end
end

