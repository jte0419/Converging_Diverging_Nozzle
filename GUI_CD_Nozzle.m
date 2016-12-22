% CONVERGING DIVERGING NOZZLE
% Written by: JoshTheEngineer
% Started: 10/07/16
% Updated: 10/07/16 - Started code
%          10/07/16 - Works as intended
%          10/09/16 - Made some functions to make code easier to read
%          10/09/16 - Added M and P plots, and shock nozzle plot

function varargout = GUI_CD_Nozzle(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_CD_Nozzle_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_CD_Nozzle_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% --- Executes just before GUI_CD_Nozzle is made visible.
function GUI_CD_Nozzle_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_CD_Nozzle_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% ON STARTUP
set(handles.figureCDNozzle,'Resize','On');
set(handles.figureCDNozzle,'Units','Normalized');
set(handles.plotP,'XTickLabel',[]);
set(handles.plotM,'XTickLabel',[]);
% AXES
set(handles.plotData,'Units','Normalized');
set(handles.plotP,'Units','Normalized');
set(handles.plotM,'Units','Normalized');
% CHECK BOX
set(handles.checkPlotData,'Units','Normalized');
% EDIT TEXT
set(handles.editGamma,'Units','Normalized');
set(handles.editAreaRatio,'Units','Normalized');
set(handles.editPressureRatio,'Units','Normalized');
set(handles.editMPPts,'Units','Normalized');
% FIGURE
% LIST
% PANEL
set(handles.panelInputs,'Units','Normalized');
set(handles.panelSolutions,'Units','Normalized');
set(handles.panelStatus,'Units','Normalized');
set(handles.panelPlotting,'Units','Normalized');
% POP
% PUSH
set(handles.pushExit,'Units','Normalized');
% RADIO
% TABLE
% TEXT
set(handles.textTitle,'Units','Normalized');
set(handles.textChoked,'Units','Normalized');
set(handles.textState,'Units','Normalized');
set(handles.textMSub,'Units','Normalized');
set(handles.textMSup,'Units','Normalized');
set(handles.textPePoSub,'Units','Normalized');
set(handles.textPePoSup,'Units','Normalized');
set(handles.textPePoNS,'Units','Normalized');
set(handles.textA_Astar,'Units','Normalized');
set(handles.textChokedSolution,'Units','Normalized');
set(handles.textStateSolution,'Units','Normalized');
set(handles.textMSubSolution,'Units','Normalized');
set(handles.textMSupSolution,'Units','Normalized');
set(handles.textPePoSubSolution,'Units','Normalized');
set(handles.textPePoSupSolution,'Units','Normalized');
set(handles.textPePoNSSolution,'Units','Normalized');
set(handles.textA_AstarSolution,'Units','Normalized');
set(handles.textStatus,'Units','Normalized');
set(handles.textMPPts,'Units','Normalized');

% Set tool tip strings
set(handles.textGamma,'TooltipString',...
    'Specific heat ratio');
set(handles.textAe_At,'TooltipString',...
    'Nozzle exit area to throat area');
set(handles.textPe_Po,'TooltipString',...
    'Exit-to-reservoir pressure ratio');
set(handles.textChoked,'TooltipString',...
    'Indicates whether or not the nozzle is choked');
set(handles.textState,'TooltipString',...
    'Indicates what type of process is happening through nozzle');
set(handles.textMSub,'TooltipString',...
    'Isentropic subsonic Mach number at nozzle exit');
set(handles.textMSup,'TooltipString',...
    'Isentropic supersonic Mach number at nozzle exit');
set(handles.textPePoSub,'TooltipString',...
    'Isentropic subsonic exit-to-reservoir pressure ratio');
set(handles.textPePoSup,'TooltipString',...
    'Isentropic supersonic exit-to-reservoir pressure ratio');
set(handles.textPePoNS,'TooltipString',...
    'Exit-to-reservoir pressure ratio for normal shock at nozzle exit');
set(handles.textA_Astar,'TooltipString',...
    'Normal shock location in the nozzle (if applicable)');
set(handles.textMPPts,'TooltipString',...
    'Number of nozzle points to plot on P and M plots');
set(handles.checkPlotData,'TooltipString',...
    'If not selected, only nozzle will plot (no P or M plots)');

% Call the solve function on startup to get default solution
SOLVE(handles);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% --------------------------- INITIALIZATION ---------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% EDIT ---------------- Specific Heat Ratio -------------------------------
function editGamma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT -------------------- Area Ratio ------------------------------------
function editAreaRatio_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT ------------------ Pressure Ratio ----------------------------------
function editPressureRatio_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT --------------------- Nozzle Points --------------------------------
function editNozzlePts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT ---------------------- M & P Points --------------------------------
function editMPPts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT ------------------- Reservoir Points -------------------------------
function editResPts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------- CALLBACKS ------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% CHECK ------------------- Plot Data -------------------------------------
function checkPlotData_Callback(hObject, eventdata, handles)
% Disable plot points edit text box if deselected
if (get(hObject,'Value') == 1)
    set(handles.editMPPts,'Enable','on');
elseif (get(hObject,'Value') == 0)
    set(handles.editMPPts,'Enable','off');
end

% EDIT ---------------- Specific Heat Ratio -------------------------------
function editGamma_Callback(hObject, eventdata, handles)
% Call the SOLVE function
SOLVE(handles);

% EDIT -------------------- Area Ratio ------------------------------------
function editAreaRatio_Callback(hObject, eventdata, handles)
% Call the SOLVE function
SOLVE(handles);

% EDIT ------------------ Pressure Ratio ----------------------------------
function editPressureRatio_Callback(hObject, eventdata, handles)
% Call the SOLVE function
SOLVE(handles);

% EDIT --------------------- Nozzle Points --------------------------------
function editNozzlePts_Callback(hObject, eventdata, handles)
% Do nothing, use this value in other callbacks

% EDIT ---------------------- M & P Points --------------------------------
function editMPPts_Callback(hObject, eventdata, handles)
% Do nothing, use this value in other callbacks

% EDIT ------------------- Reservoir Points -------------------------------
function editResPts_Callback(hObject, eventdata, handles)
% Do nothing, use this value in other callbacks

% FUNCTION ----------------- S O L V E ------------------------------------
function [] = SOLVE(handles)
% =========================================================================
% - Load all input parameters
% - Solve the relevant system
% - Plot the nozzle with shock system if appropriate
% - Plot the P ratio and M as a function of position
% - Display solution to screen
% =========================================================================

% Indicate status that we are solving
set(handles.textStatus,'String','Solving...',...
                       'ForegroundColor','r');
drawnow();

% Load relevant data
g     = str2double(get(handles.editGamma,'String'));                        % Ratio of specific heats
Ae_At = str2double(get(handles.editAreaRatio,'String'));                    % Exit/throat area ratio
Pe_Po = str2double(get(handles.editPressureRatio,'String'));                % Exit/reservoir pressure ratio

% Some relevant variables for simplification later
gm1   = g-1;
gp1   = g+1;
gm1o2 = gm1/2;
gogm1 = g/gm1;
gogp1 = g/gp1;

% Solve for subsonic and supersonic Mach numbers
Msub = A_M_RELATION(Ae_At,0,g,'Sub');                                       % Subsonic Mach number
Msup = A_M_RELATION(Ae_At,0,g,'Sup');                                       % Supersonic Mach number

% Pressure ratios for sub and sup
Pe_Po_Sub = (1 + gm1o2*Msub^2)^(-gogm1);                                    % Subsonic pressure ratio
Pe_Po_Sup = (1 + gm1o2*Msup^2)^(-gogm1);                                    % Supersonic pressure ratio

% Pressure ratios for normal shock at exit
Pe_P2    = 1 + ((2*gogp1)*(Msup^2-1));                                      % NS pressure ratio
Pe_Po_NS = (Pe_P2)*(Pe_Po_Sup);                                             % Back-to-reservoir pressure ratio

% Solve normal shock location if (Pe_Po < Pe_Po_sub && Pe_Po > Pe_Po_NS)
% - Use direct method (as opposed to iterative method)
% - See page 215 of Anderson
if (Pe_Po < Pe_Po_Sub && Pe_Po > Pe_Po_NS)
    % Solve for exit Mach number
    term1 = -(1/gm1);
    term2 = 1/(gm1^2);
    term3 = 2/gm1;
    term4 = (2/gp1)^(gp1/gm1);
    term5 = (Pe_Po*Ae_At)^-2;
    Me    = sqrt(term1 + sqrt(term2 + term3*term4*term5));                  % Eq. 5.28
    assignin('base','Me_NS',Me);
    
    % Exit stag-to-stat pressure ratio (isentropic)
    Poe_Pe = (1 + (gm1o2*Me^2))^(g/gm1);
    
    % Get Po2/Po1 pressure ratio
    Po2_Po1 = Poe_Pe*Pe_Po;                                                 % Eq. 5.29
    
    % Get corresponding upstream Mach number from this pressure ratio
    problem.objective = @(M) exp(-gogm1*log(((2 + gm1*M^2)/(gp1*M^2))*...
                                    (1+(2*g/gp1)*(M^2-1))) + ...
                                    log((1+(2*g/gp1)*(M^2-1)))) - Po2_Po1;
    problem.x0        = [1+1e-5 50];
    problem.solver    = 'fzero';
    problem.options   = optimset(@fzero);
    M1                = fzero(problem);
    
    num = 1 + gm1o2*M1^2;
    den = (g*M1^2) - gm1o2;
    M2  = sqrt(num/den);
    
    % Solve for the area ratios where the normal shock is located
    A1_Astar     = A_M_RELATION(0,M1,g,'');                                 % Area ratio at shock (pre-shock)
    A2_Astarstar = A_M_RELATION(0,M2,g,'');                                 % Area ratio at shock (post-shock)
    Ae_Astarstar = A_M_RELATION(0,Me,g,'');                                 % Area ratio at exit
    
    % Need this area ratio in plotting the normal shock
    assignin('base','A1_Astar',A1_Astar);
    assignin('base','A2_Astarstar',A2_Astarstar);
    assignin('base','Ae_Astar',Ae_At);
    assignin('base','Ae_Astarstar',Ae_Astarstar);
end

% Print solution to GUI
set(handles.textMSubSolution,'String',num2str(Msub));
set(handles.textMSupSolution,'String',num2str(Msup));
set(handles.textPePoSubSolution,'String',num2str(Pe_Po_Sub));
set(handles.textPePoSupSolution,'String',num2str(Pe_Po_Sup));
set(handles.textPePoNSSolution,'String',num2str(Pe_Po_NS));
set(handles.textA_AstarSolution,'String','-');

% Print choked and state conditions to GUI
if (Pe_Po >= Pe_Po_Sub)
    CDState = 'ISEN_SUB';
    set(handles.textChokedSolution,'String','NO');
    set(handles.textStateSolution,'String','Isentropic Subsonic');
elseif (Pe_Po < Pe_Po_Sub && Pe_Po > Pe_Po_NS);
    CDState = 'NS_NOZZLE';
    set(handles.textChokedSolution,'String','YES');
    set(handles.textStateSolution,'String','Normal Shock in Nozzle');
    set(handles.textA_AstarSolution,'String',num2str(A1_Astar));
elseif (Pe_Po == Pe_Po_NS);
    CDState = 'NS_EXIT';
    set(handles.textChokedSolution,'String','YES');
    set(handles.textStateSolution,'String','Normal Shock at Exit');
elseif (Pe_Po < Pe_Po_NS && Pe_Po > Pe_Po_Sup);
    CDState = 'OE';
    set(handles.textChokedSolution,'String','YES');
    set(handles.textStateSolution,'String','Overexpanded');
elseif (Pe_Po == Pe_Po_Sup)
    CDState = 'ISEN_SUP';
    set(handles.textChokedSolution,'String','YES');
    set(handles.textStateSolution,'String','Isentropic Supersonic');
elseif (Pe_Po < Pe_Po_Sup)
    CDState = 'UE';
    set(handles.textChokedSolution,'String','YES');
    set(handles.textStateSolution,'String','Underexpanded');
end

% Assign relevant variables into base workspace
assignin('base','Pe_Po_Sub',Pe_Po_Sub);
assignin('base','Pe_Po_Sup',Pe_Po_Sup);
assignin('base','Pe_Po_NS',Pe_Po_NS);
assignin('base','CDState',CDState);

% Call function to plot data and nozzle
PLOT(handles);

% FUNCTION ----------- Data and Nozzle Plotting ---------------------------
function [] = PLOT(handles)
% =========================================================================
% - Plot the nozzle shape contour
% - Plot the shock location in the nozzle if appropriate
% - Plot the Mach number vs. position
% - Plot the pressure ratio vs. position
% =========================================================================

% Indicate status that we are plotting
set(handles.textStatus,'String','Plotting...',...                           % Set status text
                       'ForegroundColor','r');                              % Change status text to red
drawnow();                                                                  % Make sure it displays immediately

% Check box values
checkPlot = get(handles.checkPlotData,'Value');                             % Get check box value
if (checkPlot == 0)
    % Clear plots so we don't get confused by previous plots
    axes(handles.plotP);    cla;                                            % Clear pressure ratio plot
    axes(handles.plotM);    cla;                                            % Clear Mach number plot
    axes(handles.plotData); cla;                                            % Clear nozzle plot
    
    % Indicate status that we are done
    set(handles.textStatus,'String','Done!',...                             % Set status text
                           'ForegroundColor','k');                          % Change status text to black
    drawnow();                                                              % Make sure it displays immediately
    return;                                                                 % Don't plot if user doesn't want to
end

% Get values from GUI edit text boxes
CDState   = evalin('base','CDState');                                       % State of the nozzle
g         = str2double(get(handles.editGamma,'String'));                    % Specific heat ratio
Ae_At     = str2double(get(handles.editAreaRatio,'String'));                % Area ratio
Pe_Po     = str2double(get(handles.editPressureRatio,'String'));            % Pressure ratio
numMPPts  = str2double(get(handles.editMPPts,'String'));                    % Number of P ratio and M points for plotting
numResPts = 25;                                                             % Number of reservoir points for plotting

% Bounds for accessing reservoir and nozzle indices
iResS = 1;                                                                  % Reservoir starting index
iResE = iResS + numResPts - 1;                                              % Reservoir ending index
iNozS = iResE + 1;                                                          % Nozzle starting index
iNozE = iNozS + numMPPts -1;                                                % Nozzle ending index

% Some relevant variables for simplification later
gm1   = g-1;
gp1   = g+1;
gm1o2 = gm1/2;
gm1og = gm1/g;
gp1o2 = gp1/2;
gogm1 = g/gm1;

% Load relevant variables for normal shock in nozzle case
if (strcmpi(CDState,'NS_NOZZLE'))
    A1_Astar     = evalin('base','A1_Astar');                               % Area ratio pre-shock
    A2_Astarstar = evalin('base','A2_Astarstar');                           % Area ratio post-shock
    Ae_Astarstar = evalin('base','Ae_Astarstar');                           % Area ratio at exit
end

% Initialize solution arrays
PPLOT = zeros(numMPPts,1);                                                  % Initialize pressure ratio array
MPLOT = zeros(numMPPts,1);                                                  % Initialize Mach number array

% SET GEOMETRY HERE - NOZZLE AND RESERVOIR
b     = 0.5;                                                                % Throat radius
minX  = -85;                                                                % Minimum X [deg]
midX  = 0;                                                                  % Middle X
maxX  = 1;                                                                  % Maximum X
XRES  = linspace(minX,midX,numResPts)';                                     % Reservoir X points
XNOZ  = linspace(midX,maxX,numMPPts)';                                      % Nozzle X points
ARES  = (-tand(XRES) + 0.5)*2;                                              % Reservoir area ratio
ANOZ  = linspace(1,Ae_At,numMPPts)';
scFc  = 0.25*(b)/abs(minX);                                                 % Scale factor
XPLOT = [(XRES).*scFc; XNOZ];                                               % Combine RES and NOZ for total X array
APLOT = [ARES; ANOZ];                                                       % Combine Res and NOZ for total area ratio array

% ==============================================
% ============= PLOTTING NOZZLE ================
% ==============================================
axes(handles.plotData);                                                     % Select appropriate axes
cla; hold on; grid on;                                                      % Get ready for plotting
xNoz = [0 1];                                                               % Nozzle X points
yNoz = [b (Ae_At/2)];                                                       % Nozzle Y points
plot(xNoz,yNoz,'k-',xNoz,-yNoz,'k-','LineWidth',3);                         % Plot nozzle
ang = linspace(-80,0,numResPts)';                                           % Angle for tangent calc
yR  = ((-tand(ang)).*(yNoz(end)/25))+b;                                     % Reservoir Y-values
xR  = ang.*(0.25*(b)/abs(minX));                                            % Reservoir X-values
plot(xR,yR,'k-',xR,-yR,'k-','LineWidth',3);                                 % Plot reservoir
plot([xNoz(end) xNoz(end)],[yNoz(end) -yNoz(end)],'k-','LineWidth',2);      % Plot line at nozzle exit
plot([0 0],[-b b],'k-','LineWidth',2);                                      % Plot line at nozzle throat
if (strcmpi(CDState,'NS_NOZZLE'))                                           % Normal shock in nozzle
    NSLoc = (A1_Astar-1)/(Ae_At-1);
    xNS   = [NSLoc NSLoc];                                                  % Construct X line
    yNS   = [0 (((Ae_At/2)-b)*NSLoc)+b];                                    % Construct Y line
    plot(xNS,yNS, 'r-',xNS,-yNS,'r-','LineWidth',3);                        % Plot shock
elseif (strcmpi(CDState,'OE'));                                             % Overexpanded (Oblique)
    xOS = [xNoz(end) xNoz(end)+0.25*xNoz(end)];                             % Construct X line
    yOS = [yNoz(end) 0.75*yNoz(end)];                                       % Construct Y line
    plot(xOS,yOS,'m-',xOS,-yOS,'m-','LineWidth',3);                         % Plot shock
elseif (strcmpi(CDState,'UE'))                                              % Underexpanded (P-M)
    xPM  = [xNoz(end) xNoz(end)+0.25*xNoz(end)];                            % Construct X line
    yPM1 = [yNoz(end) 0.65*yNoz(end)];                                      % Construct Y line 1
    yPM2 = [yNoz(end) 0.75*yNoz(end)];                                      % Construct Y line 2
    yPM3 = [yNoz(end) 0.85*yNoz(end)];                                      % Construct Y line 3
    plot(xPM,yPM1,'c-',xPM,-yPM1,'c-','LineWidth',3);                       % Plot PM wave 1
    plot(xPM,yPM2,'c-',xPM,-yPM2,'c-','LineWidth',3);                       % Plot PM wave 2
    plot(xPM,yPM3,'c-',xPM,-yPM3,'c-','LineWidth',3);                       % Plot PM wave 3
end
xlim([xR(1) max(xNoz)+0.25*max(xNoz)]);
ylim([-yNoz(end)-0.1*yNoz(end) yNoz(end)+0.1*yNoz(end)]);                   % Set Y-axis limits
% ==================================================
% ============= END PLOTTING NOZZLE ================
% ==================================================

% ===============================================
% ============= PLOTTING P AND M ================
% ===============================================
% ---- CHOKED SUBSONIC AND SUPERSONIC ISENTROPIC SOLUTIONS ----
% - Always calculate to show dashed lines on plots

% Initialize arrays
MPLOTSub = zeros(numMPPts,1);
MPLOTSup = zeros(numMPPts,1);
MPLOTRes = zeros(numResPts,1);

% Solve for Mach number of subsonic and supersonic branches
for i = 1:1:numMPPts
    if (ANOZ(i) == 1)                                                       % If we are at the throat
        MPLOTSub(i) = 1;                                                    % Set throat (choked) Mach number to unity
        MPLOTSup(i) = 1;                                                    % Set throat (choked) Mach number to unity
    else
        MPLOTSub(i) = A_M_RELATION(ANOZ(i),0,g,'Sub');                      % Subsonic Mach numbers from area ratio
        MPLOTSup(i) = A_M_RELATION(ANOZ(i),0,g,'Sup');                      % Supersonic Mach numbers from area ratio
    end
end

% Solve for P ratio
PPLOTSub = ((1 + gm1o2.*MPLOTSub.^2).^(-gogm1))';                           % Subsonic pressure ratios from Mach numbers
PPLOTSup = ((1 + gm1o2.*MPLOTSup.^2).^(-gogm1))';                           % Supersonic pressure ratios from Mach numbers

% Reservoir
for i = 1:1:numResPts
    MPLOTRes(i,1) = A_M_RELATION(APLOT(i),0,g,'Sub');
end
PPLOTRes = ((1 + gm1o2.*MPLOTRes.^2).^(-gogm1))';
% ------------------------------------------------------------

% ------------------ SUBSONIC ISENTROPIC --------------
if (strcmpi(CDState,'ISEN_SUB'))
    Me       = sqrt((2/gm1)*(Pe_Po^(-gm1og)-1));
    Ae_Astar = A_M_RELATION(0,Me,g,'Sub');
    At_Astar = (1/Ae_At)*(Ae_Astar);
    Mt       = A_M_RELATION(At_Astar,0,g,'Sub');
    Pt_Po    = (1 + gm1o2*Mt^2)^(-gogm1);
    
    for i = 1:1:length(ARES)
        Aarbres_Astar(i,1) = ARES(i)*(1/Ae_At)*Ae_Astar;
        Marbres(i,1)       = A_M_RELATION(Aarbres_Astar(i),0,g,'Sub');
    end
    Pt_Po_Arbres = (1 + gm1o2*Marbres.^2).^(-gogm1);
    
    for i = 1:1:length(ANOZ)
        Aarb_Astar(i,1) = ANOZ(i)*(1/Ae_At)*Ae_Astar;
        Marb(i,1)       = A_M_RELATION(Aarb_Astar(i),0,g,'Sub');
    end
    Pt_Po_Arb = (1 + gm1o2*Marb.^2).^(-gogm1);
    
    assignin('base','Aarb_Astar',Aarb_Astar);
    assignin('base','Marb',Marb);
    assignin('base','Pt_Po_Arb',Pt_Po_Arb);
    
    assignin('base','Me',Me);
    assignin('base','Ae_Astar',Ae_Astar);
    assignin('base','At_Astar',At_Astar);
    assignin('base','Mt',Mt);
    assignin('base','Pt_Po',Pt_Po);
end
% -----------------------------------------------------

% ----------------- NORMAL SHOCK IN NOZZLE-------------
if (strcmpi(CDState,'NS_NOZZLE'))
    pctNS        = (A1_Astar-1)/(Ae_At-1);
    midMPPts     = pctNS*maxX;
    numPreShock  = floor(pctNS*numMPPts);
    if (numPreShock == 0)
        numPreShock = 1;
    end
    numPostShock = numMPPts - numPreShock;
    XNOZ_PRE     = linspace(0,midMPPts,numPreShock)';
    ANOZ_PRE     = linspace(1,A1_Astar,numPreShock)';
    XNOZ_POST    = linspace(midMPPts,maxX,numPostShock)';
    ANOZ_POST    = linspace(A2_Astarstar,Ae_Astarstar,numPostShock)';
    
    % Initialize pre- and post-shock arrays
    MPLOT_PRE  = zeros(numPreShock,1);
    MPLOT_POST = zeros(numPostShock,1);
    PPLOT_PRE  = zeros(numPreShock,1);
    PPLOT_POST = zeros(numPostShock,1);
    
    % Pre-shock Mach number and pressure ratio
    for i = 1:1:numPreShock
        % Solve for Mach number pre-shock
        if (ANOZ_PRE(i) == 1)
            MPLOT_PRE(i) = 1;
        else
            MPLOT_PRE(i) = A_M_RELATION(ANOZ_PRE(i),0,g,'Sup');
        end
        
        % Solve for pressure ratio pre-shock
        PPLOT_PRE(i) = (1 + gm1o2*MPLOT_PRE(i)^2)^(-gogm1);
    end
    
    % Post-shock Mach number and pressure ratio
    for i = 1:1:numPostShock
        % Solve for Mach number pre-shock
        if (ANOZ_POST(i) == 1)
            MPLOT_POST(i) = 1;
        else
            MPLOT_POST(i) = A_M_RELATION(ANOZ_POST(i),0,g,'Sub');
        end
        
        % Solve for pressure ratio pre-shock
        PPLOT_POST(i) = (1 + gm1o2*MPLOT_POST(i)^2)^(-gogm1);
    end
end
% -----------------------------------------------------

% Solve for critical expansion ratio (gamma dependent)
Pcrit = (gp1o2)^(-gogm1);

% PLOT: Pressure vs. Position
axes(handles.plotP);
cla; hold on; grid on;
plot([XPLOT(1) XPLOT(end)],[Pcrit Pcrit],'k--','LineWidth',2);              % Plot throat P ratio
plot(XPLOT(iNozS:iNozE),PPLOTSub,'k--','LineWidth',2);                      % Isentropic subsonic
plot(XPLOT(iNozS:iNozE),PPLOTSup,'k--','LineWidth',2);                      % Isentropic supersonic
plot(XPLOT(iResS:iResE),PPLOTRes,'k--','LineWidth',2);                      % Reservoir
if (strcmpi(CDState,'OE') || strcmpi(CDState,'UE'))
    plot(XPLOT(iResS:iResE),PPLOTRes,'r-','LineWidth',2);
    plot(XPLOT(iNozS:iNozE),PPLOTSup,'r-','LineWidth',2);
    plot(XPLOT(iNozE),Pe_Po,'ro','MarkerFaceColor','r',...
                                 'MarkerEdgeColor','k');
elseif (strcmpi(CDState,'ISEN_SUB'))
    plot(XPLOT(iNozS:iNozE),Pt_Po_Arb,'r-','LineWidth',2);
    plot(XPLOT(iResS:iResE),Pt_Po_Arbres,'r-','LineWidth',2);
elseif (strcmpi(CDState,'NS_NOZZLE'))
    assignin('base','XNOZ_PRE',XNOZ_PRE);
    assignin('base','PPLOT_PRE',PPLOT_PRE);
    plot(XPLOT(iResS:iResE),PPLOTRes,'r-','LineWidth',2);
    plot(XNOZ_PRE,PPLOT_PRE,'r-','LineWidth',2);
    plot(XNOZ_POST,PPLOT_POST,'r-','LineWidth',2);
    plot([XNOZ_PRE(end) XNOZ_POST(1)],[PPLOT_PRE(end) PPLOT_POST(1)],...
                'r-','LineWidth',2);
else
    plot(XPLOT,PPLOT,'r-','LineWidth',2);
end
xlim([XPLOT(1) XPLOT(end)]);                                                % X-axis limits
xlabel('Position');                                                         % X label
ylabel('Pe/Po');                                                            % Y label
set(handles.plotP,'XTickLabel',[]);                                         % Get rid of X-axis tick marks

% PLOT: Mach vs. Position
axes(handles.plotM);
cla; hold on; grid on;
plot(XPLOT(iNozS:iNozE),MPLOTSub,'k--','LineWidth',2);                      % Isentropic subsonic
plot(XPLOT(iNozS:iNozE),MPLOTSup,'k--','LineWidth',2);                      % Isentropic supersonic
plot(XPLOT(iResS:iResE),MPLOTRes,'k--','LineWidth',2);                      % Reservoir
if (strcmpi(CDState,'OE') || strcmpi(CDState,'UE'))                         % For overexpanded or underexpanded solutions
    plot(XPLOT(iResS:iResE),MPLOTRes,'b-','LineWidth',2);
    plot(XPLOT(iNozS:iNozE),MPLOTSup,'b-','LineWidth',2);
elseif (strcmpi(CDState,'ISEN_SUB'))
    plot(XPLOT(iResS:iResE),Marbres,'b-','LineWidth',2);
    plot(XPLOT(iNozS:iNozE),Marb,'b-','LineWidth',2);
elseif (strcmpi(CDState,'NS_NOZZLE'))
    plot(XPLOT(iResS:iResE),MPLOTRes,'b-','LineWidth',2);
    plot(XNOZ_PRE,MPLOT_PRE,'b-','LineWidth',2);
    plot(XNOZ_POST,MPLOT_POST,'b-','LineWidth',2);
    plot([XNOZ_PRE(end) XNOZ_POST(1)],[MPLOT_PRE(end) MPLOT_POST(1)],...
                'b-','LineWidth',2);
else
    plot(XPLOT(iNozS:iNozE),MPLOT,'b-','LineWidth',2);
end
xlim([XPLOT(1) XPLOT(end)]);
xlabel('Position');
ylabel('M');
set(handles.plotM,'XTickLabel',[]);

% ===================================================
% ============= END PLOTTING P AND M ================
% ===================================================

% Indicate status that we are done
set(handles.textStatus,'String','Done!',...
                       'ForegroundColor','k');
drawnow();

% FUNCTION ----------- Area-Mach Number Relation --------------------------
function [sol] = A_M_RELATION(ARatio,M,g,SubSup)
% =========================================================================
% - Solve for either the area ratio or Mach number
%   - Depends on what the user inputs are
% =========================================================================

% Set initial value of solution so it doesn't error out
sol = inf;

% Get and set convenient variables
gp1   = g + 1;
gm1   = g - 1;
gm1o2 = gm1/2;

% Solve for area ratio or Mach number
if (ARatio == 0)                                                            % SOLVE: AREA RATIO
    sol = sqrt((1/(M^2))*(((2/gp1)*(1+gm1o2*M^2)).^(gp1/gm1)));             % Solve
elseif (M == 0)                                                             % SOLVE: MACH NUMBER
    if (strcmpi(SubSup,'Sub'))                                              % Subsonic value
        problem.objective = @(M) sqrt((1/(M^2))*(((2/gp1)*...
                                    (1+gm1o2*M^2))^(gp1/gm1))) - ARatio;    % Objective function
        problem.x0        = [1e-6 1];                                       % Solver bounds
        problem.solver    = 'fzero';                                        % Find the zero
        problem.options   = optimset(@fzero);                               % Default options
        sol               = fzero(problem);                                 % Solve
    elseif (strcmpi(SubSup,'Sup'))                                          % Supersonic value
        problem.objective = @(M) sqrt((1/(M^2))*(((2/gp1)*...
                                    (1+gm1o2*M^2))^(gp1/gm1))) - ARatio;    % Objective function
        problem.x0        = [1 50];                                         % Solver Bounds
        problem.solver    = 'fzero';                                        % Find the zero
        problem.options   = optimset(@fzero);                               % Default options
        sol               = fzero(problem);                                 % Solve
    end
end

if (sol == inf)
    fprintf('Error in A_M_RELATION function\n\n');
end

% PUSH -------------------- Exit the GUI ----------------------------------
function pushExit_Callback(hObject, eventdata, handles)
clc;
delete(handles.figureCDNozzle);
