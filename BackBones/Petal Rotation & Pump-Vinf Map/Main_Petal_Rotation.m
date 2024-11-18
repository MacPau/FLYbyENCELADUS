%% INITIALIZATION %%
clear all; close all; clc; format long g;

% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [1]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
                            % (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
                            % (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

% Define Flyby Hyperbolic Excess Velocities
pars.INPUTS.V_inf = 4.0;     %[km/s] 

% Define the crank angle to consider
pars.INPUTS.crank = [0, pi];

% Define the resonances to consider for the Petal Rotation
pars.INPUTS.Resonances = [7 1];

% Define the Petal Rotation Inputs 
pars.PetalRot.Starting_Reso  = [7 1];    % Resonance Orbit at Start of Petal Rotation
pars.PetalRot.Starting_Crank = 0;        % [rad] Crank angle at start of Petal Rotation (0 = Outbound / 180° = Inbound)
pars.PetalRot.Starting_Vinf  = 4.0;      % [km/s] Vinf at Start of Petal Rotation
pars.PetalRot.TrueLong_Obj   = deg2rad(10); % [rad] Objective of change in true longitude
pars.PetalRot.Direction      = +1;      % Direction in which to rotate (+1 = clockwise / -1 = anticlockwise)
pars.PetalRot.Num_Flybys     = 6;       % Maximum n° of flybys that the Petal Rotation can use
pars.PetalRot.Beam_Width     = 1000;    % N° of solutions to store at each Petal Rotation expansion step

% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 15e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 60;   %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 
pars.INPUTS.maxDV_Defect   = 0.15;  %[km/s] Maximum DV defect allowed for a flyby

% Define starting epoch
pars.INPUTS.epoch0 = 0;       % Days passed since MJD2000 

% Retrieve Central Body (Planet) Parameters
[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral); %[km3/s2],[km],[km] & [km]

% Retrieve Desired Moon Parameters 
if pars.INPUTS.idCentral == 3
    pars.Moon.OrbRad = 384748; pars.Moon.mu  = getAstroConstants('Moon','mu'); %[km],[km3/s2]
    pars.Moon.EquRad = getAstroConstants('Moon','Radius'); pars.Moon.hmin = 50;  %[km], [km]
elseif pars.INPUTS.idCentral == 5
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = jupMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
elseif pars.INPUTS.idCentral == 6
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
else
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = uranusMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
end

for idx = 1:size(pars.INPUTS.idMoon,2)
    pars.Moon.Vel(idx)    = sqrt(pars.Planet.mu/pars.Moon.OrbRad(idx));           %[km/s] Moon Orbital velocity
    pars.Moon.Period(idx) = 2*pi*sqrt(pars.Moon.OrbRad(idx)^3/pars.Planet.mu);    %[s] Moon orbital period
    pars.Moon.HillSph(idx) = pars.Moon.OrbRad(idx)*( pars.Moon.mu(idx)/(3*(pars.Moon.mu(idx) + pars.Planet.mu)))^(1/3);     %[km] Moon Hill's Sphere
end

% Define other parameters
pars.AU     = 1.49597870691*10^(8);   %[km]
pars.g0     = 9.80665;                %[m/s^2] 
pars.Day    = 86400;                  %[s]
pars.Year   = 365.25;                 %[days]
pars.JD     = 2400000.5;

% Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;    %[km]
V_infs = pars.INPUTS.V_inf;
for i = 1:length(V_infs)
    e_fly     = 1 + ((rp_flyby*V_infs(i)^2)/pars.Moon.mu); %[-] 
    delta_max(i) = 2*asin(1/e_fly);                        %[rad]
end
pars.delta_max = delta_max;
pars.rp_flyby = rp_flyby;

% Load sectors for groundtrack plots
load('sectorObj.mat');
load('sectorObj_.mat');
pars.sectorObj  = sectorObj;
pars.sectorObj_ = sectorObj_;
pars.rect       = defineRectangleMapping(); % --> save the rectangle mapping

%% COMPUTE DATABASE OF RESONANT TRANSFERS %%
Reso_TransfersDatabase = [];
for i = 1:size(pars.INPUTS.V_inf, 2)
    [Reso_Transfers, ResoStruc] = ResoTransfers_Generation(pars.INPUTS.V_inf(i), pars);
    Reso_TransfersDatabase = [Reso_TransfersDatabase; Reso_Transfers];
end

%% COMPUTE DATABASE OF PSEUDO-RESONANT TRANSFERS %%
% Define which 81 pseudo-resonance to use
pars.INPUTS.remove81 = 1;  %[0] = add 1 rev to the outbound-inbound short case
                           %[1] = do not add one rev to the outbound-inbound short case

PseudoReso_TransfersDatabase = [];
for i = 1:size(pars.INPUTS.V_inf, 2)
    [PseudoReso_Transfer, PseudoResoStruct] = PseudoResoTransfers_Generation(pars.INPUTS.V_inf(i), pars);
    PseudoReso_TransfersDatabase = [PseudoReso_TransfersDatabase; PseudoReso_Transfer];
end

%% CONSTRUCT THE PETAL ROTATION SEQUENCES %%

[PetalRots_Seqs, PetalRots_FlybysData] = PetalRotation_BuildUp(Reso_TransfersDatabase, PseudoReso_TransfersDatabase, pars);


%% ANALYZE RESULTS %%
%%%% IF YOU WANT TO SEARCH FOR A SPECIFIC SEQUENCE, UNCOMMENT THIS %%%%%
% % Search for a petal rotation doing a specific sequence of transfers
% Specific_Seq = [pars.PetalRot.Starting_Reso; 6 1; pars.PetalRot.Starting_Reso; 6 1; pars.PetalRot.Starting_Reso; 6 1;...
%     pars.PetalRot.Starting_Reso];
% 
% idx_sol = [];
% 
% % Loop through each cell in PetalRots_Seqs
% for i = 1:length(PetalRots_Seqs)
% 
%     Current_seq = PetalRots_Seqs{i};
%     Current_seq = Current_seq(:,1:2);
%     
%     % Check if the sequence matches Specific_Seq
%     if isequal(Current_seq, Specific_Seq)
%         idx_sol = [idx_sol; i];
%         break;
%     end
% end
% idx = idx_sol;


%%% IF NOT, USE THE FOLLOWING CODE TO SELECT WHICH PETAL ROTATION TO PLOT 
% Select the index of a solution to analyze
idx = 1;
PetalRotation_Sequence = cell2mat(PetalRots_Seqs(idx));
PetalRotation_FlybyData = cell2mat(PetalRots_FlybysData(idx));

% Plot the Petal Rotation Planet-Centric Trajectory
[Fig_PetalRot] = Plot_PetalRotation_Trajectory(PetalRotation_Sequence, PetalRotation_FlybyData, pars);

% Plot the Petal Rotation Ground-Track
for i = 2:size(PetalRotation_Sequence, 1)
    node_in = PetalRotation_Sequence(i-1,7:9);
    node_out = PetalRotation_Sequence(i, 3:5);
    Flyby_Data(i) = Flyby_BuildUp(node_in, node_out, pars);
end

Flyby_Data(1) = [];
colors = cool(length(Flyby_Data));
fig1 = figure( 'Color', [1 1 1] );
hold on;
plotTextureLatLong(pars.INPUTS.idMoon , pars.INPUTS.idCentral , 1);
plotSquares(pars, 1);  
axis normal;
for i = 1:size(Flyby_Data, 2)
    Plot_Flyby_GT(Flyby_Data(i), colors(i,:));
end

% Plot Trade-Off Change in True Longitude VS DV for all solutions
for i = 1:size(PetalRots_FlybysData, 2)
    Data_aux = cell2mat(PetalRots_FlybysData(i));

    DV_plot(i)      = sum([Data_aux.DV_Cost]); %[km/s]
    TrueLong_Change = Data_aux(end).Change_TrueLong';
    TrueLong_ChangePlot(i) = rad2deg(abs(TrueLong_Change(end))); %[degrees] 
end

FigComparison = figure('Color', [1 1 1]);clf;set(FigComparison,'defaulttextinterpreter','latex') ;
hold on; grid on; 
plot(DV_plot, TrueLong_ChangePlot, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
xlabel('$\Delta V$-cost [km/s]', 'FontSize', 18);
ylabel('$\Delta$ True Longitude of the Flyby Position wrt Initial [degrees]', 'FontSize', 18);

%% PETAL ROTATION CHARACTERISTICS %%
TOF_total = sum(PetalRotation_Sequence(2:end,6)/86400); %[days]

DTrue_Long = PetalRotation_FlybyData.Change_TrueLong;
Tot_DTrueChange = rad2deg(abs(DTrue_Long(end)));

Ratio_Change = Tot_DTrueChange/TOF_total;  %[degrees/day]
