%% MAIN CODE FOR CRANK OVER THE TOP ANALYSIS %%

%% INITIALIZATION %%
clear all; close all; clc; format long g;
% addpath(genpath([pwd '/functions']));

% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [1]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
                            % (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
                            % (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

% Define Flyby Hyperbolic Excess Velocities
pars.INPUTS.V_inf = 4;     %[km/s] 

% Define the resonances to consider
% pars.INPUTS.Resonances = [5 1; 6 1; 7 1; 8 1];
pars.INPUTS.Resonances = [7 1];

% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 5;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 

% Define starting epoch
pars.INPUTS.epoch0 = 0;    % Days passed since MJD2000

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

for i = 1:size(pars.INPUTS.idMoon,2)
    pars.Moon.Vel(i)    = sqrt(pars.Planet.mu/pars.Moon.OrbRad(i));           %[km/s] Moon Orbital velocity
    pars.Moon.Period(i) = 2*pi*sqrt(pars.Moon.OrbRad(i)^3/pars.Planet.mu);    %[s] Moon orbital period
    pars.Moon.HillSph(i) = pars.Moon.OrbRad(i)*( pars.Moon.mu(i)/(3*(pars.Moon.mu(i) + pars.Planet.mu)))^(1/3);     %[km] Moon Hill's Sphere
end

% Define other parameters
pars.AU     = 1.49597870691*10^(8);   %[km]
pars.g0     = 9.80665;                %[m/s^2] 
pars.Day    = 86400;                  %[s]
pars.Year   = 365.25;                 %[days]
pars.JD     = 2400000.5;

% Load sectors for groundtrack plots
load('sectorObj.mat');
load('sectorObj_.mat');
pars.sectorObj  = sectorObj;
pars.sectorObj_ = sectorObj_;
pars.rect       = defineRectangleMapping(); % --> save the rectangle mapping

%% DEFINE NODES BASED ON CRANK DISCRETIZATION %%
% Compute data of the resonance
[rp_res, ra_res, E_res, alfa_res, i_max, ~] = Tisserand_Resonances(pars.INPUTS.Resonances, pars.INPUTS.V_inf, pars.Moon.OrbRad, pars.Planet.mu, pars);                                    

% Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;           %[km
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);    %[-] 
delta_max = 2*asin(1/e_fly);                                      %[rad]
pars.delta_max = delta_max;

% Determine maximum crank change (since only cranking is done)
alfa_in = alfa_res; alfa_out = alfa_res;
Delta_Crank_max = acos((cos(pars.delta_max) - cos(alfa_in)*cos(alfa_out))/(sin(alfa_in)*sin(alfa_out)));

%% BUILD THE CRANK OVER THE TOP SEQUENCE %%
% Define COT strategy
pars.INPUTS.COT.Start_Crank         = 180; % Starting Crank (INBOUND: 180°, OUTBOUND: 0°)
pars.INPUTS.COT.Latitude_Peri_Limit = -70; % Latitude limit to consider South Pole flybys [°]
pars.INPUTS.COT.Crank_Direction     = +1;  % Cranking Direction (-1 for negative, +1 for positive)
pars.INPUTS.COT.Pump_Angle          = alfa_res; % Pump angle of the selected resonance for the COT [rad]
pars.INPUTS.COT.Max_Crank_Change    = Delta_Crank_max; %[rad] Max. crank change

% Build-Up the COT
tic
[COT_Data] = COT_BuildUp(pars);
toc

%% POST-PROCESS %%
dim = size(pars.INPUTS.idMoon,1);
Extract_Flybys =  struct('nodein', cell(1,size(dim,1)), ...
    'nodeout', cell(1,size(dim,1)), ...
    'vvinfin', cell(1,size(dim,1)), ...
    'vvinfou', cell(1,size(dim,1)), ...
    'lats', cell(1,size(dim,1)),...
    'longs', cell(1,size(dim,1)),...
    'rp_lat', cell(1,size(dim,1)), ...
    'rp_long', cell(1,size(dim,1)),...
    'fly_States', cell(1,size(dim,1)),...
    'fly_tts', cell(1,size(dim,1)),...
    'altitudes', cell(1,size(dim,1)),...
    'State_In', cell(1,size(dim,1)),...
    'State_Out', cell(1,size(dim,1)),...
    'State_planet', cell(1,size(dim,1)),...
    'OE_In', cell(1,size(dim,1)),...
    'OE_Out',cell(1,size(dim,1)));

% Find flybys mapping the South-Pole
idx_extract = [];
for i = 1:size(COT_Data, 2)
    rp_lat = COT_Data(i).rp_lat;

    if rp_lat <= pars.INPUTS.COT.Latitude_Peri_Limit
        Extract_Flybys(end+1) = COT_Data(i);
        idx_extract(end+1)    = i; 
    end
end

Extract_Flybys(1) = [];

% Plot the flyby groundtracks
colors = cool(length(Extract_Flybys));
fig1 = figure( 'Color', [1 1 1] );
hold on;
plotTextureLatLong(pars.INPUTS.idMoon , pars.INPUTS.idCentral , 1);
plotSquares(pars, 1);  
axis normal;
for i = 1:size(Extract_Flybys, 2)
    Plot_Flyby_GT(Extract_Flybys(i), colors(i,:));
end

% Plot of the SC True Anomaly at Flyby Moment
TA = rad2deg(arrayfun(@(x) x.OE_In(end), Extract_Flybys));
fig2 = figure( 'Color', [1 1 1] );
hold on; grid on;
Leg = linspace(1, size(Extract_Flybys,2), size(Extract_Flybys,2));
plot(Leg, TA, '.', 'Color', 'blue', 'MarkerSize',18 ); 
xlabel('Flyby Number','FontSize',20);
ylabel('SC True Anomaly (degrees)','FontSize',20)

TA_all = rad2deg(arrayfun(@(x) x.OE_In(end), COT_Data));
fig2 = figure( 'Color', [1 1 1] );
hold on; grid on;
Leg_all = linspace(1, size(COT_Data,2), size(COT_Data,2));
plot(Leg_all, TA_all, '.', 'Color', 'blue', 'MarkerSize',18 ); 
plot(Leg_all(idx_extract), TA_all(idx_extract), '.', 'Color', 'red', 'MarkerSize',18 ); 
xlabel('Flyby Number','FontSize',20);
ylabel('SC True Anomaly (degrees)','FontSize',20)
legend('All COT flybys', 'Flybys over the South Pole')

% Plot the flybys around the Moon sphere
fig3 = Plot_Flybys_MoonOrb(Extract_Flybys, pars);

% Plot the planeto-centric orbit of the SC
fig4 = Plot_Flybys_CentralBOrb(Extract_Flybys, pars);

%% COST OF PERFORMING A MANEUVER TO GET BACK TO EQUATORIAL AFTER 5 FLYBYS %%
% Extract outgoing node of the 5th flyby
node_last = Extract_Flybys(1).nodeout;

% Define the inclination of the orbit we want to achieve
incli_target = 0;
[DV_cost] = DV_Incl_Change(node_last, incli_target, pars);

