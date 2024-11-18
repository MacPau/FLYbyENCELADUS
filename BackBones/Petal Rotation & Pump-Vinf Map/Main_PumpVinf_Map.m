%% CODE TO PLOT THE PUMP - Vinf MAP %%

%% INITIALIZATION %%
clear all; close all; clc; format long g;
% addpath(genpath([pwd '/functions']));

% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [1]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
                            % (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
                            % (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

% Define Flyby Hyperbolic Excess Velocities
pars.INPUTS.V_inf{1} = 3.8:0.01:4.2;     %[km/s] 

% Define the crank angle to consider
pars.INPUTS.crank = 0;

% Define the resonances to consider for the Petal Rotation
% pars.INPUTS.Resonances = [5 1; 11 2; 6 1; 13 2; 19 3; 20 3; 7 1; 22 3; 23 3; 15 2; 8 1;];
pars.INPUTS.Resonances = [5 1; 11 2; 6 1; 13 2; 7 1; 15 2; 8 1];

% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 15e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 10;   %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 
pars.INPUTS.maxDV_Defect   = 0;    %[km/s] Maximum DV defect allowed for a flyby

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
V_infs = cell2mat(pars.INPUTS.V_inf);
for i = 1:length(V_infs)
    e_fly     = 1 + ((rp_flyby*V_infs(i)^2)/pars.Moon.mu); %[-] 
    delta_max(i) = 2*asin(1/e_fly);                        %[rad]
end
pars.delta_max = delta_max;
pars.rp_flyby = rp_flyby;


%% COMPUTE PARAMETERS OF THE RESONANCES %%
for idx = 1:size(pars.INPUTS.idMoon,1)
    [Resos_3D] = Tisserand_Resonances_3D(idx, pars);
end

%% COMPUTE PARAMETERS OF THE PSEUDO-RESONANCES %%
pars.INPUTS.remove81 = 1;

[Pseudo_Resos] = PseudoResonances(pars);

%% PLOT THE PUMP V-INF MAP %%
[PumpVinf_Map] = Plot_PumpVinf_Map(Resos_3D, Pseudo_Resos, pars);


