function [Output] = COT_BuildUp(pars)

% This function builds up the Crank over the Top sequence for a given moon,
% hyperbolic excess velocity and resonance ratio.

% Author: Jose Carlos Garcia Mateas
% Last revision: 23/08/2024

%% INPUTS %%
% - pars: Structure containing parameters defining the problem.

%% OUTPUTS %%
% - Output: Structure containing the flyby parameters defined as output of
%           the function.

%% CHANGELOG %%
% - 23/08/2024, J.C Garcia Mateas: updated the fields in "Output" to
%               include new outputs from the Flyby_BuildUp function.

%% FUNCTION %%
% Retrive Inputs
max_DCrank = pars.INPUTS.COT.Max_Crank_Change;     % [rad] Maximum change in crank angle per flyby
kin        = deg2rad(pars.INPUTS.COT.Start_Crank); % [°] COT Starting crank angle
signDk     = pars.INPUTS.COT.Crank_Direction;      % Cranking Direction (-1 for negative, +1 for positive)
alfain     = pars.INPUTS.COT.Pump_Angle;           % [rad] Pump Angle of the COT (considered constant since only cranking)
vinfin     = pars.INPUTS.V_inf;                    % [km/s] Hyperbolic excess velocity

% Initialize function Output
dim = size(pars.INPUTS.idMoon,1);
Output = struct('nodein', cell(1,size(dim,1)), ...
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

% Determine Crank Bounds 
if kin == 0
    if signDk == +1
        kmax = pi;
    elseif signDk == -1
        kmax = -pi;
    end
elseif kin == pi
    if signDk == +1
        kmax = 2*pi;
    elseif signDk == -1
        kmax = 0;
    end
end

k = unique(kin:signDk*max_DCrank:kmax, 'rows', 'stable');
k = linspace(kin, kmax, length(k)+1);

% Loop through all the flybys in the COT
for indk = 1:length(k)-1

    kin = k(indk);
    kou = k(indk+1);

    nodein  = [vinfin, alfain, kin];
    nodeout = [vinfin, alfain, kou];

    % Compute flyby parameters
    [Flyby] = Flyby_BuildUp(nodein, nodeout, pars);

    % Store the results
    Output(indk)  = Flyby;

end

end