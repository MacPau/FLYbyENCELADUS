function [Output] = Flyby_BuildUp(nodein, nodeout, pars)

% This function builds-up a flyby, determining the properties of the flyby
% and the ground-track associated to it.

% Author: J.C Garcia Mateas
% Last revision: 23/08/2024

%% INPUTS %%
% - nodein: Vector of size (1x3) containing the Vinf, pump and crank angles
%           of the incoming node, expressed in [km/s, rad, rad].
%
% - nodeout: Vector of size (1x3) containing the Vinf, pump and crank angles
%           of the outgoing node, expressed in [km/s, rad, rad].
%
% - pars: structure containing problem parameters & constants.

%% OUTPUTS %%
% - Output: Structure with 12 fields containing all the information
%           computed for the flyby.

%% CHANGE LOG %%
% - 22/05/2024 - J.C Garcia Mateas: switched to using the function car2kep
%                to compute the SC orbital elements to avoid problems of
%                singularities.
%
% - 23/08/2024 - J. P Sanchez: included more fields for more data from the
%                flyby in the function output.

%% FUNCTION %%
% Initialize the output
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

% Retrieve initial epoch
epoch0    = pars.INPUTS.epoch0;

% Incoming node
vinfin  = nodein(1);
alphain = nodein(2);
kin     = nodein(3);

% Outgoing node (post-flyby)
vinfou  = nodeout(1);
alphaou = nodeout(2);
kou     = nodeout(3);

% Converting Incoming Node into Planeto Centric Cartesian state vector in TCN frame
[vvinfin, rr_in, vv_in, ~] = vinfAlphaCrank_to_VinfCART(vinfin, alphain, kin, epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);

% Converting Outgoing Node into Planeto-Centric Cartesian state vector in TCN frame
[vvinfou, rr_ou, vv_ou, vvga]  = vinfAlphaCrank_to_VinfCART(vinfou, alphaou, kou, epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);

% Compute the flyby ground-track parameters
[states, tt4states, Alts, ~, lats, longs, ~, ~, rp_lat, rp_long] = groundTrackFlyByMoon(vvinfin, vvinfou, pars.delta_max, rr_ou, vvga, pars.Moon.mu, pars);

% Find Saturn centric Orbital Elements of Incoming & Outgoing
OE_in  = car2kep([rr_in, vv_in], pars.Planet.mu);
OE_out = car2kep([rr_ou, vv_ou], pars.Planet.mu);

% OE_in  = Cartesian_to_OE([rr_in, vv_in], pars.Planet.mu);
% OE_out = Cartesian_to_OE([rr_ou, vv_ou], pars.Planet.mu);

% Store Flyby outputs
Output.nodein    = nodein;
Output.nodeout   = nodeout;
Output.vvinfin    = vvinfin;
Output.vvinfou   = vvinfou;
Output.lats      = lats';
Output.longs     = longs';
Output.rp_lat    = rp_lat;
Output.rp_long   = rp_long;
Output.fly_States = states;
Output.fly_tts = tt4states;
Output.altitudes = Alts';
Output.State_In  = [rr_in, vv_in];
Output.State_Out = [rr_ou, vv_ou];
Output.State_planet  = [rr_in, vvga];
Output.OE_In     = OE_in;
Output.OE_Out    = OE_out;

end