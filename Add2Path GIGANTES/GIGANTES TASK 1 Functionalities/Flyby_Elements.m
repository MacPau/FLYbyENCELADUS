function [Output] = Flyby_Elements(nodein, nodeout, pars)

% This function builds-up a flyby, determining the properties of the flyby
% and the ground-track associated to it.

% Author: J.C Garcia Mateas
% Last revision: 04/06/2024

%% INPUTS %%
% - nodein: Vector of size (1x3) containing the Vinf, pump and crank angles
%           of the incoming node, expressed in [km/s, rad, rad].
%
% - nodeout: Vector of size (1x3) containing the Vinf, pump and crank angles
%           of the outgoing node, expressed in [km/s, rad, rad].
%
% - pars: structure containing problem parameters & constants.

%% OUTPUTS %%
% - Output: Structure with 6 fields containing saturn centric information
%           computed for the flyby.

%% CHANGE LOG %%
% 04/06/2024, J.C Garcia Mateas: added new output fields to the structure,
%             being the DV cost of the flyby (for DV defects), the True
%             Longitude of the SC and the Change in true longitude for the
%             computation of the petal rotations.


%% FUNCTION %%
% Initialize the output
dim = size(pars.INPUTS.idMoon,1);
Output = struct('nodein', cell(1,size(dim,1)), ...
    'nodeout', cell(1,size(dim,1)), ...
    'State_In', cell(1,size(dim,1)),...
    'State_Out', cell(1,size(dim,1)),...
    'OE_In', cell(1,size(dim,1)),...
    'OE_Out', cell(1,size(dim,1)),...
    'DV_Cost', cell(1,size(dim,1)),...
    'TrueLongitudes', cell(1,size(dim,1)),...
    'Change_TrueLong',cell(1,size(dim,1)));

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

% Find Saturn centric Orbital Elements of Incoming & Outgoing
OE_in  = car2kep([rr_in, vv_in], pars.Planet.mu);
OE_out = car2kep([rr_ou, vv_ou], pars.Planet.mu);

% OE_in  = Cartesian_to_OE([rr_in, vv_in], pars.Planet.mu);
% OE_out = Cartesian_to_OE([rr_ou, vv_ou], pars.Planet.mu);

% Store Flyby outputs
Output.nodein    = nodein;
Output.nodeout   = nodeout;
Output.State_In  = [rr_in, vv_in];
Output.State_Out = [rr_ou, vv_ou];
Output.OE_In     = OE_in;
Output.OE_Out    = OE_out;
Output.DV_Cost   = [];
Output.TrueLongitudes = [];
Output.Change_TrueLong = [];

end