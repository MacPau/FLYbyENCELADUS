function [DV_cost] = DV_Incl_Change(node_last, incli_target, pars)

% This function computes the cost of performing an inclination change for
% an elliptical orbit. This maneuver is carried out at the opposite node.

%% INPUTS %%

%% OUTPUTS %%

%% FUNCTION %%
% Extract the Vinf, pump & crank from the input node
vinf = node_last(1); alpha = node_last(2); k = node_last(3);

% Compute Saturn Centric State vector
[vvinfin, rr_sc, vv_sc, ~] = vinfAlphaCrank_to_VinfCART(vinf, alpha, k, pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);

% Compute associated Orbital Elements
OE_sc = car2kep([rr_sc, vv_sc], pars.Planet.mu);
incli = OE_sc(3); %[rad]

% "Propagate" to the next node (change true anomaly by 180°)
OE_sc_nextNode = OE_sc;
OE_sc_nextNode(end) = OE_sc(end) + pi; %[rad]

% Compute Saturn Centric State Vector at the next node
state_nextNode = kep2car(OE_sc_nextNode, pars.Planet.mu); %[km] & [km/s]
vv_sc_nextNode = state_nextNode(4:6); %[km/s] (small because we are near apoapsis)

% Compute the SC flight path angle 
A = OE_sc_nextNode(2)*sin(OE_sc_nextNode(end));
B = 1 + OE_sc_nextNode(2)*cos(OE_sc_nextNode(end));
gamma = atan(A/B);

% other option is Eqt. 2.56 of Strange PhD thesis, sign of gamma would be the opposite @ opposite node)
% gamma = asin( (vinf/norm(vv_sc_nextNode))*sin(alpha)*cos(k) ); %[rad]

% Determine the inclination change
Incli_change = abs(incli - incli_target); %[rad]

% Compute the DV required 
DV_cost = 2*norm(vv_sc_nextNode)*cos(gamma)*sin(Incli_change/2);  %[km/s]



end