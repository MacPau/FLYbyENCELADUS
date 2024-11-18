function [vinfCAR, rr, vv] = vinfAlphaCrank2car(NODE, SVMoon, muPlanet)

% DESCRIPTION
% This function computes the v-infinity vector in cartesian coordinates
% given a set of v-infinity magnitude, pump and crank angles.
%
% INPUT
% - NODE      : (v-infinity magnitude [L/T], pump angle [rad], crank angle [rad])
% - SVMoon    : Full state vector of the GA body [L, L/T]
% - muPlanet  : Gravitational constant of the planet [L^3/T^2]
%
% OUTPUT
% - vinfCAR : 1x3 vector with v-infinity in cartesian coordinates (vinfx,
%             vinfy, vinfx) [L/T]
% - rr      : 1x3 vector of spacecraft position at epoch [L]
% - vv      : 1x3 vector of spacecraft velocity at epoch [L/S]
%
%--------------------------------------------------------------------------
% unwrap the node
vinf_norm=NODE(1); % Hyperbolic excess velocity at ga
alpha=NODE(2); % Pump angle [RAD]
k=NODE(3); % Crank angle [RAD]

%--------------------------------------------------------------------------
% unwrap the ga object state vector
rrga=SVMoon(1:3); 
vvga=SVMoon(4:6); 

kepga=car2kep(SVMoon,muPlanet);

vinfTCN = vinf_norm.*[ cos(alpha), -sin(alpha)*sin(k), sin(alpha)*cos(k) ];

thga    = kepga(end);
gamma   = pi/2 - sign( rrga(1)*vvga(2) - rrga(2)*vvga(1) )*acos(dot( rrga./norm(rrga), vvga./norm(vvga) )); 

s = sin( thga - gamma );
c = cos( thga - gamma );

rot_mat = [ -s, 0, c; c, 0, s; 0, 1, 0 ];
vinfCAR = (rot_mat*vinfTCN')';

rr = rrga;
vv = vvga + vinfCAR;

end