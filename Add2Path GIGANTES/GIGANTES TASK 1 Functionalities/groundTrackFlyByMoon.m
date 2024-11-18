function [yy, tt, RR, VV, lats, longs, vvinfouBM, delta, lat, long] = ...
    groundTrackFlyByMoon(vvinfin, vvinfou, deltaMax, rr, vv, muPL, pars)

% This function computes the ground-track of the flyby hyperbola.

% Author: A. Bellome
% Last revision: 15/11/2024

%% INPUTS %%
% - vvinfin: Incoming V_infinity vector expressed in Main Body centric
%            coordinates and in [km/s]. Vector of size [1 x 3].
%
% - vvinfou: Outgoing V_infinity vector expressed in Main Body centric
%            coordinates and in [km/s]. Vector of size [1 x 3].
%
% -deltaMax: Maximum bending angle for the incoming infinity velocity,
%            expressed in [rad].
%
% - rr: Flyby body (Moon) position vector expressed in Main Body centric
%       coordinates and in [km]. Vector of size [1 x 3].
%
% - vv: Flyby body (Moon) velocity vector expressed in Main Body centric
%       coordinates and in [km/s]. Vector of size [1 x 3].
%
% - muPL: Gravitational parameter of the flyby body (the Moon) in [km3/s2].
%
% - pars: Structure containing problem parameters.

%% OUTPUTS %%
% - yy: Matrix of size [2*npoints x 6] containing the state vector
%       (position & velocity) resulting from the forward & backward propagation
%       from the flyby periapsis point. Units [km] & [km/s].
%
% - tt: Matrix of size [2*npoints x 1] containing the time steps of the 
%       propagation of the flyby hyperbola.
%
% - RR: Matrix of size [2*npoints x 1] containing the norm of the position 
%       vector of the flyby hyperbola. Units [km].
%
% - VV: Matrix of size [2*npoints x 1] containing the norm of the velocity 
%       vector of the flyby hyperbola. Units [km/s].
%
% - lats: Matrix of size [2*npoints x 1] containing the Latitude of the 
%         ground-track points. Units [rad].
%
% - longs: Matrix of size [2*npoints x 1] containing the Longitude of the 
%          ground-track points. Units [rad].
%
% - vvinfouBM: Outgoing V_infinity vector expressed in Body-Fixed (Moon
%              centric) coordinates and in [km/s]. Vector of size [1 x 3].
%
% - delta: Flyby bending angle in [rad].
%
% - lat: Latitude of the flyby periapsis in [rad].
%
% - long: Longitude of the flyby periapsis in [rad].

%% CHANGES %%
% - 16/05/2024, J.C Garcia Mateas: added input & outputs description; 
%               added "pars" as function input & changed propagation steps.
%
% - 26/08/2024, J.P Sanchez: changed the propagation of the flyby hyperbola
%               to an ODE to overcome errors when doing analytically.

%% FUNCTION %%
if nargin == 6
    npoints = 10e3;
else
    npoints = pars.GroundTr.npoints;
    t_prop  = pars.GroundTr.t_prop; % Time of flyby hyperbola propagation [minutes]
end

% --> outgoing cartesian elements (before the DV-defect)
[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfou, deltaMax, vv);

% --> convert in body-fixed reference frame (this is GTOC6 ref. frame)
b1        = -rr./norm(rr);
b3        = cross(rr, vv)./norm(cross(rr, vv));
b2        = cross(b3,b1);

% --> convert in body-fixed reference frame (this is similar to MIDAS ref.
% frame TCN) this is more similar to Campagnola but still the crank angle k
% is in the opposite direction
% b1 = rr./norm(rr);                             % N
% b2 = vv./norm(vv);                             % T
% b3 = cross( b1, b2 )./norm( cross( b1, b2 ) ); % C

Rm        = [ b1' b2' b3' ]';
vvinfin   = [Rm*vvinfin']';
vvinfouBM = [Rm*vvinfouBM']';

% --> retrieve all the info. on the flyby hyperbola
Energy     = 0.5*norm(vvinfin)^2;
sma        = -muPL/(2*Energy);
ecc        = 1/(sin(delta/2));
rp         = sma*(1 - ecc);
hhat       = cross( vvinfin, vvinfouBM )./norm(cross( vvinfin, vvinfouBM ));
vp         = sqrt( norm(vvinfin)^2 + 2*muPL/rp );
rrp        = rp.*( vvinfin - vvinfouBM )./norm( vvinfin - vvinfouBM );
vvp        = vp.*cross( hhat, rrp./rp );
[tt1, yy1] = propagateKeplerODE(rrp, vvp, linspace(0, -t_prop*60, npoints), muPL);

tt1        = flip(tt1);
yy1        = flip(yy1);

[tt2, yy2] = propagateKeplerODE(rrp, vvp, linspace(0, t_prop*60, npoints), muPL);
yy         = [yy1; yy2];
tt         = [tt1; tt2];
RR         = vecnorm( yy(:,1:3)' )';  % [km] Norm of the SC position vector
VV         = vecnorm( yy(:,4:6)' )';  % [km/s] Norm of the SC velocity vector

% --> only save until it reaches vinf
differenceVV = abs(VV - norm(vvinfin));
indxs        = find( differenceVV >=1e-6 );

if ~isempty(indxs)
    tt = tt(indxs);
    yy = yy(indxs,:);
    VV = VV(indxs);
    RR = RR(indxs);
end

% --> compute the groundtracks
lats         = asin( yy(:,3)./RR );         % [rad] Latitude of the ground-track points
longs         = atan2( yy(:,2), yy(:,1) );   % [rad] Longitude of the ground-track points

% Compute the flyby periapsis point
rp_vec = rrp;
if rp_vec(3) == 0
    lat = 0;
else
    x   = asin(rp_vec(3)/norm(rp_vec)) * 180/pi;  % Latitude
    lat = x;
end

if lat == 90 || lat == -90
    long = 0;
else
    % Adjust the longitude to be within -180 to 180 degrees
    if rp_vec(1)>=0
        y   = atan(rp_vec(2)/rp_vec(1)) * 180/pi;     % Longitude
    elseif rp_vec(1)<0 && rp_vec(2)>=0
        y = atan(rp_vec(2)/rp_vec(1)) * 180/pi + 180;
    else
        y = atan(rp_vec(2)/rp_vec(1)) * 180/pi - 180;
    end
    long = y;
end

% --> make it comparable with Campagnola --> shift to east-longitude
longs   = wrapTo2Pi(longs);
long = wrapTo360(long);

end


%% auxiliary function
function [tt, yy] = propagateKeplerODE(rvec, vvec, timevector, muPL)

%Equation of motion ¨r+mu*r/R^3 = T/m
F=@(t,x)   [x(4); %dx/dt=Vx
    x(5); %dy/dt=Vy
    x(6); %dz/dt=Vz
    -muPL*x(1)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3);
    -muPL*x(2)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3);
    -muPL*x(3)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)];

options=odeset('RelTol',1e-6,'AbsTol',1e-7,'Refine',50);
[tt,yy]=ode45(F,timevector,[rvec vvec],options);

end

%% NOTES
% A function propagateKEPLER is also found within the repository, which
% containts a bug in the newton-raphson for hyperbolic transfers. It is
% currently preferible to use the ODE propagation, but ideally this should
% be substituted soon. 