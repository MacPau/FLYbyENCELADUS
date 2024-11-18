function [vpip, eip, Eip, aip, delta, thinf, tpip, thsoi] = Vinf2HyperbolaSOI(vinf, rpip, mu, rsoi)

% computes hyberbola geometry from vinf and periapsis radius, w.r.t. the
% central body (mu). If requested in input, the time of flight towards the
% planet sphere of influence is also computed. 

Eip   = 0.5*vinf^2;
aip   = -mu/(2*Eip);
vpip  = sqrt(2*(Eip + mu/rpip));
eip   = 1 + rpip*vinf^2/mu;
delta = 2*asin(1/eip);
mu_rp = mu/rpip;
thinf = acos(-mu_rp/(norm(vinf)^2+mu_rp));

% time from pericentre to SOI (tpip)
if nargin > 3
    thsoi = wrapToPi(acos((aip*(1-eip^2) - rsoi)/(eip*rsoi)));
    M0    = 0;
    M1    = theta2M(thsoi, eip);
    nh    = sqrt(mu/(-aip^3));
    tpip  = ((M1 - M0))/nh; % time to reach the pericentre of incoming hyperbola (seconds)
else
    tpip  = [];
    thsoi = [];
end
    
end