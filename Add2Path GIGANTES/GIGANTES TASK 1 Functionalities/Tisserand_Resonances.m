function [rp_res,ra_res,E_res, alfa_res, i_max, trA_moonX] = Tisserand_Resonances(Resonances,V_inf, a_moon, mu_planet, pars)

% This function computes the orbital energy and periapsis radius of the
% orbits corresponding to Moon-Spacecraft resonances in the Earth-Moon
% system, so as to plot the corresponding points in the Tisserand graph.

% Author: José Carlos García
% Last revision: 22/06/2022
% Modified: Joan Pau Sánchez (5/7/2024) removed unecessary inputs

%% COMPUTING ORBITAL ENERGY & PERIAPSIS RADIUS OF RESONANT ORBITS %%

N = Resonances(1);   % Moon nº of revolutions                                                   
M = Resonances(2);   % Spacecraft nº of revolutions

T_moon = 2*pi*sqrt(a_moon^3/mu_planet); % [seconds] Period of the Moon

Tsc = (N/M)*T_moon;                             %[seconds] Spacecraft orbital period 

asc = (mu_planet*(Tsc/(2*pi))^2)^(1/3);         %[km] SC semi-major axis 

V_moon = sqrt(mu_planet/a_moon);           %[km/s] Moon Orbital velocity

e_res  = sqrt(1 - (a_moon/asc)*(0.5*(3 - (V_inf/V_moon)^2 - (a_moon/asc)))^2);   %[-] Spacecraft eccentricity
rp_res = asc*(1-e_res)/pars.Planet.EquRad;                                       %[km] Periapsis radius of Spacecraft orbit
ra_res = asc*(1+e_res)/pars.Planet.EquRad;                                       %[km] Apoapsis radius of Spacecraft orbit
E_res  = -mu_planet/(2*asc);                                                     %[km^2/s^2] Orbital energy

a_res    = asc;
alfa_res = acos((V_moon/(2*V_inf))*(1-(V_inf/V_moon)^2 - (a_moon/a_res)));  %[rad] Pump angle associated to the resonance

if isreal(alfa_res)
    i_max    = atan( (V_inf*sin(alfa_res))/(V_moon + V_inf*cos(alfa_res))); %[rad] Maximum inclination that can be achieved with this resonance
else
    i_max = []; rp_res = []; ra_res = []; E_res = []; alfa_res = [];
end

trA_moonX=acos(a_res/e_res*(1-e_res^2)/a_moon-1/e_res); % [rad] True anomaly of intersection


end