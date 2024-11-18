function [Pos_ECI, Vel_ECI] = OE_to_Cartesian(OE_0,mu)

% This function has been developed with the objective of converting 
% Keplerian elements into a state vector in Earth-Centered Inertial 
%(ECI) coordinates. 

%Inputs:     - 1x6 vector containing orbital elements 
%            - mu , the gravitational parameter in [km^3/s2]

%Outputs:    - Pos_ECI: 3x1 vector with position components X,Y,Z in [km]
%            - Vel_ECI:3x1 vector with velocity components X,Y,Z in [km/s]

%% STEP 0: Identification of orbital elements
a = OE_0(1); %Semi-major axis [km]
e = OE_0(2); %Eccentricity [-]
i = OE_0(3); %Inclination [rad]
w = OE_0(4); %Argument of the perigee [rad]
raan = OE_0(5); %Righ ascension of the ascending node [rad]
theta = OE_0(6); %True Anomaly [rad]

%% STEP 1: Position and velocity vectors in the PERIFOCAL FRAME %%
%Position vector
Rx_pf = ((a*(1-e^2))/(1+e*cos(theta)))*cos(theta);
Ry_pf = ((a*(1-e^2))/(1+e*cos(theta)))*sin(theta);

Pos_pf = [Rx_pf; Ry_pf; 0]; %[km]

%Velocity vector
Vx_pf = sqrt(mu/(a*(1-e^2)))*(-sin(theta));
Vy_pf = sqrt(mu/(a*(1-e^2)))*(e+cos(theta));

Vel_pf = [Vx_pf; Vy_pf; 0]; %[km/s]

%% STEP 2: Rotation from PERIFOCAL FRAME to ECI %%
Matrix = [-sin(raan)*cos(i)*sin(w)+cos(raan)*cos(w) -sin(raan)*cos(i)*cos(w)-cos(raan)*sin(w) sin(raan)*sin(i);...
    cos(raan)*cos(i)*sin(w)+sin(raan)*cos(w) cos(raan)*cos(i)*cos(w)-sin(raan)*sin(w) -cos(raan)*sin(i);...
    sin(i)*sin(w) sin(i)*cos(w) cos(i)];

Pos_ECI = Matrix*Pos_pf; %[km]
Vel_ECI = Matrix*Vel_pf; %[km/s]

end 
