function [OE_f] = Cartesian_to_OE(state_vector_f, mu)

% This function has been developed in order to convert vectors expressed in
% Cartesian coordinates to Keplerian orbital elements.

% Inputs:   - 1x6 vector containing ECI components for position and
%             velocity [X Y Z Vx Vy Vz] in km and km/s
%           - mu, the gravitational parameter in [km^3/s2]

%Outputs:   - 1x6 vector containing the orbital elements, arranged as 
%             [a e i w raan theta]

% The function has been verified by comparing the outputs with those
% retrieved using the similar function from the CURTIS book MATLAB files &
% with the function car2kep of the ASTRA Toolbox.

% Author: José Carlos García
% Last revision: 31/07/2022

%% FUNCTION CODE %%
position = [state_vector_f(1) state_vector_f(2) state_vector_f(3)]; %[km]
velocity = [state_vector_f(4) state_vector_f(5) state_vector_f(6)]; %[km/s]


%STEP 1: Calculating the distance and the speed
r = norm(position);
v = norm(velocity);

%STEP 2: Calculating radial velocity and specific angular momentum
v_r = (position(1)*velocity(1) + position(2)*velocity(2) + position(3)*velocity(3))/r; %radial velocity [km/s]
h = cross(position,velocity); %specific ang. momentum
h_mag = norm(h); %magnitude of the specific ang. momentum

%STEP 3: Calculating the inclination
i = acosd(h(3)/h_mag);

%STEP 4: Define the node line
K_vector = [0 0 1];
N_vector = cross(K_vector,h);
N_mag = norm(N_vector);

%STEP 5: Calculate the right ascension of the ascending node
if N_vector(2) >= 0
    raan = acosd(N_vector(1)/N_mag);
else
    raan = 360 - acosd(N_vector(1)/N_mag);
end

%STEP 6: Calculate the eccentricity
e_vector = (1/mu)*((v^2-mu/r)*position - r*v_r*velocity); %Vector
e = sqrt(1+(h_mag^2/mu^2)*(v^2 - ((2*mu)/r))); %Magnitude [-]

%STEP 7: Calculating the argument of the perigee
if e_vector(3) >= 0
    w = acosd(dot(N_vector,e_vector)/(N_mag*e));
else
    w = 360 - acosd(dot(N_vector,e_vector)/(N_mag*e));
end

%STEP 8: Calculating the true anomaly
if v_r >= 0
    theta = acosd(dot(e_vector/e,position/r));
else
    theta = 360 - acosd(dot(e_vector/e,position/r));
end

%STEP 9: Calculating the semi-major axis
r_p = (h_mag^2/mu)*(1/(1+e*cosd(0))); %[km]
r_a = (h_mag^2/mu)*(1/(1+e*cosd(180))); %[km]

a = 0.5*(r_p + r_a);

%Creating vector with orbital elements
OE_f = [a e i w raan theta]; %[- / km / degrees]

end








