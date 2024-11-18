function [Resos_Output] = Tisserand_Resonances_3D(idx, pars)

% This function computes the variables associated to a certain resonance in
% order to plot the 3D Tisserand Graph.

% Author: J.C. Garcia Mateas
% Last revision: 17/04/2024

%% INPUTS %%
%
% idx: Index corresponding to the gravity assist body being considered from
%      the ones defined in pars.INPUTS.idMoon.
%
% pars: Structure containing problem data and parameters.

%% OUTPUTS %%
%
% Resos_Output: Structure containing the output parameters associated to
%               the resonances. Size is [1 x j], where j is the number of
%               resonances being computed, and it has 7 fields. Each field
%               is a 3D matrix [j,i,m], where "j" is for each resonance, 
%               "i" is for each V_infinity and "m" for each crank angle.

%% FUNCTION %%
% Define some parameters
V_c   = pars.Moon.Vel(idx);         %[km/s] Orbital Velocity of the GA body 
r_enc = pars.Moon.OrbRad(idx);      %[km] Radius of encounter @ flyby

V_inf = cell2mat(pars.INPUTS.V_inf(idx)); %[km/s] Hyperbolic excess velocity

% Loop through each resonance
for j = 1:size(pars.INPUTS.Resonances,1)

    N = pars.INPUTS.Resonances(j,1);   % Moon nº of revolutions                                                   
    M = pars.INPUTS.Resonances(j,2);   % Spacecraft nº of revolutions
    
    Tsc = (N/M)*pars.Moon.Period(idx);               %[seconds] Spacecraft orbital period 
    
    a_sc = (pars.Planet.mu*(Tsc/(2*pi))^2)^(1/3);   %[km] SC/Resonance semi-major axis 

    for i = 1:size(V_inf,2)
    
        alfa_r(j,i) = acos((V_c/(2*V_inf(i)))*(1-(V_inf(i)/V_c)^2 - (r_enc/a_sc)));  % [rad] Pump angle associated to the resonance

        % If the resonance exists
        if isreal(alfa_r(j,i))          

            % Store the pump angle associated to the resonance
            alfa_res(j,i) = alfa_r(j,i);
            
            % Loop through each crank angle
            for m = 1:length(pars.INPUTS.crank) 

                % Store crank for the function Output
                crank_res(j,i,m) = pars.INPUTS.crank(m);

                % Compute spacecraft inclination
                i_res(j,i,m) = atan(sin(pars.INPUTS.crank(m))*(sin(alfa_r(j,i))/(V_c/V_inf(i) + cos(alfa_r(j,i)))));  %[rad] Inclination
                
                % Compute spacecraft eccentricity
                e_res(j,i,m) = sqrt( 1 - (r_enc/a_sc)*( (3 - (V_inf(i)/V_c)^2  - (r_enc/a_sc))/(2*cos(i_res(j,i,m))))^2); %[-] Eccentricity
                
                % Compute orbit periapsis & apoapsis
                r_p_res(j,i,m) = a_sc*(1-e_res(j,i,m));       %[km] Perigee radius 
                r_a_res(j,i,m) = a_sc*(1+e_res(j,i,m));       %[km] Apogee radius

                % Compute argument of periapsis
                w_res(j,i,m) = acos((1/e_res(j,i,m))*((a_sc/r_enc)*(1 - e_res(j,i,m)^2) - 1));  %[rad] Arg. Periapsis

                Energy_res(j,i,m) = -pars.Planet.mu/(2*a_sc);  %[km^2/s^2]
            end
        end
    end

end


% Create Output structure
Resos_Output = struct('r_p', cell(1,size(V_c,1)), ...
    'r_a', cell(1,size(V_c,1)), ...
    'Energy', cell(1,size(V_c,1)),...
    'e_sc', cell(1,size(V_c,1)),...
    'i_sc', cell(1,size(V_c,1)),...
    'w', cell(1,size(V_c,1)),...
    'alfa', cell(1,size(V_c,1)),...
    'crank',cell(1,size(V_c,1)));

Resos_Output.r_p = r_p_res;
Resos_Output.r_a = r_a_res;
Resos_Output.Energy = Energy_res;
Resos_Output.e_sc = e_res;
Resos_Output.i_sc = i_res;
Resos_Output.w = w_res;
Resos_Output.alfa = alfa_res;
Resos_Output.crank = crank_res;


end
