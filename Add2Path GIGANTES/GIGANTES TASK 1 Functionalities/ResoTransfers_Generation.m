function [Resonances, ResonanceStruct] = ResoTransfers_Generation(V_inf, pars)

% This function computes and stores all feasible resonances given a
% user-defined hyperbolic excess velocity, a set of resonance ratios and a
% body of interest.

% Author: J.C. Garcia Mateas
% Last revision: 10/06/2024

%% INPUTS %%
% - V_inf: Hyperbolic excess velocity of the flyby in [km/s]
%
% - pars: Structure containing different constants & parameters

%% OUTPUTS %%
% - Resonances: Matrix of size [N x 10], where each row is a resonance, and
%               the columns have information regarding the resonant
%               transfer as follows: col1 = N° moon revolutions; col2 = n°
%               spacecraft revs; col3 = Vinf [km/s], col4 = pump angle
%               [rad], col5 = crank angle [rad], col6 = TOF [seconds], col
%               7 = col3; col8 = col4; col9 = col6; col10 = transfer type 
%               where 11 = Inbound/Inbound and 88 = Outbound/Outbound.
%               The first node (col 3-5) is the node departing the moon,
%               while the second node (colums 7 to 9) is the node that
%               arrives again at the moon after "flying" for the TOF.
% 
% - ResonanceStruct: Structure containing the same information as the
%                    previous output "Resonances" but in structure form.

%% CHANGES %%
% 11/03/2024 - J.C. Garcia Mateas: Added the pump angle as the 4th column in the output
%              variable "Resonances" and the spacecraft TOF in the 6th column.
%
% 22/03/2024 - J.C. Garcia Mateas: Added the possibility of computing the
%              nodes from a list of user-defined resonance ratios.
%
% 28/05/2024 - J.C. Garcia Mateas: Changed the output format of Resonances
%              & ResonanceStruct to include the columns regarding the
%              arrival node (node2) to the next Moon encounter after the
%              TOF "propagation".
%
% 10/06/2024 - J.C Garcia Mateas: revised the function, added comments on
%              inputs/outputs.

%% FUNCTION %%
% Initialize Resonances Matrix & Structure
Resonances = [];

% Define Moon Circular Orbital Velocity
V_moon = pars.Moon.Vel;

% Generate all non-repeated combinations (both increasing and decreasing)
if isempty(pars.INPUTS.Resonances)
    for i = pars.INPUTS.Min_Reso:pars.INPUTS.Max_Reso
        for j = pars.INPUTS.Min_Reso:pars.INPUTS.Max_Reso
            if i ~= j
                % Flyby Pump angle
                alfa = acos((V_moon^2*(2 - (i/j)^(2/3) ) - V_inf^2 - V_moon^2)/(2*V_inf*V_moon)); %[rad]
    
                % Spacecraft TOF
                TOF = (j/i)*pars.Moon.Period; %[seconds]
    
                if isreal(alfa)
                    for k = 1:size(pars.INPUTS.Crank,2)
                        % Introduce the Crank Discretization
                        Resonances = [Resonances; [i, j, V_inf, alfa, pars.INPUTS.crank(k), TOF, V_inf, alfa, pars.INPUTS.crank(k), 0]];
                    end
                end
            end
        end
    end

else % If a resonance list has been introduced

    for i = 1:size(pars.INPUTS.Resonances, 1)

        % Define n° of Moon revolutions
        N = pars.INPUTS.Resonances(i,1); 

        % Define n° of Spacecraft revolutions
        M = pars.INPUTS.Resonances(i,2);

        % Flyby Pump angle 
        alfa = acos( (V_moon^2*(2 - (M/N)^(2/3) ) - V_inf^2 - V_moon^2)/(2*V_inf*V_moon) ); %[rad]

        % Spacecraft TOF
        TOF = (N/M)*pars.Moon.Period; %[seconds]

        if isreal(alfa)
            % Introduce the Crank Discretization
            for k = 1:size(pars.INPUTS.crank,2)

                % If resonant transfer is Outbound/Outbound
                if pars.INPUTS.crank(k) >= 0 && pars.INPUTS.crank(k) < pi || pars.INPUTS.crank(k) == 2*pi   
                    Resonances = [Resonances; [N, M, V_inf, alfa, pars.INPUTS.crank(k), TOF, V_inf, alfa, pars.INPUTS.crank(k), 88]];

                % If resonant transfer is Inbound/Inbound
                elseif pars.INPUTS.crank(k) >= pi && pars.INPUTS.crank(k) < 2*pi                           
                    Resonances = [Resonances; [N, M, V_inf, alfa, pars.INPUTS.crank(k), TOF, V_inf, alfa, pars.INPUTS.crank(k), 11]];
                end

            end
        end

    end

end

% Assigning values to Output
ResonanceStruct = struct('N', cell(size(Resonances,1), 1), 'M', cell(size(Resonances,1), 1),...
    'V_inf1', cell(size(Resonances,1), 1), 'Alfa1', cell(size(Resonances,1), 1), 'Crank1', cell(size(Resonances,1), 1), 'TOF_SC', cell(size(Resonances,1), 1),...
     'V_inf2', cell(size(Resonances,1), 1), 'Alfa2', cell(size(Resonances,1), 1), 'Crank2', cell(size(Resonances,1), 1), 'Type', cell(size(Resonances,1), 1));

for i = 1:size(Resonances,1)
    ResonanceStruct(i).N = Resonances(i,1);        % Number of Gravity Assist Body revolutions
    ResonanceStruct(i).M = Resonances(i,2);        % Number of Spacecraft revolutions
    ResonanceStruct(i).V_inf1 = Resonances(i,3);   % Hyperbolic Excess velocity [km/s]
    ResonanceStruct(i).Alfa1 = Resonances(i,4);    % Pump Angle [rad]
    ResonanceStruct(i).Crank1 = Resonances(i,5);   % Crank Angle [rad]
    ResonanceStruct(i).TOF_SC = Resonances(i,6);   % SC Time of Flight [seconds]
    ResonanceStruct(i).V_inf2 = Resonances(i,7);   % Hyperbolic Excess velocity [km/s]
    ResonanceStruct(i).Alfa2 = Resonances(i,8);    % Pump Angle [rad]
    ResonanceStruct(i).Crank2 = Resonances(i,9);   % Crank Angle [rad]
    ResonanceStruct(i).Type = Resonances(i,10);    % Transfer Type [0 = unkown]
end

end