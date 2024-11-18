function [PseudoReso_Transfer, PseudoResoStruct] = PseudoResoTransfers_Generation(V_inf, pars)

% This function computes and stores all feasible pseudoresonances given a
% user-defined hyperbolic excess velocity, a set of resonance ratios and a
% body of interest.

% Author: J.C. Garcia Mateas
% Last revision: 28/05/2024

%% INPUTS %%
% - V_inf: Hyperbolic excess velocity of the flyby in [km/s]
%
% - pars: Structure containing different constants & parameters

%% OUTPUTS %%
% - PseudoReso_Transfer: Matrix of size [N x 10], where each row is a 
%                        pseudo-resonant transfer, with  column 1 being 
%                        the n° of Gravity Assist Body revolutions,
%                        column 2 the n° of spacecraft revolutions, col 3 the
%                        V_inf [km/s], cols 4  & 5 the Pump and Crank angle 
%                        respectively in [rad] of the incoming node, columns 
%                        6, 7, 8 the Vinf, pump & crank of the outgoing node,
%                        column 9 the TOF in [seconds] and column 10 the type
%                        of transfer [18] = Inbound-/Outbound; 
%                        [81] = Outbound/Inbound.
% 
% - PseudoResoStruct: Structure containing the same information as the
%                     previous output but in another form.

%% CHANGES %%
% - 28/05/2024 - J.C. Garcia Mateas: Rearranged outputs to make them compatible
%                with those of the ResoTransfers_Generation function.

%% FUNCTION %%
% Initialize Pseudo-Resonances Matrix & Structure
PseudoReso_Transfer = [];

% Define type of Pseudo-Resonant Transfers
type  = [18, 81]; % [18] = Inbound-/Outbound; [81] = Outbound/Inbound

for i = 1:size(pars.INPUTS.Resonances,1)
    
    % Define n° of Moon revolutions 
    N = pars.INPUTS.Resonances(i,1); 

    % Define n° of Spacecraft revolutions
    M = pars.INPUTS.Resonances(i,2);
    
    % Compute long (+) and short (-) pseudo-resonant transfers
    for j = 1:size(type, 2)

        % Solve the pseudo-resonant transfer
        [vinf1, alpha1, crank1, vinf2, alpha2, crank2, TOF] = pseudoResTransf(type(j), N, M, V_inf, pars.INPUTS.idMoon, pars.INPUTS.idCentral, pars.INPUTS.remove81);

        if ~isnan(vinf1)
            % Store the results
            PseudoReso_Transfer = [PseudoReso_Transfer; [N, M, vinf1, alpha1, crank1, TOF, vinf2, alpha2, crank2, type(j)]];
        end
    end
end

% Assigning values to Output
PseudoResoStruct = struct('N', cell(size(PseudoReso_Transfer,1), 1), 'M', cell(size(PseudoReso_Transfer,1), 1),...
    'V_inf1', cell(size(PseudoReso_Transfer,1), 1), 'Alfa1', cell(size(PseudoReso_Transfer,1), 1), 'Crank1', cell(size(PseudoReso_Transfer,1), 1),...
    'TOF_SC', cell(size(PseudoReso_Transfer,1), 1),'V_inf2', cell(size(PseudoReso_Transfer,1), 1), 'Alfa2', cell(size(PseudoReso_Transfer,1), 1),...
     'Crank2', cell(size(PseudoReso_Transfer,1), 1),'Type', cell(size(PseudoReso_Transfer,1), 1));

for i = 1:size(PseudoReso_Transfer,1)
    PseudoResoStruct(i).N = PseudoReso_Transfer(i,1);       % Number of Gravity Assist Body revolutions
    PseudoResoStruct(i).M = PseudoReso_Transfer(i,2);       % Number of Spacecraft revolutions
    PseudoResoStruct(i).V_inf1 = PseudoReso_Transfer(i,3);  % Hyperbolic Excess velocity [km/s]
    PseudoResoStruct(i).Alfa1 = PseudoReso_Transfer(i,4);   % Pump Angle [rad]
    PseudoResoStruct(i).Crank1= PseudoReso_Transfer(i,5);   % Crank Angle [rad]
    PseudoResoStruct(i).TOF_SC = PseudoReso_Transfer(i,6);  % SC Time of Flight [seconds]
    PseudoResoStruct(i).V_inf2 = PseudoReso_Transfer(i,7);  % Hyperbolic Excess velocity [km/s]
    PseudoResoStruct(i).Alfa2 = PseudoReso_Transfer(i,8);   % Pump Angle [rad]
    PseudoResoStruct(i).Crank2 = PseudoReso_Transfer(i,9);  % Crank Angle [rad]
    PseudoResoStruct(i).Type = PseudoReso_Transfer(i,10);   % Type of transfer [18] = Inbound-/Outbound; [81] = Outbound/Inbound
end

end