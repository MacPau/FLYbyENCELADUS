function [Output_PseudoRes] = PseudoResonances(pars)

% This function computes all the possible pseudo-resonant transfers for a
% given moon, V_inf and set of resonance ratios.

% Author: J.C Garcia Mateas
% Last revision: 22/05/2024

%% INPUTS %%
% - pars: Structure containing the problem constant & parameters

%% OUTPUTS %%
% - Output_PseudoRes: Structure with 15 fields and size N x 1 (where N = 
%                     n° of resonance ratios being considered x 2)

%% FUNCTION %%
% Retrieve Vinf levels
V_inf = cell2mat(pars.INPUTS.V_inf); %[km/s] Hyperbolic excess velocity

% Define type of Pseudo-Resonant Transfers
type  = [18, 81]; % [18] = Inbound-/Outbound; [81] = Outbound/Inbound

for i = 1:size(pars.INPUTS.Resonances,1)
    
    N = pars.INPUTS.Resonances(i,1); M = pars.INPUTS.Resonances(i,2);

    for m = 1:size(V_inf,2)

        for j = 1:size(type, 2)
            
            % Solve the pseudo-resonant transfer
            [vinf1, alpha1, crank1, vinf2, alpha2, crank2, tof1] = pseudoResTransf(type(j), N, M, V_inf(m), pars.INPUTS.idMoon, pars.INPUTS.idCentral, pars.INPUTS.remove81);    
            nodein  = [vinf1, alpha1, crank1]; 
            nodeout = [vinf2, alpha2, crank2];
    
            if ~isnan(nodein)
                % Extract maximum bending for this Vinf
                delta_max = pars.delta_max;        % Original values
                delta_max_aux = pars.delta_max(m); % Bending for this Vinf
                pars.delta_max = delta_max_aux;

                % Compute the parameters associated to the flyby
                Flyby_Data = Flyby_BuildUp(nodein, nodeout, pars);

                pars.delta_max = delta_max; % Rewrite original values
        
                % Add pseudo-resonant transfer parameters to the structure
                Flyby_Data.TOF  = tof1;    %[seconds] Time of Flight
                Flyby_Data.Reso = [N, M];  % Resonance Ratio
                Flyby_Data.Type = type(j); % Transfer type
        
                % Store parameters associated to the pseudo-resonant transfer & flyby
                Output(i,m,j) = Flyby_Data;
            end
        end
    end
end

Output_PseudoRes = Output;

end