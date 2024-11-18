function [PumpVinf_Map] = Plot_PumpVinf_Map(Resos_3D, Pseudo_Resos, pars)

% This function plots the Pump - Vinf Map for a set of resonances and
% pseudo-resonances.

% Author: J.C Garcia Mateas
% Last revision: 24/05/2024

%% INPUTS %%
% - Resos_3D: Structure with 7 fields containing the data regarding
%             resonances at a Moon. The fields will be matrices of size A x
%             B, where A = n° of resonance ratios & B = n° of Vinfs at
%             which the values have been computed.
%
% - Pseudo_Resos: Structure of size I x M x J, where I = n° of resonance
%             ratios considered, M = n° of discretized Vinfs & J = 2 (types
%             of pseudo-resonant transfer IO / OI). 
%
% - pars: Structure with 15 fields containing the data regarding
%         pseudo-resonances with a Moon. The dimensions will be A x B x 2,
%         where the 2 comes from the two types of pseudoresonances that
%         exist, ie, long transfer (+) and short (-).

%% OUTPUTS %%
% - PumpVinf_Map: Figure of the Pump-Vinf map.

%% FUNCTION %%

V_infs = cell2mat(pars.INPUTS.V_inf);
Pumps_Resos = Resos_3D.alfa;

% Define colormap for each orbit
colors = jet(length(pars.INPUTS.Resonances));

PumpVinf_Map = figure('Color', [1 1 1]);clf;set(PumpVinf_Map,'defaulttextinterpreter','latex') ;
hold on; grid on;
for i = 1:length(pars.INPUTS.Resonances)
    % Plot the Resonances
    Pumps_aux = Pumps_Resos(i, :);
    Pump_Resos_plot = Pumps_aux(Pumps_aux ~= 0);
    V_infs_plot = V_infs(Pumps_aux ~= 0);
     plot(V_infs_plot, rad2deg(Pump_Resos_plot), 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName',[num2str(pars.INPUTS.Resonances(i,1)) ':' num2str(pars.INPUTS.Resonances(i,2))]);
end
legend show; hold on;
xlabel('Hyperbolic Excess Velocity [km/s]', 'FontSize', 18)
ylabel('Pump Angle [degrees]', 'FontSize', 18)

% Extract Pump for Pseudo-Resonances
for i = 1:length(pars.INPUTS.Resonances)
    for j = 1:2

        % Results for one pseudo-resonance & one type of transfer
        Pseudo_Res_aux = Pseudo_Resos(i,:,j);

        % Clean elements which are empty due to pseudo-resonance unfeasible
        non_empty_idx = false(1, length(Pseudo_Res_aux));
        for k = 1:length(Pseudo_Res_aux)
            non_empty_idx(k) = ~isempty(Pseudo_Res_aux(k).nodein);
        end
        Pseudo_Res_plot = Pseudo_Res_aux(non_empty_idx);

        % Extract the pump & Vinf values for the plot
        for m = 1:size(Pseudo_Res_plot,2)
            V_inf_plot(m) = Pseudo_Res_plot(m).nodein(1);
            Alfa_plot(m) = Pseudo_Res_plot(m).nodein(2);
        end

        if j == 1 % Type 1 of pseudo-resonance, long case (+)
            plot(V_inf_plot, rad2deg(Alfa_plot), '--' ,'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', [num2str(pars.INPUTS.Resonances(i,1)) ':' num2str(pars.INPUTS.Resonances(i,2)) '+'] );
            clear V_inf_plot; clear Alfa_plot;
        elseif j == 2
            plot(V_inf_plot, rad2deg(Alfa_plot), ':' ,'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', [num2str(pars.INPUTS.Resonances(i,1)) ':' num2str(pars.INPUTS.Resonances(i,2)) '-'] );
            clear V_inf_plot; clear Alfa_plot;
        end

        clear Pseudo_Res_plot;

    end

end

end