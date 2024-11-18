function [Fig] = Plot_Flybys_CentralBOrb(Flyby_Data, pars)

% This function plots the orbit view of a flyby around a the main central 
% body (i.e: Jupiter or Saturn).

% Author: J.C Garcia Mateas
% Last revision: 15/11/2024

%% INPUTS %%
% - Flyby_Data: Structure with 12 fields containing the data regarding the
%               flybys to be plotted. The size of the structure (1 x A)
%               corresponds to the n° of flybys (A) to be plotted.
%
% - pars: Structure containing different parameters & constants.

%% OUTPUTS %%
% - Fig: Figure of the plotted planeto-centric trajectories.

%% FUNCTION %%
% Extract info of the first flyby
node1 = Flyby_Data(1).nodein;

% Compute the planeto-centric orbit
[~, rrin, vvin] = vinfAlphaCrank_to_VinfCART(node1(1), node1(2), node1(3), pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
kepin           = car2kep([ rrin, vvin ], pars.Planet.mu);
tof             = 2*pi*sqrt(kepin(1)^3/pars.Planet.mu);                 %[seconds]
[~, yy]         = propagateKepler_tof(rrin, vvin, tof, pars.Planet.mu); %[km] & [km/s]

% Compute the periapsis point of the 1st planeto-centric orbit
for i = 1:size(yy,1)
    yy_norm(i) = norm(yy(i,1:3));
end
[~, idx] = min(yy_norm);

% kep_peri             = kepin;
% kep_peri(end)        = 0;    % TA @ periapsis
% [peri_pos, peri_vel] = OE_to_Cartesian(kep_peri ,pars.Planet.mu); %[km] & [km/s]

% Define colormap for each orbit
colors = cool(length(Flyby_Data));

% Plot the figure
Fig = figure('Color', [1 1 1]);clf;set(Fig,'defaulttextinterpreter','latex') ;
hold on; grid on; axis equal; 
hold on; grid on;
plotMoons(pars.INPUTS.idMoon, pars.INPUTS.idCentral, 1);
hold on;
plot3(yy(:,1), yy(:,2), yy(:,3), 'LineWidth', 2, 'Color', colors(1,:), 'DisplayName', ['First flyby orbit, k_{in}  = ' num2str(round(rad2deg(node1(3)),1)) char(176)]);  % 1st Orbit
plot3(yy(1,1), yy(1,2), yy(1,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'DisplayName',' Flyby point');                         % Point where flybys occur
plot3(yy(idx,1), yy(idx,2), yy(idx,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'DisplayName','Periapsis of 1^{st} SC orbit'); % Periapsis point of Spacecraft Central Body 1st Orbit
% plot3(peri_pos(1), peri_pos(2), peri_pos(3), 'o', 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green', 'DisplayName','Periapsis of SC Orbit: '); % Periapsis point of Spacecraft Central Body Orbit

for i = 2:length(Flyby_Data)

    node = Flyby_Data(i).nodein;

    % Retrieve Cartesian state & propagate orbit
    [~, rr, vv] = vinfAlphaCrank_to_VinfCART(node(1), node(2), node(3), pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
    [~, yy]     = propagateKepler_tof(rr, vv, tof, pars.Planet.mu);

    if i == length(Flyby_Data)

        hold on;
        plot3( yy(:,1), yy(:,2), yy(:,3), 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', ['Last flyby orbit, k_{in}  = ' num2str(round(rad2deg(node(3)),1)) char(176)]);

    else

        hold on;
        plot3( yy(:,1), yy(:,2), yy(:,3), 'LineWidth', 2, 'Color', colors(i,:), 'HandleVisibility', 'off');
        
    end

end

legend show;

end