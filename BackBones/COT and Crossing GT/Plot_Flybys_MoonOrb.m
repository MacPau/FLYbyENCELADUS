function [Fig] = Plot_Flybys_MoonOrb(Flyby_Data, pars)

% This function plots the orbit view of a flyby around a body.

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
% Plot the figure
Fig = figure('Color', [1 1 1]);clf;set(Fig,'defaulttextinterpreter','latex') ;
hold on; grid on; axis equal; 
[x, y, z] = sphere(200);
h = surfl(x.*pars.Moon.EquRad, y.*pars.Moon.EquRad, z.*pars.Moon.EquRad);  % Plot of the Moon
set(h, 'FaceAlpha', 0.3, 'EdgeColor', 'None', 'FaceColor', [0.4660 0.6740 0.1880], 'handlevisibility', 'off');

% Define colormap for flyby trajectories
colors = cool(length(Flyby_Data));

for i = 1:length(Flyby_Data)

    states    = Flyby_Data(i).fly_States;
    altitudes = Flyby_Data(i).altitudes;

    pltotlim = pars.Moon.EquRad + 3500;
    indxs    = find(altitudes <= pltotlim);
    states   = states(indxs,:);

    if i == 1

        hold on;
        node = Flyby_Data(i).nodein;
        kin  = node(end);
        plot3(states(:,1), states(:,2), states(:,3), 'Color', colors(i,:), 'Linewidth', 2, 'DisplayName', ['First flyby, k_{in}: ' num2str(round(rad2deg(kin),1)) char(176)]);

    elseif i == length(Flyby_Data)

        hold on;
        node = Flyby_Data(i).nodein;
        kou  = node(end);
        plot3(states(:,1), states(:,2), states(:,3), 'Color', colors(i,:), 'Linewidth', 2, 'DisplayName', ['Last flyby, k_{in}: ' num2str(round(rad2deg(kou),1)) char(176)]);

    else
        hold on;
        plot3(states(:,1), states(:,2), states(:,3), 'Color', colors(i,:), 'Linewidth', 2, 'HandleVisibility', 'off' );
       
    end

    plotSingleAxis(states(1,1:3), states(10,1:3), colors(i,:));
    plotSingleAxis(states(end-10,1:3), states(end,1:3), colors(i,:));

end

view([45 5]);
xlabel('x [km]'); ylabel('y [km]'); zlabel( 'z [km]' );

plotSingleAxis([0 0 0], [1 0 0]*5e3);

legend( 'Location', 'best' );


end
