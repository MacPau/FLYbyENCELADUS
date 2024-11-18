function [] = Plot_Flyby_GT(Flyby_Data, color)

% This function plots the Ground-track of a Flyby.

% Author: Jose Carlos garcia Mateas
% Last revision: 17/05/2024

%% INPUTS %%
%
% - Flyby_Data: Structure with 16 fields containing all the information and
%               parameters of the flyby linking the Vinf, pump and crank
%               defined in the fields "nodein" and "nodeout".
%
% - color: Double of size 1 x 3 defining the color of the ground-track to
%          be plotted.

%% FUNCTION %%
% Extract Flyby Data
lats    = Flyby_Data.lats;
longs   = Flyby_Data.longs;
rp_lat  = Flyby_Data.rp_lat;
rp_long = Flyby_Data.rp_long;

% Find if there is a jump in longitude (passing from near 0 to 360)
stop_index = find(diff(longs) > 0, 1);

% If there's a jump, split the vectors into two parts
if ~isempty(stop_index) && ~isequal(stop_index, 1)

    % Split the vector at the negative difference indices
    Long_1 = longs(1:stop_index);
    Lat_1  = lats(1:stop_index);
    Long_2 = longs(stop_index+1:end);
    Lat_2  = lats(stop_index+1:end);
    
    plot(rad2deg(Long_1), rad2deg(Lat_1), '-', 'Color', color, 'LineWidth', 2);
    plot(rad2deg(Long_2), rad2deg(Lat_2), '-', 'Color', color, 'LineWidth', 2);
    plot(rad2deg(Long_1(1)), rad2deg(Lat_1(1)), 'o', 'MarkerSize', 24,'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black', 'MarkerSize', 5 ); % Entry point of the flyby hyperbola
    plot(rad2deg(Long_2(end)), rad2deg(Lat_2(end)), 'o', 'MarkerSize', 24,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerSize', 5 ); % Exit point of the flyby hyperbola
    plot(rp_long, rp_lat, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 5 );

    lineHandles = findobj(gca, 'Type', 'line');
    lineHandles = flipud(lineHandles);
    set(lineHandles(1), 'HandleVisibility', 'off');

    legend('Groundtrack', 'Flyby Entry Point', 'Flyby Exit Point' ,'Flyby Periapsis');

else
%     plot(rad2deg(longs), rad2deg(lats), '-', 'Color', color, 'LineWidth', 2, 'DisplayName', 'Groundtrack');
%     plot(rad2deg(longs(1)), rad2deg(lats(1)), 'o', 'MarkerSize', 24,'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue', 'MarkerSize', 5, 'DisplayName', 'Flyby Entry Point' ); % Entry point of the flyby hyperbola
%     plot(rad2deg(longs(end)), rad2deg(lats(end)), 'o', 'MarkerSize', 24,'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerSize', 5, 'DisplayName', 'Flyby Exit Point' ); % Exit point of the flyby hyperbola
%     plot(rp_long, rp_lat, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 5, 'DisplayName', 'Flyby Periapsis Point');

    plot(rad2deg(longs), rad2deg(lats), '-', 'Color', color, 'LineWidth', 2);
    plot(rad2deg(longs(1)), rad2deg(lats(1)), 'o', 'MarkerSize', 24,'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black', 'MarkerSize', 5 ); % Entry point of the flyby hyperbola
    plot(rad2deg(longs(end)), rad2deg(lats(end)), 'o', 'MarkerSize', 24,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerSize', 5); % Exit point of the flyby hyperbola
    plot(rp_long, rp_lat, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 5);

    lineHandles = findobj(gca, 'Type', 'line');
    lineHandles = flipud(lineHandles);
    legend('Groundtrack', 'Flyby Entry Point', 'Flyby Exit Point' ,'Flyby Periapsis');

end

% legend('Groundtrack below 300km', 'Flyby Entry Point', 'Flyby Exit Point' ,'Flyby Periapsis');
legend show;



end