function [fig] = plotTextureLatLong(pl, idcentral, holdon)

% This function plots a planet/moon texture in a longitude-latitude map.

% Author: Andrea Bellome
% Last revision: 01/05/2024

%%  INPUTS %%
% - pl        : ID of the planet or moon. Depends upon the idcentral.
% - idcentral : ID of the central body (1. Sun, 5. Jupiter, 6. Saturn)
% - holdon    : optional input. If 1, then overwrites the current figure,
%               else, opens a new figure. If not specified, or empty, it
%               will open a new figure.

%% OUTPUTS %%
% - fig : figure structure

%% CHANGE LOG %%
%  01/05/2024, J.C. Garcia Mateas: added the 180° shift in longitude to
%              align surface features with lat/lon of the mercator. Added
%              also new xticks for that.

%% FUNCTION %%

if nargin == 2
    holdon = 0;
elseif nargin == 3
    if isempty(holdon)
        holdon = 0;
    end
end

if holdon == 0
    figure( 'Color', [1 1 1] );
end

hold on; axis equal;
xlabel( 'Longitude - deg' ); ylabel( 'Latitude - deg' );

name = tissPlId2Name(idcentral, pl);
hold on;
% xlim([-180 180]); ylim([-90 90]);
xlim([0 360]); ylim([-90 90]);
xticks([0 30 60 90 120 150 180 210 240 270 300 330 360]); % Define custom x-axis tick positions
xticklabels({'0°E','30°E','60°E', '90°E','120°E','150°E', '180°E','210°E','240°E', '270°E','300°E','330°E', '360°E'}); % Customize x-axis tick labels
I = imread([name '.jpg']);
h = image(xlim, -ylim, I);
grid on;
uistack(h,'bottom');

fig = gcf;

end
