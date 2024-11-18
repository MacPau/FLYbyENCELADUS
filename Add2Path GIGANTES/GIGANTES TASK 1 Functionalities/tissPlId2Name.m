function [name] = tissPlId2Name(idcentral, pl)

% This function returns the name of the body (planet or moon) at which the
% flyby is being performed.

% Author: A. Bellome
% Last revision: 27/08/2026

%% INPUTS %%
% - idcentral: number indicating which is the central body of the system.
%
% - pl: number indicating which is the body at which the flyby is performed

%% OUTPUTS %%
% - name: string of characters containing the name of the body at which the
%         flyby ocurs.

%% CHANGELOG %%
% - 27/08/2026, J.C Garcia Mateas: added the option in the loop for Uranus
%               and its moons.

%% FUNCTION %%

if idcentral == 1 % --> central body is Sun

    if pl == 1
        name = 'Mercury';
    elseif pl == 2
        name = 'Venus';
    elseif pl == 3
        name = 'Earth';
    elseif pl == 4
        name = 'Mars';
    elseif pl == 5
        name = 'Jupiter';
    elseif pl == 6
        name = 'Saturn';
    elseif pl == 7
        name = 'Uranus';
    elseif pl == 8
        name = 'Neptune';
    elseif pl == 9
        name = 'Pluto';
    end

elseif idcentral == 5 % --> central body is Jupiter

    if pl == 1
        name = 'Io';
    elseif pl == 2
        name = 'Europa';
    elseif pl == 3
        name = 'Ganymede';
    elseif pl == 4
        name = 'Callisto';
    end

elseif idcentral == 6 % --> central body is Saturn

    if pl == 1
        name = 'Enceladus';
    elseif pl == 2
        name = 'Thetys';
    elseif pl == 3
        name = 'Dione';
    elseif pl == 4
        name = 'Rhea';
    elseif pl == 5
        name = 'Titan';
    end

elseif idcentral == 7 % --> central body is Uranus

    if pl == 1
        name = 'Miranda';
    elseif pl == 2
        name = 'Ariel';
    elseif pl == 3
        name = 'Umbriel';
    elseif pl == 4
        name = 'Titania';
    elseif pl == 5
        name = 'Oberon';
    end

end

end
