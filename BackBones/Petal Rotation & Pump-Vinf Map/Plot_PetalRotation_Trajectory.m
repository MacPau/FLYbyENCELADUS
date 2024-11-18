function [Fig] = Plot_PetalRotation_Trajectory(PetalRotation_Sequence, ~, pars)

% This function plots the trajectory and sequence of flyby locations
% associated to a petal rotation. The plot includes the SC orbits, markers
% to show the flyby locations and an arrow pointing in the Planet - Sun
% direction to provide insight on illumination conditions.

% Author: Jose Carlos Garcia Mateas
% Last revision: 27/10/2024

%% INPUTS %%
% - PetalRotation_Sequence: Matrix of size N x 10, where N = n° of
%                           transfers associated to the petal rotation. The
%                           information in the columns is: col1 = moon
%                           revs; col2 = SC revs; col3 = Vinf [km/s] ;
%                           col4 = pump angle [rad]; col5 = crank; 
%                           col6 = TOF [seconds]; cols 7,8 & 9 = Vinf, pump,
%                           crank at next flyby after the TOF; col10 = type
%                           of transfer (11 = Inbound - Inbound / 18 =
%                           Inbound - Outbound / 81 = Out-In / 88 = Out-Out)
%
% - pars: Structure containing different parameters.

%% OUTPUTS %%
% - Fig: Figure of the petal rotation trajectory

%% FUNCTION %%
colors = cool(size(PetalRotation_Sequence,1));

% Plot the figure
Fig = plotMoons( pars.INPUTS.idMoon, pars.INPUTS.idCentral);
axis equal; hold on;
for i = 1:size(PetalRotation_Sequence, 1)

    PseudoRes = PetalRotation_Sequence(i,:);
    TOF(i)    = PseudoRes(6);
    if i == 1
        [~, rr1, vv1, ~] = vinfAlphaCrank_to_VinfCART(PseudoRes(3), PseudoRes(4), PseudoRes(5) , pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
        [~, yy] = propagateKeplerODE(rr1, vv1, linspace(0, TOF(i), 500), pars.Planet.mu);

        marker_start(i,:) = [yy(1,1), yy(1,2), yy(1,3)];
        marker_end(i,:) = [yy(end,1), yy(end,2), yy(end,3)];

        plot3( yy(:,1), yy(:,2), yy(:,3), 'LineWidth', 1.5,  'Color', colors(i,:));
        plot3(yy(1,1), yy(1,2), yy(1,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black', 'MarkerSize', 8 );
        plot3(yy(end,1), yy(end,2), yy(end,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black', 'MarkerSize', 8 );

    elseif i == 2
        [~, rr1, vv1, ~] = vinfAlphaCrank_to_VinfCART(PseudoRes(3), PseudoRes(4), PseudoRes(5), pars.INPUTS.epoch0 + sum(TOF(1:i-1))/86400, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
        [~, yy] = propagateKeplerODE(rr1, vv1, linspace(0, TOF(i), 500), pars.Planet.mu);

        marker_start(i,:) = [yy(1,1), yy(1,2), yy(1,3)];
        marker_end(i,:) = [yy(end,1), yy(end,2), yy(end,3)];

        plot3( yy(:,1), yy(:,2), yy(:,3), 'LineWidth', 1.5,  'Color', colors(i,:));
        plot3(yy(end,1), yy(end,2), yy(end,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 8 );

    else
        [~, rr1, vv1, ~] = vinfAlphaCrank_to_VinfCART(PseudoRes(3), PseudoRes(4), PseudoRes(5), pars.INPUTS.epoch0 + sum(TOF(1:i-1))/86400, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
        % [~, yy]          = propagateKepler_tof(rr1, vv1, TOF(i), pars.Planet.mu);
        [~, yy] = propagateKeplerODE(rr1, vv1, linspace(0, TOF(i), 500), pars.Planet.mu);

        marker_start(i,:) = [yy(1,1), yy(1,2), yy(1,3)];
        marker_end(i,:) = [yy(end,1), yy(end,2), yy(end,3)];

        plot3( yy(:,1), yy(:,2), yy(:,3), 'LineWidth', 1.5,  'Color', colors(i,:));
        plot3(yy(1,1), yy(1,2), yy(1,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 8 );
        plot3(yy(end,1), yy(end,2), yy(end,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 8 );
    end

    % Convert last state vector to OE & Store
    OE_SC_end(i,:)  = car2kep(yy(end,:),pars.Planet.mu);
    clear yy;
end


% Retrieving Sun - Planet line to plot Illumination direction 
[r,v, kep] = EphSS_car(pars.INPUTS.idCentral , pars.INPUTS.epoch0);
unit_vect = [-r(1)/norm(r), -r(2)/norm(r), r(3)/norm(r)]*5*pars.Planet.EquRad; % Rotate 180° for line to be Planet - Sun

% Plotting the unit vector as an arrow
start_point = [0, 0, 0];
quiver3(start_point(1), start_point(2), start_point(3), unit_vect(1), unit_vect(2), unit_vect(3),'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Create markers for the legend
h1 = plot3(NaN, NaN, NaN, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black', 'MarkerSize', 8);
h2 = plot3(NaN, NaN, NaN, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 8);
h3 = plot3(NaN, NaN, NaN, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 8);
h4 = plot3(NaN, NaN, NaN, 'r', 'LineWidth', 2); 

legend([h1, h2, h3], {'Initial Flyby Position', 'Intermediate Flybys Position', 'Final Flyby Position'}, 'Location', 'best');
% legend([h1, h2, h3, h4], {'Initial Flyby Position', 'Intermediate Flybys Position', 'Final Flyby Position', 'Planet - Sun direction'}, 'Location', 'best');


% Retrieve SC True Longitude at the initial flyby
TrueLong_0 = mod(sum(OE_SC_end(1, 4:6)), 2*pi);

% Compute change in True Longitude achieved by each transfer
for i = 2:size(OE_SC_end,1)
        TrueLong(i)         = mod(sum(OE_SC_end(i,4:6)), 2*pi); % [rad]
        TrueLong_Diff(i, :) = TrueLong(i) - TrueLong_0;         % [rad]
end

end


%% Auxiliary function
function [tt, yy] = propagateKeplerODE(rvec, vvec, timevector, muPL)

%Equation of motion ¨r+mu*r/R^3 = T/m
F=@(t,x)   [x(4); %dx/dt=Vx
    x(5); %dy/dt=Vy
    x(6); %dz/dt=Vz
    -muPL*x(1)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3);
    -muPL*x(2)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3);
    -muPL*x(3)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)];

options=odeset('RelTol',1e-6,'AbsTol',1e-7,'Refine',50);
[tt,yy]=ode45(F,timevector,[rvec vvec],options);

end