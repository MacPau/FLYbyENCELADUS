%% DEMO SCRIPT TO DESIGN & VISUALIZE  GROUNDTRACKS %%




%% CHOOSE YOUR FLYBY CHARACTERISTICS %%
% this part of the script can allow to chose the change of pumb and crank
% based on a menu of periapsis positions.

[rrga, vvga] = approxEphem_CC(parsPLOT.INPUTS.idMoon, parsPLOT.INPUTS.epoch0, parsPLOT.INPUTS.idCentral); % epoch0 provides only a dummy state vector for the moon.
%--------------------------------------------------------------------------
% STEP 1: Chose your starting node - NODEIN DEFINITION

NODEIN  = NODEIN_Matrix(indexNode,:); %[km/s, rad, rad]
% display info in screen
disp(['NODE IN : ',num2str(NODEIN(1)),' km/s ', num2str(NODEIN(2)),' rads ',num2str(NODEIN(3)),' rads '])


% Determine maximum bending due to flyby
rp_flyby  = parsPLOT.INPUTS.Flyby.min_h + parsPLOT.Moon.EquRad;           %[km
e_fly     = 1 + ((rp_flyby*NODEIN(1)^2)/parsPLOT.Moon.mu);    %[-]
delta_max = 2*asin(1/e_fly);                                      %[rad]
parsPLOT.delta_max = delta_max;

SetofAvailableOptions=plot_vinfinity_sphere(NODEIN,[],parsPLOT);
% Select the adequate positions of the periapsis among the available set,
% open the matrix SetofAvailableOptions and chose among the longitude and
% latitude availables
% SetofAvailableOptions(1) = longitude of the periapsis position
% SetofAvailableOptions(2) = latitude of the periapsis position
% SetofAvailableOptions(3) = Delta Pumb
% SetofAvailableOptions(4) = Delta Crank

% plot all possible flybys

for iPF=1:length(SetofAvailableOptions(:,1))

NODEOUT=[NODEIN(1), NODEIN(2)+SetofAvailableOptions(iPF,3), NODEIN(3)+SetofAvailableOptions(iPF,4)];

CARNODEIN=vinfAlphaCrank2car(NODEIN, [rrga vvga],  parsPLOT.Moon.mu);
CARNODEOUT=vinfAlphaCrank2car(NODEOUT, [rrga vvga], parsPLOT.Moon.mu);
delta = acos( dot(CARNODEIN, CARNODEOUT)./(norm(CARNODEIN)*norm(CARNODEOUT)));

% display info in screen
disp(['NODE OUT : ',num2str(NODEOUT(1)),' km/s ', num2str(NODEOUT(2)),' rads ',num2str(NODEOUT(3)),' rads '])
disp(['Deflection angle is ',num2str(delta),' rads'])
disp(['Maximum Deflection angle is ',num2str(delta_max),' rads'])

% Compute flyby parameters
Flyby(iPF) = Flyby_BuildUp(NODEIN, NODEOUT, parsPLOT);

end

%% PLOT GROUDNTRACK %%
colors = cool(length(Flyby));
fig1 = figure( 'Color', [1 1 1] );
hold on;
plotTextureLatLong(parsPLOT.INPUTS.idMoon , parsPLOT.INPUTS.idCentral , 1);
plotSquares(parsPLOT, 1);
axis normal;
for i = 1:size(Flyby, 2)
    Plot_Flyby_GT(Flyby(i), colors(i,:));
end

%% Plot the v_infity_sphere of the flyby
function [SetofAvailableOptions]=plot_vinfinity_sphere(NODEIN,NODEOUT,pars)


INComing_direction=[cos(NODEIN(2)), sin(NODEIN(2))*cos(NODEIN(3)), -sin(NODEIN(2))*sin(NODEIN(3))];


% plot all possible positions assuming maximum deflection
RotationDirection=cross(INComing_direction,[0 0 1])/norm(cross(INComing_direction,[0 0 1]));
v1 = eulerAxisAngle(INComing_direction,RotationDirection,pars.delta_max);

thetaRot=linspace(0,2*pi);
SetofAvailableOptions=zeros(length(thetaRot),3);
AvailablePeriapis_alphas_kas=zeros(length(thetaRot),7);
for i=1:length(thetaRot)
    SetofAvailableOptions(i,:)=eulerAxisAngle(v1,INComing_direction,thetaRot(i));
end


PossibleDelta_alphas=linspace(-pars.delta_max,pars.delta_max);
SetofAvailableOptions=zeros(2*length(PossibleDelta_alphas),4);
for j=1:length(PossibleDelta_alphas)
    cos_Delta_Crank=(cos(pars.delta_max)-(cos(NODEIN(2))^2*cos(PossibleDelta_alphas(j))-cos(NODEIN(2))*sin(NODEIN(2))*sin(PossibleDelta_alphas(j))))/(sin(NODEIN(2))^2*cos(PossibleDelta_alphas(j))+sin(NODEIN(2))*cos(NODEIN(2))*sin(PossibleDelta_alphas(j)));
    if (cos_Delta_Crank-1)<1e-6 && cos_Delta_Crank>1
        cos_Delta_Crank=1;
    end
    Delta_Crank_2max=[acos(cos_Delta_Crank) -acos(cos_Delta_Crank)];
    % NODEOut Definition
    AvailableNODEOUT1=[NODEIN(1), NODEIN(2)+PossibleDelta_alphas(j), NODEIN(3)+Delta_Crank_2max(1)];
    AvailableNODEOUT2=[NODEIN(1), NODEIN(2)+PossibleDelta_alphas(j), NODEIN(3)+Delta_Crank_2max(2)];

    OUTgoing_direction1=[cos(AvailableNODEOUT1(2)), sin(AvailableNODEOUT1(2))*cos(AvailableNODEOUT1(3)), -sin(AvailableNODEOUT1(2))*sin(AvailableNODEOUT1(3))];
    OUTgoing_direction2=[cos(AvailableNODEOUT2(2)), sin(AvailableNODEOUT2(2))*cos(AvailableNODEOUT2(3)), -sin(AvailableNODEOUT2(2))*sin(AvailableNODEOUT2(3))];

    PossiblePeri1=(INComing_direction-OUTgoing_direction1)/norm(INComing_direction-OUTgoing_direction1);
    LAT1         = asin(PossiblePeri1(3));         % [rad] Latitude of the ground-track points
    LONG1        = atan2( PossiblePeri1(2),PossiblePeri1(1) )+pi/2;   % [rad] Longitude of the ground-track points (added rotation of the x axis from the moon velocity to the radial direction)
    PossiblePeri2=(INComing_direction-OUTgoing_direction2)/norm(INComing_direction-OUTgoing_direction2);
    LAT2         = asin(PossiblePeri2(3));         % [rad] Latitude of the ground-track points
    LONG2        = atan2( PossiblePeri2(2),PossiblePeri2(1) )+pi/2;   % [rad] Longitude of the ground-track points

    Summary=[rad2deg(LONG1) rad2deg(LAT1)  PossibleDelta_alphas(j) Delta_Crank_2max(1); rad2deg(LONG2) rad2deg(LAT2) PossibleDelta_alphas(j) Delta_Crank_2max(2)];

    SetofAvailableOptions(2*j-1:2*j,:)=Summary;
end

end
