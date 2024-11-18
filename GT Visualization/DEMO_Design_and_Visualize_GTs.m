%% DEMO SCRIPT TO PLOT GROUNDTRACKS %%
clear all; close all; clc; warning off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Name: [DEMO_Design_and_Visualize_GTs.m]
% Author: [Jose Carlos García & Joan Pau Sánchez]
% Date: 15/11/2024
% Version: 1.0
% 
% Description:
% This script computes and plots a ground track. In particular, it plots
% the v-Infinity_sphere of the fly-by, it plots the ground trach in
% mercator projections and it plots the 3D fly-by with the b-plane. The
% script purpose is to serve as a base to prepare visualizations of fly-by
% sequences.  To design a ground track both the entrance assymptote
% (NODEIN) and the exit assymptote (NODEOUT) need to be  defined.
%
% Usage:
% To design a ground track both the entrance assymptote (NODEIN) and the
% exit assymptote (NODEOUT) need to be  defined.
% Example:
% >> run('DEMO_Design_and_Visualize_GTs.m')
%
% Inputs:
%
% Outputs:
%
% Dependencies:
% - MATLAB R2021b or later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%% DEFINE PLANET, MOON PARAMETERS & CONSTANTS %%
% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = 1; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
                            % (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
                            % (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral); %[km3/s2],[km],[km] & [km]

% Retrieve Desired Moon Parameters 
if pars.INPUTS.idCentral == 3
    pars.Moon.OrbRad = 384748; pars.Moon.mu  = getAstroConstants('Moon','mu'); %[km],[km3/s2]
    pars.Moon.EquRad = getAstroConstants('Moon','Radius'); pars.Moon.hmin = 50;  %[km], [km]
elseif pars.INPUTS.idCentral == 5
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = jupMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
elseif pars.INPUTS.idCentral == 6
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
else
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = uranusMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
end

for i = 1:size(pars.INPUTS.idMoon,2)
    pars.Moon.Vel(i)    = sqrt(pars.Planet.mu/pars.Moon.OrbRad(i));           %[km/s] Moon Orbital velocity
    pars.Moon.Period(i) = 2*pi*sqrt(pars.Moon.OrbRad(i)^3/pars.Planet.mu);    %[s] Moon orbital period
    pars.Moon.HillSph(i) = pars.Moon.OrbRad(i)*( pars.Moon.mu(i)/(3*(pars.Moon.mu(i) + pars.Planet.mu)))^(1/3);     %[km] Moon Hill's Sphere
end
%--------------------------------------------------------------------------
%% DEFINE FLYBY PARAMETERS & CONSTANTS %%

% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 15;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 

% Define starting epoch
pars.INPUTS.epoch0 = 0;

% Load sectirs for groundtrack plots
load('sectorObj.mat');
load('sectorObj_.mat');
pars.sectorObj  = sectorObj;
pars.sectorObj_ = sectorObj_;
pars.rect       = defineRectangleMapping(); % --> save the rectangle mapping

%% COMPUTE RESONANT & PSEUDO-RESONANT TRANSFERS OF INTEREST %%
pars.INPUTS.Resonances = [7 1];
pars.INPUTS.V_inf      = 4;       %[km/s] 
pars.INPUTS.crank      = [0, pi]; %[rad]

% COMPUTE DATABASE OF RESONANT TRANSFERS %
Reso_TransfersDatabase = [];
for i = 1:size(pars.INPUTS.V_inf, 2)
    [Reso_Transfers, ResoStruc] = ResoTransfers_Generation(pars.INPUTS.V_inf(i), pars);
    Reso_TransfersDatabase = [Reso_TransfersDatabase; Reso_Transfers];
end

% PSEUDO-RESONANT TRANSFERS %
% Define which 81 pseudo-resonance to use
pars.INPUTS.remove81 = 1;  %[0] = add 1 rev to the outbound-inbound short case
                           %[1] = do not add one rev to the outbound-inbound short case

PseudoReso_TransfersDatabase = [];
for i = 1:size(pars.INPUTS.V_inf, 2)
    [PseudoReso_Transfer, PseudoResoStruct] = PseudoResoTransfers_Generation(pars.INPUTS.V_inf(i), pars);
    PseudoReso_TransfersDatabase = [PseudoReso_TransfersDatabase; PseudoReso_Transfer];
end

%% COMPUTE FLYBY CHARACTERISTICS %%
% Define incoming & outgoing nodes (Custom values can be defined if desired)
vinfin  = pars.INPUTS.V_inf;
alfain  = Reso_Transfers(1, 4);     
alfaou  = alfain;
kin     = pi/2 -  0.021490214656051; 
kou     = pi/2;
nodein  = [vinfin, alfain, kin]; %[km/s, rad, rad]
nodeout = [vinfin, alfaou, kou];

% Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;           %[km
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);    %[-] 
delta_max = 2*asin(1/e_fly);                                      %[rad]
pars.delta_max = delta_max;

% Compute flyby parameters
[Flyby] = Flyby_BuildUp(nodein, nodeout, pars);

%% PLOT GROUDNTRACK %%
colors = cool(length(Flyby));
fig1 = figure( 'Color', [1 1 1] );
hold on;
plotTextureLatLong(pars.INPUTS.idMoon , pars.INPUTS.idCentral , 1);
plotSquares(pars, 1);  
axis normal;
for i = 1:size(Flyby, 2)
    Plot_Flyby_GT(Flyby(i), colors(i,:));
end

%% PLOT FLYBY DISTANCE EVOLUTION %%
MaxDistance = 1000; % [km] Distance under which flyby Ground-Track is ploted

States   = Flyby.fly_States;
Distance = sqrt(States(:,1).^2+States(:,2).^2+States(:,3).^2);

Fig = figure('Color', [1 1 1]);clf;set(Fig,'defaulttextinterpreter','latex') ;
hold on; grid on; 
semilogy(Flyby.fly_tts, Distance, 'LineWidth', 2);
xlabel('Time to / from flyby periapsis below Max. Distance Level [seconds]', 'FontSize', 18);
ylabel('Distance to Moon center [km]', 'FontSize', 18);

%% PLOT VINF SPHERE %%

SetofAvailableOptions = plot_vinfinity_sphere(nodein,nodeout,pars);

%% PLOT BPLANE DISPLAY OF FLYBY %%

plot_bplane_flyby(Flyby, pars);

%% FUNCTION PLOT v_infity_sphere of the flyby %%
function [SetofAvailableOptions]=plot_vinfinity_sphere(NODEIN,NODEOUT,pars)


INComing_direction=[cos(NODEIN(2)), sin(NODEIN(2))*cos(NODEIN(3)), -sin(NODEIN(2))*sin(NODEIN(3))];

% Create figure
figure1 = figure('WindowState','maximized');
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% plot incoming direction
h1 = quiver3(0,0,0,INComing_direction(1),INComing_direction(2),INComing_direction(3),'LineWidth',2,'Color',[0 0.498039215803146 0],...
    'AutoScaleFactor',0.889999985694885, 'DisplayName','Incoming Direction');
axis equal
% plot vinf sphere
     nfaces = 20;
    [X,Y,Z] = sphere(nfaces);
% 4.2. ploting
surf (X, Y, Z,'Parent',axes1,'FaceColor','none','EdgeColor',[0.800000011920929 0.800000011920929 0.800000011920929]);

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


if ~isempty(NODEOUT)
OUTgoing_direction=[cos(NODEOUT(2)), sin(NODEOUT(2))*cos(NODEOUT(3)), -sin(NODEOUT(2))*sin(NODEOUT(3))]; %Outgoing directin = ORANGE color

h2 = quiver3(0,0,0,OUTgoing_direction(1),OUTgoing_direction(2),OUTgoing_direction(3),'LineWidth',2,...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'AutoScaleFactor',0.889999985694885,'DisplayName','Outgoing Direction');

PeriapsisPoint=(INComing_direction-OUTgoing_direction)/norm(INComing_direction-OUTgoing_direction); % Periapsis point = RED color
h3 = quiver3(0,0,0,PeriapsisPoint(1),PeriapsisPoint(2),PeriapsisPoint(3),'LineWidth',2,...
     'Color',[0.6350 0.0780 0.1840],...
    'AutoScaleFactor',0.889999985694885,'DisplayName','Periapsis Position');
end

% % PLOT THE CIRCLE OF CONSTANT ALFA
% Create zlabel
zlabel('z (Negative out of plane direction)');
% Create ylabel
ylabel('y (radial outward direction)');
% Create xlabel
xlabel('x (direction of planetary motion)');
view(axes1,[-10.8869329523003 10.4574443141852]);
grid(axes1,'on');
axis(axes1,'tight');
% hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1]);
% Create legend
legend([h1, h2, h3]);
set(gca, 'Clipping', 'off');
end

%% FUNCTION PLOT THE BPLANE AND 3D ORBIT OF THE FLYBY %%
function plot_bplane_flyby(Flyby, pars)

States=Flyby.fly_States;
Distance=sqrt(States(:,1).^2+States(:,2).^2+States(:,3).^2);

MaxDistance=max(Distance); % Define Distance under which the trajectory will be ploted.

[~,indexBpoint]=min(Distance);
Bpoint=States(indexBpoint,1:3);
EntranceHyperbola=States(1,:);
ExitHyperbola=States(end,:);

%--------------------------------------------------------------------------
% Define The B Plane Components
% The intention is to define a b-plane on which the coordinate Y shows a
% vector within the equator of the moon and perpendicular to the direction
% of motion of the spacecraft
T=EntranceHyperbola(4:6)/norm(EntranceHyperbola(4:6)); % Tangent Direction (i.e., perpendicular to the B Plane)
Z=[0 0 1]; % polar direction for the moon
S= crossFast(Z,T)/norm(crossFast(Z,T));
N =crossFast(T,S);
%--------------------------------------------------------------------------
% Draw Planet: Enceladus
UnitDistance = pars.Moon.OrbRad; %[km]
scale = 1;
figure1 = figure('Color',[1 1 1]);
% Create axes
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YColor',[1 1 1],...
    'XColor',[1 1 1]);
hold(axes1,'all');
MoonFB='Enceladus';
[hfig] = drawPlanet(MoonFB,[0 0 0],figure1,scale);
axis equal
view(0,90)
%--------------------------------------------------------------------------
% Draw the BPlane directions
CMPos=[0 0 0];
p1=MaxDistance*T;
Vinf = line([ CMPos(1) p1(1)],[ CMPos(2) p1(2)],[ CMPos(3) p1(3)],'Color',[0 0 0],'DisplayName','Relative Velocity Direction');
% Arrow head
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);
r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);
H=surf(x+p1(1),y+p1(2),z*h+p1(3));
set(H,'FaceColor','k','EdgeColor','k')
P=p1-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p1)
%--------------------------------------------------------------------------
CMPos=[0 0 0];
p2=MaxDistance*S;
line([ CMPos(1) p2(1)],[ CMPos(2) p2(2)],[ CMPos(3) p2(3)],'Color',[0 0 0],'DisplayName','b-plane Equatorial Direction')
% Arrow head
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);
r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);
H=surf(x+p2(1),y+p2(2),z*h+p2(3));
set(H,'FaceColor','k','EdgeColor','k')
P=p2-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p2)
%--------------------------------------------------------------------------
CMPos=[0 0 0];
p3=MaxDistance*N;
line([ CMPos(1) p3(1)],[ CMPos(2) p3(2)],[ CMPos(3) p3(3)],'Color',[0 0 0],'DisplayName','b-plane polar')
% Arrow head
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);
r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);
H=surf(x+p3(1),y+p3(2),z*h+p3(3));
set(H,'FaceColor','k','EdgeColor','k')
P=p3-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p3)
%--------------------------------------------------------------------------
% Quiver -Direction Planet Velocity (Leading Edge)
% VInfinDir=Flyby.vvinfin/norm(Flyby.vvinfin); %relative direction of motion
CMPos=[0 0 0];
State_planet=Flyby.State_planet;
rr=State_planet(1:3); vv=State_planet(4:6);
b1        = -rr./norm(rr);
b3        = cross(rr, vv)./norm(cross(rr, vv));
b2        = cross(b3,b1);
Rm        = [ b1' b2' b3' ]';
LEdge=State_planet(1,4:6)/norm(State_planet(1,4:6));
LEdge=(Rm*LEdge')'; % In general this should be the same as [0 1 0]

p1=MaxDistance*LEdge;
%    p2=[xEarth-0.3 0 0];
Moon_Velocity = line([ CMPos(1) p1(1)],[ CMPos(2) p1(2)],[ CMPos(3) p1(3)],'Color',[0.4660 0.6740 0.1880], 'LineWidth', 2, 'DisplayName','Moon Velocity Direction');
% Arrow head
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);
r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);
H=surf(x+p1(1),y+p1(2),z*h+p1(3));
set(H,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor',[0.4660 0.6740 0.1880])
P=p1-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p1)

% Plot b-plane as a transparent patch
PatchPoints=MaxDistance*[(N+S)' (N-S)' -(N+S)' (S-N)'];

 Xpoints_bplane=PatchPoints(1,:);
 Ypoints_bplane=PatchPoints(2,:);
 Zpoints_bplane=PatchPoints(3,:);
 patch(Xpoints_bplane,Ypoints_bplane,Zpoints_bplane, [0 0.749019622802734 0.749019622802734],'FaceAlpha',0.2);

% Plot Leading-trailing edge and sub-planet-antiplanet
[z,y,x]=cylinder(pars.Moon.EquRad(1),50);
  x=x(1,:);y=y(1,:);z=z(1,:);
 GreatCircle1=plot3(x,y,z,'LineWidth',2,'Color', [0 0 0],'DisplayName', 'Leading-Trailing Edge Circle' );
[x,z,y]=cylinder(pars.Moon.EquRad(1),50);
 x=x(1,:);y=y(1,:);z=z(1,:);

% Plot flyby
SC_Trajectory = plot3(States(pars.GroundTr.npoints/2:pars.GroundTr.npoints +pars.GroundTr.npoints/2 ,1),...
    States(pars.GroundTr.npoints/2:pars.GroundTr.npoints + pars.GroundTr.npoints/2,2), States(pars.GroundTr.npoints/2:pars.GroundTr.npoints + pars.GroundTr.npoints/2,3),...
    'LineWidth',2.5, 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'SC Trajectory');
 plot3(Bpoint(1),Bpoint(2),Bpoint(3),'MarkerFaceColor',[1 0 0],...
     'MarkerEdgeColor',[1 0 0],...
     'MarkerSize',4,...
     'Marker','o',...
     'LineStyle','none');
Bpoint_GT=pars.Moon.EquRad(1)*Bpoint/norm(Bpoint);

 Flyby_Peri = plot3(Bpoint_GT(1),Bpoint_GT(2),Bpoint_GT(3),'MarkerFaceColor',[1 0 0],...
     'MarkerEdgeColor',[1 0 0],...
     'MarkerSize',4,...
     'Marker','o',...
     'LineStyle','none', 'DisplayName', 'Periapsis Point');

% Plot Ground-track on the Sphere
States_GT=zeros(length(States(:,1)),3);
for iGT=1:length(States_GT(:,1)) 
    States_GT(iGT,:)=pars.Moon.EquRad(1)*States(iGT,1:3)/norm(States(iGT,1:3));
end
 GT = plot3(States_GT(:,1), States_GT(:,2), States_GT(:,3),'LineWidth',2.5,'Color',[0 1 1], 'DisplayName', 'Ground-Track');
camva('auto');
legend([Moon_Velocity, GreatCircle1, SC_Trajectory, GT, Flyby_Peri]);
set(gca, 'Clipping', 'off');
end