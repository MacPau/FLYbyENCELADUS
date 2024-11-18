function plot_bplane_flyby(Flyby, pars, usecurrentFigure)

%--------------------------------------------------------------------------
% Preraring Figure Object
if usecurrentFigure~=1
    figure1 = figure('Color',[1 1 1]);
    % Create axes
    axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YColor',[1 1 1],...
        'XColor',[1 1 1]);
    hold(axes1,'all');
else
    figure1=gcf;
    axes1=gca;
end
%--------------------------------------------------------------------------


States=Flyby.fly_States;
Distance=sqrt(States(:,1).^2+States(:,2).^2+States(:,3).^2);

MaxDistance=max(Distance); % Define Distance under which the trajectory will be ploted.

%
% indexDel=find(Distance>MaxDistance);
% States(indexDel,:)=[]; Distance(indexDel)=[]; % deleting states that are above the maximum distance.

[~,indexBpoint]=min(Distance);
Bpoint=States(indexBpoint,1:3);
EntranceHyperbola=States(1,:);
ExitHyperbola=States(end,:);

%--------------------------------------------------------------------------
% should not the velocity at the beging of the fly by, be directed on the
% same direction as Flyby.vvinfin?


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
UnitDistance = pars.Moon.OrbRad; % Km, from takubo et al.
scale = 1;
MoonFB='Enceladus';
[hfig] = drawPlanet(MoonFB,[0 0 0],figure1,scale);
axis equal
view(0,90)
%--------------------------------------------------------------------------
% Draw the BPlane directions
CMPos=[0 0 0];
p1=MaxDistance*T;
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) p1(1)],[ CMPos(2) p1(2)],[ CMPos(3) p1(3)],'Color',[0 0 0],'DisplayName','Relative Velocity Direction')
% Arrow head
% scale factor
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
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) p2(1)],[ CMPos(2) p2(2)],[ CMPos(3) p2(3)],'Color',[0 0 0],'DisplayName','b-plane Equatorial Direction')
% Arrow head
% scale factor
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
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) p3(1)],[ CMPos(2) p3(2)],[ CMPos(3) p3(3)],'Color',[0 0 0],'DisplayName','b-plane polar')
% Arrow head
% scale factor
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
line([ CMPos(1) p1(1)],[ CMPos(2) p1(2)],[ CMPos(3) p1(3)],'Color',[0 0.749019622802734 0.749019622802734])
% Arrow head
% scale factor
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
set(H,'FaceColor',[0 0.749019622802734 0.749019622802734],'EdgeColor',[0 0.749019622802734 0.749019622802734])
P=p1-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p1)
%--------------------------------------------------------------------------
% Plot b-plane as a transparent patch

PatchPoints=MaxDistance*[(N+S)' (N-S)' -(N+S)' (S-N)'];

Xpoints_bplane=PatchPoints(1,:);
Ypoints_bplane=PatchPoints(2,:);
Zpoints_bplane=PatchPoints(3,:);
patch(Xpoints_bplane,Ypoints_bplane,Zpoints_bplane, [0 0.749019622802734 0.749019622802734],'FaceAlpha',0.2);


%--------------------------------------------------------------------------
% Plot Leading-trailing edge and sub-planet-antiplanet

[z,y,x]=cylinder(pars.Moon.EquRad(1),50);
x=x(1,:);y=y(1,:);z=z(1,:);
GreatCircle1=plot3(x,y,z,'LineWidth',2,'Color',[0 1 1],'DisplayName', 'Leading-Trailing Edge Circle' )

[x,z,y]=cylinder(pars.Moon.EquRad(1),50);
x=x(1,:);y=y(1,:);z=z(1,:);
GreatCircle2=plot3(x,y,z,'LineWidth',2,'Color',[1 0 1],'DisplayName', 'Antiplanet-Subplanet Edge Circle' )

%---------------------------------------------------------------------------
% % plot fly_by
plot3(States(:,1), States(:,2), States(:,3),'LineWidth',2.5)
% % How periapsis point

plot3(Bpoint(1),Bpoint(2),Bpoint(3),'MarkerFaceColor',[1 1 0],...
    'MarkerEdgeColor',[0.200000002980232 1 0],...
    'MarkerSize',4,...
    'Marker','o',...
    'LineStyle','none');
Bpoint_GT=pars.Moon.EquRad(1)*Bpoint/norm(Bpoint);

plot3(Bpoint_GT(1),Bpoint_GT(2),Bpoint_GT(3),'MarkerFaceColor',[1 1 0],...
    'MarkerEdgeColor',[0.200000002980232 1 0],...
    'MarkerSize',4,...
    'Marker','o',...
    'LineStyle','none');

%  % Show entrance of the assymptote
%  line([0 EntranceHyperbola(1)],[0 EntranceHyperbola(2)], [0 EntranceHyperbola(3)])
% % % Show exit of the assymptote
%  line([0 ExitHyperbola(1)],[0 ExitHyperbola(2)], [0 ExitHyperbola(3)])

%-------------------------------------------------------------------------
% Plot Ground track on the Sphere
States_GT=zeros(length(States(:,1)),3);
for iGT=1:length(States_GT(:,1))
    States_GT(iGT,:)=pars.Moon.EquRad(1)*States(iGT,1:3)/norm(States(iGT,1:3));
end
plot3(States_GT(:,1), States_GT(:,2), States_GT(:,3),'LineWidth',2.5,'Color',[0 0 0])

view(65, 0)
end
