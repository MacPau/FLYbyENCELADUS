function plotfullSaturnTransfer(Full_route,pars)
%--------------------------------------------------------------------------
% preliminaries
DAYS2SECS=3600*24;
%--------------------------------------------------------------------------
% Epremeris Moons Set up
KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
SVTitan_func = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
SVEnceladus_func = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
%--------------------------------------------------------------------------
nTitan=sqrt(pars.Planet.mu/pars.Moon.OrbRad(2)^3);
nEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);


M0Titan=Full_route(1,5);
M0Enceladus=Full_route(1,6);

nlegs=length(Full_route(:,1));
%--------------------------------------------------------------------------
% for each leg do the following

% Plot reference frame
Plot_Saturn_Equatorial_frame_withEnceladus_Titan_Orbits(pars)



FigPlot_saturnTrajectory=gcf;
AxesPlot_saturnTrajectory=gca;

for iLeg=1:nlegs

    TOF_seconds = Full_route(iLeg,7)*3600*24; %[seconds]

    % Departure Position & Velocity
    DepFlag=Full_route(iLeg,1);

    if DepFlag==5
        % Titan position & velocity
        MAnomaly_MoonDeparture=Full_route(iLeg,5);
        SVatTitan=SVTitan_func(MAnomaly_MoonDeparture);
        r_moon_departure=SVatTitan(1:3);
        v_moon_departure=SVatTitan(4:6);
    elseif DepFlag==1
        % Enceladus position & velocity
        MAnomaly_MoonDeparture=Full_route(iLeg,6);
        SVatEnceladus=SVEnceladus_func(MAnomaly_MoonDeparture);
        r_moon_departure=SVatEnceladus(1:3);
        v_moon_departure=SVatEnceladus(4:6);
    end


    % Arrival Position & Velocity
    ArrFlag=Full_route(iLeg,8);
    if ArrFlag==5
        % Titan position & velocity
        MAnomaly_MoonArrival=Full_route(iLeg,12);
        SVatTitan=SVTitan_func(MAnomaly_MoonArrival);
        r_moon_arrival=SVatTitan(1:3);
        v_moon_arrival=SVatTitan(4:6);
    elseif ArrFlag==1
        % Enceladus position & velocity
        MAnomaly_MoonArrival=Full_route(iLeg,13);
        SVatEnceladus=SVEnceladus_func(MAnomaly_MoonArrival);
        r_moon_arrival=SVatEnceladus(1:3);
        v_moon_arrival=SVatEnceladus(4:6);
    end

    %--------------------------------------------------------------------------
    vinf_plus_departure= vinfAlphaCrank2car(Full_route(iLeg,2:4), [r_moon_departure v_moon_departure], pars.Planet.mu);
    v_sc_departure=vinf_plus_departure+v_moon_departure;

    vinf_minus_arrival= vinfAlphaCrank2car(Full_route(iLeg,9:11), [r_moon_arrival v_moon_arrival], pars.Planet.mu);
    v_sc_arrival=vinf_minus_arrival+v_moon_arrival;
    % Compute distance with Titan
    %--------------------------------------------------------------------------
    % Propagate from departure from Titan and up to 30 days
    % I should generate the state vector from the Lambert arc and from the
    % resonance calculation to ensure that I obtain the same results.
    % in the first instance, I will simply use the Lambert arc assuming that
    % these are the same (which should be if everything was correct).

    SV_Spacecraft=[r_moon_departure v_sc_departure];
    OE_Spacecraft=car2kep(SV_Spacecraft,pars.Planet.mu);
    % Propagation time
    tf=TOF_seconds; % [sec] A trial and error check will show that you need at

    %Equation of motion ¨r+mu*r/R^3 = 0
    F2BDyn=@(t,x)   [x(4); %dx/dt=Vx
        x(5); %dy/dt=Vy
        x(6); %dz/dt=Vz
        -pars.Planet.mu*x(1)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3); %dVx/dt
        -pars.Planet.mu*x(2)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3); %dVx/dt
        -pars.Planet.mu*x(3)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)]; %dVx/dt

    %Propagation
    options=odeset('RelTol',1e-6,'AbsTol',1e-6);
    [Tstep,SVspacecraft]=ode45(F2BDyn,[0:3600:tf],SV_Spacecraft,options);   % Any ODE solver should


    %Propagation of Departure Moom
    options=odeset('RelTol',1e-6,'AbsTol',1e-6);
    [~,SVoftitan]=ode45(F2BDyn,Tstep,[r_moon_departure v_moon_departure],options);   % Any ODE solver should


    %Propagation of Arrival Moon
    options=odeset('RelTol',1e-6,'AbsTol',1e-6);
    [~,SVBackwardProp]=ode45(F2BDyn,[0 -tf],[r_moon_arrival v_moon_arrival],options);
    [~,SVofEnceladus]=ode45(F2BDyn,Tstep,SVBackwardProp(end,:),options);   % Any ODE solver should


    plot3(AxesPlot_saturnTrajectory,r_moon_departure(1),r_moon_departure(2),r_moon_departure(3),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',6,...
        'Marker','o',...
        'LineStyle','none');

    plot3(AxesPlot_saturnTrajectory,SVspacecraft(:,1),SVspacecraft(:,2),SVspacecraft(:,3))
    view(0,90)
    plot3(AxesPlot_saturnTrajectory,r_moon_arrival(1),r_moon_arrival(2),r_moon_arrival(3),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',6,...
        'Marker','o',...
        'LineStyle','none');

    % distance_titan=sqrt((SVspacecraft(:,1)-SVoftitan(:,1)).^2+(SVspacecraft(:,2)-SVoftitan(:,2)).^2+(SVspacecraft(:,3)-SVoftitan(:,3)).^2)/pars.Moon.HillSph(2);
    % distance_Enceladus=sqrt((SVspacecraft(:,1)-SVofEnceladus(:,1)).^2+(SVspacecraft(:,2)-SVofEnceladus(:,2)).^2+(SVspacecraft(:,3)-SVofEnceladus(:,3)).^2)/pars.Moon.HillSph(2);

    %
    % % Create figure
    % figure1 = figure;
    % % Create axes
    % axes1 = axes('Parent',figure1);
    % hold(axes1,'on');
    % % Create multiple line objects using matrix input to semilogy
    % semilogy1 = semilogy(Tstep/3600/24,[distance_titan distance_Enceladus],'LineWidth',2);
    % set(semilogy1(1),'Color',[0 0 0], 'DisplayName','Titan Distance');
    % set(semilogy1(2),'Color',[1 0 0], 'DisplayName','Enceladus Distance');
    % % Create line
    % line([0 Tstep(end)/3600/24],[1 1],'Parent',axes1,'LineWidth',2,'LineStyle','--',...
    %     'Color',[0.313725501298904 0.313725501298904 0.313725501298904],'DisplayName','1 SOI');
    % % Uncomment the following line to preserve the Y-limits of the axes
    % ylim(axes1,[0.1 60]);
    % box(axes1,'on');
    % hold(axes1,'off');
    % % Set the remaining axes properties
    % set(axes1,'YMinorTick','on','YScale','log');
    % % Create ylabel
    % ylabel('Distance [Titan SOIs]');
    % % Create xlabel
    % xlabel('ToF [Days]');
    % % Create legend
    % legend1 = legend(axes1,'show');
    % set(legend1,...
    %     'Position',[0.612738090708142 0.124603171878391 0.281071433101382 0.12619047891526]);



    %--------------------------------------------------------------------------

end

end

%% Auxiliary Function
%--------------------------------------------------------------------------
% PlotSaturn_Equatorial_frame_withEnceladus_Titan_Orbits.m
function Plot_Saturn_Equatorial_frame_withEnceladus_Titan_Orbits(pars)
RS=getAstroConstants('Saturn','radius');
% 2. Plot Axes
scale = 1;
figure1 = figure('Color',[1 1 1]);
% Create axes
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YColor',[1 1 1],...
    'XColor',[1 1 1]);
hold(axes1,'all');
%--------------------------------------------------------------------------
[hfig] = drawPlanet('Saturn',[0 0 0],figure1,scale);
%--------------------------------------------------------------------------
axis equal
view(0,90)
%--------------------------------------------------------------------------
% Plot Enceladus Orbit as a circle
rEnceladus=pars.Moon.OrbRad(1);
theta=linspace(0,2*pi);
xE=rEnceladus*cos(theta);
yE=rEnceladus*sin(theta);
plot3(xE,yE,zeros(1,length(yE)),'LineWidth',1,...
    'Color',[0.380392163991928 0.380392163991928 0.380392163991928]);

rTitan=pars.Moon.OrbRad(2);
theta=linspace(0,2*pi);
xT=rTitan*cos(theta);
yT=rTitan*sin(theta);
plot3(xT,yT,zeros(1,length(yE)),'LineWidth',1,...
    'Color',[0.380392163991928 0.380392163991928 0.380392163991928]);
%--------------------------------------------------------------------------
% Plot Axis
%--------------------------------------------------------------------------
% X Axis
CMPos=[0 0 0];
p1=[10*RS 0 0];
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) p1(1)],[ CMPos(2) p1(2)],[ CMPos(3) p1(3)],'Color',[0 0 0])
%--------------------------------------------------------------------------
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
% Y axis
CMPos=[0 0 0];
py=[0 10*RS 0];
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) py(1)],[ CMPos(2) py(2)],[ CMPos(3) py(3)],'Color',[0 0 0])
% Arrow head Y
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

H=surf(x+py(1),y+py(2),z*h+py(3));
set(H,'FaceColor','k','EdgeColor','k')
P=py-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,py)
%--------------------------------------------------------------------------
% Z axis
CMPos=[0 0 0];
pz=[0 0  7*RS];
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) pz(1)],[ CMPos(2) pz(2)],[ CMPos(3) pz(3)],'Color',[0 0 0])
% Arrow head Y
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

H=surf(x+pz(1),y+pz(2),z*h+pz(3));
set(H,'FaceColor','k','EdgeColor','k')
P=pz-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,pz)
%--------------------------------------------------------------------------
end