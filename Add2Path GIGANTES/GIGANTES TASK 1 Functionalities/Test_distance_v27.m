function [MinDistance]=Test_distance_v27(NextResonance,FullNode, pars)
% Tests
%--------------------------------------------------------------------------
% Epremeris Moons Set up
KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
SVTitan_func = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
SVEnceladus_func = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
%--------------------------------------------------------------------------
PeriodEnceladus=2*pi*sqrt(pars.Moon.OrbRad(1)^3/pars.Planet.mu);
nTitan=sqrt(pars.Planet.mu/pars.Moon.OrbRad(2)^3);
nEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);

% ToF_leg_indays=FullNode(7);
% TOF_seconds = ToF_leg_indays*3600*24; %[seconds]

TrA_Titan=FullNode(12); % end position
TrA_Enceladus=FullNode(13);

% Departure Position & Velocity
DepFlag=FullNode(1);
if DepFlag==5
    % Titan position & velocity
    SVatTitan=SVTitan_func(TrA_Titan);
    r_moon_departure=SVatTitan(1:3);
    v_moon_departure=SVatTitan(4:6);
elseif DepFlag==1
    % Enceladus position & velocity
    SVatEnceladus=SVEnceladus_func(TrA_Enceladus);
    r_moon_departure=SVatEnceladus(1:3);
    v_moon_departure=SVatEnceladus(4:6);
end

% Enceladus position & velocity
SVatEnceladus=SVEnceladus_func(TrA_Enceladus);
r_moon_arrival=SVatEnceladus(1:3);
v_moon_arrival=SVatEnceladus(4:6);

% Lambert Arc boundary conditions and time of flight
r1=r_moon_departure;
r2=r_moon_arrival;

%--------------------------------------------------------------------------

[vinfCAR_departure] = vinfAlphaCrank2car(FullNode(2:4), [r_moon_departure v_moon_departure], pars.Planet.mu);
v_sc_departure=vinfCAR_departure+v_moon_departure;

[vinfCAR_arrival] = vinfAlphaCrank2car(FullNode(9:11), [r_moon_arrival v_moon_arrival], pars.Planet.mu);
v_sc_arrival=vinfCAR_arrival+v_moon_arrival;

% norm(v_sc_arrival-v_sc_departure) % should be zero.

%--------------------------------------------------------------------------

% vinfTitan_temp=norm(v_sc_departure-v_moon_departure);
% vinfEnceladus_temp=norm(v_sc_arrival-v_moon_arrival);

% Compute distance with Titan
%--------------------------------------------------------------------------
% Propagate from departure from Titan and up to 30 days
% I should generate the state vector from the Lambert arc and from the
% resonance calculation to ensure that I obtain the same results.
% in the first instance, I will simply use the Lambert arc assuming that
% these are the same (which should be if everything was correct).

SV_Spacecraft=[r2 v_sc_arrival];
% OE_Spacecraft=car2kep(SV_Spacecraft,pars.Planet.mu);
% Propagation time


tf=NextResonance(1)*PeriodEnceladus; % [sec] A trial and error check will show that you need at

%Equation of motion ¨r+mu*r/R^3 = 0
F2BDyn=@(t,x)   [x(4); %dx/dt=Vx
    x(5); %dy/dt=Vy
    x(6); %dz/dt=Vz
    -pars.Planet.mu*x(1)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3); %dVx/dt
    -pars.Planet.mu*x(2)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3); %dVx/dt
    -pars.Planet.mu*x(3)/(sqrt(x(1)^2+x(2)^2+x(3)^2)^3)]; %dVx/dt

%Propagation
options=odeset('RelTol',1e-6,'AbsTol',1e-6);
[Tstep,SVspacecraft]=ode45(F2BDyn,[0:3600/2:tf],SV_Spacecraft,options);   % Any ODE solver should


%Propagation of Titan
options=odeset('RelTol',1e-6,'AbsTol',1e-6);
SV_Titan=SVTitan_func(TrA_Titan);
[~,SVoftitan]=ode45(F2BDyn,Tstep,SV_Titan,options);   % Any ODE solver should


%Propagation of Enceladus
options=odeset('RelTol',1e-6,'AbsTol',1e-6);
n_Enceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);
% Enceladus position & velocity
SV_Enceladus=SVEnceladus_func(TrA_Enceladus);
[~,SVofEnceladus]=ode45(F2BDyn,Tstep,SV_Enceladus,options);   % Any ODE solver should

% 
% % Plot reference frame
% run PlotSaturn_Equatorial_frame_withEnceladus_Titan_Orbits.m
% plot3(SVspacecraft(:,1),SVspacecraft(:,2),SVspacecraft(:,3))
% view(0,90)


distance_titan=sqrt((SVspacecraft(:,1)-SVoftitan(:,1)).^2+(SVspacecraft(:,2)-SVoftitan(:,2)).^2+(SVspacecraft(:,3)-SVoftitan(:,3)).^2)/pars.Moon.HillSph(2);
distance_Enceladus=sqrt((SVspacecraft(:,1)-SVofEnceladus(:,1)).^2+(SVspacecraft(:,2)-SVofEnceladus(:,2)).^2+(SVspacecraft(:,3)-SVofEnceladus(:,3)).^2)/pars.Moon.HillSph(2);


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
[MinDistance,indexMin]=min(distance_titan);
TimeMinDistance=Tstep(indexMin);
if min(distance_Enceladus)~=0
    error
end

end