%% BLOCK3_Enceladus_Titan_LambertArc.m
%--------------------------------------------------------------------------
% Constants & Parameters
RAD=pi/180;
Days2Sec=3600*24;
nEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);
nTitan=sqrt(pars.Planet.mu/pars.Moon.OrbRad(2)^3);
%--------------------------------------------------------------------------
% Epremeris Moons Set up
KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
SVTitan_func = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
SVEnceladus_func = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
%--------------------------------------------------------------------------
% Unwrap the FullNode information.
% MTitan=Mtitan_atstart; % initial true anomaly of Titan.
% 
% ToF_leg_indays=FullNode(7)+ToF2LastNode; %TBR

Mtitan_atstartofleg=FullNode(12);
MEnceladus=FullNode(13);

EnceladusNode_info=FullNode(9:11);
%----------------------------------------------------------------------
% Detection of Close approaches in the next E-E path
% Is the Departure Node I or O? if kepSCatENceladus(6)<pi is O, if >pi
% is I.
SVatEnceladus=SVEnceladus_func(MEnceladus);
rMoonDeparture=SVatEnceladus(1:3);
vMoonDeparture=SVatEnceladus(4:6);

vinfEnceladus_vec=vinfAlphaCrank2car(EnceladusNode_info, SVatEnceladus, pars.Planet.mu);

vSC=vinfEnceladus_vec+vMoonDeparture;
[kepSCatENceladus]=car2kep([rMoonDeparture vSC],pars.Planet.mu);

%----------------------------------------------------------------------
% Check position of Titan in the following two closes approaches. There
% must be one option for inward departures and another one for outward
% departures.
ecc_sc=kepSCatENceladus(2);
sma_sc=kepSCatENceladus(1);
n_sc=sqrt(pars.Planet.mu/sma_sc^3);

E1=acos(1/ecc_sc-pars.Moon.OrbRad(1)/sma_sc/ecc_sc);
M1=E1-ecc_sc*sin(E1);
t1=M1/n_sc;

E2=acos(1/ecc_sc-pars.Moon.OrbRad(2)/sma_sc/ecc_sc);
M2=E2-ecc_sc*sin(E2);
t2=M2/n_sc;

Tsc=2*pi*sqrt(sma_sc^3/pars.Planet.mu);

if kepSCatENceladus(6)>=pi
    %--------------------------------------------------------------------------
    % Inward Departure Option
    % ToF_nextCrossing=
    TwoTimeofFlight=[t2+t1 Tsc-t2+t1];
elseif kepSCatENceladus(6)<pi
    %--------------------------------------------------------------------------
    % Outward Departure Option
    % ToF_nextCrossing=
    TwoTimeofFlight=[t2-t1 Tsc-t2-t1];
end

% TargetingTimeVector=[-3:0.01:3]*3600*24;

ArrivalDate=[TwoTimeofFlight(1)-Days2Sec:3600/2:TwoTimeofFlight(2)+Days2Sec];

%--------------------------------------------------------------------------
% Maximum Deflection for Initial Fly-by
rp_min=pars.Moon.EquRad(1)+pars.Moon.hmin(1);
ecchyperbolic=1+rp_min*norm(vinfEnceladus_vec)^2/pars.Moon.mu(1);
DeltaMax_Enceladus=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby
%--------------------------------------------------------------------------
rDE=rMoonDeparture;
vDE=vMoonDeparture;
dvVector=zeros(1,length(ArrivalDate));
LambertParam_rev0=zeros(length(ArrivalDate),2);
for iL=1:length(ArrivalDate)

    % This position must be the same as r1 if resonant.
    ArrivalDate_corrected4Lambert=ArrivalDate(iL);
    ArrivalMeanNomaly=mod(Mtitan_atstartofleg+ArrivalDate_corrected4Lambert*nTitan,2*pi);
    SVatTitan_LambertArrival=SVTitan_func(ArrivalMeanNomaly);
    rAT=SVatTitan_LambertArrival(1:3);
    vAT=SVatTitan_LambertArrival(4:6);

    ToF_Leg_inSeconds=ArrivalDate_corrected4Lambert;
    % Titan position & velocity
    orbitType=0;
    Nrev=0;
    Ncase=0;
    optionsLMR=0;
    [smaLambert,pLambert,eccLambert,ERROR,v1,v2] = lambertMR(rDE,rAT,ToF_Leg_inSeconds,pars.Planet.mu,orbitType,Nrev,Ncase,optionsLMR);
% OK, there is an issue with the Lambert arc, out of plane transfers will
% need a somewhat more sofisticated lambert arc to be used, otherwise the
% lambert arc will always give the inplane transfer, since this is a 180
% degrees transfer. 
    %-------------------------------------------------------------------------
    % DV Defect Calculation for the departure Assymptote
    vinf_minus_departure=vinfEnceladus_vec;
    vinf_plus_departure=v1-vDE;
    % necessary deflection
    delta_at_Enceladus=acos(dot(vinf_minus_departure,vinf_plus_departure)/norm(vinf_minus_departure)/norm(vinf_plus_departure));
    if delta_at_Enceladus<=DeltaMax_Enceladus
        DV_at_Enceladus=abs(norm(vinf_plus_departure)-norm(vinf_minus_departure));
    else
        DV_at_Enceladus=sqrt(norm(vinf_plus_departure)^2+norm(vinf_minus_departure)^2-2*norm(vinf_plus_departure)*norm(vinf_minus_departure)*cos(DeltaMax_Enceladus-delta_at_Enceladus));
    end

    % %--------------------------------------------------------------------------
    % %  DV Defect for Titan
    % vinf_minus_Titan=v2-vAT;
    % %--------------------------------------------------------------------------
    % % Maximum Deflection for final Fly-by
    % rp_min=pars.Moon.EquRad(2)+pars.Moon.hmin(2);
    % ecchyperbolic=1+rp_min*norm(vinf_minus_Titan)^2/pars.Moon.mu(2);
    % DeltaMax_Titan=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby

    % how it would be the vector for the same resonance?
    % [rp_res,ra_res,E_res, alfa_res, i_max, trA_moonX] = Resonance_matching(Resonances,V_inf, a_moon, mu_planet, pars)
    dvVector(iL)=DV_at_Enceladus;

end

%--------------------------------------------------------------------------
% One complete revolution for those transfers with two revolutions of the
% spacecraft
dvVector_rev1=[];
ArrivalDate2=[];
LambertParam_rev1=[];
if NextResonance(2)==2
    ArrivalDate2=[ArrivalDate+Tsc];

dvVector_rev1=zeros(1,length(ArrivalDate2));
LambertParam_rev1=zeros(length(ArrivalDate2),2);
for iL=1:length(ArrivalDate2)

    % This position must be the same as r1 if resonant.
    ArrivalDate_corrected4Lambert=ArrivalDate2(iL);
    ArrivalMeanNomaly=mod(Mtitan_atstartofleg+ArrivalDate_corrected4Lambert*nTitan,2*pi);
    SVatTitan_LambertArrival=SVTitan_func(ArrivalMeanNomaly);
    rAT=SVatTitan_LambertArrival(1:3);
    vAT=SVatTitan_LambertArrival(4:6);

    ToF_Leg_inSeconds=ArrivalDate_corrected4Lambert;
%--------------------------------------------------------------------------
% Two cases Lambert
    % Titan position & velocity
    orbitType=0;
    Nrev=1;
    Ncase=0;
    optionsLMR=0;
    [smaLambert,pLambert,eccLambert,ERROR,v1,v2] = lambertMR(rDE,rAT,ToF_Leg_inSeconds,pars.Planet.mu,orbitType,Nrev,Ncase,optionsLMR);
    %-------------------------------------------------------------------------
    % DV Defect Calculation for the departure Assymptote
    vinf_minus_departure=vinfEnceladus_vec;
    vinf_plus_departure=v1-vDE;
    % necessary deflection
    delta_at_Enceladus=acos(dot(vinf_minus_departure,vinf_plus_departure)/norm(vinf_minus_departure)/norm(vinf_plus_departure));
    if delta_at_Enceladus<=DeltaMax_Enceladus
        DV_at_Enceladus1=abs(norm(vinf_plus_departure)-norm(vinf_minus_departure));
    else
        DV_at_Enceladus1=sqrt(norm(vinf_plus_departure)^2+norm(vinf_minus_departure)^2-2*norm(vinf_plus_departure)*norm(vinf_minus_departure)*cos(DeltaMax_Enceladus-delta_at_Enceladus));
    end
    %--------------------------------------------------------------------------
    %  DV Defect for Titan
    vinf_minus_Titan1=v2-vAT;

% Two cases Lambert
    % Titan position & velocity
    orbitType=0;
    Nrev=1;
    Ncase=1;
    optionsLMR=0;
    [smaLambert,pLambert,eccLambert,ERROR,v1,v2] = lambertMR(rDE,rAT,ToF_Leg_inSeconds,pars.Planet.mu,orbitType,Nrev,Ncase,optionsLMR);
    %-------------------------------------------------------------------------
    % DV Defect Calculation for the departure Assymptote
    vinf_minus_departure=vinfEnceladus_vec;
    vinf_plus_departure=v1-vDE;
    % necessary deflection
    delta_at_Enceladus=acos(dot(vinf_minus_departure,vinf_plus_departure)/norm(vinf_minus_departure)/norm(vinf_plus_departure));
    if delta_at_Enceladus<=DeltaMax_Enceladus
        DV_at_Enceladus2=abs(norm(vinf_plus_departure)-norm(vinf_minus_departure));
    else
        DV_at_Enceladus2=sqrt(norm(vinf_plus_departure)^2+norm(vinf_minus_departure)^2-2*norm(vinf_plus_departure)*norm(vinf_minus_departure)*cos(DeltaMax_Enceladus-delta_at_Enceladus));
    end
    % %--------------------------------------------------------------------------
    % %  DV Defect for Titan
    % vinf_minus_Titan2=v2-vAT;
    % 
    % 
    % %--------------------------------------------------------------------------
    % % Maximum Deflection for final Fly-by
    % rp_min=pars.Moon.EquRad(2)+pars.Moon.hmin(2);
    % ecchyperbolic=1+rp_min*norm(vinf_minus_Titan)^2/pars.Moon.mu(2);
    % DeltaMax_Titan=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby

    % how it would be the vector for the same resonance?
    % [rp_res,ra_res,E_res, alfa_res, i_max, trA_moonX] = Resonance_matching(Resonances,V_inf, a_moon, mu_planet, pars)
    [Value, index]=min([DV_at_Enceladus1 DV_at_Enceladus2]);
    dvVector_rev1(iL)=min([DV_at_Enceladus1 DV_at_Enceladus2]);
    LambertParam_rev1(iL,:)=[Nrev index-1];
end

end

ArrivalDate=[ArrivalDate ArrivalDate2];
dvVector=[dvVector dvVector_rev1];
LambertParam=[LambertParam_rev0; LambertParam_rev1];

% 
% Create figure
% figure1 = figure;
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% % Create plot
% plot(ArrivalDate/3600/24,dvVector);
% % Create ylabel
% ylabel('DV [km/s]');
% % Create xlabel
% xlabel('TimeofFlight [days]');
% line([TwoTimeofFlight(1)/3600/24 TwoTimeofFlight(1)/3600/24],[0 2])
% line([TwoTimeofFlight(2)/3600/24 TwoTimeofFlight(2)/3600/24],[0 2])
% line(xlim,[AcceptableDV_threshold AcceptableDV_threshold])
% ylim(axes1,[0 2]);
% box(axes1,'on');
% hold(axes1,'off');


indexFeasible=find(dvVector<AcceptableDV_threshold);
