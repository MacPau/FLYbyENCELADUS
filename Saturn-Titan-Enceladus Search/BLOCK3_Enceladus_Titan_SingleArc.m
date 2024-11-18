% BLOCK3_Enceladus_Titan_SingleArc.m
%--------------------------------------------------------------------------
% Store and plot Lambert arc

ToF_Leg_inSeconds_minimum=ArrivalDate(indexFeasible(iLA));

ArrivalMeanNomaly=mod(Mtitan_atstartofleg+ToF_Leg_inSeconds_minimum*nTitan,2*pi);
SVatTitan_LambertArrival=SVTitan_func(ArrivalMeanNomaly);
rAT=SVatTitan_LambertArrival(1:3);
vAT=SVatTitan_LambertArrival(4:6);

%--------------------------------------------------------------------------
% Multi-Rev Lambert now
LambertCases=LambertParam(indexFeasible(iLA),:);
%--------------------------------------------------------------------------
orbitType=0;
Nrev=LambertCases(1);
Ncase=LambertCases(2);
optionsLMR=0;
[smaLambert,pLambert,eccLambert,ERROR,v1,v2] = lambertMR(rDE,rAT,ToF_Leg_inSeconds_minimum,pars.Planet.mu,orbitType,Nrev,Ncase,optionsLMR);
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
vinf_minus_Titan=v2-vAT;
% %--------------------------------------------------------------------------
% % Maximum Deflection for final Fly-by
% rp_min=pars.Moon.EquRad(2)+pars.Moon.hmin(2);
% ecchyperbolic=1+rp_min*norm(vinf_minus_Titan)^2/pars.Moon.mu(2);
% DeltaMax_Titan=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby
% %---------------------------------------------------------------------------
% Anatomy of Transfers
[KepSCDeparture]=car2kep([rDE v1],pars.Planet.mu);
[KepSCArrival]=car2kep([rAT v2],pars.Planet.mu);
AnatomyE2T=[1 1];
if KepSCDeparture(6)>pi
    AnatomyE2T(1)=-1;
end
if KepSCArrival(6)>pi
    AnatomyE2T(2)=-1;
end

%% Store New Lambert Arc Leg

%--------------------------------------------------------------------------
% Here an Enceladus Titan Route is inserted
NewLeg(1)= 1; % Departure Enceladus node

EncNODE_Departure = car2vinfAlphaCrank(vinf_plus_departure, [rDE vDE]);
NewLeg(2:4)= EncNODE_Departure; % Node Departure at Enceladus
NewLeg(5)= Mtitan_atstartofleg; % mean anomaly Titan.
NewLeg(6)= MEnceladus; % mean anomaly Enceladus.
NewLeg(7)= ToF_Leg_inSeconds_minimum/3600/24; % Full period of one spacecraft revolution [days]

NewLeg(8)= 5; % Arrival Titan node
TitanNODE_Arrival = car2vinfAlphaCrank(vinf_minus_Titan, SVatTitan_LambertArrival);
NewLeg(9:11)= TitanNODE_Arrival; % Node Arrival Titan after Conexion Transfer

NewLeg(12)= mod(NewLeg(5)+nTitan*NewLeg(7)*Days2Sec,2*pi); % mean anomaly Titan.
NewLeg(13)= mod(NewLeg(6)+nEnceladus*NewLeg(7)*Days2Sec,2*pi);% mean anomaly Enceladus.
