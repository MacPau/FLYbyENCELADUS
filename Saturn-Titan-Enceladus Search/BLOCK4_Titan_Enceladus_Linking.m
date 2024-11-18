% Block4_Titan_Enceladus_linking.m
% The search attemps to find new links with the set of Enceladus Nodes in
% the database. Note however that the true latitude of the nodes is
% now different of the original. 
for i1=1:length(BLOCK1nodes)

    nodenumber=i1;
    NodeTarget=BLOCK1nodes{nodenumber}.NodeInit_standard;

    CurrentLeg=FullNode;

    % Main Functions to Encapsulate
    %--------------------------------------------------------------------------
    % Constants & Parameters Initialization
    RAD=pi/180;
    Days2Sec=3600*24;
    nEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);
    nTitan=sqrt(pars.Planet.mu/pars.Moon.OrbRad(2)^3);
    %--------------------------------------------------------------------------
    % Epremeris Moons Set up (Circular Co-planar)
    KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
    KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
    SVTitan_func = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
    SVEnceladus_func = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
    %--------------------------------------------------------------------------

    toF_transfer_seconds=NodeTarget(7)*Days2Sec; % I should allow for a small variability in ToF, to be decided how much, but provably a few hours up and down, probably using fminbound.

    % toF_transfer_bnds= [toF_transfer_seconds-Days2Sec/2 toF_transfer_seconds+Days2Sec/2];
    fun_ToFAdjustement=@(x)Adjusting_TOF_forLAMBERT(x,CurrentLeg,NodeTarget,pars);

    FixedTOF_DV=fun_ToFAdjustement(toF_transfer_seconds);

    DVLimit_threshold_to_attempt_reoptimization=0.5;
    if FixedTOF_DV<DVLimit_threshold_to_attempt_reoptimization
        %[ToF_adjustement,DVAdjusted,exitflag] = fminbnd(fun_ToFAdjustement,toF_transfer_bnds(1),toF_transfer_bnds(2));
        [ToF_adjustement,DVAdjusted,exitflag] = fminsearch(fun_ToFAdjustement,toF_transfer_seconds);
    else
        ToF_adjustement=toF_transfer_seconds;
    end

    % Departure point
    TitanNODE_minus=CurrentLeg(9:11);
    MTitan_now=CurrentLeg(12);
    SVDeparture=SVTitan_func(MTitan_now);
    r_moon_Departure=SVDeparture(1:3);
    v_moon_Departure=SVDeparture(4:6);
    v_inf_minus_titan = vinfAlphaCrank2car(TitanNODE_minus, SVDeparture, pars.Planet.mu);

    % Arrival Point
    EnceladusNODE_plus=NodeTarget(9:11);
    MEnceladus_now=CurrentLeg(13);
    MEnceladus_Arrival=mod(MEnceladus_now+nEnceladus*ToF_adjustement,2*pi);
    SVArrival=SVEnceladus_func(MEnceladus_Arrival);
    r_moon_Arrival=SVArrival(1:3);
    v_moon_Arrival=SVArrival(4:6);
    v_inf_plus_enceladus=vinfAlphaCrank2car(EnceladusNODE_plus, SVArrival, pars.Planet.mu);


    % Lambert Arc Linkage Generation
    orbitType=0;
    Nrev=0;
    Ncase=0;
    optionsLMR=0;
    [smaLambert,pLambert,eccLambert,ERROR,v1,v2] = lambertMR(r_moon_Departure,r_moon_Arrival,ToF_adjustement,pars.Planet.mu,orbitType,Nrev,Ncase,optionsLMR);

    % Calculation of the DV Defect at Departure

    %--------------------------------------------------------------------------
    % Maximum Deflection for Departure Fly-by
    rp_min=pars.Moon.EquRad(2)+pars.Moon.hmin(2);
    ecchyperbolic=1+rp_min*norm(TitanNODE_minus(1))^2/pars.Moon.mu(2);
    DeltaMax_Titan=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby

    %--------------------------------------------------------------------------
    % Calculation of dv_defect for closest match possible (NOTE it may not be necessarily feasible)
    v_inf_plus_titan=v1-v_moon_Departure;

    % necessary deflection
    delta_at_Titan=acos(dot(v_inf_minus_titan,v_inf_plus_titan)/norm(v_inf_minus_titan)/norm(v_inf_plus_titan));
    if delta_at_Titan<=DeltaMax_Titan
        DV_at_Titan=abs(norm(v_inf_plus_titan)-norm(v_inf_minus_titan));
    else
        DV_at_Titan=sqrt(norm(v_inf_plus_titan)^2+norm(v_inf_minus_titan)^2-2*norm(v_inf_plus_titan)*norm(v_inf_minus_titan)*cos(DeltaMax_Titan-delta_at_Titan));
    end

    % Calculation of the DV Defect at Arrival
    %--------------------------------------------------------------------------
    % Maximum Deflection for Departure Fly-by
    rp_min=pars.Moon.EquRad(1)+pars.Moon.hmin(1);
    ecchyperbolic=1+rp_min*norm(EnceladusNODE_plus(1))^2/pars.Moon.mu(1);
    DeltaMax_Enceladus=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby

    %--------------------------------------------------------------------------
    % Calculation of dv_defect for closest match possible (NOTE it may not be necessarily feasible)
    v_inf_minus_enceladus=v2-v_moon_Arrival;

    % necessary deflection
    delta_at_Enceladus=acos(dot(v_inf_minus_enceladus,v_inf_plus_enceladus)/norm(v_inf_minus_enceladus)/norm(v_inf_plus_enceladus));
    if delta_at_Enceladus<=DeltaMax_Enceladus
        DV_at_Enceladus=abs(norm(v_inf_plus_enceladus)-norm(v_inf_minus_enceladus));
    else
        DV_at_Enceladus=sqrt(norm(v_inf_plus_enceladus)^2+norm(v_inf_minus_enceladus)^2-2*norm(v_inf_plus_enceladus)*norm(v_inf_minus_enceladus)*cos(DeltaMax_Enceladus-delta_at_Enceladus));
    end

    %--------------------------------------------------------------------------
    TotalDV=DV_at_Titan+DV_at_Enceladus;

    %--------------------------------------------------------------------------
    if TotalDV<AcceptableDV_threshold
        % disp(['Fixed TOF DV is ',num2str(FixedTOF_DV),' km/s, while variable TOF DV is ', num2str(DVAdjusted),' km/s'])
        % Define NewLeg
        %--------------------------------------------------------------------------
        % Here an Titan Enceladus Route is inserted
        NewLeg(1)= 5; % Departure Titan node

        TitanNODE_Departure = car2vinfAlphaCrank(v_inf_plus_titan, [r_moon_Departure v_moon_Departure]);
        NewLeg(2:4)= TitanNODE_Departure; % Node Departure at Titan

        NewLeg(5)=MTitan_now; % mean anomaly of Titan
        NewLeg(6)=MEnceladus_now; % mean anomaly of Enceladus

        NewLeg(7)= ToF_adjustement/3600/24; % Time of Transfer to next intersection, with Enceladus.

        NewLeg(8)= 1; % Arrival Enceladus node

        EncelNODE_Arrival = car2vinfAlphaCrank(v_inf_minus_enceladus, [r_moon_Arrival v_moon_Arrival]);
        NewLeg(9:11)= EncelNODE_Arrival; % Enceladus Arrival Node
        NewLeg(12)= mod(NewLeg(5)+nTitan*NewLeg(7)*Days2Sec,2*pi); % mean anomaly Titan.
        NewLeg(13)= mod(NewLeg(6)+nEnceladus*NewLeg(7)*Days2Sec,2*pi);% mean anomaly Enceladus.
        AnatomyT2E=[1 1];
        if (mod(TitanNODE_Departure(3),2*pi)>pi/2)&(mod(TitanNODE_Departure(3),2*pi)<3*pi/2)
            AnatomyT2E(1)=-1;
        end

        if (mod(EncelNODE_Arrival(3),2*pi)>pi/2)&(mod(EncelNODE_Arrival(3),2*pi)<3*pi/2)
            AnatomyT2E(2)=-1;
        end

        Expanded_route_anatomy=[treeExploration{iN}.routeAnatomy; 0 0 AnatomyT2E TitanNODE_Departure(1); BLOCK1nodes{nodenumber}.ResAnatomyTargeted];
        Expanded_route=[treeExploration{iN}.route; NewLeg];
        Expanded_routeDVCosts=[treeExploration{iN}.routeDVCosts; TotalDV];
        treeExploration_temp{NewLayerCounter}.routeAnatomy=Expanded_route_anatomy;
        treeExploration_temp{NewLayerCounter}.route=Expanded_route;
        treeExploration_temp{NewLayerCounter}.routeDVCosts=Expanded_routeDVCosts;
        % Overall DV Cost of the Route
        treeExploration_temp{NewLayerCounter}.TotalDVCost=sum(Expanded_routeDVCosts);

        NewLayerCounter=NewLayerCounter+1;
    end

end

% [minval,minIndex]=min(TotalDV)



%% AUXILIARY FUNCTIONS

function [TotalDV]=Adjusting_TOF_forLAMBERT(toF_transfer_seconds,CurrentLeg,NodeTarget,pars)

% Main Functions to Encapsulate
%--------------------------------------------------------------------------
% Constants & Parameters Initialization
% RAD=pi/180;
% Days2Sec=3600*24;
nEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);
% nTitan=sqrt(pars.Planet.mu/pars.Moon.OrbRad(2)^3);
%--------------------------------------------------------------------------
% Epremeris Moons Set up (Circular Co-planar)
KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
SVTitan_func = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
SVEnceladus_func = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
%--------------------------------------------------------------------------
% Departure point
TitanNODE_minus=CurrentLeg(9:11);
MTitan_now=CurrentLeg(12);
SVDeparture=SVTitan_func(MTitan_now);
r_moon_Departure=SVDeparture(1:3);
v_moon_Departure=SVDeparture(4:6);
v_inf_minus_titan = vinfAlphaCrank2car(TitanNODE_minus, SVDeparture, pars.Planet.mu);

% Arrival Point
EnceladusNODE_plus=NodeTarget(9:11);
MEnceladus_now=CurrentLeg(13);
MEnceladus_Arrival=mod(MEnceladus_now+nEnceladus*toF_transfer_seconds,2*pi);
SVArrival=SVEnceladus_func(MEnceladus_Arrival);
r_moon_Arrival=SVArrival(1:3);
v_moon_Arrival=SVArrival(4:6);
v_inf_plus_enceladus=vinfAlphaCrank2car(EnceladusNODE_plus, SVArrival, pars.Planet.mu);
% Lambert Arc Linkage Generation
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;
[~,~,~,~,v1,v2] = lambertMR(r_moon_Departure,r_moon_Arrival,toF_transfer_seconds,pars.Planet.mu,orbitType,Nrev,Ncase,optionsLMR);

% Calculation of the DV Defect at Departure
%--------------------------------------------------------------------------
% Maximum Deflection for Departure Fly-by
rp_min=pars.Moon.EquRad(2)+pars.Moon.hmin(2);
ecchyperbolic=1+rp_min*norm(TitanNODE_minus(1))^2/pars.Moon.mu(2);
DeltaMax_Titan=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby

%--------------------------------------------------------------------------
% Calculation of dv_defect for closest match possible (NOTE it may not be necessarily feasible)
v_inf_plus_titan=v1-v_moon_Departure;

% necessary deflection
delta_at_Titan=acos(dot(v_inf_minus_titan,v_inf_plus_titan)/norm(v_inf_minus_titan)/norm(v_inf_plus_titan));
if delta_at_Titan<=DeltaMax_Titan
    DV_at_Titan=abs(norm(v_inf_plus_titan)-norm(v_inf_minus_titan));
else
    DV_at_Titan=sqrt(norm(v_inf_plus_titan)^2+norm(v_inf_minus_titan)^2-2*norm(v_inf_plus_titan)*norm(v_inf_minus_titan)*cos(DeltaMax_Titan-delta_at_Titan));
end

% Calculation of the DV Defect at Arrival

%--------------------------------------------------------------------------
% Maximum Deflection for Departure Fly-by
rp_min=pars.Moon.EquRad(1)+pars.Moon.hmin(1);
ecchyperbolic=1+rp_min*norm(EnceladusNODE_plus(1))^2/pars.Moon.mu(1);
DeltaMax_Enceladus=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby

%--------------------------------------------------------------------------
% Calculation of dv_defect for closest match possible (NOTE it may not be necessarily feasible)
v_inf_minus_enceladus=v2-v_moon_Arrival;

% necessary deflection
delta_at_Enceladus=acos(dot(v_inf_minus_enceladus,v_inf_plus_enceladus)/norm(v_inf_minus_enceladus)/norm(v_inf_plus_enceladus));
if delta_at_Enceladus<=DeltaMax_Enceladus
    DV_at_Enceladus=abs(norm(v_inf_plus_enceladus)-norm(v_inf_minus_enceladus));
else
    DV_at_Enceladus=sqrt(norm(v_inf_plus_enceladus)^2+norm(v_inf_minus_enceladus)^2-2*norm(v_inf_plus_enceladus)*norm(v_inf_minus_enceladus)*cos(DeltaMax_Enceladus-delta_at_Enceladus));
end

%--------------------------------------------------------------------------
TotalDV=DV_at_Titan+DV_at_Enceladus;

end


