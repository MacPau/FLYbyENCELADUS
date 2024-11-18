% function []=Res2Node(Resonance_data, )

% DESCRIPTION
% This function computes the fly-by node description and the two potential
% state vectors of any given resonance or pseudo-resonance.
%
% INPUT
% Resonance_data



% OUTPUT

CurrentLeg=FullNode;
% Departure point
TitanNODE_Arrival=CurrentLeg(9:11);
MTitan_now=CurrentLeg(12);
SVDeparture=SVTitan_func(MTitan_now);
r_moon_Departure=SVDeparture(1:3);
v_moon_Departure=SVDeparture(4:6);
% v_inf_minus_titan = vinfAlphaCrank2car(TitanNODE_minus, SVDeparture, pars.Planet.mu);

%--------------------------------------------------------------------------
% Preliminaries
Period_titan = 2*pi*sqrt(pars.Moon.OrbRad(2)^3/pars.Planet.mu);
[~, rr, vv] = vinfAlphaCrank2car(TitanNODE_Arrival, SVDeparture, pars.Planet.mu);
[kep_sc] = car2kep([rr vv], pars.Planet.mu);
Period_spacecraft = 2*pi*sqrt(kep_sc(1)^3/pars.Planet.mu);
%--------------------------------------------------------------------------
% are we inbound or outbound?
TypeFlag=81; % Outbound-Inbound
if mod(kep_sc(6),2*pi)>pi
    TypeFlag=18; % Inbound-outbound
end


% what is my current semimajor axis ratio?
Period_titan/Period_spacecraft;


% so, it seems as if 1:2 resonance with Titan is very close by.
% are we inbound or outbound? I am going Inbound

if TypeFlag==18

    for indexRes=1:length(PseudoRes_Titan_List(:,1))

        N=PseudoRes_Titan_List(indexRes,1);
        M=PseudoRes_Titan_List(indexRes,2);

        vinf=TitanNODE_Arrival(1);
        idmoon=5;
        idcentral=6;

        [vinf1, alpha1, crank1, vinf2, alpha2, crank2, tofsc] = wrap_pseudoResTransf(TypeFlag, N, M, vinf, idmoon, idcentral );

        [~, rr, vv] = vinfAlphaCrank2car([vinf1 alpha1 crank1], SVDeparture, pars.Planet.mu);
        [kep_sc_out] = car2kep([rr vv], pars.Planet.mu);
        Period_spacecraft_out = 2*pi*sqrt(kep_sc_out(1)^3/pars.Planet.mu);

        % tofsc/Period_titan;
        % tofsc/Period_spacecraft_out;

        % Period_titan/Period_spacecraft_out;
        % Can I calculate the DV for such a transition?

        if (~isnan(tofsc))&&(tofsc~=0)

            RelVel_in = vinfAlphaCrank2car(TitanNODE_Arrival, SVatTitan_LambertArrival, pars.Planet.mu);
            RelVel_out = vinfAlphaCrank2car([vinf1 alpha1 crank1], SVatTitan_LambertArrival, pars.Planet.mu);
            rpmin_flyby=pars.Moon.EquRad(2)+pars.Moon.hmin(2);
            [dv, alpha, alpha_A] = findDV(RelVel_in, RelVel_out, pars.Moon.mu(2), rpmin_flyby);
        else
            dv=1e99;
        end

        if dv<AcceptableDV_threshold

            previousLeg=CurrentLeg;
            %--------------------------------------------------------------------------
            % Here an Enceladus Titan Route is inserted
            NewLeg(1)= 5; % Departure Enceladus node

            TitanNODE_Departure = [vinf1 alpha1 crank1];
            NewLeg(2:4)= TitanNODE_Departure; % Node Departure at Enceladus
            NewLeg(5)= previousLeg(12); % mean anomaly Titan.
            NewLeg(6)= previousLeg(13); % mean anomaly Enceladus.
            NewLeg(7)= tofsc/3600/24; % Full period of one spacecraft revolution [days]

            NewLeg(8)= 5; % Arrival Titan node
            TitanNODE_Arrival_final = [vinf2, alpha2, crank2];

            NewLeg(9:11)= TitanNODE_Arrival_final; % Node Arrival Titan after Conexion Transfer
            NewLeg(12)= mod(NewLeg(5)+nTitan*NewLeg(7)*Days2Sec,2*pi); % mean anomaly Titan.
            NewLeg(13)= mod(NewLeg(6)+nEnceladus*NewLeg(7)*Days2Sec,2*pi);% mean anomaly Enceladus.

            AnatomyT2T=[1 1];
            if (mod(TitanNODE_Departure(3),2*pi)>pi/2)&(mod(TitanNODE_Departure(3),2*pi)<3*pi/2)
                AnatomyT2T(1)=-1;
            end

            if (mod(TitanNODE_Arrival_final(3),2*pi)>pi/2)&(mod(TitanNODE_Arrival_final(3),2*pi)<3*pi/2)
                AnatomyT2T(2)=-1;
            end

            Expanded_route_anatomy=[treeExploration{iN}.routeAnatomy; N M AnatomyT2T TitanNODE_Departure(1)];
            Expanded_route=[treeExploration{iN}.route; NewLeg];
            Expanded_routeDVCosts=[treeExploration{iN}.routeDVCosts; dv];
            treeExploration_temp{NewLayerCounter}.routeAnatomy=Expanded_route_anatomy;
            treeExploration_temp{NewLayerCounter}.route=Expanded_route;
            treeExploration_temp{NewLayerCounter}.routeDVCosts=Expanded_routeDVCosts;
            % Overall DV Cost of the Route
            treeExploration_temp{NewLayerCounter}.TotalDVCost=sum(Expanded_routeDVCosts);

            NewLayerCounter=NewLayerCounter+1;
        end
    end
elseif TypeFlag==81

    for indexRes=1:length(PseudoRes_Titan_List(:,1))

        N=PseudoRes_Titan_List(indexRes,1);
        M=PseudoRes_Titan_List(indexRes,2);

        vinf=TitanNODE_Arrival(1);
        idmoon=5;
        idcentral=6;

        [vinf1, alpha1, crank1, vinf2, alpha2, crank2, tofsc] = wrap_pseudoResTransf(TypeFlag, N, M, vinf, idmoon, idcentral );

        [~, rr, vv] = vinfAlphaCrank2car([vinf1 alpha1 crank1], SVDeparture, pars.Planet.mu);
        [kep_sc_out] = car2kep([rr vv], pars.Planet.mu);
        Period_spacecraft_out = 2*pi*sqrt(kep_sc_out(1)^3/pars.Planet.mu);


        % tofsc/Period_titan;
        % tofsc/Period_spacecraft_out;


        % Can I calculate the DV for such a transition?

        if (~isnan(tofsc))&&(tofsc~=0)

            RelVel_in = vinfAlphaCrank2car(TitanNODE_Arrival, SVatTitan_LambertArrival, pars.Planet.mu);
            RelVel_out = vinfAlphaCrank2car([vinf1 alpha1 crank1], SVatTitan_LambertArrival, pars.Planet.mu);
            rpmin_flyby=pars.Moon.EquRad(2)+pars.Moon.hmin(2);
            [dv, alpha, alpha_A] = findDV(RelVel_in, RelVel_out, pars.Moon.mu(2), rpmin_flyby);
        else
            dv=1e99;
        end


        if dv<AcceptableDV_threshold

            previousLeg=CurrentLeg;
            %--------------------------------------------------------------------------
            % Here an Enceladus Titan Route is inserted
            NewLeg(1)= 5; % Departure Enceladus node

            TitanNODE_Departure = [vinf1 alpha1 crank1];
            NewLeg(2:4)= TitanNODE_Departure; % Node Departure at Enceladus
            NewLeg(5)= previousLeg(12); % mean anomaly Titan.
            NewLeg(6)= previousLeg(13); % mean anomaly Enceladus.
            NewLeg(7)= tofsc/3600/24; % Full period of one spacecraft revolution [days]

            NewLeg(8)= 5; % Arrival Titan node
            TitanNODE_Arrival_final = [vinf2, alpha2, crank2];

            NewLeg(9:11)= TitanNODE_Arrival_final; % Node Arrival Titan after Conexion Transfer
            NewLeg(12)= mod(NewLeg(5)+nTitan*NewLeg(7)*Days2Sec,2*pi); % mean anomaly Titan.
            NewLeg(13)= mod(NewLeg(6)+nEnceladus*NewLeg(7)*Days2Sec,2*pi);% mean anomaly Enceladus.

            AnatomyT2T=[1 1];
            if (mod(TitanNODE_Departure(3),2*pi)>pi/2)&(mod(TitanNODE_Departure(3),2*pi)<3*pi/2)
                AnatomyT2T(1)=-1;
            end

            if (mod(TitanNODE_Arrival_final(3),2*pi)>pi/2)&(mod(TitanNODE_Arrival_final(3),2*pi)<3*pi/2)
                AnatomyT2T(2)=-1;
            end

            Expanded_route_anatomy=[treeExploration{iN}.routeAnatomy; N M AnatomyT2T TitanNODE_Departure(1)];
            Expanded_route=[treeExploration{iN}.route; NewLeg];
            Expanded_routeDVCosts=[treeExploration{iN}.routeDVCosts; dv];
            treeExploration_temp{NewLayerCounter}.routeAnatomy=Expanded_route_anatomy;
            treeExploration_temp{NewLayerCounter}.route=Expanded_route;
            treeExploration_temp{NewLayerCounter}.routeDVCosts=Expanded_routeDVCosts;
            % Overall DV Cost of the Route
            treeExploration_temp{NewLayerCounter}.TotalDVCost=sum(Expanded_routeDVCosts);
            NewLayerCounter=NewLayerCounter+1;
        end
    end
end