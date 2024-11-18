%--------------------------------------------------------------------------
%% Block1_Generation_NODE1.m

% What is this block doing? 
% this block must compute all possible Enceladus nodes used in the search
% and the linked Titan nodes that define the phase between both moons,
% considering the first an encounter on its first intersection. 

nodenumber=1;
for jR=1:length(Resonances_Enceladus(:,1)) % loop over defined resonances

    for i=1:length(vInfEnceladus_list)  % loop over the v_inf List

        vinfEnceladus=vInfEnceladus_list(i);
        VelEnceladus   = sqrt(pars.Planet.mu/pars.Moon.OrbRad(1));           %[km/s] Moon Orbital velocity
        [rp_res,ra_res,E_res, alfa_res, ~, trA_moonCrossing] = Tisserand_Resonances(Resonances_Enceladus(jR,:),vinfEnceladus, pars.Moon.OrbRad(1),pars.Planet.mu, pars);

% Since this is block 1, we seek nodes that intersect with enceladus,
% departing from Titan. 

% generating Titan Nodes with crossing at enceladus
        if (rp_res<pars.Moon.OrbRad(1)/pars.Planet.EquRad)&(ra_res>pars.Moon.OrbRad(2)/pars.Planet.EquRad) % is there an intersection to ensure a Moon to Moon transfer? if so, follow on. 

            trA_Enceladus=[pi-trA_moonCrossing pi+trA_moonCrossing]; % [RAD]

            % Four Possible Time of Flights
            ecc_res=(ra_res-rp_res)/(ra_res+rp_res);
            sma_res=(ra_res+rp_res)/2*pars.Planet.EquRad;
            n_res=sqrt(pars.Planet.mu/sma_res^3);

            % crossing point at lower altitude
            theta1=acos(sma_res/pars.Moon.OrbRad(1)/ecc_res*(1-ecc_res^2)-1/ecc_res);
            E1=acos(1/ecc_res-pars.Moon.OrbRad(1)/sma_res/ecc_res);
            M1=E1-ecc_res*sin(E1);
            t1=M1/n_res;

            % crossing point at higher altitude
            theta2=acos(sma_res/pars.Moon.OrbRad(2)/ecc_res*(1-ecc_res^2)-1/ecc_res);
            E2=acos(1/ecc_res-pars.Moon.OrbRad(2)/sma_res/ecc_res);
            M2=E2-ecc_res*sin(E2);
            t2=M2/n_res;

            Tsc=2*pi*sqrt(sma_res^3/pars.Planet.mu);

            % FourTimeofFlight=[t2-t1 t2+t1 Tsc-t2-t1 Tsc-t2+t1];

            % Finish the porkchop plot with fixed resonances as in the tisserand with
            % k=0.
            %--------------------------------------------------------------
            % Compute analytically state vectors at the beginning and at
            % the end of the transfer. Reference frame is defined such that
            % perifocal and inertial reference frame correspond at t=0.
            %--------------------------------------------------------------
            % Beginning Transfer - Outward Position
            Theta_departure=theta2; % watch out here, it is either I or O and then this changes.
            smallOmega_res=-Theta_departure; % To keep true latitude fixed.

            h=sqrt(pars.Planet.mu*sma_res*(1-ecc_res^2));

            r1_norm_analytic=h^2/pars.Planet.mu/(1+ecc_res*cos(Theta_departure));

            r1_analytic(1)=r1_norm_analytic*cos(Theta_departure);
            r1_analytic(2)=r1_norm_analytic*sin(Theta_departure);
            r1_analytic(3)=0;

            v1_analytic(1)=-pars.Planet.mu/h*sin(Theta_departure);
            v1_analytic(2)=pars.Planet.mu/h*(ecc_res+cos(Theta_departure));
            v1_analytic(3)=0;

            A2_smallOmega=[cos(-smallOmega_res) sin(-smallOmega_res) 0;...
                           -sin(-smallOmega_res) cos(-smallOmega_res) 0;
                           0 0 1];

            r1_analytic_O=A2_smallOmega*r1_analytic';
            v1_analytic_O=A2_smallOmega*v1_analytic';
            %--------------------------------------------------------------
            % Beginning Transfer - Inward Position
            trLatitude_departure=0;
            Theta_departure=2*pi-theta2; % watch out here, it is either I or O and then this changes.
            smallOmega_res=-Theta_departure; % To keep true latitude fixed.

            h=sqrt(pars.Planet.mu*sma_res*(1-ecc_res^2));

            r1_norm_analytic=h^2/pars.Planet.mu/(1+ecc_res*cos(Theta_departure));

            r1_analytic(1)=r1_norm_analytic*cos(Theta_departure);
            r1_analytic(2)=r1_norm_analytic*sin(Theta_departure);
            r1_analytic(3)=0;

            v1_analytic(1)=-pars.Planet.mu/h*sin(Theta_departure);
            v1_analytic(2)=pars.Planet.mu/h*(ecc_res+cos(Theta_departure));
            v1_analytic(3)=0;

            A2_smallOmega=[cos(-smallOmega_res) sin(-smallOmega_res) 0;...
                           -sin(-smallOmega_res) cos(-smallOmega_res) 0;
                           0 0 1];

            r1_analytic_I=A2_smallOmega*r1_analytic';
            v1_analytic_I=A2_smallOmega*v1_analytic';            
%--------------------------------------------------------------------------
% Summary of Transfer solutions [p q delta_theta delta_time state_vector_departure]
            FourOptions=[-1 -1 (2*pi-theta1)-(2*pi-theta2)  t2-t1 r1_analytic_I' v1_analytic_I';... % I-I Downtransfer
                          1 -1 (2*pi-theta1)-theta2  Tsc-t2-t1 r1_analytic_O' v1_analytic_O';...    % O-I Downtransfer
                         -1 1 theta2+theta1  t2+t1 r1_analytic_I' v1_analytic_I';...               % I-O Downtransfer
                          1 1 2*pi-(theta2-theta1)  Tsc-t2+t1 r1_analytic_O' v1_analytic_O'];      % O-O Downtransfer
%--------------------------------------------------------------------------
% Lambert arc to compute the v1 and v2 of the trajectory. This should be
% substituted by other than the Lambert arc.
%--------------------------------------------------------------------------
% Epremeris Moons Set up
            KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
            KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
            SVTitan_func = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
            SVEnceladus_func = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
%--------------------------------------------------------------------------
% Lambert Arc
            TrA_Titan=0; % Titan is initialized as at zero true anomaly at Epoch 0

            for iL=1:4
                p=FourOptions(iL,2);
                TrA_Enceladus=FourOptions(iL,3);
                TOF_seconds = FourOptions(iL,4); %[seconds]
                % Titan position & velocity
                SVatTitan=SVTitan_func(TrA_Titan);
                rTitan=SVatTitan(1:3);
                vTitan=SVatTitan(4:6);

                % Enceladus position & velocity
                SVatEnceladus=SVEnceladus_func(TrA_Enceladus);
                rEnceladus=SVatEnceladus(1:3);
                vEnceladus=SVatEnceladus(4:6);

                % Lambert Arc boundary conditions and time of flight
                r1=rTitan;
                r2=rEnceladus;

                v1=FourOptions(iL,8:10);
                [sv2,error] = kepPro([r1 v1],TOF_seconds,pars.Planet.mu);
                % check
                if norm(sv2(1:3)-r2)>1
                    msg='INTERSECTION AT ENCELADUS IS NOT CORRECT';
                    error(msg)
                else
                    v2=sv2(4:6);
                end
                

                vinfTitan_temp=norm(v1-vTitan);
                vinfEnceladus_temp=norm(v2-vEnceladus);

                % Standard Node Description
                vinfCAR_atTitan=v1-vTitan;
                SVTitan=[rTitan vTitan];
                [TitanNode] = car2vinfAlphaCrank(vinfCAR_atTitan, SVTitan);

                vinfCAR_atEnceladus=v2-vEnceladus;
                SVEnceladus=[rEnceladus vEnceladus];
                [EncelNode] = car2vinfAlphaCrank(vinfCAR_atEnceladus, SVEnceladus);
               
                % NodeInitialization_temporalStore = [MATitan trAEnceladus ToF_days LambertTypes VinfTitan vInfEnceladus]
                % Departure point Calculations
                MAEnceladus=mod(TrA_Enceladus-nEnceladus*TOF_seconds,2*pi); % unsure if this does what I want it to do, TBR.
                ToF_days=TOF_seconds/3600/24;
                % Arrival point Calculations
                MATitan=mod(0+nTitan*TOF_seconds,2*pi);
                % NodeInitialization_temporalStore= [NDeparture(1) MTitan(1) MEnceladus(1) ToF(1) LambertType(2) VinfDeparture(3) nArrival(1) Mtitan(1) MEnceladus(1) VinfArrival(3)]
                BLOCK1nodes{nodenumber}.ResAnatomyTransfer=[0 0 FourOptions(iL,1:2) norm(v1-vTitan)];
                BLOCK1nodes{nodenumber}.ResAnatomyTargeted=[Resonances_Enceladus(jR,:) p p vinfEnceladus];
                BLOCK1nodes{nodenumber}.NodeInitialization=[5 0 ToF_days vinfCAR_atTitan 1 TrA_Enceladus v2-vEnceladus];
                % [DepNode_Moon_Index NODECODE MeanTitan MeanEnceladus ToF_transfer ArrNode_Moon_Index NODECODE MeanTitan MeanEnceladus]
                Menceladus_atstart=mod(TrA_Enceladus-nEnceladus*ToF_days*Days2Sec,2*pi);
                MTitan_atend=mod(TrA_Titan+nTitan*ToF_days*Days2Sec,2*pi);
                BLOCK1nodes{nodenumber}.NodeInit_standard=[5 TitanNode TrA_Titan Menceladus_atstart ToF_days 1 EncelNode MTitan_atend TrA_Enceladus];
                % BLOCK1nodes{nodenumber}.NodeInit_standard=[5 0 ToF_days TitanNode 1 TrA_Enceladus EncelNode];
                BLOCK1nodes{nodenumber}.State_vector_departure_analitic=FourOptions(iL,5:10);
                nodenumber=nodenumber+1;
                % Anatomy of the resonance [M N p q]
                % in our cases p and q are suppose to be the same
                % [6 1 -1 -1 Vinf   

            end
        end
    end
end
%--------------------------------------------------------------------------
% clean zeros
% indexDEL=find(NodeInitialization_temporalStore(:,3)==0);
% NodeInitialization_temporalStore(indexDEL,:)=[];