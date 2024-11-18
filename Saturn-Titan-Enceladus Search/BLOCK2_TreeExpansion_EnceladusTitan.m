%% TreeExpansion_EnceladusTitan_GIGANTES.m
%--------------------------------------------------------------------------
% expansion of the tree
NewLayerCounter=1;
NumofCompletePaths=0;
clear treeExploration_temp
for iN=1:length(treeExploration)
    if sum(treeExploration{iN}.route(:,8)==1)<TargetEnceladusFB % stop condition if maximum number of Enceladus fly-bys is reached
        %--------------------------------------------------------------------------
        % Expand Tree Level by branching each current branch
        if treeExploration{iN}.route(end,8)==1 % End Node in Enceladus

            Mtitan_atstart=treeExploration{iN}.route(end,5); % mean anomaly of Titan at the start of the leg
            Menceladus_atstart=treeExploration{iN}.route(end,6);% mean anomaly of Enceladus at the start of the leg
            FullNode=treeExploration{iN}.route(end,:); % provide only the latest leg in the route

            %--------------------------------------------------------------------------
            % Distance Check with Titan
            %--------------------------------------------------------------------------
            % Test(FullNode,Mtitan_atstart,Menceladus_atstart,ToF2LastNode, pars)
            % [MissDistance, ToF_res]=Intersection_Titan_check(FullNode,Mtitan_atstart, Menceladus_atstart, ToF2LastNode, pars,0);
            % Change the distance threshold function for a numerical one at the
            % moment, the analytical one is rather tricky to contamplate all
            % potential cases.
            NextResonance=treeExploration{iN}.routeAnatomy(end,1:2);
            MissDistance=Test_distance_v27(NextResonance,FullNode, pars);

            % MissDistance
            if MissDistance>ThresholdDistance2routingTitan_Obligation
                %--------------------------------------------------------------------------
                % Here an Enceladus Enceladus Route is inserted
                NewLeg(1)= 1; % Departure Enceladus node
                NewLeg(2:4)= FullNode(9:11); % Node Departure at Enceladus
                NewLeg(5)= FullNode(12); % mean anomaly Titan.
                NewLeg(6)= FullNode(13); % mean anomaly Enceladus.
                NewLeg(7)= NextResonance(1)*PeriodEnceladus/Days2Sec; % ToF of next leg, which corresponds to full period of one spacecraft revolution [days]
                NewLeg(8)= 1; % Arrival Enceladus node
                NewLeg(9:11)= FullNode(9:11); % Node Arrival Enceladus after one resonance
                NewLeg(12)= mod(FullNode(12)+nTitan*NewLeg(7)*Days2Sec,2*pi); % mean anomaly Titan.
                NewLeg(13)=FullNode(13); % mean anomaly Enceladus.

                Expanded_route_anatomy=[treeExploration{iN}.routeAnatomy; treeExploration{iN}.routeAnatomy(end,:)];
                Expanded_route=[treeExploration{iN}.route; NewLeg];
                Expanded_routeDVCosts=[treeExploration{iN}.routeDVCosts; 0];
                treeExploration_temp{NewLayerCounter}.routeAnatomy=Expanded_route_anatomy;
                treeExploration_temp{NewLayerCounter}.route=Expanded_route;
                treeExploration_temp{NewLayerCounter}.routeDVCosts=Expanded_routeDVCosts;

                % Overall DV Cost of the Route
                % routeDVCosts(find(isnan(Expanded_routeDVCosts)))=[];
                treeExploration_temp{NewLayerCounter}.TotalDVCost=sum(Expanded_routeDVCosts);

                NewLayerCounter=NewLayerCounter+1;
            end

            if (MissDistance<=ThresholdDistance2routingTitan_Optional)
                %--------------------------------------------------------------
                % Here an Enceladus Titan Lambert arc is inserted
                % disp('first time here')
                % compute and store the Lambert arc that Minimizes DV to fly-by
                run BLOCK3_Enceladus_Titan_PorkChop.m
                if ~isempty(indexFeasible)
                    for iLA=1:length(indexFeasible)
                        % here we now know there is a feasible leg that reaches
                        % Titan.
                        run BLOCK3_Enceladus_Titan_SingleArc.m
                        DV_at_Enceladus_previous=DV_at_Enceladus;
                        NewLeg_previous=NewLeg;

                        Anatomy=[1 1];
                        if (mod(NewLeg_previous(4),2*pi)>pi/2)&(mod(NewLeg_previous(4),2*pi)<3*pi/2)
                            Anatomy(1)=-1;
                        end

                        if (mod(NewLeg_previous(11),2*pi)>pi/2)&(mod(NewLeg_previous(11),2*pi)<3*pi/2)
                            Anatomy(2)=-1;
                        end

                        Expanded_route_anatomy=[treeExploration{iN}.routeAnatomy; 0 0 Anatomy NewLeg_previous(2)];
                        Expanded_route=[treeExploration{iN}.route; NewLeg_previous];
                        Expanded_routeDVCosts=[treeExploration{iN}.routeDVCosts; DV_at_Enceladus_previous];
                        treeExploration_temp{NewLayerCounter}.routeAnatomy=Expanded_route_anatomy;
                        treeExploration_temp{NewLayerCounter}.route=Expanded_route;
                        treeExploration_temp{NewLayerCounter}.routeDVCosts=Expanded_routeDVCosts;
                        % Overall DV Cost of the Route
                        % routeDVCosts(find(isnan(Expanded_routeDVCosts)))=[];
                        treeExploration_temp{NewLayerCounter}.TotalDVCost=sum(Expanded_routeDVCosts);

                        NewLayerCounter=NewLayerCounter+1;

                    end
                    %----------------------------------------------------------
                end

            end

        elseif treeExploration{iN}.route(end,8)==5 % Last node in Titan


            Mtitan_atstart=treeExploration{iN}.route(end,5); % mean anomaly of Titan at the start of the leg
            Menceladus_atstart=treeExploration{iN}.route(end,6);% mean anomaly of Enceladus at the start of the leg
            FullNode=treeExploration{iN}.route(end,:); % provide only the latest leg in the route
            %------------------------------------------------------
            % Expand the following leg from Titan
            % Lambert arc between Node Titan and Node Enceladus
            % allowing two defects at enceladus and Titan.
            % run CheckAllResonances.m
            PreviousCounter=NewLayerCounter;
            run BLOCK4_Titan_Enceladus_Linking.m

            % if PreviousCounter<NewLayerCounter
            %     disp('here')
            % end
            % Here now check pseudoresonances with Titan
            % Pseu
            if (sum(treeExploration{iN}.route(:,8)==5)<nConsecutiveTitanFBs&&PreviousCounter==NewLayerCounter)
                Mtitan_atstart=treeExploration{iN}.route(end,5); % mean anomaly of Titan at the start of the leg
                Menceladus_atstart=treeExploration{iN}.route(end,6);% mean anomaly of Enceladus at the start of the leg
                FullNode=treeExploration{iN}.route(end,:); % provide only the latest leg in the route
                run BLOCK4_Titan_Titan_PseudoRes.m
            end
        end
    else
        treeExploration_temp{NewLayerCounter}=treeExploration{iN};
        NewLayerCounter=NewLayerCounter+1;
        NumofCompletePaths=NumofCompletePaths+1;
    end

end

% treeExploration = treeExploration(~cellfun(@isempty, treeExploration));
treeExploration=treeExploration_temp;