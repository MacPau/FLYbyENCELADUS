clear all; close all; clc; format short g; warning off;
%--------------------------------------------------------------------------
%% Initialization of Algorithm (Constants, Constraints and other variables)
%--------------------------------------------------------------------------
% Constants & Parameters
RAD=pi/180;
Days2Sec=3600*24;
% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [1, 5]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
% (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
% (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)
% Retrieve Central Body (Planet) Parameters
[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral); %[km3/s2],[km],[km] & [km]
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
PeriodEnceladus=2*pi*sqrt(pars.Moon.OrbRad(1)^3/pars.Planet.mu);
nEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);
nTitan=sqrt(pars.Planet.mu/pars.Moon.OrbRad(2)^3);
for i = 1:size(pars.INPUTS.idMoon,2)
    pars.Moon.Vel(i)    = sqrt(pars.Planet.mu/pars.Moon.OrbRad(i));           %[km/s] Moon Orbital velocity
    pars.Moon.Period(i) = 2*pi*sqrt(pars.Moon.OrbRad(i)^3/pars.Planet.mu);    %[s] Moon orbital period
    pars.Moon.HillSph(i) = pars.Moon.OrbRad(i)*( pars.Moon.mu(i)/(3*(pars.Moon.mu(i) + pars.Planet.mu)))^(1/3);     %[km] Moon Hill's Sphere
end
%--------------------------------------------------------------------------
% Epremeris Moons Set up (Circular Co-planar)
KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
SVTitan = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
SVEnceladus = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
%--------------------------------------------------------------------------
% v infinity bound for enceladus 
vInfEnceladus_bounds=[3 6];
% v infinity bound for titan 
vInfTitan_bounds=[2 7];
%--------------------------------------------------------------------------
ThresholdDistance2routingTitan_Obligation=1.5; % Titan close approach at which a Lambert arc to Titan is inserted
ThresholdDistance2routingTitan_Optional=3; % Titan close approach at which a Lambert arc to Titan is inserted
AcceptableDV_threshold=0.05; MaxDV_Allowed=0.05; % km/s
TargetEnceladusFB=8;
nConsecutiveTitanFBs=5;
numMaxOfLoops=TargetEnceladusFB+nConsecutiveTitanFBs;
FileNameStore='Saturn_Enceladus_Titan_exploration_results_4GIGANTES.mat';
% Pseudo Resonances Considered at Titan
PseudoRes_Titan_List=[1 0; 0 1; 1 1; 2 1; 1 2; 2 2; 3 2; 2 3];
tic
%% EXPAND Step 1 - Initialization of Phase
vInfStep=0.05;
vInfEnceladus_list=[vInfEnceladus_bounds(1):vInfStep:vInfEnceladus_bounds(2)];
% Resonances_Enceladus = [11 2; 6 1; 13 2; 7 1; 15 2; 8 1; 17 2];
Resonances_Enceladus = [6 1; 7 1; 8 1];
% Generate all possible Nodes to chose from the vector Resonances_Enceladus
run BLOCK1_Generation_NODE1.m
%--------------------------------------------------------------------------
% Initial Tree expansion
% Expand tree to all possible nodes
% treeExploration.route=[];
numberofNodes=length(BLOCK1nodes);
for iN=1:numberofNodes
    % route=[Node1 Node2]
    treeExploration{iN}.routeAnatomy=[BLOCK1nodes{iN}.ResAnatomyTransfer; BLOCK1nodes{iN}.ResAnatomyTargeted];
    treeExploration{iN}.route=BLOCK1nodes{iN}.NodeInit_standard;
    treeExploration{iN}.routeDVCosts=0;
end
disp(['Loop ', num2str(1),' out of ', num2str(numMaxOfLoops),': ', num2str(length(treeExploration)),' paths being tracked'])
%--------------------------------------------------------------------------
% OPTIONAL
% By recalling the generation of nodes, we construct a tree search using a
% lower number of starting nodes, but allow possible changes to larger
% period resonances, if these linkages exist. 
% Resonances_Enceladus = [6 1; 7 1; 8 1];
% Resonances_Enceladus = [11 2; 6 1; 13 2; 7 1; 15 2; 8 1; 17 2];
% Resonances_Enceladus = [7 1];
% run BLOCK1_Generation_NODE1.m
%% EXPAND Level 2
% Step 2.1 - Define Current Node
% expand for numMaxOfLoops loops
treeExploration_DeadBranches=[];
DeadEndCount=1;
TreeGrowth=[1 length(treeExploration)];
for iE=2:numMaxOfLoops
    run BLOCK2_TreeExpansion_EnceladusTitan.m
    %--------------------------------------------------------------------------
    % Filter out branches not satisfying DV constraint
    Delta_V=zeros(1,length(treeExploration));
    for iN=1:length(treeExploration)
        Delta_V(iN)=treeExploration{iN}.TotalDVCost;
    end
    % Remove Cases with high DB
    undexMaxDV=find(Delta_V<MaxDV_Allowed);
    treeExploration = treeExploration(undexMaxDV);
%--------------------------------------------------------------------------
   disp(['Loop ', num2str(iE),' out of ', num2str(numMaxOfLoops),': ', num2str(length(treeExploration)),' paths being tracked'])
   disp(['Completed Paths ', num2str(NumofCompletePaths)])
   clear treeExploration_temp
   TreeGrowth=[TreeGrowth; iE length(treeExploration)];
   % Saving at the end of the loop
   save(FileNameStore);
end
%% Postprocessing and Pruning
Delta_V=zeros(1,length(treeExploration));
Score_TidalSpread=zeros(1,length(treeExploration));
Score_TidalSpread2=zeros(1,length(treeExploration));
Score_ScientificInteres=zeros(1,length(treeExploration));
Delta_longitude=zeros(1,length(treeExploration));
Total_time_of_flight=zeros(1,length(treeExploration));
for iN=1:length(treeExploration)
    % Analyse each individual route
    Entire_route=treeExploration{iN}.route;

    % Enceladus FlyBys
    indexRelevantFlyBys=find(Entire_route(:,8)==1);
    numEnceladusFB(iN)=length(indexRelevantFlyBys);
    % Enhance TreeExploration Structure
    treeExploration{iN}.numEnceladusFB=numEnceladusFB(iN);

    Delta_V(iN)=treeExploration{iN}.TotalDVCost;

    % Diversity on Flybys
    RelevantNodes=[Entire_route(indexRelevantFlyBys,9) Entire_route(indexRelevantFlyBys,10) Entire_route(indexRelevantFlyBys,11)/RAD];
    Unique_nodes = table(round(RelevantNodes(:,1)/0.500),RelevantNodes(:,2),round(RelevantNodes(:,3)/90));
    Score_ScientificInteres(iN)=height(unique(Unique_nodes));


    % True Longitude of FlyBys
    LongitudeofallEnceladusNodes=Entire_route(indexRelevantFlyBys,13);
    Unique_nodes = table(round(LongitudeofallEnceladusNodes/RAD/15));
    Score_TidalSpread(iN)=height(unique(Unique_nodes));

    % Longitude of furthest Enceladus FBs - relevant to illumination
    % Conditions
    Initial_Longitude=Entire_route(indexRelevantFlyBys(1),13);
    LongitudeofallEnceladusNodes=Entire_route(indexRelevantFlyBys,13);
    Delta_longitude(iN)=max(abs((Initial_Longitude-LongitudeofallEnceladusNodes)))/RAD;

    Travel_matrix=zeros(length(LongitudeofallEnceladusNodes));
    for itm=1:length(LongitudeofallEnceladusNodes)
        for jtm=itm:length(LongitudeofallEnceladusNodes)
            Travel_matrix(itm,jtm)=abs(LongitudeofallEnceladusNodes(itm)-LongitudeofallEnceladusNodes(jtm));
        end
    end
    Score_TidalSpread2(iN)=height(unique(Unique_nodes))*sum(Travel_matrix(:));

    %Overall Time of Flight
    Total_time_of_flight(iN)=sum(Entire_route(:,7)); % days
    %--------------------------------------------------------------------------
    % Enhance TreeExploration Structure
    treeExploration{iN}.Score_ScientificInteres= Score_ScientificInteres(iN);
    treeExploration{iN}.Delta_longitude_excursion=Delta_longitude(iN);
    treeExploration{iN}.TotalTimeofFlight=Total_time_of_flight;
end
  save(FileNameStore);

  TT=toc/3600;
  disp(['Simulation Completed in ',num2str(TT),' h'])

%% Final Plotting and Analysing

% Create figure 
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create plot
plot(Delta_longitude,Delta_V,'MarkerFaceColor',[0 0 0],'Marker','o','LineStyle','none');
% Create ylabel
ylabel('\Deltav [km/s]');
% Create xlabel
xlabel('\DeltaL [degrees]');
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',12,'FontWeight','bold',...
    'XGrid','on','YGrid','on');

% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create plot
plot3(Delta_longitude,Total_time_of_flight,Delta_V ,'MarkerFaceColor',[0 0 0],'Marker','o','LineStyle','none');
zlabel('\Deltav [km/s]');
ylabel('Total Time of Flight [days]');
xlabel('\DeltaL [degrees]');
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',12,'FontWeight','bold',...
    'XGrid','on','YGrid','on');
view(60, 45)

%% bin 
% % ploting the maximum longitude tour
% [~,indexMax]=max(Delta_longitude);
% Full_route=treeExploration{indexMax}.route;
% totalDV=sum(treeExploration{indexMax}.routeDVCosts)
% routeAnatomy=treeExploration{indexMax}.routeAnatomy;
% plotfullSaturnTransfer(Full_route,pars)
% 
% % Create figure
% figure1 = figure;
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% % Create plot
% plot(Score_ScientificInteres,Delta_V,'MarkerFaceColor',[0 0 0],'Marker','o','LineStyle','none');
% % Create ylabel
% ylabel('\Deltav [km/s]');
% % Create xlabel
% xlabel('Diversity of Flybys');
% box(axes1,'on');
% hold(axes1,'off');
% % Set the remaining axes properties
% set(axes1,'FontName','Times New Roman','FontSize',12,'FontWeight','bold',...
%     'XGrid','on','YGrid','on');
% 
% 
% 
% % ploting the maximum Diversity  tour
% [~,indexMax_Score]=max(Score_ScientificInteres);
% Full_route=treeExploration{indexMax_Score}.route;
% totalDV=sum(treeExploration{indexMax_Score}.routeDVCosts)
% routeAnatomy=treeExploration{indexMax_Score}.routeAnatomy;
% plotfullSaturnTransfer(Full_route,pars)
% 
% 
% % Create figure
% figure1 = figure;
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% % Create plot
% plot(Score_ScientificInteres,Total_time_of_flight,'MarkerFaceColor',[0 0 0],'Marker','o','LineStyle','none');
% % Create ylabel
% ylabel('Total Time of Flight [days]');
% % Create xlabel
% xlabel('Diversity of Flybys');
% box(axes1,'on');
% hold(axes1,'off');
% % Set the remaining axes properties
% set(axes1,'FontName','Times New Roman','FontSize',12,'FontWeight','bold',...
%     'XGrid','on','YGrid','on');
% 
% indexMax_Score=find(Score_ScientificInteres==max(Score_ScientificInteres));
% 
% treeExploration_SubSet=treeExploration(indexMax_Score);
% indexValue=[1:1:length(indexMax_Score)];
% InfoSubSet=[indexValue' Delta_V(indexMax_Score)' Total_time_of_flight(indexMax_Score)' Delta_longitude(indexMax_Score)']
% 
% % ploting the maximum Diversity  tour
% IndexCheck=3010
% Full_route=treeExploration_SubSet{IndexCheck}.route;
% totalDV=sum(treeExploration_SubSet{IndexCheck}.routeDVCosts)
% routeAnatomy=treeExploration_SubSet{IndexCheck}.routeAnatomy;
% plotfullSaturnTransfer(Full_route,pars)
% 
% 
% 