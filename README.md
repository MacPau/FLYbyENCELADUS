# Flyby Missions to Enceladus using MATLAB
This folder contains FLYbyENCELADUS toolbox (100% MATLAB-based) for designing flyby missions to Enceladus. FLYbyENCELADUS is meant as a preliminary design tool for the pseudo-orbiter phase of ESA's next Large-class mission (L4), a mission to Saturn which was recommended by the Expert's Committee as part of ESA's Voyage 2050 plan. There are four main application in this tool:

1. GT Visualization - Ground Track Visualization containts demo scripts to plot the ground tracks of a two Node link at Enceladus. Please, refer to GIGANTES Flyby Mission Technical Report for the definition of a Moon encounter Node. A flyby is represented by two Enounters Node at a given Moon and the representation of the Ground Track is given as a mercator projection, as well as in a 3D moon centric with b-plane representation.
2. Backbones - This contains the necessary tools (Demo scripts and functions) to create the main types of short flyby sequences for science applications at Enceladus. These includes Crank-over-the-top Sequences (COT), Crossing Ground tracks, Petal rotation's and Leading and Trailing flybys.
3. Saturn-Titan-Enceladus Search - This contains the demoscripts and main algorithms to run a tree search exploration of Enceladus sequences that may be encountering Titan as a gravitational leverage to change encounter conditions at Enceladus.
4. Saturn-Titan-Enceladus TestCase BuildUp - This containts the final demoscripts that, given the sequence database provided by the Saturn-Titan-Enceladus Search, they build a series of specific ground track characteristics to explore Enceladus. 

These algorithms are the result of the work carried out at ISAE-SUPAERO under the ESA project GIGANTES in 2024. 

The corresponding Python implementation can be found in ESA's MIDAS tool at:

```bash
https://midas.io.esa.int/midas/api_reference/generated/????
```

This readme.md shows how to install and run a preliminary design of fly-by missions to Enceladus with test scripts provided.

## Installation & Requirements

To work with FLYbyENCELADUS, one can simply clone the repository in the local machine:

```bash
git clone "https://github.com/MacPau/FLYbyENCELADUS.git"
```

Only invited developers can contribute to the folder, and each should create a separate branch or fork, otherwise push requests will not be accepted on main branch modifications. This work is under [European Space Agency Public License (ESA-PL) Permissive (Type 3) - v2.4 licence](https://essr.esa.int/license/european-space-agency-public-license-v2-4-permissive-type-3). See [LICENSE](https://github.com/andreabellome/saturn_moon_tours/blob/main/LICENSE) file for details.

To cite this software, please cite Bellome [[4]](#4) and Bellome et al. [[5]](#5).

To run a full exploration of Saturn system, the following system requirements are recommended:
+ CPU six-core from 2.6 GHz to 3.6 GHz
+ RAM minimum 16 GB
+ Any version of [MATLAB](https://it.mathworks.com/products/matlab.html)>2021b 

## Usage and test cases

To use the repository, one finds four main folders, one for each of the functionalities described above (GT Visualization, Backbones, Saturn-Titan-Enceladus Search, Saturn-Titan-Enceladus TestCase BuildUp). To run the Demo_scripts within these folders, certain libraries need to be in the Matlab path. These libreries are included under the name Add2Path XXXX. Note Add2Path THIRDPARTY contains certain typical standard functionalities (conversions and Lambert Arcs) that can be easily substituted by others. These functions are not supported by the license provided. 

### 1. GT Visualization
This folder containts two main DEMO_scripts.m that shall be used as base to prepare specific visualizations of the ground tracks achieved with sequences of flybys. 

The reference script is [DEMO_Plot_GroundTracks.m](https://github.com/update_update_update.m). This provided the minimum base to plot a mercator projection of a flyby given the flyby entry Node description and flyby departure Node. 

The script [DEMO_Design_and_Visualize_GTs.m](https://github.com/update_update_update.m) adds the possibility to use the v_infinity sphere to understand the available flyby configurations and also plots the fly-by in a 3D representation with b-plane, to help understanding the physical behaviour observed in the mercator projection. 

### 2. Backbones
This folder contains the codes required to construct typical sequences of science driven flybys. In particular, one sub-folder contains the codes for achieving crank-over-the-top (COT) sequences and crossing ground-tracks, while another sub-folder contains the codes to compute petal rotations and the pump-Vinf map.

### 2.1 COT sequences
This script enables the user to compute and plot the ground-tracks and spacecraft orbit around the central body for COT sequences. The reference script to run the full construction is [Main_COT.m](https://github.com/MacPau/GIGANTES-Repo/tree/main/BACKBONES/COT%20and%20Crossing%20GT)

The script starts by defining several inputs, such as the central body and moon of interest at which the COT will be built, the hyperbolic excess velocity at the moon  (which can be introduced as a vector with several values), resonance ratios of interest to consider and the epoch at which to obtain the COT. 

```matlab
%% INITIALIZATION %%
clear all; close all; clc; format long g;

% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [1]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
                            % (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
                            % (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

% Define Flyby Hyperbolic Excess Velocities
pars.INPUTS.V_inf = 4;     %[km/s] 

% Define the resonances to consider
pars.INPUTS.Resonances = [7 1];

% Define starting epoch
pars.INPUTS.epoch0 = 0;    % Days passed since MJD2000
```
The initialization part also includes a section to define the parameters relative to the ground-track computation. In particular, the user can define the minimum flyby altitude, the maximum altitude above which the ground-track will not be plotted, as well as the number of time-steps and time since periapsis to be used for the ground-track propagation. It must be noted that all these parameters should be adapted depending on the moon of interest.

```matlab
% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 5;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 
```
The algorithm then computes the parameters associated to the resonances being considered for the COT, and determines the maximum change in crank angle that can be achieved with one flyby at the moon of interest. Note that the current codes only support pure-cranking COT's, and not the build-up of sequences also achieving pump-angle changes. 

Before building up the COT, the user can define further characteristics, including the starting crank-angle, the cranking direction (which leads to North to South or South to North coverage) and, if desired, a bound in the latitude of the periapsis, which can be used to construct partial-COTs that do not cover the full hemisphere.

```matlab
%% BUILD THE CRANK OVER THE TOP SEQUENCE %%
% Define COT strategy
pars.INPUTS.COT.Start_Crank         = 180; % Starting Crank (INBOUND: 180°, OUTBOUND: 0°)
pars.INPUTS.COT.Latitude_Peri_Limit = -70; % Latitude limit to consider South Pole flybys [°]
pars.INPUTS.COT.Crank_Direction     = +1;  % Cranking Direction (-1 for negative, +1 for positive)
pars.INPUTS.COT.Pump_Angle          = alfa_res; % Pump angle of the selected resonance for the COT [rad]
pars.INPUTS.COT.Max_Crank_Change    = Delta_Crank_max; %[rad] Max. crank change

% Build-Up the COT
tic
[COT_Data] = COT_BuildUp(pars);
toc
```
The resulting COT sequence data can then be post-processed to retain only the flybys respecting the pars.INPUTS.COT.Latitude_Peri_Limit constraint. The resulting ground-tracks, evolution of the spacecraft true anomaly (to visualize the switch in the type of flyby Inbound to Outbound or viceversa when a crank of 90° is reached) and the spacecraft orbit around the Moon and the central body can be obtained.

```matlab
% Plot the flyby groundtracks
fig1 = figure( 'Color', [1 1 1] );
hold on;
plotTextureLatLong(pars.INPUTS.idMoon , pars.INPUTS.idCentral , 1);
plotSquares(pars, 1);  
axis normal;
for i = 1:size(Extract_Flybys, 2)
    Plot_Flyby_GT(Extract_Flybys(i), colors(i,:));
end

% Plot the flybys around the Moon sphere
fig3 = Plot_Flybys_MoonOrb(Extract_Flybys, pars);

% Plot the planeto-centric orbit of the SC
fig4 = Plot_Flybys_CentralBOrb(Extract_Flybys, pars);
```
### 2.2 Crossing Ground-tracks
This script enables the user to compute and plot the ground-tracks and spacecraft orbit around the central body for sequences leading to ground-tracks crossing each other. The reference script to run the full construction is [Main_CrossingGT.m](https://github.com/MacPau/GIGANTES-Repo/tree/main/BACKBONES/COT%20and%20Crossing%20GT). 

To obtain these ground-tracks at Enceladus, where the bending achievable at high relative velocities is very small and Titan encounters occur regularly, the approach followed has been to envision two partial-COTs, that is, to construct two sequences of a few flybys following a pure-cranking strategy but in which the cranking direction is switched between each other, thus using the principles of the COT codes described in the previous section.  

```matlab
% Define COT strategy N°1
pars.INPUTS.COT.Crank_Direction     = -1;  % Cranking Direction (-1 for negative, +1 for positive)

% Build-Up the COT N°1
tic
[COT_Data1] = COT_BuildUp(pars);
toc

% Define COT strategy N°2
pars.INPUTS.COT.Crank_Direction     = +1;  % Cranking Direction (-1 for negative, +1 for positive)

% Build-Up the COT N°2
[COT_Data2] = COT_BuildUp(pars);
```
In the post-process section, the user can then defined the regions of latitude at which the periapsis of the crossing-ground tracks are desired, which will lead to extracting from the COT sequences the flyby nodes and corresponding information leading to achieving these sequences. A plot of one example at Enceladus is provided below.
```matlab
%% POST-PROCESS %%
% Define bounds for Latitudes of Interest
pars.INPUTS.COT.Min_Lat = -45; %[degrees]
pars.INPUTS.COT.Max_Lat = -35; %[degrees]
```
![image](https://github.com/user-attachments/assets/38e4b3c6-e7fe-4420-a8f2-6e48b4139595)

### 2.3 Petal Rotations 
To construct a petal rotation, that is, a sequence of flybys and transfers which is used to rotate the line of apsides of the spacecraft orbit and thus change the location of the flybys in the moon's orbit, the script to run is [Main_Petal_Rotation.m](https://github.com/MacPau/GIGANTES-Repo/blob/main/BACKBONES/Petal%20Rotation%20%26%20Pump-Vinf%20Map/Main_Petal_Rotation.m)

The main script starts by allowing the user to define the input parameters in the Initialization section as was done previously for the COT sequence. Here, it is worth noting a couple of extra inputs. Firstly, the change in true longitude of the encounter (with respect to the initial flyby location) that is desired, which is defined as an objective such that pairs of pseudo-resonant transfers are concatenated until the rotation objective defined by this parameter is achieved. Secondly, the rotation direction, which is coupled with the starting crank and enables to specify wether the line of apsides will be rotated clockwise or counterclockwise. Finally, the user can define a maximum number of flybys to consider, after which the petal rotation is terminated even if the change in true longitude objective has not been achieved.

```matlab
% Define the Petal Rotation Inputs 
pars.PetalRot.Starting_Reso  = [7 1];    % Resonance Orbit at Start of Petal Rotation
pars.PetalRot.Starting_Crank = 0;        % [rad] Crank angle at start of Petal Rotation (0 = Outbound / 180° = Inbound)
pars.PetalRot.Starting_Vinf  = 4.0;      % [km/s] Vinf at Start of Petal Rotation
pars.PetalRot.TrueLong_Obj   = deg2rad(10); % [rad] Objective of change in true longitude
pars.PetalRot.Direction      = +1;      % Direction in which to rotate (+1 = clockwise / -1 = anticlockwise)
pars.PetalRot.Num_Flybys     = 6;       % Maximum n° of flybys that the Petal Rotation can use
pars.PetalRot.Beam_Width     = 1000;    % N° of solutions to store at each Petal Rotation expansion step
```
the script computes the database of resonant and pseudo-resonant transfers for the resonance ratios and hyperbolic excess velocities specified using the following function calls:

```matlab
%% COMPUTE DATABASE OF RESONANT TRANSFERS %%
Reso_TransfersDatabase = [];
for i = 1:size(pars.INPUTS.V_inf, 2)
    [Reso_Transfers, ResoStruc] = ResoTransfers_Generation(pars.INPUTS.V_inf(i), pars);
    Reso_TransfersDatabase = [Reso_TransfersDatabase; Reso_Transfers];
end

%% COMPUTE DATABASE OF PSEUDO-RESONANT TRANSFERS %%
% Define which 81 pseudo-resonance to use
pars.INPUTS.remove81 = 1;  %[0] = add 1 rev to the outbound-inbound short case
                           %[1] = do not add one rev to the outbound-inbound short case

PseudoReso_TransfersDatabase = [];
for i = 1:size(pars.INPUTS.V_inf, 2)
    [PseudoReso_Transfer, PseudoResoStruct] = PseudoResoTransfers_Generation(pars.INPUTS.V_inf(i), pars);
    PseudoReso_TransfersDatabase = [PseudoReso_TransfersDatabase; PseudoReso_Transfer];
end
```
Once the database of transfer options has been defined, the function to construct the petal rotation sequences is called.
```matlab
[PetalRots_Seqs, PetalRots_FlybysData] = PetalRotation_BuildUp(Reso_TransfersDatabase, PseudoReso_TransfersDatabase, pars);
```
The user can then select which of the solutions obtained he wants to post-process by defining the parameter "idx" and then running the rest of the script, which yields values on the rate of change in true longitude, the cost in terms of $\Delta$ V (measured as a $\Delta$-V defect between the arriving and outgoing hyperbolic excess velocities at the flyby), and plots of the orbit around the central body and the associated ground-tracks. An example of the Saturn-centric orbit view for a petal rotation sequence at 4 km/s at Enceladus using the 7:1 pseudo-resonances is presented as an example below.

![image](https://github.com/user-attachments/assets/e9594307-4eae-4646-8c15-a8ba92fdfdf3)


#### 2.3.1 Pump - $v_\infty$ Map
A particularly useful tool when designing petal rotation sequences and, more generally, pseudo-orbiter phases, is the Pump - $v_\infty$ map, a graphical tool where contours for different resonant and pseudo-resonant transfers are plotted as function of pump angle ($\alpha$) and hyperbolic excess velocity $v_\infty$. This map can be plotted for any moon using the script [Main_PumpVinf_Map.m](https://github.com/MacPau/GIGANTES-Repo/blob/main/BACKBONES/Petal%20Rotation%20%26%20Pump-Vinf%20Map/Main_PumpVinf_Map.m). An example of this graph for Enceladus at the $v_\infty$ of interest for the L4 mission is shown below.

![Pump_VinfMap_Large_withZoom](https://github.com/user-attachments/assets/8b6bd42c-7cec-421d-bb6d-a5f5d41752b1)

### 2.4 Leading and Trailing Edge Flybys
Finally, achieving flybys of the leading and trailing edge's of the mooon's can be done by employing pump-changing flybys, since a pump-up flyby will increase the orbit period by reducing the pump angle and thus fly over the trailing edge (longitude of 90° East) while a pump-down flyby will decrease orbit energy and thus fly over the leading edge (located at a longitude of 270° East). Thus, to see if ground-tracks can be achieved covering this longitudes the same codes as for the petal rotation can be used, and the user only has to define the number of flybys input as equal to 1 as follows:
```matlab
pars.PetalRot.Num_Flybys     = 1;       % Maximum n° of flybys that the Petal Rotation can use
```
Doing so enables to obtain the parameters of the flyby linking the starting resonant orbit and the post-flyby pseudo-resonance, with the ground-track being plotted using the [Plot_Flyby_GT.m] function. An example of a 4 km/s leading-edge flyby at Enceladus is shown below.

![image](https://github.com/user-attachments/assets/6499df45-0761-415b-8709-da51a0f09578)

### 3. Saturn-Titan-Enceladus Search

As described in GIGANTES Flyby Mission Technical Report, Titan is the main leverage to achieve large changes of fly-by conditions at Enceladus. This folder then presents two main Demo scripts to explore the search space of Enceladus Fly-by missions in orbits that intersect Titan’s orbit. 

[DEMO_Saturn_Enceladus_Titan_PseudoOrb_SearchEngine.m]() provides a simple set up that should run within 1 to 2h in a system as that recommended in section Installation & Requirements. 

Instead, [DEMO_Saturn_Enceladus_Titan_PseudoOrb_BeamSearch.m]() triggers an exploration that may take on the order of 24h to be completed. The search space expands very rapidly, so a basic beam search algorithm is inserted such that the exploration can be carried out in a reasonable time. 

Both Demo scripts require the definition of a set of problem parameters. These include: The relative velocity bounds for Enceladus and Titan, the distance threshold to force a fly-by to Titan and the distance threshold to attempt an encounter to Titan, the DV threshold for every single event (fly-by) and the sum, the maximum number of Enceladus flybys sought, the maximum number of Titan flybys allowed in the sequence, the name of the file where the database is stored and the pseudo-resonances allowed at Titan (if feasible from the DV point of view).
```matlab
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
```
Step 1 of the expansion generates all the Titan-Enceladus Nodes that respect the relative velocities constraints with both Enceladus and Titan and belonging to the Enceladus resonances defined in variable Resonances_Enceladus. The generation of nodes is discretized not only by the variable Resonances_Enceladus, but also by the relative velocity step applied at Enceladus relative velocity (as by variable vInfStep)
```matlab
vInfStep=0.05;
vInfEnceladus_list=[vInfEnceladus_bounds(1):vInfStep:vInfEnceladus_bounds(2)];
% Resonances_Enceladus = [11 2; 6 1; 13 2; 7 1; 15 2; 8 1; 17 2];
Resonances_Enceladus = [6 1; 7 1; 8 1];
% Generate all possible Nodes to chose from the vector Resonances_Enceladus
run BLOCK1_Generation_NODE1.m
```
Finally, the Expand Level 2, expands the tree in a breath first process as described by the algorithm presented in GIGANTES Flyby Mission Technical Report. 
```matlab
%% EXPAND Level 2
% Step 2.1 - Define Current Node
% expand for numMaxOfLoops loops
for iE=2:numMaxOfLoops
    run BLOCK2_TreeExpansion_EnceladusTitan.m
    save(FileNameStore);
end
```

### 4. Saturn-Titan-Enceladus TestCase BuildUp
Folder Saturn-Titan-Enceladus TestCase Buildup proposes three examples of construction of Enceladus Fly-by Missions. The basic process is to search for the desired type of sequence of fly-bys within the database computed with the DEMO_Saturn_Enceladus_Titan_PseudoOrb_SearchEngine.m algorithms. 
It should be noted that the output of the search carried out in DEMO_Saturn_Enceladus_Titan_PseudoOrb_SearchEngine.m contains only “void” fly-bys at Enceladus. This is, fly-bys that have no deflection or change between the arrival and final conditions. Hence, the purpose of the Demo_scripts in this section is to reconstruct the sequences by adding specific backbones or types of flybys at Enceladus. 


