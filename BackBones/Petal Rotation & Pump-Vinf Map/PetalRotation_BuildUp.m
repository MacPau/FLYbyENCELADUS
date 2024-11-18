function [PetalRots_Seqs, PetalRots_FlybysData] = PetalRotation_BuildUp(Reso_TransfersDatabase, PseudoReso_TransfersDatabase, pars)

% This function computes a petal rotation sequence given a set of flyby
% options and the desired change in argument of periapsis.

% Author: J.C. Garcia Mateas
% Last revision: 02/11/2024

%% INPUTS %%
% - Reso_TransfersDatabase: Matrix of size A x 10, contaning the data for
%            the resonant transfers being considered, where A is the n° of
%            options. Col1 = N (N) of moon revolutions; col2 = n° of
%            spacecraft revolutions; col3, col4 & col5; outgoing Vinf [km/s] 
%            pump angle [rad] & crank angle [rad]; col6: TOF [seconds]; 
%            col7,8 & 9: incoming Vinf, pump & crank after the transfer.
%            Col10 = type of transfer (88 = Outbound/Outbound / 11 =
%            Inbound-Inbound / 81 = Outbound-Inbound / 18 = Inbound - Outbound)
%
% - PseudoReso_TransferDatabase: Matrix of size B x 10 containing the data
%           for all the pseudo-resonant transfers being considered. The
%           columns contain the same information as the one described above
%           for Reso_TransfersDatabase.
%
% - pars: Structure containing problem constants and parameters

%% OUTPUTS %%
%
% - PetalRotations_Seqs: Cell of size 1 x N containing the matrices
%               associated to the different petal rotation options. Each
%               solution will be a matrix of size C x 10, where C is the
%               number of transfers in the sequence and the columns of the
%               matrix are the same as described in INPUTS Reso_TransfersDatabase
%
% - PetalRots_FlybysData: Cell of size 1 x N containing structures with N-1
%               flybys of the petal rotation sequence and the data
%               associated to each flyby defined in several fields.

%% CHANGELOG %%


%% BUILD THE PETAL ROTATION SEQUENCES %%
% Initialize Flag to indicate when desired petal rotation is achieved
Flag_PetalRot_BuildUP = 0;

% Initialize Flag regarding first segment of the Petal Rotation
Flag_Flyby_1 = 1;

% Initialize Flyby Counter
Flyby_Count = 0;


while Flag_PetalRot_BuildUP == 0

    if Flag_Flyby_1 == 1 % If the first flyby must be computed

        % Retrieve Resonant & Pseudo-Resonant Transfers to be considered
        [Reso_Options, PseudoReso_Opts] = Identify_Transfers(Reso_TransfersDatabase, PseudoReso_TransfersDatabase, pars);

        % Initialize variables
        Buffer_FlybyData = [];

        tic
        % Explore possible flybys & construct petal rotation sequences
        [Buffer_Transfers, Buffer_Transfers_IDX, Buffer_FlybyData] = Transfer_Combinations(Reso_Options, PseudoReso_Opts, Buffer_FlybyData, pars);

        % Rank solutions according to change in true longitude
        [Buffer_Transfers, Buffer_FlybyData] = Ranking_Solutions(Buffer_Transfers, Buffer_FlybyData, pars);

        toc

        % Change Flag of First join after accommplished
        Flag_Flyby_1 = 0;

        Flyby_Count = Flyby_Count + 1;

         % Check if maximum number of flybys has been reached
         if Flyby_Count == pars.PetalRot.Num_Flybys 
             Flag_PetalRot_BuildUP = 1;
         end

    else

        % Convert the Transfers Buffer into 3D Matrix Form 
        Buff_Transfers_Matrix = cat(3, Buffer_Transfers{:});

        tic

        % Explore possible flybys & construct petal rotation sequences
        [Buffer_Transfers, Buffer_Transfers_IDX, Buffer_FlybyData] = Transfer_Combinations(Buff_Transfers_Matrix, PseudoReso_Opts, Buffer_FlybyData, pars);

        % Rank solutions according to change in true longitude
        [Buffer_Transfers, Buffer_FlybyData] = Ranking_Solutions(Buffer_Transfers, Buffer_FlybyData, pars);

        toc

         Flyby_Count = Flyby_Count + 1;

         % Check if maximum number of flybys has been reached
         if Flyby_Count == pars.PetalRot.Num_Flybys 
             Flag_PetalRot_BuildUP = 1;
         end

    end

end


%% PREPARE THE FUNCTION OUTPUT %%
% Define function outputs
PetalRots_Seqs       = Buffer_Transfers;
PetalRots_FlybysData = Buffer_FlybyData;

end