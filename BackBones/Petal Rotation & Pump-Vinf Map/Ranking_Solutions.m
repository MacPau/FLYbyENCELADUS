function [Buffer_PetalRots_Ranked, Buffer_FlybyData_Ranked] = Ranking_Solutions(Buffer_PetalRotSeqs, Buffer_FlybyData, pars)

% This function ranks a set of transfer sequences according to a
% user-defined criteria. Currently, it ranks according to change in true
% longitude of the encounter (to be changed & improved in the future).

% Author: J.C. Garcia Mateas
% Last revision: 11/06/2024

%% INPUTS %%
% Buffer_PetalRotSeqs: Cell of size 1 x N, where N = n° of petal rotation 
%                      sequences stored up to now. Each cell is a matrix of 
%                      size A x 10, where A = n° of transfers involved in
%                      the petal rotation and the 10 columns contain the
%                      information characterizing each transfer, such that:
%                      col1 = N° moon revolutions; col2 = n° spacecraft revs; 
%                      col3 = Vinf [km/s], col4 = pump angle [rad], col5 = crank angle [rad], 
%                      col6 = TOF [seconds], col 7 = col3; col8 = col4;
%                      col9 = col6; col10 = transfer type.
%
% - Buffer_FlybyData: Cell of size [1xC], where C is the n° of petal rotation 
%                        sequences stored. Each cell is a structure with
%                        several fields containing the information
%                        regarding the flybys involved in the petal
%                        rotation.
%
% - pars: Structure containing problem constants and parameters

%% OUTPUTS %%
% - Buffer_PetalRots_Ranked: Cell of size 1 x M, containing M petal rotations
%                      ranked according to user-defined criteria. Each cell 
%                      is a matrix with the same characteristics as described 
%                      for the input Buffer_PetalRotSeqs. A beam width is
%                      applied to retain only a maximun n° of petal
%                      rotation sequences.
%
% - Buffer_FlybyData_Ranked: Cell of size 1 x M containing the flyby data
%                           of the petal rotations sequences stored in the
%                           output Buffer_PetalRots_Ranked.

%% CHANGELOG %%
% 11/06/2024, J.C. Garcia Mateas: added comments to the function.
%
% 01/07/2024, J.C. Garcia Mateas: eliminated propagation to compute state
%             vector at the end of the transfer; added DV as an objective
%             to the ranking of the solutions; updated the computation of 
%             the change in true longitude so that now is every 2 flybys.
%
% 08/10/2024, J.C Garcia Mateas: corrected bug in the computation of the
%             change in true longitude because it must be computed every
%             two flybys but this corresponds to the "odd" rows in
%             PetalRotation_Sequence because row1 is just the initial
%             condition and not a flyby. Thus, mod(Num_Flys, 2) ~= 0 in
%             line 101 to compute true longitude change.

%% FUNCTION %%

for k = 1:size(Buffer_PetalRotSeqs, 2)

    % Extract Petal Rotation Sequence & Data
    PetalRotation_Sequence = cell2mat(Buffer_PetalRotSeqs(k));
    PetalRotation_Data     = cell2mat(Buffer_FlybyData(k));

    for i = 1:size(PetalRotation_Sequence, 1)
    
        PseudoRes = PetalRotation_Sequence(i,:);
        TOF(i)    = PseudoRes(6);
    
        [~, rr_end, vv_end, ~] = vinfAlphaCrank_to_VinfCART(PseudoRes(7), PseudoRes(8), PseudoRes(9) , sum(TOF)/86400, pars.INPUTS.idMoon, pars.INPUTS.idCentral);

        OE_SC_end(i,:)  = car2kep([rr_end, vv_end], pars.Planet.mu);

        % Retrieve SC True Longitudes
        TrueLong_0  = mod(sum(OE_SC_end(1,4:6)), 2*pi);
        TrueLong(i) = mod(sum(OE_SC_end(i,4:6)), 2*pi);
        TrueLong_Diff(i, :) = TrueLong(i) - TrueLong_0;         % [rad]
    end

    % Store True Longitude in the Flyby Data structure
    if isempty(PetalRotation_Data(end).TrueLongitudes)
        PetalRotation_Data(end).TrueLongitudes  = TrueLong;
        PetalRotation_Data(end).Change_TrueLong = TrueLong_Diff;
    else
        PetalRotation_Data.TrueLongitudes(end+1)  = TrueLong;
        PetalRotation_Data.Change_TrueLong(end+1) = TrueLong_Diff;
    end

    % Update the Buffer
    Buffer_FlybyData(k) = {PetalRotation_Data};

    % Clear variables for each petal rotation sequence
    clear TOF; clear OE_SC_end; clear TrueLong; clear TrueLong_Diff;
end

% Extract the accumulated change in True Longitude for Ranking purpose
for i = 1:size(Buffer_FlybyData, 2)
    PetalRotation_Data = cell2mat(Buffer_FlybyData(i));
    
    Num_Flys = size(PetalRotation_Data(end).Change_TrueLong, 1);

    if mod(Num_Flys, 2) ~= 0 % When flybys are located in the desired rotation direction (occurs every 2 flybys, but since the first row is the IC, we look at odd numbers)

        % Retrieve the change in the longitude of the flyby (location in moons orbit) 
        DTrue_Long(i) = PetalRotation_Data(end).Change_TrueLong(end);

        % Accumulated DV cost
        DV_cost(i)    = sum(cell2mat({PetalRotation_Data.DV_Cost}));   %[km/s]

    else
        DTrue_Long(i) = 0; % Assign the same value to all options since no comparison to be made in this flyby

        % Accumulated DV cost
        DV_cost(i)    = sum(cell2mat({PetalRotation_Data.DV_Cost}));   %[km/s]
    end
end


% Rank according to change in true longitude (primary objective) & DV defect cost the ones remaining
[B, idx] = sortrows([abs(DTrue_Long); DV_cost]', [1, 2]); 


% Re-arrange the Buffer of sequences according to the ranking
Buffer_PetalRots_Ranked = Buffer_PetalRotSeqs(idx);
Buffer_FlybyData_Ranked = Buffer_FlybyData(idx);

% Apply Beam Width 
Beam_Width = pars.PetalRot.Beam_Width;

if Beam_Width > numel(Buffer_PetalRots_Ranked )
    Buffer_PetalRots_Ranked  = Buffer_PetalRots_Ranked (1,1:end);
    Buffer_FlybyData_Ranked = Buffer_FlybyData_Ranked(1,1:end);
else
    Buffer_PetalRots_Ranked = Buffer_PetalRots_Ranked (1,1:Beam_Width);
    Buffer_FlybyData_Ranked = Buffer_FlybyData_Ranked(1,1:Beam_Width);
end

end