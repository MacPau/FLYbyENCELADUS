function [Buffer_PetalRot_Seqs, Buffer_PetalRot_idx, Buffer_FlybyData] = Transfer_Combinations(PetalRotSeqs_In, PseudoResos_Options, Buffer_FlybyData_in, pars)

% This function computes the different combinations of transfers that can
% be used to perform a petal rotation.

% Author: J.C Garcia Mateas
% Last revision: 02/11/2024

%% INPUTS %%
% - PetalRotSeqs_In: 3D Matrix of size [A x B x C], where A is the n° of
%                    transfers in to the petal rotation sequence, B is the
%                    number of fields containing the data (Vinf, alfa,
%                    crank, TOF etc) of the transfer and C is the n° of
%                    petal rotation sequences stored.
%
% - PseudoResos_Options: Matrix of size [N x 10], where each row is a 
%                        pseudo-resonant transfer, with  column 1 being 
%                        the n° of Gravity Assist Body revolutions,
%                        column 2 the n° of spacecraft revolutions, col 3 the
%                        V_inf [km/s], cols 4  & 5 the Pump and Crank angle 
%                        respectively in [rad] of the incoming node, columns 
%                        6, 7, 8 the Vinf, pump & crank of the outgoing node,
%                        column 9 the TOF in [seconds] and column 10 the type
%                        of transfer [18] = Inbound-/Outbound; 
%                        [81] = Outbound/Inbound.
%
% - Buffer_FlybyData_in: Cell of size [1xC], where C is the n° of petal rotation 
%                        sequences stored. Each cell is a structure with
%                        several fields containing the information
%                        regarding the flybys involved in the petal
%                        rotation.
%
% - pars: Structure containing problem constants and parameters

%% OUTPUTS %%
% - Buffer_PetalRot_Seqs: Cell of size 1 x N, where N = n° of petal rotation
%                         sequences stored. Each cell is a matrix of size Y
%                         x 10, where Y = n° of transfers in the petal
%                         rotation and the 10 columns are the parameters
%                         characterizing the transfer, ordered as detailed
%                         in the "PseudoResos_Options" input description.
%
% - Buffer_PetalRot_idx: Cell of size 1 x N. Each cell is a row vector of
%                        size Y containing the indices from the
%                        PseudoResos_Options matrix which are used to
%                        construct each sequence (useless for now).
%
% - Buffer_FlybyData: Cell of size [1xN]. Each cell is a structure with
%                        several fields containing the information
%                        regarding the flybys involved in the petal
%                        rotation. Each row in each structure corresponds
%                        to a transfer of the petal rotation.

%% FUNCTION %%
% Initialize Output
Buffer_PetalRot_Seqs = [];
Buffer_PetalRot_idx = [];
Buffer_FlybyData = [];

% parfor i = 1:size(Seqs_Y2, 1)
for i = 1:size(PseudoResos_Options,1)

    % Retrieve Pseudo-resonant Transfer option
    transfer_out = PseudoResos_Options(i,:);
    node_out     = transfer_out(3:5);

    % Retrieve Pseudo-resonant Transfer outgoing node info (vinf, alpha, k)
    vinfou  = node_out(1);
    alphaou = node_out(2);
    kou     = node_out(3);

    for j = 1:size(PetalRotSeqs_In, 3)

        % Retrieve the full existing Petal Rotation sequence
        Full_transfer = PetalRotSeqs_In(:,:,j);

        % Retrieve the last transfer from the Petal Rotation sequence
        transfer_in = PetalRotSeqs_In(end,:,j);

        % Check that this combination involves different pseudo-resonances
        if ~isequal(transfer_in, transfer_out)
            node_in     = transfer_in(7:9);
    
            % Retrieve incoming node info (vinf, alpha, k)
            vinfin  = node_in(1);
            alphain = node_in(2);
            kin     = node_in(3);
    
            % Change in Crank
            Delta_Crank = abs(kou - kin); 
    
            % Compute required bending change for flyby to be feasible
            d_req = acos(cos(alphain)*cos(alphaou) + sin(alphain)*sin(alphaou)*cos(Delta_Crank)); %[rad]
    
            % Retrieve max. bending at defined Vinf 
            id = find(pars.INPUTS.V_inf == vinfin);
    
            % Compute DV defect
            vvrel_A = vinfAlphaCrank_to_VinfTCN(vinfin, alphain, kin); % [km/s] Vinf velocity vector at arrival to the Moon (before flyby)
            vvrel_D = vinfAlphaCrank_to_VinfTCN(vinfou, alphaou, kou); % [km/s] Vinf velocity vector at departure from the Moon (after flyby)
            dv      = findDV(vvrel_A, vvrel_D, pars.Moon.mu, pars.rp_flyby); % [km/s] DV defect
    
            % Check flyby feasibility & associated DV cost
            if d_req <= pars.delta_max(id)
    
                % Re-write delta_max for correct usage in Flyby_BuildUp function
                delta_max_aux = pars.delta_max;
                pars.delta_max = pars.delta_max(id);
    
                % Compute the parameters associated to the flyby
    %             Flyby_Data = Flyby_BuildUp(node_in, node_out, pars);
                Flyby_Data = Flyby_Elements(node_in, node_out, pars);
    
                % Assign DV cost
                Flyby_Data.DV_Cost = 0;
    
                if isempty(Buffer_FlybyData_in)
    
                    % Store the outputs
                    Buffer_PetalRot_Seqs{end+1} = [Full_transfer; transfer_out];
                    Buffer_PetalRot_idx{end+1} = [j, i];
                    Buffer_FlybyData{end+1} = Flyby_Data;
    
                else
    
                    % Retrieve the Flyby Data for this sequence
                    Flyby_Data_prev = cell2mat(Buffer_FlybyData_in(j));
    
                    % Store the outputs
                    Buffer_PetalRot_Seqs{end+1} = [Full_transfer; transfer_out];
                    Buffer_PetalRot_idx{end+1} = [j, i];       
                    Buffer_FlybyData{end+1} = [Flyby_Data_prev,  Flyby_Data];
    
                end
    
                % Rewrite delta_max after usage in Flyby_BuildUp function
                pars.delta_max = delta_max_aux;
    
    
            %%% IF DV IS RESPECTED %%     
            elseif dv <= pars.INPUTS.maxDV_Defect
                % Re-write delta_max for correct usage in Flyby_BuildUp function
                delta_max_aux = pars.delta_max;
                pars.delta_max = pars.delta_max(id);
    
                % Compute the parameters associated to the flyby
    %             Flyby_Data = Flyby_BuildUp(node_in, node_out, pars);
                Flyby_Data = Flyby_Elements(node_in, node_out, pars);
    
                % Assign DV cost
                Flyby_Data.DV_Cost = dv;
    
                if isempty(Buffer_FlybyData_in)
    
                    % Store the outputs
                    Buffer_PetalRot_Seqs{end+1} = [Full_transfer; transfer_out];
                    Buffer_PetalRot_idx{end+1} = [j, i];
                    Buffer_FlybyData{end+1} = Flyby_Data;
    
                else
    
                    % Retrieve the Flyby Data for this sequence
                    Flyby_Data_prev = cell2mat(Buffer_FlybyData_in(j));
    
                    % Store the outputs
                    Buffer_PetalRot_Seqs{end+1} = [Full_transfer; transfer_out];
                    Buffer_PetalRot_idx{end+1} = [j, i];       
                    Buffer_FlybyData{end+1} = [Flyby_Data_prev,  Flyby_Data];
    
                end
    
                % Rewrite delta_max after usage in Flyby_BuildUp function
                pars.delta_max = delta_max_aux;
            end
        end
    end
end

end