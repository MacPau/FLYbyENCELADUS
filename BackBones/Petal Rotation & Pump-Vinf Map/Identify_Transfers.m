function [Reso_Options, PseudoReso_Opts] = Identify_Transfers(Reso_TransfersDatabase, PseudoReso_TransfersDatabase, pars)

% This function identifies the transfers involved in a petal rotation
% provided a user-defined starting pure resonance, Vinf, crank angle and/or
% rotation direction (clockwise/counterclockwise).

% Author: J.C. Garcia Mateas
% Last revision: 02/11/2024

%% INPUTS %%
% - Reso_TransfersDatabase: Matrix of size Nx10 containing the database of
%                           resonant transfers. The columns are: col1 = n°
%                           of moon revolutions; col2 = n° of spacecraft
%                           revs; col3 = Vinf [km/s]; col4 = pump angle [rad] 
%                           after the flyby; col 5: crank angle [rad] after
%                           the flyby; col 6: Time of flight [seconds]
%                           between the two flybys of the resonant
%                           transfer; col 7 = Vinf [km/s] at next flyby
%                           entrance; col 8 = pump angle [rad] at next
%                           flyby entrance; col 9 = crank angle [rad] at
%                           next flyby entrance; col 10 = type of transfer
%                           (O = Outbound/Outbound / 1 = Inbound/Inbound).
%
% - PseudoReso_TransferDatabase: Matrix of size B x 10 containing the data
%           for all the pseudo-resonant transfers being considered. The
%           columns contain the same information as the one described above
%           for Reso_TransfersDatabase.
%
% - pars: Structure containing problem constants and parameters

%% OUTPUTS %%
%
% - Reso_Options: Matrix of size A x 10, where A is the number of resonant
%                 transfers options to be used in the petal rotation. The
%                 columns correspond to the same information as described
%                 above for Reso_TransfersDatabase.
%
% - PseudoReso_Opts: Matrix of size B x 10, where B is the n° of
%       pseudo-resonant transfers which can be used in the petal rotation.
%       The columns correspond to the same information as described above.

%% FUNCTION %%
% Starting Resonance, Vinf and Crank
Reso  = pars.PetalRot.Starting_Reso;
Vinf  = pars.PetalRot.Starting_Vinf;
Crank = pars.PetalRot.Starting_Crank;   % [rad] Crank angle at start of Petal Rotation (0 = Outbound / 180° = Inbound)
Direction = pars.PetalRot.Direction;    % Direction in which to rotate (+1 = clockwise / -1 = anticlockwise)


% Find resonant transfers which have the user defined ratio, Vinf & crank
idx = Reso_TransfersDatabase(:,1) == Reso(1) & Reso_TransfersDatabase(:,2) == Reso(2) & Reso_TransfersDatabase(:,3) == Vinf & Reso_TransfersDatabase(:,5) == Crank;
Reso_Options = Reso_TransfersDatabase(idx,:);

% Find pseudo-resonances having the same Vinf
idx2 = PseudoReso_TransfersDatabase(:,3) == Reso_Options(1,3);
PseudoReso_Opts = PseudoReso_TransfersDatabase(idx2,:);

end