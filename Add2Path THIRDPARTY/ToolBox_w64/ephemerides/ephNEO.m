% ------------------------ MEX GATEWAY FUNCTION ---------------------------
%
% ephNEO.m - Ephemerides of NEOs read from a text file database.
%
% PROTOTYPE:
%   [kep, mass, M] = ephNEO(time, id)
%   ephNEO('load', filename)
%   ephNEO('unload')
%
% DESCRIPTION:
%   This function returns the orbital parameters, the mass, the mean
%   anomaly and the name of NEOs. Each NEO is identified by an id.
%   NEO data is read from a database of ephemerides.
%   The database is read the first time the function is invoked and then
%   stored into persistent memory.
%   If no database is loaded explicitly, the function tries to load the
%   ephNEO_data.txt.
%   The ID of the first NEO in the database is conventionally 12.
%
%   There are 3 ways of calling this function:
%       [kep, mass, M] = ephNEO(time, id)
%           If no database is loaded, loads the default database and
%           returns the ephemerides.
%       ephNEO('load', filename)
%           Loads the database specified in filename. The input filename is
%           optional, if not specified, tries to read the default database.
%       ephNEO('unload')
%           Unloads the database from memory, freeing the persistent
%           memory allocation.
%
%   The data file should be:
%       numberOfNEOs
%       n a e i Omega omega meanAnom epoch mag mass refDate name comment
%       n a e i Omega omega meanAnom epoch mag mass refDate name comment
%       n a e i Omega omega meanAnom epoch mag mass refDate name comment
%   Where:
%       numberOfNEOs    Number of objects in database.
%       n               Id of NEO. The first ID should be 12. Not actually
%                       used.
%       a               Semimajor axis, AU.
%       e               Eccentricity.
%       i               Inclination, deg.
%       Omega           Right ascension of ascending node, deg.
%       omega           Anomaly of pericentre, deg.
%       meanAnom        Mean anomaly at epoch, deg.
%       epoch           Epoch of ephemeris, MJD.
%       mag             Magnitude. Used to compute the mass, if the one in
%                       the database is 0.
%       mass            Mass, kg. If 0, then mag is used to estimate the mass.
%       refDate         Date at which the ephemeris was added, or any other
%                       reference date.
%       name            Name (max 80 characters, no spaces allowed)
%       comment         Any comment, even with spaces.
%   Note:
%       All floating point numbers, except numberOfNEOs and n, must contain
%       the decimal separator ".": so zero must be "0.0", etc. Exponential 
%       notation is possible.
%
% INPUT:
%   time    Epoch, in MJD2000 [d].
%   id      NEO identifier. It is a number starting from 12 (because the
%           identifiers from 1 to 11 are reserved for the Solar System).
%
% OUTPUT:
%   kep     Keplerian parameters. It is a 6 entry vector:
%               [a e i Om om wom]
%           where:
%               a is the semimajor axis [km];
%               e is the eccentricity;
%               i is the inclination [rad];
%               Om is the anomaly of the ascending node [rad];
%               om is the anomaly of the pericentre [rad];
%               wom is the true anomaly (from the pericentre) [rad].
%   mass    Mass of the NEO [kg]. It can be read from the database, or, if
%           not available, estimated by an approximate equation.
%   M       Mean anomaly at time [rad].
%   
% NON-STANDARD LIBRARIES:
%   mathUtils.h, astroConstants
%
% ORIGINAL VERSION:
%   Paolo de Pascale, 01/11/2004, MATLAB, NeoEphemeris.m
%   
% AUTHOR:
%   Matteo Ceriotti, 09/12/2009
%
% CHANGELOG:
%
% -------------------------------------------------------------------------
