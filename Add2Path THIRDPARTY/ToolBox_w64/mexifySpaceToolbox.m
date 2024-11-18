function mexifySpaceToolbox()
% mexifySpaceToolbox - Compiles and creates MEX files for the Space
%   Toolbox.
%
% AUTHOR:
%   Matteo Ceriotti, 18/12/2009

tic;

%% Script to mexify the Space Toolbox
disp(' ');
disp(' ');
disp('----------------------------------------------');
disp('Creation of the mex files of the Space Toolbox');
disp('----------------------------------------------');

%% Extension and thus folder where the mex file shall be saved:
if isunix
    currentdir = '.';
    slash = '/';
    obj = '.o';
else
    currentdir = '.';
    slash = '\';
    obj = '.obj';
end

%% Libraries
fprintf('\n');
disp('CREATING LIBRARIES');

str = ['mex ' ...
    'src' slash 'transfer_highThrust.c ' ...
    'src' slash 'mathUtils.c ' ...
    'src' slash 'astroConstants.c ' ...
    'src' slash 'dynamics.c ' ...
    'src' slash 'keplerianMotion.c ' ...
    'src' slash 'anomConv.c ' ...
    'src' slash 'conversion.c ' ...
    'src' slash 'ephNEO.c ' ...
    'src' slash 'ephemerides.c ' ...
    '-c -DMEXCOMPILE'];
    %'src' slash 'relative_motion.c ' ...
    
eval(str);
fprintf('done\n');

%% transfer_highThrust
outputDir = 'transfer_highThrust';
if ~exist([currentdir slash outputDir],'dir'), mkdir(outputDir); end
funName = {'lambertMR'};
libs = {'transfer_highThrust'};

callMex(funName, libs, outputDir, slash, obj);

%% Ephemerides
outputDir = 'ephemerides';
if ~exist([currentdir slash outputDir],'dir'), mkdir(outputDir); end
funName = {'ephNEO', 'astroConstants', 'uplanet', 'ephSS_car', 'ephSS_kep'};
libs = {'ephemerides', 'conversion', 'astroConstants', 'mathUtils', 'ephNEO', 'anomConv'};

callMex(funName, libs, outputDir, slash, obj);

%% Conversion
outputDir = 'conversion';
if ~exist([currentdir slash outputDir],'dir'), mkdir(outputDir); end
funName = {'car2kep', 'kep2car'};
libs = {'conversion', 'mathUtils'};

callMex(funName, libs, outputDir, slash, obj);

%% Conversion/Time
outputDir = ['conversion' slash 'time'];
if ~exist([currentdir slash outputDir],'dir'), mkdir(outputDir); end
funName = {'date2jd', 'date2mjd', 'date2mjd2000', 'fracday2hms', 'jd2date', ...
    'jd2mjd', 'jd2mjd2000', 'mjd2date', 'mjd2jd', 'mjd2mjd2000', 'mjd20002date', ...
    'mjd20002jd', 'mjd20002mjd'};
libs = {'conversion', 'mathUtils'};

callMex(funName, libs, outputDir, slash, obj);

%% Keplerian Motion
outputDir = 'keplerianMotion';
if ~exist([currentdir slash outputDir],'dir'), mkdir(outputDir); end
funName = {'kepEq_f', 'kepEq_t', 'kepPro'};
libs = {'keplerianMotion', 'mathUtils', 'conversion', 'ephemerides', 'ephNEO', 'anomConv'};

callMex(funName, libs, outputDir, slash, obj);

%% Tools
outputDir = 'tools';
if ~exist([currentdir slash outputDir],'dir'), mkdir(outputDir); end
funName = {'qck','cartProd'};
libs = {'mathUtils'};

callMex(funName, libs, outputDir, slash, obj);

%% Delete libraries
fprintf('\n');
disp('DELETING LIBRARIES');
fprintf('done\n');
delete(['*' obj]);

%%
disp(' ');
disp(['All mex files created and saved in ' num2str(ceil(toc)) ' seconds.']);
disp('----------------------------------------------');
disp(' ');
disp(' ');

%%
function callMex(funName, libs, outputDir, slash, obj)

%%
mexName = cell(1,length(funName));
for i = 1:length(funName)
    mexName{i} = ['mex' upper(funName{i}(1)) funName{i}(2:end)];
end

%%
fprintf('\n');
for i = 1:length(funName)
    fprintf([funName{i} ': ']);
    str = ['mex src' slash mexName{i} '.c '];
    for j = 1:length(libs)
        str = [str libs{j} obj ' '];
    end
    str = [str '-outdir ' outputDir ' -output ' funName{i} ' -DMEXCOMPILE'];
    eval(str)
    fprintf('done\n');
end