tic;

%% Script to mexify the spaceART toolbox
disp(' ');
disp(' ');
disp('-------------------------------------------');
disp('Deleting the mex files of the Space Toolbox');
disp('-------------------------------------------');

confirm = input('Are you sure (y/n)? ', 's');
if confirm ~= 'y' && confirm ~= 'Y'
    return
end
    
%% Extension and thus folder where the mex file shall be saved:
ext = mexext;
if isunix
    slash = '/';
    obj = '.o';
else
    slash = '\';
    obj = '.obj';
end

%% Delete
fprintf('\n');
outputDir = {'transfer_highThrust', 'conversion', ['conversion' slash 'time'], 'dynamics', 'ephemerides', 'orbital_transfers', 'relative_motion', 'swingby', 'tools'};

for i = 1:length(outputDir)
    disp(['Deleting ' outputDir{i} ' mex files...']);
    delete([outputDir{i} slash '*.mex*']);
end

disp(['Deleting object libraries...']);
delete(['*' obj]);

%%
disp(' ');
disp('All mex files and libraries deleted.');
disp('------------------------------------');
