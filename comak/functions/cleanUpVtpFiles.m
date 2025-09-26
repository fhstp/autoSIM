function cleanUpVtpFiles(vtp2keep, rootDirectory, prefix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cleaning  *.vtp files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used to find and deleted all user selected *.vtp files in
% the jam subfolders only to save harddisc space.
%
% Define which *.vtp files I want to keep, the rest will be deleted to save hard disc space (if selected in settins)
% vtp2keep = 		{'femur_cartilage', 'tibia_cartilage', 'patella_cartilage', ...
%				     'femur_bone', 'tibia_bone', 'patella_bone', 'ACL', 'PCL'};
%
% rootDirectory = 'D:\BriansSimulationGround\WalkingInVR_SEBT\';
% prefix = 		  'standard';
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
%
%
% Last changed:         04/2023
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Delete *.vtp files which I don`t need to save hard disc space

%% Prepare input
rootDirectory = fullfile(rootDirectory);

% Start timer
tStartVtp = tic;
disp('**********************************************************');
disp('Starting to clean *.vtp files ...');
disp('**********************************************************');

% Find dirs
dir2check = dir(fullfile(rootDirectory, '\**\jam\*.vtp'));

% Delete files
vtpFiles = dir2check; % Find all *.vtp files
delCnt = 0;
for i_vtp = 1:length(vtpFiles)
    % Define file to delte
    file2delete = fullfile(vtpFiles(i_vtp).folder,vtpFiles(i_vtp).name);
    % Check if file should be kept or delted
    if ~contains(file2delete, vtp2keep)
        delete(file2delete)
        delCnt = delCnt + 1;
    end
end

disp(' '); % New line
disp(char(strcat('>>>>>', {' '}, num2str(delCnt), {' '}, '*.vtp files deleted by user.')));


%% FINAL %%
% End timer
tEndVtp = toc(tStartVtp);
disp('  ');
disp('**********************************************************');
fprintf('>>>>> *.vtp CleanUp duration: %d minutes and %0.f seconds.\n', floor(tEndVtp/60), rem(tEndVtp,60));
disp('**********************************************************');