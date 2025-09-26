function path2extLoadFile = createExtLoadsFile(path2enf, path2motGRF, trial_name, firstContactL, firstContactR)
% -------------------------------------------------------------------------
% This function creates an external load file for the COMAK workflow. It
% reads the force plate contacts from an Vicon Nexus specific database
% *.enf file (a simple text file) and creates the external load file using
% the open sim % API. This file assumes three force plates. This can be
% easily adjusted in the code below.

% INPUT:
% path2enf:     path to the *.enf file containing the force plate contacts
% path2motGRF:  path to the ground reaction force file (*.mot)
% trial_name:   name of the trial,(string)
% firstContact: indicates the side with which body segment was first in contact
%               with the force plate, e.g. calcn_l
% OUTPUT:       ext. loads file

% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Last changed:         04/2020
% -------------------------------------------------------------------------

% Import OpenSim API
import org.opensim.modeling.*

% Optional Inputs
if ~exist('firstContactL','var') && ~exist('firstContactR','var')
    firstContactL = 'calcn_l';
    firstContactR = 'calcn_r';
end
% Create folder information
[rootFolder,~]=fileparts(path2enf);

% Create external loads

% Get *.enf file containing force plate contact information
enf_filename = path2enf;

% Read the enf file to a cell array
fid = fopen(enf_filename,'rt');
enf_cell = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
enf_cell = enf_cell{1,1};

%% Find information about force plates in file
% NOTE: if you change this you need to change it also in prepareInputData!

enf_idx = find(not(cellfun('isempty',strfind(enf_cell,'FP'))));

% Several force plates are here assumed, add lines for addtional forces
% plates. You will have to add them in the next block "Harmonize ..." as well

for j = 1: length(enf_idx)
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP1'); FP1 = enf_cell{enf_idx(j)}(5:end); end
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP2'); FP2 = enf_cell{enf_idx(j)}(5:end); end
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP3'); FP3 = enf_cell{enf_idx(j)}(5:end); end
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP4'); FP4 = enf_cell{enf_idx(j)}(5:end); end
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP5'); FP5 = enf_cell{enf_idx(j)}(5:end); end
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP6'); FP6 = enf_cell{enf_idx(j)}(5:end); end
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP7'); FP7 = enf_cell{enf_idx(j)}(5:end); end
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP8'); FP8 = enf_cell{enf_idx(j)}(5:end); end
    if strcmp(enf_cell{enf_idx(j)}(1:3), 'FP9'); FP9 = enf_cell{enf_idx(j)}(5:end); end
end

%% Harmonize different enf file outputs and create FB cell for loop
if exist('FP1','var') && exist('FP2','var') && exist('FP3','var') && exist('FP4','var') && exist('FP5','var') && exist('FP6','var') && exist('FP7','var') && exist('FP8','var') && exist('FP9','var')
    FPs = {FP1, FP2, FP3, FP4, FP5, FP6, FP7, FP8, FP9};
elseif exist('FP1','var') && exist('FP2','var') && exist('FP3','var') && exist('FP4','var') && exist('FP5','var') && exist('FP6','var') && exist('FP7','var') && exist('FP8','var')
    FPs = {FP1, FP2, FP3, FP4, FP5, FP6, FP7, FP8};
elseif exist('FP1','var') && exist('FP2','var') && exist('FP3','var') && exist('FP4','var') && exist('FP5','var') && exist('FP6','var') && exist('FP7','var')
    FPs = {FP1, FP2, FP3, FP4, FP5, FP6, FP7};
elseif exist('FP1','var') && exist('FP2','var') && exist('FP3','var') && exist('FP4','var') && exist('FP5','var') && exist('FP6','var')
    FPs = {FP1, FP2, FP3, FP4, FP5, FP6};
elseif exist('FP1','var') && exist('FP2','var') && exist('FP3','var') && exist('FP4','var') && exist('FP5','var')
    FPs = {FP1, FP2, FP3, FP4, FP5};
elseif exist('FP1','var') && exist('FP2','var') && exist('FP3','var') && exist('FP4','var')
    FPs = {FP1, FP2, FP3, FP4};
elseif exist('FP1','var') && exist('FP2','var') && exist('FP3','var')
    FPs = {FP1, FP2, FP3};
elseif exist('FP1','var') && exist('FP2','var')
    FPs = {FP1, FP2};
elseif exist('FP1','var')
    FPs = {FP1};
end

% Replace other version of enf file info with the standard used here
FPs = strrep(FPs, 'FPL','Left');
FPs = strrep(FPs, 'FPR','Right');
FPs = strrep(FPs, 'FPI','Invalid');

%% Create external loads file
Loads = ExternalLoads();

% Loads.setName('Loads');
Loads.setDataFileName(path2motGRF); %grf-file

% Loads.setExternalLoadsModelKinematicsFileName(path2IK); %motion file
Loads.setLowpassCutoffFrequencyForLoadKinematics(6);

for i = 1:length(FPs)
    FPnr = num2str(i);
    switch FPs{i}
        case 'Invalid'
            % do nothing
        case 'Right'
            stor = Storage(10000, '');
            ForceR = ExternalForce(stor);
            ForceR.setName(strcat('FP',FPnr));
            ForceR.setAppliedToBodyName(firstContactR);
            ForceR.setForceIdentifier(strcat('ground_force_',FPnr,'_v'));
            ForceR.setPointIdentifier(strcat('ground_force_',FPnr,'_p'));
            ForceR.setTorqueIdentifier(strcat('ground_torque_',FPnr,'_'));
            ForceR.setPointExpressedInBodyName('ground');
            ForceR.setForceExpressedInBodyName('ground');
            Loads.cloneAndAppend(ForceR);
        case 'Left'
            stor = Storage(10000, '');
            ForceL = ExternalForce(stor);
            ForceL.setName(strcat('FP',FPnr));
            ForceL.setAppliedToBodyName(firstContactL);
            ForceL.setForceIdentifier(strcat('ground_force_',FPnr,'_v'));
            ForceL.setPointIdentifier(strcat('ground_force_',FPnr,'_p'));
            ForceL.setTorqueIdentifier(strcat('ground_torque_',FPnr,'_'));
            ForceL.setPointExpressedInBodyName('ground');
            ForceL.setForceExpressedInBodyName('ground');
            Loads.cloneAndAppend(ForceL);
    end
    % Create xml node
    names.ExternalLoads_setting = strcat(trial_name,'_ExternalLoads.xml');
end
% Write external loads xml file to root folder
cd(rootFolder);
Loads.print(names.ExternalLoads_setting);

% Output
path2extLoadFile = strcat(rootFolder,'\',names.ExternalLoads_setting);


%% Clear variables except output to prevet memory leak.
clearvars -except path2extLoadFile
end