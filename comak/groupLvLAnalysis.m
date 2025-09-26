clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run group level analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is intended to support the group level analysis of the output
% from the comak workflow. It gets all data from a set of data
% (condition/prefix) and collects all data in a single struct. For this
% purpose it reduces the data on a per-subject level by eg. taking the
% average of all available trials.
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Last changed:         04/2022


% TodDos / Notes:
% - ...
% - ...


%% User Settings ---------------------------------------------------------
% Root directory containing all working directories to process.
rootDirectory = 'D:\ResearchStayTestingGround\TestDataNorm\DFN01yyyy-mm-dd\2020-07-23\comak-groupData\';

% Decide which methods to use to reduce trials of one individual, e.g. by
% using Sangeux's "most representative Method".
aggregateData = 'avg-outlierReduced'; % mostRepTrial, firstTrialOnly, avg-all-omitnan, avg-outlierReduced
thresholdS = 2;

% Note: this functions does not work correctly if NO prefix was used. It might then use all files, including the ones with a prefix
condition = 'Dynamic'; %,'mSEBT_PM','mSEBT_PL'};     % default = 'Dynamic' or one of these: 'mSEBT_AT','mSEBT_PM','mSEBT_PL'
prefix = 'standardone';         % default = 'standard',
files2skip = '';                % default = {''}; or {'Abc', 'Cde', 'Efg'}

% Set thresholds for IK errors
maxIKcrit = 2.8;
rmsIKcrit = 2;

% Decide which side to plot
sides2Plot = {'Right', 'Left'};

useYLim = false;
myxLabel = '% gait cycle';
%% Shared plotting settings

% Line width
line_width_consistency = 1;
line_width_avg  = 1.8;
line_width_std  = 1.2;

% Line style
lineStyle_avg = '-';
lineStyle_std = '--';
lineStyle_consistency = '-';

% Line color
line_alpha_avg = 0.8;
line_alpha_consistency = 0.3;
ColAvg =    [0, 0, 0, line_alpha_avg]; % black
ColData =   [1, 0, 0, line_alpha_consistency; ... % red
             0, 0, 1, line_alpha_consistency]; % blue

shadings = true;
consistency = true;

%% Add repo path ---------------------------------------------------------
% Add repo path to workflow
% Get path where start file is located on disc and set that as the repoPath
% Note that the start file needs to to located in its original location,
% within the repo top lvl folder.
tmpPath = mfilename('fullpath');
[repoPath, ~, ~] = fileparts(tmpPath);
repoPath = strcat(repoPath,'\');

% Add all repo subfolders to the Matlab path
addpath(genpath(fullfile(repoPath)));

%% PREPARE some things ----------------------------------------------------
% Toggle debug mode
%dbstop if error
%dbclear if error

% Go to folder
cd(rootDirectory);
Folder = cd;

% Find all relevat files
tmpList = dir(fullfile(Folder, strcat('*', prefix,'*', condition, '*-JAM-Results-All-Trials.mat')));

% Make big overall group lvl file with all data
[groupLvl_Data, usedFileList] = getAllVariables(tmpList, aggregateData, files2skip, thresholdS, maxIKcrit, rmsIKcrit);

% Rename fieldnames for plotting ...
groupLvl_Data = renameStructFields(groupLvl_Data);

%% Plot the data of the left and right side
close all

for i = 1 : length(sides2Plot)

    % Necessary if you rerun this block
    vars2delete = who('ax*');
    if ~isempty(vars2delete)
        clear(vars2delete{:});
    end

    % Set side
    sideLong = lower(sides2Plot{i});
    side = sideLong(1);
    
    % Set color
    ColTmp = ColData(i,:);
    
    % Set figure
    hf = figure;
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 40, 25]);
    N = num2str(width(groupLvl_Data.(side).ContactForces.tf_contact_casting_total_contact_force_y)); % This needs a rethink!
    sgtitle(char(strcat('<', sideLong ,'>', 'side only. N = ', {' '}, N)));
    rowN = 8; % number of rows of subplots
    colN = 3; % number of cols of subplots
    idx_sp = 1;

    % Patello Femoral Contact Force fromressure
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactPressure.PFcontactForceTot, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([-0.5 0.5]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'Total Patello-femoral', 'Contact Force', '(from Pressure)', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Shear force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactPressure.PFcontactForceMed, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 0.1]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'PF Medial Contact Force', '(from Pressure)', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Shear force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactPressure.PFcontactForceLat, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 0.1]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({' PF Lateral Contact Force', '(from Pressure)', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;
    
    % Patello Femoral Contact Pressure Med/Lat
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    dat = (groupLvl_Data.(side).ContactPressure.PFpeakPressureTot.Variables)/1000000; % to Mega Pascal
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 10]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'patello-femoral', 'Total peak pressure', '(from Pressure)', '[Mega Pascal]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % medial peak pressure
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    dat = (groupLvl_Data.(side).ContactPressure.PFpeakPressureMed.Variables)/1000000; % to Mega Pascal;
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 10]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'patello-femoral', 'Medial peak pressure', '(from Pressure)', '[Mega Pascal]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Lateral peak pressure
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    dat = (groupLvl_Data.(side).ContactPressure.PFpeakPressureLat.Variables)/1000000; % to Mega Pascal;
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 10]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'patello-femoral', 'Lateral peak pressure', '(from Pressure)', '[Mega Pascal]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;    
    
    % Tibio Femoral Contact Pressure
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactPressure.TFvForceTot, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 6]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'tibio-femoral', 'Vertcial Contact Force', '(from Pressure)', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Shear force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactPressure.TFvForceMed, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 6]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'Medial Contact Force', '(from Pressure)', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Shear force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactPressure.TFvForceLat, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 6]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'Lateral Contact Force', '(from Pressure)', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Tibio Femoral Contact Force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactForces.tf_contact_casting_total_contact_force_y, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81))*-1;
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 6]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'tibio-femoral', 'Vertcial Contact Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Shear force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactForces.tf_contact_casting_total_contact_force_x, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 6]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'x Contact Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Shear force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactForces.tf_contact_casting_total_contact_force_z, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([-0.5 0.5]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'z Contact Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Patella Femoral Contact Force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactForces.pf_contact_casting_total_contact_force_y, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([-0.3 0.3]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'patella-femoral', 'Inf-superior PF Shear Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Shear force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactForces.pf_contact_casting_total_contact_force_x, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 2]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'PF Contact Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Shear force
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).ContactForces.pf_contact_casting_total_contact_force_z, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81))*-1;
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([-0.5 0.5]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'Medio-Lateral Shear Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Inverse Kinematicss
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    dat = groupLvl_Data.(side).KinematicsComak.knee_flex_ipsilateral.Variables;
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    yline(0,'--');
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([-20 70]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'KNEE', 'Flex.-Ext. Angle', '[deg.]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Inverse Dynamics
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    [var1, var2] = getHeaderAlignedDataSet(groupLvl_Data.(side).InverseDynamics.knee_flex_moment_ipsilateral, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables));
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    yline(0,'--');
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([-1 1]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'Flex.-Ext. Moment', '[N/kg]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Inverse Kinematics
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    dat = groupLvl_Data.(side).KinematicsComak.hip_flex_ipsilateral.Variables;
    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    yline(0,'--');
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([-40 40]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'HIP', ' Flex.-Ext. Angle', '[deg.]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % ACL Ligament Load
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    data =  groupLvl_Data.(side).ContactForces.ACLam1_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLam2_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLam3_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLam4_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLam5_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLam6_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLpl1_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLpl2_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLpl3_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLpl4_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLpl5_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.ACLpl6_force_total.Variables;

    data = array2table(data,'VariableNames',groupLvl_Data.(side).ContactForces.ACLam1_force_total.Properties.VariableNames);
    [var1, var2] = getHeaderAlignedDataSet(data, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));

    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 1]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'ACL (am + pl)', 'Total Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % PCL Ligament Load
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    data =  groupLvl_Data.(side).ContactForces.PCLal1_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLal2_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLal3_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLal4_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLal5_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLpm1_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLpm2_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLpm3_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLpm4_force_total.Variables + ...
        groupLvl_Data.(side).ContactForces.PCLpm5_force_total.Variables;

    data = array2table(data,'VariableNames',groupLvl_Data.(side).ContactForces.PCLal1_force_total.Properties.VariableNames);
    [var1, var2] = getHeaderAlignedDataSet(data, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));

    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 1]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'PCL (al + pl)', 'Total Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 2;

    % Muscle Forces (extensors)
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    data =   groupLvl_Data.(side).MuscleForces.forceset_tfl.Variables + ...
             groupLvl_Data.(side).MuscleForces.forceset_vaslat.Variables + ...
             groupLvl_Data.(side).MuscleForces.forceset_vasmed.Variables + ...
             groupLvl_Data.(side).MuscleForces.forceset_vasint.Variables;

    data = array2table(data,'VariableNames',groupLvl_Data.(side).MuscleForces.forceset_tfl.Properties.VariableNames);
    [var1, var2] = getHeaderAlignedDataSet(data, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));

    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    if useYLim; ylim([0 3]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'KNEE Extensors', 'Muscle Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;

    % Muscle Forces (flexors)
    if exist(strcat('ax',num2str(idx_sp)),'var'); eval(strcat('subplot(ax',num2str(idx_sp),')')); else; eval(strcat('ax',num2str(idx_sp),'= subplot(rowN,colN,idx_sp);')); end
    data =  groupLvl_Data.(side).MuscleForces.forceset_bflh.Variables + ...
        groupLvl_Data.(side).MuscleForces.forceset_bfsh.Variables + ...
        groupLvl_Data.(side).MuscleForces.forceset_semimem.Variables + ...
        groupLvl_Data.(side).MuscleForces.forceset_semiten.Variables + ...
        groupLvl_Data.(side).MuscleForces.forceset_gaslat.Variables + ...
        groupLvl_Data.(side).MuscleForces.forceset_gasmed.Variables + ...
        groupLvl_Data.(side).MuscleForces.forceset_sart.Variables;

    data = array2table(data,'VariableNames',groupLvl_Data.(side).MuscleForces.forceset_bflh.Properties.VariableNames);
    [var1, var2] = getHeaderAlignedDataSet(data, groupLvl_Data.(side).MetaData.Bodymass);
    dat = (var1.Variables ./ (var2.Variables .*9.81));

    [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat);
    hold on;
    if consistency; plot(dat, 'LineStyle', lineStyle_consistency, 'LineWidth', line_width_consistency, 'Color', ColTmp); end
    plot(avgCurve, 'LineStyle', lineStyle_avg, 'LineWidth', line_width_avg, 'Color', ColAvg);
    plot(avgPlusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    plot(avgMinusSDCurve, 'LineStyle', lineStyle_std, 'LineWidth', line_width_std, 'Color', ColAvg);
    if shadings; shade(1:length(dat),avgMinusSDCurve,1:length(dat),avgPlusSDCurve,'FillType',[2 1],'FillColor', ColTmp(1:3), 'FillAlpha', ColTmp(4),'LineStyle', 'none'); end
    %if useYLim; ylim([0 5]); end
    xlim([0 100]);
    xlabel(myxLabel,'FontWeight','bold');
    ylabel({'KNEE Flexors', 'Muscle Force', '[N/(kg * g)]'},'FontWeight','bold');
    idx_sp = idx_sp + 1;



    %% Save plot
    %exportgraphics(hf, [rootDirectory, sideLong,'_','Overview_MeanConsistency.png'],'Resolution','600');

end