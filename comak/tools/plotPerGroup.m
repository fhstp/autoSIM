close all
path2NormKCF = 'C:\LokaleDaten\CloudStation\MyPapers\2020_GuidedGrowth\code\osimProcessing-GitLab\setupFiles\Norm\OrthoLoad\OrthoLoad_KCFnorm.mat';
condition = 'Dynamic';

%% Effect of contact energy
workingDirectory = fullfile('D:\TestingEffektSecConstrainSim\OSS\');
prefixesCE = {'CE200-002','CE200-004','CE200-008','CE200-009','CE200-010','CE200-011'};
myTitle = 'Effects of different ContactEnergy Settings';
plotKCFPerGroup(workingDirectory, prefixesCE, path2NormKCF, myTitle);

%% Effect of frontal knee alignment
workingDirectory = fullfile('C:\LokaleDaten\CloudStation\MyPapers\2020_GuidedGrowth\data\AK-data\VarValg-experiments_single\');
prefixesVarusValgus = {'standard', 'var5', 'var3', 'var2', 'valg2', 'valg3', 'valg5',};
myTitle = 'Effects of different frontal knee alignments';
plotKCFPerGroup(workingDirectory, prefixesVarusValgus, path2NormKCF, myTitle);

%% Effect of different Splines
workingDirectory = fullfile('D:\TestingEffektSecConstrainSim\OSS\');
prefixesCE = {'CE200-orig','CE200-002','CE200-004','CE200-008','CE200-009','CE200-010','CE200-011'};
myTitle = 'Effects of different Splines';
plotKCFPerGroup(workingDirectory, prefixesCE, path2NormKCF, myTitle);


%% Plot
function [h] = plotKCFPerGroup(workingDirectory, prefixes, path2NormKCF, myTitle)
% This function looks for all comak output files in the working directory
% and plots all of the data group-wise. It does not take care about left or
% right steps.

% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
side = 'r';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Part one - plotting based on contact force

% Define body side to plot
side2plot = strcat('_', side, '_');

% Find all results files in working directory
files = dir(fullfile(strcat(workingDirectory, 'JAM\**\', '*FinalResultsStruct.mat')));

% Make list form struct
for i = 1 : length(files)
    filesTmp{i} = files(i).name;
    folderTmp{i} = files(i).folder;
end

k = 1;
for i = 1 : length(prefixes)
    
    % Get files matching the prefix & the side
    dat = filesTmp(contains(filesTmp, prefixes(i)));
    dat = dat(contains(dat, side2plot));
    datFolder = folderTmp(contains(filesTmp, prefixes(i)));
    datFolder = datFolder(contains(datFolder, side2plot(1:end-1))); % folders end with _r or _l
    
    % Load file
    plottingDat = load(char(fullfile(datFolder, dat)));
    
    % Load data - in some cases the array is only 100 frames long - for
    % visualization purpose I interpolated to 101;
    cf_v(:,k) = interpft(plottingDat.results.ContactForces.('tf_contact.casting.total.contact_force_y')*-1,101)/plottingDat.results.metadata.BW;
    cf_ap(:,k) = interpft(plottingDat.results.ContactForces.('tf_contact.casting.total.contact_force_x'),101)/plottingDat.results.metadata.BW;
    cf_ml(:,k) = interpft(plottingDat.results.ContactForces.('tf_contact.casting.total.contact_force_z'),101)/plottingDat.results.metadata.BW;
    legendNames{k} = prefixes{i};
    TO(k) = (plottingDat.results.metadata.TO - plottingDat.results.metadata.IC)/(plottingDat.results.metadata.ICi - plottingDat.results.metadata.IC);
    k = k+1;
    
end

%% PLot data
hfig = figure;
set(hfig,'units','centimeters','position',[0,0,28,16]);
sgtitle(myTitle,'FontWeight','bold');

% PLot norm first
% Load NormData
load(path2NormKCF);
KFCfn = {'Fz','Fy','Fx'};
TOidx = round(mean(TO)*100);
% Select contact forces (3)
tf_contactF = {'tf_contact.casting.total.contact_force_y','tf_contact.casting.total.contact_force_x','tf_contact.casting.total.contact_force_z'};

for i = 1:3

% Plot OrthoLoad norm data
    % Adjust sign for norm data
    curveKFC_mean_stance = interpft(mean(KFCnorm.NormBW.(KFCfn{i}),2),TOidx);
    
    if contains(tf_contactF{i},'force_y') || contains(tf_contactF{i},'force_x')
        curveKFC_mean_stance = curveKFC_mean_stance*-1;
    end
    
    % Prepare data
    curveKFC_mean_stance = padarray(curveKFC_mean_stance,101-TOidx,'post')/100; % padd zeros to end, change to multipel of BW
    KFC_sd_stance = interpft(std(KFCnorm.NormBW.(KFCfn{i}),0,2),TOidx)/100; % change to multipel of BW
    KFC_sd_stance = padarray(KFC_sd_stance,101-TOidx,'post'); % padd zeros to end
    
    % Plot norm
    forces_time = 1:101;
    subplot(1,3,i)
    hold on
    plot(forces_time,curveKFC_mean_stance,'k','LineWidth',1);
    plot(forces_time,curveKFC_mean_stance+2*KFC_sd_stance,'k--','LineWidth',0.5);
    plot(forces_time,curveKFC_mean_stance-2*KFC_sd_stance,'k--','LineWidth',0.5);
end

% PLot data
subplot(2,3,1)
n = size(cf_v,2);
for i = 1:1:n
    lh = plot(cf_v(:, i), 'LineWidth', 2);
    hold on
    lh.Color = [lh.Color 1 - (i - 1) / n];
end
xlim([0 100]);
ylim([0 5]);
ylabel({'- FORCE -','Vertical TF contact force','[N/body mass * g]'}, 'FontWeight', 'bold');

subplot(2,3,2)
for i = 1:1:n
    lh = plot(cf_ap(:, i), 'LineWidth', 2);
    hold on
    lh.Color = [lh.Color 1 - (i - 1) / n];
end
xlim([0 100]);
ylim([-1 1]);
ylabel({'Ant-Post TF contact force' ,'[N/body mass * g]'}, 'FontWeight', 'bold');

subplot(2,3,3)
for i = 1:1:n
    lh = plot(cf_ml(:, i), 'LineWidth', 2);
    hold on
    lh.Color = [lh.Color 1 - (i - 1) / n];
    hleg(i) = lh;
end
xlim([0 100]);
ylim([-1 1]);
ylabel({'Med-Lat TF contact force', '[N/body mass * g]'}', 'FontWeight', 'bold');

% Last tweaks
lgd = legend(hleg, legendNames);
lgd.Box = 'off';
lgd.Location = 'northeast';
lgd.FontSize = 5;

% %% Part 2 plotting based on pressure
% 
% % Define body side to plot
% side2plot = strcat('_', side);
% 
% % Find all results files in working directory
% files = dir(fullfile(strcat(workingDirectory, 'JAM\**\', '*-JAM-Results-All-Trials.mat')));
% 
% plottingDat = load(char(fullfile(files.folder, files.name)));
% filesTmp = fieldnames(plottingDat.subResults.perTrial);
% 
% k = 1;
% for i = 1 : length(prefixes)
%     
%     % Get files matching the prefix & the side
%     idx_files = contains(filesTmp, prefixes(i));
%     idx_side = contains(filesTmp, side2plot);
%     idx = and(idx_files, idx_side); % combine both logical vectors
% 
%     % Load data - in some cases the array is only 100 frames long - for
%     % visualization purpose I interpolated to 101;
%     cf_tot(:,k) = plottingDat.subResults.perTrial.(filesTmp{find(idx)}).PressureForceVtp.vForceTot_BW;
%     cf_med(:,k) = plottingDat.subResults.perTrial.(filesTmp{find(idx)}).PressureForceVtp.vForceMed_BW;
%     cf_lat(:,k) = plottingDat.subResults.perTrial.(filesTmp{find(idx)}).PressureForceVtp.vForceLat_BW;
%     legendNames{k} = prefixes{i};
%     k = k+1;
%     
% end

%% PLot data
% hfig2 = figure;
% set(hfig2,'units','centimeters','position',[0,0,28,8]);
% sgtitle(myTitle,'FontWeight','bold');

% PLot data
subplot(2,3,4)
n = size(cf_tot,2);
for i = 1:1:n
    lh = plot(cf_tot(:, i), 'LineWidth', 2);
    hold on
    lh.Color = [lh.Color 1 - (i - 1) / n];
end
xlim([0 100]);
ylim([0 5]);
xlabel('% gait cycle', 'FontWeight', 'bold', 'linewidth', 2);
ylabel({'- PRESSURE-based -', 'Total TF contact force','[N/body mass * g]'}, 'FontWeight', 'bold');

subplot(2,3,5)
for i = 1:1:n
    lh = plot(cf_med(:, i), 'LineWidth', 2);
    hold on
    lh.Color = [lh.Color 1 - (i - 1) / n];
end
xlim([0 100]);
ylim([0 5]);
xlabel('% gait cycle', 'FontWeight', 'bold');
ylabel({'Medial TF contact force' ,'[N/body mass * g]'}, 'FontWeight', 'bold');

subplot(2,3,6)
for i = 1:1:n
    lh = plot(cf_lat(:, i), 'LineWidth', 2);
    hold on
    lh.Color = [lh.Color 1 - (i - 1) / n];
    hleg(i) = lh;
end
xlim([0 100]);
ylim([0 5]);
xlabel('% gait cycle', 'FontWeight', 'bold');
ylabel({'Lateral TF contact force', '[N/body mass * g]'}', 'FontWeight', 'bold');

% Last tweaks
lgd = legend(hleg, legendNames);
lgd.Box = 'off';
lgd.Location = 'northeast';
lgd.FontSize = 5;

h = hfig;
end

