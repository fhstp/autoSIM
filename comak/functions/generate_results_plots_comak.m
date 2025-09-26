function [results, hfigKinem, hfigMA, hfigCPM_TF, hfigCPM_PF] = generate_results_plots_comak(trialType, workingDirectory, path2KneeKinem, path2ContactF, path2MuscleAct, path2setupFiles, path2IDReport, path2MuscleForceReport, path2normData, path2origGeometry, statesFileName, side, cTO, TO, IC, cIC, BW, processVtp)
%% Plot Walking Simulation Results
%==========================================================================

%% Some house keeping variables
line_width = 2;
side = lower(side(1));
path.KneeKinematics = path2KneeKinem;
path.ContactForces = path2ContactF;
path.MuscleActivation = path2MuscleAct;
path.MuscleForces = path2MuscleForceReport;
path.NormEMG = char(fullfile(strcat(path2normData,'\Lencioni_2019\Lencioni2019_EMGnorm.mat')));
path.NormKCF = char(fullfile(strcat(path2normData,'\OrthoLoad\OrthoLoad_KCFnorm.mat')));
path.ID = path2IDReport;


% Select muscles to display
msls = {['soleus_',side],['gasmed_',side],['recfem_',side], ...
    ['glmax1_',side], ['bfsh_',side], ['bflh_',side], ...
    ['perlong_',side], ['vasmed_',side], ['tibant_',side], ...
    ['semimem_',side],['glmed1_',side],['glmin1_',side]};

% Norm muscles
msls_norm = {'Soleus','Gastrocnemius_Medialis','Rectus_Femoris', ...
    'Gluteus_Maximus', 'Biceps_Femoris', 'Biceps_Femoris',  ...
    'Peroneus_Longus','Vastus_Medialis','Tibialis_Anterior', ...
    'nan'   , 'nan'   , 'nan'};

% Select tibio-femoral kinematics (6)
tf_coords = {...
    {['knee_flex_',side];['knee\_flex\_',side];'Flex.';'Angle[°]'};...
    {['knee_add_',side];['knee\_add\_',side];'Add.';'Angle[°]'};...
    {['knee_rot_',side];['knee\_rot\_',side];'Int. Rotation';'Angle[°]'};...
    {['knee_tx_',side];['knee\_tx\_',side];'Ant. Translation';'Translation [mm]'};...
    {['knee_ty_',side];['knee\_ty\_',side];'Super. Translation';'Translation [mm]'};...
    {['knee_tz_',side];['knee\_tz\_',side];'Lat. Translation';'Translation [mm]'}};

% Select patello-femoral kinematics (6)
pf_coords = {...
    {['pf_flex_',side];['pf\_flex\_',side];'Flex.';'Angle[°]'};...
    {['pf_rot_',side];['pf\_rot\_',side];'Rot.';'Angle[°]'};...
    {['pf_tilt_',side];['pf\_tilt\_',side];'Tilt';'Angle[°]'};...
    {['pf_tx_',side];['pf\_tx\_',side];'Ant. Translation';'Translation [mm]'};...
    {['pf_ty_',side];['pf\_ty\_',side];'Super. Translation';'Translation [mm]'};...
    {['pf_tz_',side];['pf\_tz\_',side];'Lat. Translation';'Translation [mm]'}};

% Select contact forces (3)
tf_contactF = {'tf_contact.casting.total.contact_force_x','tf_contact.casting.total.contact_force_y','tf_contact.casting.total.contact_force_z'};


%% Initialiez struct for results
results = struct();

%% Set figure 1
hfigKinem = figure;
set(hfigKinem,'units','centimeters','position',[0,0,29,21]);
orient(hfigKinem,'landscape');
sgtitle(strcat('Results Report for:',{' '}, strrep(statesFileName,'_',' ')),'FontWeight','bold');

%% Get Inverse Dynamics path.ID
[values_data, values_labels, ~] = read_opensim_mot(path.ID); %('../results/inverse-dynamics/Dynamic11_l_inverse-dynamics');

% Store data 2 table
results.InverseDynamicscs = array2table(values_data, 'VariableNames',values_labels);

%% Get & Plot Kinematics
[values_data, values_labels, ~] = read_opensim_mot(path.KneeKinematics); %('../results/comak/walking_values.sto');
kinem_time = values_data(:,1);

% Store data 2 table
results.Kinematics = array2table(values_data, 'VariableNames',values_labels);

% Plot tf kinematics
for i = 1:length(tf_coords)
    % Get tf translation
    ind = find(contains(values_labels,tf_coords{i}{1}));
    if(contains(tf_coords{i}{4},'Translation'))
        data = values_data(:,ind)*1000;
    else
        data = values_data(:,ind);
    end
    
    % Get TO and IC index
    [~, ICidx] = min(abs(kinem_time-IC));
    [~, cTOidx] = min(abs(kinem_time-cTO));
    
    % cut data
    data = data(ICidx:end);
    
    % Plot the data
    subplot(5,3,i);
    plot(kinem_time(ICidx:end),data,'LineWidth',line_width, 'color', [43,140,190]/255)
    title([tf_coords{i}{3} ' (' tf_coords{i}{2} ')'])
    xlim([kinem_time(ICidx),kinem_time(end)]);
    
    % Add cTO, TO
    hold on;
    y_min = (min(min(data)));
    y_max = (max(max(data)));
    line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    line([cIC cIC],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    ylim([y_min y_max]);
    
    if i == 1 || i == 4
        ylabel(tf_coords{i}{4})
    end
end

% Plot pf kinematics
for i = 1:length(pf_coords)
    % Get pf translation
    ind = find(contains(values_labels,pf_coords{i}{1}));
    if(contains(pf_coords{i}{4},'Translation'))
        data = values_data(:,ind)*1000;
    else
        data = values_data(:,ind);
    end
    
     % Get TO and IC index
    [~, ICidx] = min(abs(kinem_time-IC));
    [~, cTOidx] = min(abs(kinem_time-cTO));
    
    % cut data
    data = data(ICidx:end);
    
    % Plot the data
    subplot(5,3,i+6);
    plot(kinem_time(ICidx:end),data,'LineWidth',line_width,'color', [166,189,219]/255)
    title([pf_coords{i}{3} ' (' pf_coords{i}{2} ')']);
    xlim([kinem_time(ICidx),kinem_time(end)])
    
    % Add cTO, TO
    hold on;
    y_min = (min(min(data)));
    y_max = (max(max(data)));
    line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    line([cIC cIC],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    ylim([y_min y_max]);
    
    if i == 1 || i == 4
        ylabel(pf_coords{i}{4})
    end
end

%% Get & Plot Contact Forces
[forces_data, forces_labels, forces_header] = read_opensim_mot(path.ContactForces);%'../results/joint-mechanics/walking_ForceReporter_forces.sto');

% Store data 2 table BUT delete muscles as they do not show correct values.
% I am only intersted in everything except the muscles!
msl2exclude = {'addbrev', 'addlong', 'addmagProx', 'addmagMid', 'addmagDist', 'addmagIsch', ...
               'bflh', 'bfsh', 'edl', 'ehl', 'fdl', 'fhl', 'gaslat', 'gasmed', 'gem', 'glmax1', ...
               'glmax2', 'glmax3', 'glmed1', 'glmed2', 'glmed3', 'glmin1', 'glmin2', 'glmin3', 'grac', ...
               'iliacus', 'pect', 'perbrev', 'perlong', 'pertert', 'piri', 'psoas', 'quadfem', 'recfem', ...
               'sart', 'semimem', 'semiten', 'soleus', 'tfl', 'tibant', 'tibpost', 'vasint', 'vaslat', 'vasmed'};

% add side to muscles cell array
msl2exclude = cellfun(@(c)[c '_' side],msl2exclude,'uni',false);

% Write cleaned data to results
dataTmp = array2table(forces_data, 'VariableNames',forces_labels);
varnames = dataTmp.Properties.VariableNames;
col2exclude = not(ismember(varnames, msl2exclude));
cleanedDataTmp = dataTmp(:,col2exclude);
results.ContactForces = cleanedDataTmp;

% Load NormData
load(path.NormKCF);
KFCfn = {'Fy','F','Fx'};

for i = 1:length(tf_contactF)
    ind = find(contains(forces_labels,tf_contactF{i}));
    
    forces_time = forces_data(:,1);
    tmp_data = forces_data(:,ind)/BW; % normalize to BW
    
    % Get TO and IC index
    [~, TOidx] = min(abs(forces_time-TO));
    %[~, cTOidx] = min(abs(forces_time-cTO));
    %[~, ICidx] = min(abs(forces_time-IC));
    
    subplot(5,3,i+12);
    hDat = plot(forces_time,tmp_data,'LineWidth',line_width, 'Color', [221,28,119]/255);
    title(strrep(tf_contactF{i},'_',' '));
    xlim([forces_time(1),forces_time(end)]);
    xlabel('Time [s]');
    if i == 1; ylabel('Body Weight'); end    
    hold on
    
    % Plot OrthoLoad norm data
    % Adjust sign for norm data
    if contains(tf_contactF{i},'force_x')
        curveKFC_mean_stance = interpft(mean(KFCnorm.NormBW.(KFCfn{i}),2),TOidx)*-1;
    else
        curveKFC_mean_stance = interpft(mean(KFCnorm.NormBW.(KFCfn{i}),2),TOidx);
    end
    
    % Prepare data
    curveKFC_mean_stance = padarray(curveKFC_mean_stance,length(forces_time)-TOidx,'post')/100; % padd zeros to end, change to multipel of BW
    KFC_sd_stance = interpft(std(KFCnorm.NormBW.(KFCfn{i}),0,2),TOidx)/100; % change to multipel of BW
    KFC_sd_stance = padarray(KFC_sd_stance,length(forces_time)-TOidx,'post'); % padd zeros to end
    
    % Plot Norm Values only if walking trial
    if strcmp(trialType, 'walking')     
    % Plot walking norm
    hNorm = plot(forces_time,curveKFC_mean_stance,'k','LineWidth',line_width);
    plot(forces_time,curveKFC_mean_stance+2*KFC_sd_stance,'k--','LineWidth',0.5);
    plot(forces_time,curveKFC_mean_stance-2*KFC_sd_stance,'k--','LineWidth',0.5);
    end
    
    % Add cTO, TO
    hold on;
    y_min = (min([min(curveKFC_mean_stance-2*KFC_sd_stance); min(tmp_data)]));
    y_max = (max([max(curveKFC_mean_stance+2*KFC_sd_stance); max(tmp_data)]));
    hcTO = line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    line([cIC cIC],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    hTO = line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    ylim([y_min y_max]);
    xlim([forces_time(1) forces_time(end)]);
    hold off
    
    if i == 3
        h_tmp = findobj(gca,'Type','line');
        if strcmp(trialType, 'walking')
            hleg = legend([hcTO hTO hNorm hDat],'cont. lat. Toe-off/IC','Toe-off','OrthoLoad (2SD) [%BW])', 'TF KCF', 'Location','southwest');
        else
            hleg = legend(hDat, 'TF KCF', 'Location','southwest');
        end
        legend boxoff
        set(hleg,'FontSize',6);
    end
    
end

%% Get & Plot Muscle Activation
% Get Muscle Forces
[values_data, values_labels, ~] = read_opensim_mot(path.MuscleForces); %('../results/inverse-dynamics/Dynamic11_l_inverse-dynamics');

% Store data 2 table
results.MuscleForces = array2table(values_data, 'VariableNames',values_labels);

% Set figure 2
hfigMA = figure;
set(hfigMA,'units','centimeters','position',[0,0,29,21]);
orient(hfigMA,'landscape');
sgtitle(char(strcat('Results Report (Muscle Activations) for:',{' '}, strrep(statesFileName,'_',' '))),'FontWeight','bold');

% Plot Muscle Activation
[act_data, act_labels, ~] = read_opensim_mot(path.MuscleActivation);%'../results/comak/walking_activation.sto');
act_time = act_data(:,1);
msls_names = strrep(msls,'_', ' ');

% Store data 2 table
results.MuscleActivation = array2table(act_data, 'VariableNames',act_labels);

% Load NormData
load(path.NormEMG);

for i = 1:length(msls)
    
    % Get IC index
    [~, ICidx] = min(abs(act_time-IC));
    act_time_cropped = act_time(ICidx:end);
    
    subplot(4,3,i);
    % Find index for act_labels
    ind = find(contains(act_labels,msls{i}));
    
    hold on;
    % Plot data
    %ind = find(contains(act_labels,msls{i}));
    hDat = plot(act_time_cropped,act_data(ICidx:end,ind),'LineWidth',line_width, 'Color', [233,163,201]/255);
    box on;
    ylim([0 1]);
    xlim([act_time_cropped(1),act_time_cropped(end)]);
    ylabel(msls_names{i});

    % Plot Norm Values only if walking trial
    if strcmp(trialType, 'walking')

        if ~strcmp(msls_norm{i},'nan')
            % Amplitude normalization
            for m = 1:size(EMGnorm.(msls_norm{i}),2)
                tmp_signal = EMGnorm.(msls_norm{i})(:,m);
                emg_dat(:,m) = (tmp_signal - (min(tmp_signal))) / (((max(tmp_signal))) - (min(tmp_signal))); % Min-Max normalization
            end
            norm_mean = mean(emg_dat,2);
            sf = max(act_data(:,ind))/max(norm_mean); % calculate scaling factor to have similar peak value between norm and muscle activations
            norm_std = std(emg_dat,0,2);
            norm_std = interpft(norm_std,length(act_time_cropped)); % apply scaling factor and normalize to data length (IC-ICi)

            norm_curve = interpft(norm_mean*sf,length(act_time_cropped)); % apply scaling factor
            hNorm = plot(act_time_cropped,norm_curve,'k','LineWidth',1.5);
            plot(act_time_cropped,norm_curve+1*norm_std,'k--','LineWidth',0.5);
            plot(act_time_cropped,norm_curve-1*norm_std,'k--','LineWidth',0.5);

            % Set plotting thresholds for walking
            y_min = round(min([min(norm_curve-2*norm_std); min(act_data(ICidx:end,ind))]));
            y_max = round(max([max(norm_curve+2*norm_std); max(act_data(ICidx:end,ind))]));
        end
    else
        % Set plotting thresholds for non-walking trials
        y_min = 0;
        y_max = 1;

    end
     
    
    % Add cTO, TO
    hold on;
    line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    hcTO = line([cIC cIC],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    hTO = line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    ylim([0 1]);
    %ylim([y_min y_max]);
    xlim([act_time_cropped(1) act_time_cropped(end)]);
    hold off
    
    if i > 9
        xlim([act_time_cropped(1),act_time_cropped(end)]);
        xlabel('Time [s]');
    end

    if i == 9
        h_tmp = findobj(gca,'Type','line');
        if strcmp(trialType, 'walking')
            hleg = legend([hNorm hDat hcTO hTO],'Lencioni 2019 (1SD), scaled to peak MA', 'COMAK Muscle Activation', 'cont. lat. Toe-off/IC','Toe-off', 'Location','northeast');
        else
            hleg = legend(hDat, 'COMAK Muscle Activation', 'Location','northeast');
        end
        legend boxoff
        set(hleg,'FontSize',6);
    end
end

%% Plot contact pressure map
if processVtp
    % Get jam folder
    partsPath = regexp(fileparts(path.KneeKinematics),'\','split');
    folder = fullfile(strjoin(partsPath(1:end-2),'\'), statesFileName(1:end-1),'jam'); %, 'jam-ascii');
    geometryPath = fullfile(path2origGeometry); % use the none-scaled files instead of scaled files %fullfile(workingDirectory,'Geometry');
    [hfigCPM_TF, hfigCPM_PF] = vtp_time_plotter(folder, statesFileName, geometryPath, [0, 15, 30, 45, 60, 80, 99], BW, 'false', side, 'gait cycle');
else
    hfigCPM_TF = NaN;
    hfigCPM_PF = NaN;
end
end