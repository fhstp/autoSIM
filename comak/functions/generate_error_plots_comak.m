function [error, hfig] = generate_error_plots_comak(path2IKerrorReport,path2ActivationErrorReport,path2ForceErrorReport,path2IDErrorReport,statesFileName, cTO, TO)
% Initialize error struct to collect data
error = struct();

%% Inverse kinematics
% Add to error struct
IK = struct();
[IK.values_data, IK.values_labels, ~] = read_opensim_mot(path2IKerrorReport);

% Find all relevant variables
fn = {'time', 'total_squared_error','marker_error_RMS', 'marker_error_max'};
idx = find(contains(IK.values_labels,fn));
IKerror_labels = IK.values_labels(idx);
IKerror_labels = strrep(IKerror_labels,'_',' ');
IKerror_values = IK.values_data(:,idx);

% Add to error struct
error.IKerror = array2table(IKerror_values, 'VariableNames',IKerror_labels);

%% COMAK__activation.sto
STO = struct();
[STO.values_data, STO.values_labels, ~] = read_opensim_mot(path2ActivationErrorReport);

% Find all reserve activation variables
idx_r = find(contains(STO.values_labels,'reserve'));
idx_t = find(contains(STO.values_labels,'time'));
idx = [idx_t; idx_r];
reserve_actuators_lables = STO.values_labels(idx);
reserve_actuators_values = STO.values_data(:,idx);

% Clean labels
for j = 2: length(reserve_actuators_lables) % exclude time
    idx = strfind(reserve_actuators_lables{j},'/');
    reserve_actuators_lables{j} = reserve_actuators_lables{j}(idx(end)+1:end);
end

% Add to error struct
error.ReserveActivation = array2table(reserve_actuators_values, 'VariableNames',reserve_actuators_lables);

%% COMAK__force.sto
FORCE = struct();
[FORCE.values_data, FORCE.values_labels, ~] = read_opensim_mot(path2ForceErrorReport);

% Find all reserve activation variables
idx_r = find(contains(FORCE.values_labels,'reserve'));
idx_t = find(contains(FORCE.values_labels,'time'));
idx = [idx_t; idx_r];
force_reserve_lables = FORCE.values_labels(idx);
force_reserve_values = FORCE.values_data(:,idx);

% Clean labels
for j = 2: length(force_reserve_lables) % exclude time
    idx = strfind(force_reserve_lables{j},'/');
    force_reserve_lables{j} = force_reserve_lables{j}(idx(end)+1:end);
end

% Add to error struct
error.ForceReserve = array2table(force_reserve_values, 'VariableNames',force_reserve_lables);

%% Inverse Dynamics - Residuals forces on pelvis
ID = struct();
[ID.values_data, ID.values_labels, ~] = read_opensim_mot(path2IDErrorReport);

% Find all pelvis variables
fn = {'time','pelvis_tx_force','pelvis_ty_force','pelvis_tz_force'};
idx = find(contains(ID.values_labels,fn));
pelvis_residuals_lables = ID.values_labels(idx);
pelvis_residuals_values = ID.values_data(:,idx);
pelvis_residual_labels = strrep(fn,'_',' ');

% Add to error struct
error.ResidualForcesPelvis = array2table(pelvis_residuals_values, 'VariableNames',pelvis_residuals_lables);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot everything to error-report %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hfig = figure;
set(hfig,'units','centimeters','position',[0,0,21,29]);
%set(fig,'units','normalized', 'position', [0 0 0.3 0.8])
sgtitle(char(strcat('Error-Report for:',{' '}, strrep(statesFileName,'_',' '))),'FontWeight','bold');

%% Subplot 1 - Inverse Kinematics
markersize = 10;
jit = 0.2;
boxwidth = 1;
col = [224,243,219; 168,221,181; 67,162,202]/255;
fn = IKerror_labels;
dat = (table2array(error.IKerror));

for j = 1 : 3
    k = j+1; % exclude time
    subplot(4,3,j)
    x = dat(:,k); % exclude time
    x = x * 100; % express in cm
    raincloudErr(x,col(j,:),1,boxwidth,markersize,jit,hfig);
    hold on
    xticks(1);
    title(strrep(fn{k},'_',' ')); % exclude time
    if j == 1; ylabel('IK Error (Centimeter)','FontWeight','bold'); end
    ax = gca;                   % gca = get current axis
    ax.YAxis.Visible = 'on';   % remove y-axis
    ax.XAxis.Visible = 'off';   % remove x-axis
    ylim([0 ceil(max(x))*1.05]);
    xlim([0 2]);
    
    % Plot threshold
    cond = fn{k};
    switch cond
        case 'marker error RMS'
            plot([0,2],[2,2],'r--');
            ylim([0 ceil(max(x))*1.05]);
            h_tmp = findobj(gca,'Type','line');
            hleg = legend([h_tmp(1) h_tmp(6)],'Limit', 'Median', 'Location','southeast');
            legend boxoff
            set(hleg,'FontSize',5);
        case 'marker error max'
            plot([0,2],[4,4],'r--');
            ylim([0 ceil(max(x))*1.05]);
            
            % Plot name of marker wiht max. error
            
    end
end

%% Subplot 2 - IK against time
subplot(4,3,[4 5 6])
dat = (table2array(error.IKerror));
dat = dat * 100; % express in centimeter

% Plot wiht different colors
for j = 1:3
    jj = j+1;
    plot(dat(:,jj:jj),'color',col(j,:), 'LineWidth', 1.5);
    hold on;
end
hold off;
xlabel('Frames/Time steps','Fontweight','bold')
ylabel('IK Errors (Centimeter)','Fontweight','bold')
legend(strrep(IK.values_labels(2:end),'_',' '),'box', 'on','location', 'northeastoutside');
xlim([0 height(error.IKerror)]);

%% Subplot 3 - Reserve Actuators
subplot(4,3,[7 8 9])
dat = (table2array(error.ReserveActivation));
plot(dat(:,1),dat(:,2:end),'LineWidth', 1.5);
%xlabel('Frames/Time steps','Fontweight','bold');

% Add IC, ICi, TOc, TO
hold on;
y_min = min(min(dat));
y_max = max(max(dat));
line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '--', 'LineWidth', 0.5);

ylabel('Reserve Actuators','Fontweight','bold');
ylim([y_min y_max]);
xlim([(dat(1,1)) dat(end,1)]);


%% Subplot 4 - Residual Forces Pelvis
subplot(4,3,[10 11 12])
dat = table2array(error.ResidualForcesPelvis);
plot(dat(:,1),dat(:,2:4), 'LineWidth', 1.5);

% Add IC, ICi, TOc, TO
hold on;
y_min = -250;%min(min(dat));
y_max = max(max(dat));
line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);

ylabel('Residual Forces Pelvis [N]','Fontweight','bold');
xlabel('Frames/Time steps','Fontweight','bold');
legend([pelvis_residual_labels(2:end),'cont. lat. Toe-off','Toe-off']);
ylim([y_min y_max]);
xlim([(dat(1,1)) dat(end,1)]);
yline(50,'r--','+/- 50N Threshold','LineWidth',1.5);
yline(-50,'r--','LineWidth',1.5);
end