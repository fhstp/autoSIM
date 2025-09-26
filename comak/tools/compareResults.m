% Compare two force reporter files, e.g. walking example with my
% workflow output

path2example = 'C:\LokaleDaten\CloudStation\MyPapers\2020_GuidedGrowth\code\opensim-jam-master\opensim-jam-release\examples\walking\results\joint-mechanics\walking_ForceReporter_forces.sto';
path2mydata = 'C:\LokaleDaten\CloudStation\MyPapers\2020_GuidedGrowth\data\AK-data\Data2Verify_WalkingExample\JAM\autoExtLoad_overground_17_r\jam\comak_autoExtLoad_overground_17_r__ForceReporter_forces.sto';
line_width = 1.5;
side = 'r';

% Read walking example
[example_data, example_labels, example_header] = read_opensim_mot(path2example);

% Read my data
[my_data, my_labels, my_header] = read_opensim_mot(path2mydata);

% Select muscles to display
msls = {['soleus_',side],['gasmed_',side],['recfem_',side], ...
    ['glmax1_',side], ['bfsh_',side], ['bflh_',side], ...
    ['perlong_',side], ['vasmed_',side], ['tibant_',side], ...
    ['semimem_',side],['glmed1_',side],['glmin1_',side]};

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

%% Set figure 1
hfigKinem = figure;
set(hfigKinem,'units','centimeters','position',[0,0,29,8]);
orient(hfigKinem,'landscape');
sgtitle('Comparison Walking Example vs. My Data','FontWeight','bold');

%% Plot Knee Kinematics
time = example_data(:,1);

for i = 1:length(tf_contactF)
    % Example data
    ind = find(contains(example_labels,tf_contactF{i}));
    forces_time = example_data(:,1);
    example_data_tmp = example_data(:,ind);
    
    % My data
    ind = find(contains(my_labels,tf_contactF{i}));
    my_data_tmp = my_data(:,ind);
    
    subplot(1,3,i);
    plot(forces_time,example_data_tmp,'k-','LineWidth',line_width);
    hold on
    plot(forces_time,my_data_tmp,'r--','LineWidth',line_width);
    title(strrep(tf_contactF{i},'_',' '));
    xlim([time(1),time(end)]);
    xlabel('Time [s]');
    if i == 1; 
        ylabel('Force [N]');
    legend('Walking example','My data','location','southwest', 'box', 'off');
    end
    
    
end

% Plot tf kinematics
% for i = 1:length(tf_coords)
%     % Get tf translation
%     ind = find(contains(example_labels,tf_coords{i}{1}));
%     if(contains(tf_coords{i}{4},'Translation'))
%         ex_data = example_data(:,ind)*1000;
%         my_data = my_data(:,ind)*1000;
%     else
%         ex_data = example_data(:,ind);
%         my_data = my_data(:,ind);
%     end
%     
%     % Plot the data
%     subplot(5,3,i);
%     plot(time,ex_data,'LineWidth',line_width, 'color', [43,140,190]/255);
%     hold on
%     plot(time,my_data,'LineWidth',line_width, '--r');
%     
%     title([tf_coords{i}{3} ' (' tf_coords{i}{2} ')'])
%     xlim([time(1),time(end)]);
%     set(gca,'xticklabel',[]);
%     
%     if i == 1 || i == 4
%         ylabel(tf_coords{i}{4})
%     end
% end
% 
% % Plot pf kinematics
% for i = 1:length(pf_coords)
%     % Get pf translation
%     ind = find(contains(example_labels,pf_coords{i}{1}));
%     if(contains(pf_coords{i}{4},'Translation'))
%         ex_data = example_data(:,ind)*1000;
%         my_data = my_data(:,ind)*1000;
%     else
%         ex_data = example_data(:,ind);
%         my_data = mydata(:,ind);
%     end
%     
%     % Plot the data
%     subplot(5,3,i+6);
%     plot(time,ex_data,'LineWidth',line_width,'color', [166,189,219]/255)
%     title([pf_coords{i}{3} ' (' pf_coords{i}{2} ')']);
%     xlim([time(1),time(end)])
%     set(gca,'xticklabel',[]);
%     if i == 1 || i == 4
%         ylabel(pf_coords{i}{4})
%     end
% end