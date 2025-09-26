function [hfig1,hfig2] = vtp_time_plotter(folder, fileName, geometryPath, timeStamp, bodyWeight, normBW, side, x_Label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting *.vtp files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used to plot the tibio-femoral contact pressure. For this
% purpose it reads the vtp-ascii files, gets the contact pressure, and
% plots the pressure onto the the *.stl files from the geometry folder.
% Pressure is normalized to body weight. It uses a specific colormap which
% stored as a separate file. this should be be stored on the same level as
% this function (cmap-vtpPlotter.mat).
%
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Acknowledgements:     Many thx to Bernhard Dumphart who developed parts 
%                       of the code.
%
% Example: h = vtp_time_plotter(folder, fileName, geometryPath, ...
%                               timeStamp, body_weight, side, x_Label);
% Input: 
%         folder =          path to the specific jam-ascii folder for one 
%                           subject
%         fileName =        hardcoded file name, e.g. 'Dynamic11_l'
%         geometryPath =    path to the geometry folder of the COMAK 
%                           pipeline. Use unscaled geometries, e.g. from
%                           the setupFiles folder!
%         timeStamp =       vector containing the time stamps to be 
%                           displayed, e.g. [10, 20, 30, 40, 50, 60]
%         bodyWeight =      boday mass * g
%         normBW =          'true', 'false', define if pressure should be
%                           bodyweight normalized or not
%         side =            'l', or 'r' to indicate if it is a left or 
%                           right trial
%         x_label =         label of the x-axis, e.g. % stance or % gait 
%                           cycle
%
% Last changed:         04/2021
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myFont = 'Arial';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare input vars
side = lower(side); % make lower case in case input is upper case
normBW = lower(normBW); % make lower case in case input is upper case
if strcmp(normBW, 'false'); bodyWeight = 1; end % set norm factor to 1 if no normalization is wanted 
fileName = strrep(fileName, '_', ' '); % replace underscore

% Define if right or left side will be plotted (hardcoded file names for geometries)
% I need to use the unscaled version here because the read stl function has
% a problem with the scaled file of the tibia and fibula!
if strcmp(side, 'r')
    % Right side
    stl_fem_cart_file = 'lenhart2015-R-femur-cartilage.stl';
    stl_tib_cart_file = 'lenhart2015-R-tibia-cartilage.stl';
    stl_pat_cart_file = 'lenhart2015-R-patella-cartilage.stl';
    stl_fem_bone_file = 'lenhart2015-R-femur-bone.stl';
    stl_tib_bone_file = 'lenhart2015-R-tibia-bone.stl';
    stl_fib_bone_file = 'lenhart2015-R-fibula-bone.stl';
    stl_pat_bone_file = 'lenhart2015-R-patella-bone.stl';
    
elseif strcmp(side, 'l')
    % Left side
    stl_fem_cart_file = 'lenhart2015-R-femur-cartilage_mirror.stl';
    stl_tib_cart_file = 'lenhart2015-R-tibia-cartilage_mirror.stl';
    stl_pat_cart_file = 'lenhart2015-R-patella-cartilage_mirror.stl';
    stl_fem_bone_file = 'lenhart2015-R-femur-bone_mirror.stl';
    stl_tib_bone_file = 'lenhart2015-R-tibia-bone_mirror.stl';
    stl_fib_bone_file = 'lenhart2015-R-fibula-bone_mirror.stl';
    stl_pat_bone_file = 'lenhart2015-R-patella-bone_mirror.stl';
end

% Set number of subplots
nTime = length(timeStamp);

% Define the global file list
file_list = struct2table(dir(folder));

% Get stl file paths
stl_fem_cart = stlread(fullfile(geometryPath, stl_fem_cart_file));
stl_tib_cart = stlread(fullfile(geometryPath, stl_tib_cart_file));
stl_pat_cart = stlread(fullfile(geometryPath, stl_pat_cart_file));
stl_fem_bone = stlread(fullfile(geometryPath, stl_fem_bone_file));
stl_tib_bone = stlread(fullfile(geometryPath, stl_tib_bone_file)); % I need to use the unscaled version here because the read stl function has a problem with the scaled file.
stl_fib_bone = stlread(fullfile(geometryPath, stl_fib_bone_file)); % I need to use the unscaled version here because the read stl function has a problem with the scaled file.
stl_pat_bone = stlread(fullfile(geometryPath, stl_pat_bone_file));

% Create vtp file list (fem, tib, pat) and only contact files not the mesh files
filesVtpFemContact = file_list.name(contains(file_list.name, 'contact_femur'));
filesVtpTibContact = file_list.name(contains(file_list.name, 'contact_tibia'));
filesVtpPatContact = file_list.name(contains(file_list.name, 'contact_patella'));

% Check if files were found, otherwise rais error.
if isempty(filesVtpFemContact)
    warning('No *.vtp files found - please check settings. Did you set <Yes> to process *.vtp files even though none were processed during simulation?');
    pause();
end

%% Prepare plotting
hfig1 = figure;
set(hfig1,'units','centimeters','position',[0,0,20,5]);
hold on

% Set colormap
load(fullfile(fileparts(mfilename('fullpath')), 'cmap-vtpPlotter.mat'));
colormap(cmap);

%% Femur
for i = 1:nTime
    % Define vtp file for the femur contact pressure
%     vtp_File = xml2struct(char(fullfile(folder,filesVtpFemContact(find(contains(filesVtpFemContact,strcat('_', num2str(timeStamp(i)), '.vtp')))))));

    % Read vtk file with toolbox:
    currentFile = fullfile(folder,filesVtpFemContact(find(contains(filesVtpFemContact,strcat('_', num2str(timeStamp(i)), '.vtp')))));
    outTmp = vtkRead(char(currentFile));
    
    % Convert pascal to mega pascal and normalize by BW of patient (max range to 0.2)
    % 'target_triangle_pressure' == index 4 from 1-9
    % 'target_tirangle_pressure_tf_conact' == index 3
%     pressure = split(vtp_File.VTKFile.PolyData.Piece.CellData.DataArray{1, 3}.Text );
%     pressure = str2num(char(pressure));
%     pressure = pressure / 1000000 / bodyWeight; % Pascal to Mega Pascal (1*10^6)

    % Get tf pressure with vtk toolbox
    pressure = outTmp.cellData.target_triangle_pressure_tf_contact;
    pressure = pressure / 1000000 / bodyWeight; % Pascal to Mega Pascal (1*10^6);
    
    % Plot femur cartilage
    subaxis(2, nTime, i,'SpacingVert',0.0,'SpacingHoriz',0.0,'MR',0.20, 'ML', 0.20);
    p = patch('Vertices',stl_fem_cart.Points, 'Faces',stl_fem_cart.ConnectivityList, 'FaceVertexCData',pressure,'FaceColor','flat', 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
    axis equal off
    p.LineStyle = 'none';
    p.AmbientStrength = 0.9;
    p.DiffuseStrength = 0.8;
    p.SpecularStrength = 0.2;
    p.SpecularExponent = 25;
    hold on
    
    % Plot femur bone
    pbone = patch('Vertices',stl_fem_bone.Points, 'Faces',stl_fem_bone.ConnectivityList);
    pbone.FaceColor = [236 226 198]/255;
    pbone.LineStyle = 'none';
    pbone.AmbientStrength = 0.9;
    pbone.DiffuseStrength = 0.8;
    pbone.SpecularStrength = 0.2;
    pbone.SpecularExponent = 25;
    
    % Now change lighting etc.
    lightangle(0,0)
    lighting gouraud
    material dull
    
    % Set view
    if strcmp(side, 'l'); set(gca, 'ZDir','reverse'); end % mirror plot to have the lateral side on the same side as the tibia
    view(0,0)
    camroll(90)
    xlim([-0.05 0.05])
    %ylim([-0.105 0.105])
    zlim([-0.05 0.05])
    
end

%% Tibia & Fibula
cnt = 1;
for i = nTime+1:nTime*2
    % Tibia
    % Define vtp file for tibia contact pressure
%     vtp_File = xml2struct(char(fullfile(folder,filesVtpTibContact(find(contains(filesVtpTibContact,strcat('_', num2str(timeStamp(cnt)), '.vtp')))))));

    % Read vtk file with toolbox:
    currentFile = fullfile(folder,filesVtpTibContact(find(contains(filesVtpTibContact,strcat('_', num2str(timeStamp(cnt)), '.vtp')))));
    outTmp = vtkRead(char(currentFile));
    
    % Convert pascal to mega pascal and normalize by BW of patient (max range to 0.2)
    % 'target_triangle_pressure' == index 4 from 1-9
    % 'target_tirangle_pressure_tf_conact' == index 3
%     pressure = split(vtp_File.VTKFile.PolyData.Piece.CellData.DataArray{1, 3}.Text );
%     pressure = str2num(char(pressure));
%     pressure = pressure / 1000000 / bodyWeight; % Pascal to Mega Pascal (1*10^6)

    % Get tf pressure with vtk toolbox          
    pressure = outTmp.cellData.casting_triangle_pressure_tf_contact;
    pressure = pressure / 1000000 / bodyWeight; % Pascal to Mega Pascal (1*10^6);
    
    % Plot tibia cartilage
    subaxis(2, nTime, i,'SpacingVert',0.0,'SpacingHoriz',0.0,'MR',0.20, 'ML', 0.20);
    %p = patch('Vertices',vertices, 'Faces',faces, 'FaceVertexCData',pressure,'FaceColor','flat', 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
    p = patch('Vertices',stl_tib_cart.Points, 'Faces',stl_tib_cart.ConnectivityList, 'FaceVertexCData',pressure,'FaceColor','flat', 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
    axis equal off
    p.LineStyle = 'none';
    p.AmbientStrength = 0.9;
    p.DiffuseStrength = 0.8;
    p.SpecularStrength = 0.2;
    p.SpecularExponent = 25;
    hold on
    
    % Plot tibia bone
    pbone = patch('Vertices',stl_tib_bone.Points, 'Faces',stl_tib_bone.ConnectivityList);
    pbone.FaceColor = [236 226 198]/255;
    pbone.LineStyle = 'none';
    pbone.AmbientStrength = 0.9;
    pbone.DiffuseStrength = 0.8;
    pbone.SpecularStrength = 0.2;
    pbone.SpecularExponent = 25;
    
    % Fibula bone   
    pbone = patch('Vertices',stl_fib_bone.Points, 'Faces',stl_fib_bone.ConnectivityList);
    pbone.FaceColor = [236 226 198]/255;
    pbone.LineStyle = 'none';
    pbone.AmbientStrength = 0.9;
    pbone.DiffuseStrength = 0.8;
    pbone.SpecularStrength = 0.2;
    pbone.SpecularExponent = 25;
    
    % Now change lighting etc.
    lightangle(0,180)
    lighting gouraud
    material dull
    
    % Set view
    if strcmp(side, 'r'); set(gca, 'ZDir','reverse'); end % mirror plot to have the lateral side on the same side as the femur
    view(180,0)
    camroll(-90)
    xlim([-0.05 0.05])
    %ylim([-0.05 0.05])
    zlim([-0.05 0.05])    

    % Increase cnt
    cnt = cnt + 1;
end

%% Last plotting tweaks
% x labels, I rotated the camera, therefore I have a problem with visualizing
% the x axis. This is a dirty hack to resolve that.
xPos = 0.2:((0.7-0.2)/(nTime-1)):0.7;
for j = 1: nTime
annotation(hfig1,'textbox',...
    [xPos(j) 0.10 0.0995291005291004 0.0793650793650794],...
    'String',strcat(num2str(timeStamp(j)),'%'),...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');
end

% Anterior
annotation(hfig1,'textbox',...
    [0.452058201058201 0.88 0.0995291005291004 0.0793650793650793],...
    'String','Ant.',...
    'HorizontalAlignment','center',...
    'FontWeight','normal',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% Posterior
annotation(hfig1,'textbox',...
    [0.452058201058201 0.17 0.0995291005291004 0.0793650793650793],...
    'String','Post.',...
    'HorizontalAlignment','center',...
    'FontWeight','normal',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% Lateral
annotation(hfig1,'textbox',...
    [0.2 0.510052910052909 0.0995291005291004 0.0793650793650793],...
    'String','Lat.',...
    'HorizontalAlignment','left',...
    'FontWeight','normal',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% Medial
annotation(hfig1,'textbox',...
    [0.7 0.510052910052909 0.0995291005291004 0.0793650793650793],...
    'String','Med.',...
    'HorizontalAlignment','right',...
    'FontWeight','normal',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% xLabel
annotation(hfig1,'textbox',...
    [0.452058201058201 0.03 0.0995291005291004 0.0793650793650794],...
    'String', x_Label,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% FileName Label
annotation(hfig1,'textbox',...
    [0.20 0.93 0.404026455026455 0.0634920634920634],...
    'String', fileName,...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% Add colorbar
cb = colorbar(); 
colormap(cmap);
cb.Position = [0.82 0.17 0.02 0.74];
cb.FontWeight = 'bold';
cb.FontName = myFont;
cb.FontSize = 9;
if strcmp(normBW, 'true')
    cb.Label.String = {'Contact pressure', '[MPa/Body mass * g]'};
elseif strcmp(normBW, 'false')
    cb.Label.String = {'Contact pressure [MPa]'};
end

%%

%% Prepare plotting second figure for patella und femur contact
hfig2 = figure;
set(hfig2,'units','centimeters','position',[0,0,20,5]);
hold on

% Set colormap
load(fullfile(fileparts(mfilename('fullpath')), 'cmap-vtpPlotter.mat'));
colormap(cmap);

%% Femur
for i = 1:nTime
    % Define vtp file for the femur contact pressure
%     vtp_File = xml2struct(char(fullfile(folder,filesVtpFemContact(find(contains(filesVtpFemContact,strcat('_', num2str(timeStamp(i)), '.vtp')))))));

    % Read vtk file with toolbox:
    currentFile = fullfile(folder,filesVtpFemContact(find(contains(filesVtpFemContact,strcat('_', num2str(timeStamp(i)), '.vtp')))));
    outTmp = vtkRead(char(currentFile));
    
    % Convert pascal to mega pascal and normalize by BW of patient (max range to 0.2)
    % 'target_triangle_pressure' == index 4 from 1-9
    % 'target_tirangle_pressure_pf_conact' == index 3
%     pressure = split(vtp_File.VTKFile.PolyData.Piece.CellData.DataArray{1, 12}.Text );
%     pressure = str2num(char(pressure));
%     pressure = pressure / 1000000 / bodyWeight; % Pascal to Mega Pascal (1*10^6)

    % Get tf pressure with vtk toolbox
    pressure = outTmp.cellData.target_triangle_pressure_pf_contact;
    pressure = pressure / 1000000 / bodyWeight; % Pascal to Mega Pascal (1*10^6);
    
    % Plot femur cartilage
    subaxis(2, nTime, i,'SpacingVert',0.0,'SpacingHoriz',0.0,'MR',0.20, 'ML', 0.20);
    p = patch('Vertices',stl_fem_cart.Points, 'Faces',stl_fem_cart.ConnectivityList, 'FaceVertexCData',pressure,'FaceColor','flat', 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
    axis equal off
    p.LineStyle = 'none';
    p.AmbientStrength = 0.9;
    p.DiffuseStrength = 0.8;
    p.SpecularStrength = 0.2;
    p.SpecularExponent = 25;
    hold on
    
    % Plot femur bone
    pbone = patch('Vertices',stl_fem_bone.Points, 'Faces',stl_fem_bone.ConnectivityList);
    pbone.FaceColor = [236 226 198]/255;
    pbone.LineStyle = 'none';
    pbone.AmbientStrength = 0.9;
    pbone.DiffuseStrength = 0.8;
    pbone.SpecularStrength = 0.2;
    pbone.SpecularExponent = 25;
    
    % Now change lighting etc.
    lightangle(90,0)
    lighting gouraud
    material dull
    
    % Set view
    if strcmp(side, 'l'); set(gca, 'ZDir','reverse'); end % mirror plot to have the lateral side on the same side as the tibia
    view(90,0)
    camroll(90)
    xlim([-0.05 0.05])
    %ylim([-0.105 0.105])
    zlim([-0.05 0.05])
    
end

%% Patella
cnt = 1;
for i = nTime+1:nTime*2
    
    % Define vtp file for patella contact pressure
%     vtp_File = xml2struct(char(fullfile(folder,filesVtpPatContact(find(contains(filesVtpPatContact,strcat('_', num2str(timeStamp(cnt)), '.vtp')))))));

    % Read vtk file with toolbox:
    currentFile = fullfile(folder,filesVtpPatContact(find(contains(filesVtpPatContact,strcat('_', num2str(timeStamp(cnt)), '.vtp')))));
    outTmp = vtkRead(char(currentFile));
    
    % Convert pascal to mega pascal and normalize by BW of patient (max range to 0.2)
%     pressure = split(vtp_File.VTKFile.PolyData.Piece.CellData.DataArray{1, 3}.Text );
%     pressure = str2num(char(pressure));
%     pressure = pressure / 1000000 / bodyWeight; % Pascal to Mega Pascal (1*10^6)

    % Get pf pressure with vtk toolbox
    pressure = outTmp.cellData.casting_triangle_pressure_pf_contact;
    pressure = pressure / 1000000 / bodyWeight; % Pascal to Mega Pascal (1*10^6);
    
    % Plot patella cartilage
    subaxis(2, nTime, i,'SpacingVert',0.0,'SpacingHoriz',0.0,'MR',0.20, 'ML', 0.20);
    %p = patch('Vertices',vertices, 'Faces',faces, 'FaceVertexCData',pressure,'FaceColor','flat', 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
    p = patch('Vertices',stl_pat_cart.Points, 'Faces',stl_pat_cart.ConnectivityList, 'FaceVertexCData',pressure,'FaceColor','flat', 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
    axis equal off
    p.LineStyle = 'none';
    p.AmbientStrength = 0.9;
    p.DiffuseStrength = 0.8;
    p.SpecularStrength = 0.2;
    p.SpecularExponent = 25;
    hold on
    
    % Plot patella bone
    pbone = patch('Vertices',stl_pat_bone.Points, 'Faces',stl_pat_bone.ConnectivityList);
    pbone.FaceColor = [236 226 198]/255;
    pbone.LineStyle = 'none';
    pbone.AmbientStrength = 0.9;
    pbone.DiffuseStrength = 0.8;
    pbone.SpecularStrength = 0.2;
    pbone.SpecularExponent = 25;
    
    % Now change lighting etc.
    lightangle(90,180)
    lighting gouraud
    material dull
    
    % Set view
    if strcmp(side, 'r'); set(gca, 'ZDir','reverse'); end % mirror plot to have the lateral side on the same side as the femur
    view(90,180)
    camroll(90)
    xlim([-0.05 0.05])
    %ylim([-0.05 0.05])
    zlim([-0.05 0.05])    

    % Increase cnt
    cnt = cnt + 1;
end

%% Last plotting tweaks
% x labels, I rotated the camera, therefore I have a problem with visualizing
% the x axis. This is a dirty hack to resolve that.
xPos = 0.2:((0.7-0.2)/(nTime-1)):0.7;
for j = 1: nTime
annotation(hfig2,'textbox',...
    [xPos(j) 0.10 0.0995291005291004 0.0793650793650794],...
    'String',strcat(num2str(timeStamp(j)),'%'),...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');
end

% Anterior
annotation(hfig2,'textbox',...
    [0.452058201058201 0.88 0.0995291005291004 0.0793650793650793],...
    'String','Ant.',...
    'HorizontalAlignment','center',...
    'FontWeight','normal',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% Posterior
annotation(hfig2,'textbox',...
    [0.452058201058201 0.17 0.0995291005291004 0.0793650793650793],...
    'String','Post.',...
    'HorizontalAlignment','center',...
    'FontWeight','normal',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% Lateral
annotation(hfig2,'textbox',...
    [0.2 0.510052910052909 0.0995291005291004 0.0793650793650793],...
    'String','Lat.',...
    'HorizontalAlignment','left',...
    'FontWeight','normal',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% Medial
annotation(hfig2,'textbox',...
    [0.7 0.510052910052909 0.0995291005291004 0.0793650793650793],...
    'String','Med.',...
    'HorizontalAlignment','right',...
    'FontWeight','normal',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% xLabel
annotation(hfig2,'textbox',...
    [0.452058201058201 0.03 0.0995291005291004 0.0793650793650794],...
    'String', x_Label,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% FileName Label
annotation(hfig2,'textbox',...
    [0.20 0.93 0.404026455026455 0.0634920634920634],...
    'String', fileName,...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'FontName', myFont,...
    'FitBoxToText','off', ...
    'LineStyle','none');

% Add colorbar
cb = colorbar(); 
colormap(cmap);
cb.Position = [0.82 0.17 0.02 0.74];
cb.FontWeight = 'bold';
cb.FontName = myFont;
cb.FontSize = 9;
if strcmp(normBW, 'true')
    cb.Label.String = {'Contact pressure', '[MPa/Body mass * g]'};
elseif strcmp(normBW, 'false')
    cb.Label.String = {'Contact pressure [MPa]'};
end




end
