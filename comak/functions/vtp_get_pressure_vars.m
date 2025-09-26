function [T] = vtp_get_pressure_vars(folder, geometryPath, bodyWeight, side, idxMedPtsLeft, idxLatPtsLeft, idxMedPtsRight, idxLatPtsRight)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get *.vtp pressure variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes a set of pressure variables from the *.vtp files.
%
%
% Example: 
% folder = 'C:\Users\bhorsak\Desktop\NewFiles2Test\JAM\Dynamic11_l\jam\';
% bodyWeight = 90 * 9.81;
% geometryPath = 'C:\Users\bhorsak\Desktop\plotting\Geometry';
% side = 'l'; 
%
% T = vtp_get_pressure_vars(folder, geometryPath, bodyWeight, side);
%
% Input: 
%         folder =          path populating the vtp files
%         geometryPath =    path to the geometry folder of the COMAK 
%                           pipeline
%         bodyWeight =      boday mass * g
%         side =            'l', or 'r' to indicate if it is a left or 
%                           right trial
%
% Output: table with all vars.
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
%                       Modified by Bernhard Guggenberger - bernhard.guggenberger2@fh-joanneum.at
% Acknowledgements:     Many thx to Bryce Killian & Sam Van Rossum 
%                       (KU Leuven) and Bernhard Dumphart (UAS St. Pölten)
%                       for their great support.
%
% Last changed:         08/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare input vars
side = lower(side); % make lower case in case input is upper case

% Define if right or left side will be plotted (hardcoded file names for geometries)
% Note: to get the correct values I need to use the scaled geometries!
if strcmp(side, 'r')
    % Right side
    stl_fem_cart = 'lenhart2015-R-femur-cartilage_scaled.stl'; %_scaled
    %stl_tib_cart = 'lenhart2015-R-tibia-cartilage_scaled.stl'; %_scaled
elseif strcmp(side, 'l')
    % Left side
    stl_fem_cart = 'lenhart2015-R-femur-cartilage_mirror_scaled.stl'; %_scaled
    %stl_tib_cart = 'lenhart2015-R-tibia-cartilage_mirror_scaled.stl'; %_scaled
end

% Get stl file path
stl_fem_cart = stlread(fullfile(geometryPath, stl_fem_cart));
%stl_tib_cart = stlread(fullfile(geometryPath, stl_tib_cart));


% Define the global file list
file_list = struct2table(dir(folder));

% Create vtp file list and only contact files not the mesh files
filesVtpFemContact = file_list.name(contains(file_list.name, 'contact_femur'));

% Frame length
nFrames = length(filesVtpFemContact);

%% Initialize vars
% Pressure
peakP = nan(nFrames,1);
peakPmed = nan(nFrames,1);
peakPlat = nan(nFrames,1);
avgP = nan(nFrames,1);
avgPmed = nan(nFrames,1);
avgPlat = nan(nFrames,1);

% Total contact force from pressure
vertF_tot = nan(nFrames,1);
vertF_med = nan(nFrames,1);
vertF_lat = nan(nFrames,1);

% Contact area
areaCart = nan(nFrames,1);
contactAreaTot = nan(nFrames,1);
contactAreaMed = nan(nFrames,1);
contactAreaLat = nan(nFrames,1);

% Frames
frames2Save = nan(nFrames,1);

for i = 1 : nFrames
    
    %% Get pressure values per *.vtp file and calculate med/lat components
    frame = i-1; % *.vtp files start with 0
    
%     % Define vtp file for the femur contact pressure @ current frame 
%        vtp_File = xml2struct(char(fullfile(folder,'\jam-ascii\',filesVtpFemContact(find(contains(filesVtpFemContact,strcat('_', num2str(frame), '.vtp')))))));
    
    % Read vtk file with toolbox:
    currentFile = fullfile(folder,filesVtpFemContact(find(contains(filesVtpFemContact,strcat('_', num2str(frame), '.vtp')))));
    outTmp = vtkRead(char(currentFile));

%     % Get vertices *.vtp
%     vertices_points = split(vtp_File.VTKFile.PolyData.Piece.Points.DataArray.Text);
%     vertices_points = str2num(char(vertices_points));
%     vertices = [vertices_points(1:3:end), vertices_points(2:3:end), vertices_points(3:3:end)];
    
    % Get vertices with vtk toolbox:
    vertices = double(outTmp.points);
    
%     % Get faces *.vtp
%     face_points = split(vtp_File.VTKFile.PolyData.Piece.Polys.DataArray{1, 1}.Text);
%     face_points = str2num(char(face_points));
%     faces = [face_points(1:3:end), face_points(2:3:end), face_points(3:3:end)] + 1; % Values come from Python 0 = 1 ...
    
    % Get faces with vtk toolbox:
    faces = double(outTmp.cells);
    
%     % Get pressure *.vtp
%     % Convert pascal to mega pascal and normalize by BW of patient (max range to 0.2)
%     % 'target_triangle_pressure_tf_conact' == index 3
%      pressure = split(vtp_File.VTKFile.PolyData.Piece.CellData.DataArray{1, 3}.Text );
%      pressure = str2num(char(pressure)); % do not transform to MPa because force will we wrong then!
%      idxIsPressureTot = pressure > 0;

    % Get tf pressure with vtk toolbox
    pressure = double(outTmp.cellData.target_triangle_pressure_tf_contact);
    idxIsPressureTot = pressure > 0;
    
    %% Pressure medial/lateral
    % Get all faces for fem cart *.stl file
    x_faces = stl_fem_cart.ConnectivityList(:,1);
    y_faces = stl_fem_cart.ConnectivityList(:,2);
    z_faces = stl_fem_cart.ConnectivityList(:,3);
    
    % Only get faces which hold pressure values but for med & lat!
    facesPress = [];
    facesPress(:,1) = x_faces(idxIsPressureTot);
    facesPress(:,2) = y_faces(idxIsPressureTot);
    facesPress(:,3) = z_faces(idxIsPressureTot);
    
    % Only pressure
    onlyPressure = pressure(idxIsPressureTot);
    
    % Get index for medial
    if strcmp(side, 'l'); idxMedPts = stl_fem_cart.Points(:,3) > 0; end % see were faces belong to med. side
    if strcmp(side, 'r'); idxMedPts = stl_fem_cart.Points(:,3) < 0; end % see were faces belong to med. side
    idxMed = (1:1:length(idxMedPts))'; % create 1:x index vector in length of idxMedPts
    idxMed = idxMed(idxMedPts); %only used those idxs which belong to the med. side
    
    % Get medial faces & pressure
    medFaces = [];
    medPressure = [];
    idx = 1;
    for k = 1 : size(facesPress,1)
        if ismember(facesPress(k, :), idxMed)
            medFaces(idx, :) = facesPress(k,:);
            medPressure(idx,:) = onlyPressure(k);
            idx = idx + 1;
        end
    end
    
    % Get index for lateral
    %idxLatPts = stl_fem_cart.Points(:,3) < 0;
    if strcmp(side, 'l'); idxLatPts = stl_fem_cart.Points(:,3) < 0; end % see were faces belong to lat. side
    if strcmp(side, 'r'); idxLatPts = stl_fem_cart.Points(:,3) > 0; end % see were faces belong to lat. side
    idxLat = (1:1:length(idxLatPts))';
    idxLat = idxLat(idxLatPts);
    
    % Get lateral faces & pressure
    latFaces = [];
    latPressure = [];
    idx = 1;
    for k = 1 : size(facesPress,1)
        if ismember(facesPress(k, :), idxLat)
            latFaces(idx, :) = facesPress(k,:);
            latPressure(idx,:) = onlyPressure(k);
            idx = idx + 1;
        end
    end
    
%     p = patch('Vertices',stl_fem_cart.Points.*idxMed3, 'Faces',faces, 'FaceVertexCData', pressure,'FaceColor','flat', 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
%     zlabel('z')
%     axis equal
    
%     % Get triangle areas *.vtp
%     % 'target_triangle_area' == index 10
%       areaTri = split(vtp_File.VTKFile.PolyData.Piece.CellData.DataArray{1, 10}.Text );
%       areaTri = str2num(char(areaTri));  

    % Compute the total area: the vtp files have a bug, the varnames have '.' in
    % it, so that vtk toolbox can't read i.e. the triangle size ....
    [~, areaTri] = heronsFormula(stl_fem_cart.Points, faces); %checked, similar to vtp files

    % Compute the areas of the triangles using heron's formula where
    % pressure is active
    if isempty(medFaces)
        areaTriMed = 0;
    else
        [~,areaTriMed] = heronsFormula(stl_fem_cart.Points, medFaces);
    end
    
    if isempty(latFaces)
        areaTriLat = 0;
    else
        [~,areaTriLat] = heronsFormula(stl_fem_cart.Points, latFaces);
    end
    
    % Get triangle normal vectors % This also looks good.
    normalsTot = compute_triangle_normal(vertices, faces);
    normalsMed = compute_triangle_normal(vertices, medFaces);
    normalsLat = compute_triangle_normal(vertices, latFaces);
   
    % Compute the center of each faces triangle
    %     cenTot = [];
    %     cenTot = (vertices(faces(:,1),:)+vertices(faces(:,2),:)+vertices(faces(:,3),:))/3;
    
    %% Compute variables    
    areaCart(i) = sum(areaTri);
    
    % Total
    if ~isempty(pressure(idxIsPressureTot))
    peakP(i) = max(pressure(idxIsPressureTot));
    avgP(i) = mean(pressure(idxIsPressureTot));
    vertF_tot(i) = sum(pressure.*areaTri.*normalsTot(:,2));
    contactAreaTot(i) = sum(areaTri(idxIsPressureTot));
    end
    
    % Medial
    if ~isempty(medPressure)
    peakPmed(i) = max(medPressure);
    avgPmed(i) = mean(medPressure);
    vertF_med(i) = sum(medPressure.*areaTriMed.*normalsMed(:,2));
    contactAreaMed(i) = sum(areaTriMed);
    end
    
    % Lateral
    if ~isempty(latPressure)
    peakPlat(i) = max(latPressure);
    avgPlat(i) = mean(latPressure);
    vertF_lat(i) = sum(latPressure.*areaTriLat.*normalsLat(:,2));
    contactAreaLat(i) = sum(areaTriLat);
    end
    
    % Frames
    frames2Save(i) = frame;
    
end
%% Normalize values to BW

% Pressure
peakP_BW = peakP/bodyWeight;
peakPmed_BW = peakPmed/bodyWeight;
peakPlat_BW = peakPlat/bodyWeight;

avgP_BW = avgP/bodyWeight;
avgPmed_BW = avgPmed/bodyWeight;
avgPlat_BW = avgPlat/bodyWeight;

% Contact force
vertF_tot_BW = vertF_tot/bodyWeight;
vertF_med_BW = vertF_med/bodyWeight;
vertF_lat_BW = vertF_lat/bodyWeight;

%% Put results in table for output
T = table();
T.('FrameIndex') = frames2Save;

T.('TFpeakPressureTot') = peakP;
T.('TFpeakPressureMed') = peakPmed;
T.('TFpeakPressureLat') = peakPlat;
T.('TFpeakPressureTot_BW') = peakP_BW;
T.('TFpeakPressureMed_BW') = peakPmed_BW;
T.('TFpeakPressureLat_BW') = peakPlat_BW;


T.('TFavgPressureTot') = avgP;
T.('TFavgPressureMed') = avgPmed;
T.('TFavgPressureLat') = avgPlat;
T.('TFavgPressureTot_BW') = avgP_BW;
T.('TFavgPressureMed_BW') = avgPmed_BW;
T.('TFavgPressureLat_BW') = avgPlat_BW;

T.('TFvForceTot') = vertF_tot;
T.('TFvForceMed') = vertF_med;
T.('TFvForceLat') = vertF_lat;
T.('TFvForceTot_BW') = vertF_tot_BW;
T.('TFvForceMed_BW') = vertF_med_BW;
T.('TFvForceLat_BW') = vertF_lat_BW;

T.('TFareaCartilage') = areaCart;
T.('TFContactAreaTot') = contactAreaTot;
T.('TFContactAreaMed') = contactAreaMed;
T.('TFContactAreaLat') = contactAreaLat;

%% Calculate Pressures for the Patellofemoral joint

% Added by Bernhard Guggenberger - bernhard.guggenberger2@fh-joanneum.at
% Currently the medial and lateral components are calculated the same way
% as for the tibiofemoral pressures.

%% Prepare input vars
side = lower(side); % make lower case in case input is upper case

% Define if right or left side will be plotted (hardcoded file names for geometries)
% Note: to get the correct values I need to use the scaled geometries!
if strcmp(side, 'r')
    % Right side
    stl_pat_cart = 'lenhart2015-R-patella-cartilage_scaled.stl'; %_scaled
    %stl_tib_cart = 'lenhart2015-R-tibia-cartilage_scaled.stl'; %_scaled
elseif strcmp(side, 'l')
    % Left side
    stl_pat_cart = 'lenhart2015-R-patella-cartilage_mirror_scaled.stl'; %_scaled
    %stl_tib_cart = 'lenhart2015-R-tibia-cartilage_mirror_scaled.stl'; %_scaled
end

% Get stl file path
stl_pat_cart = stlread(fullfile(geometryPath, stl_pat_cart));

% Define the global file list
file_list = struct2table(dir(folder));

% Create vtp file list and only contact files not the mesh files
filesVtpPatContact = file_list.name(contains(file_list.name, 'contact_patella'));

% Frame length
nFrames = length(filesVtpPatContact);

%% Initialize vars
% Pressure
peakP_Pat = nan(nFrames,1);
peakPmed_Pat = nan(nFrames,1);
peakPlat_Pat = nan(nFrames,1);
avgP = nan(nFrames,1);
avgPmed = nan(nFrames,1);
avgPlat = nan(nFrames,1);

% Total contact force from pressure
contactF_tot = nan(nFrames,1);
contactF_med = nan(nFrames,1);
contactF_lat = nan(nFrames,1);

% Contact area
areaCart = nan(nFrames,1);
contactAreaTot = nan(nFrames,1);
contactAreaMed = nan(nFrames,1);
contactAreaLat = nan(nFrames,1);

% Frames
frames2Save = nan(nFrames,1);

%Calculate indices for medial and lateral side
if strcmp(side, 'l')
    idxMed = (1:1:length(idxMedPtsLeft))'; % create 1:x index vector in length of idxMedPts
    idxMed = idxMed(idxMedPtsLeft); %only used those idxs which belong to the med. side
    idxLat = (1:1:length(idxLatPtsLeft))';
    idxLat = idxLat(idxLatPtsLeft);
else
    idxMed = (1:1:length(idxMedPtsRight))'; % create 1:x index vector in length of idxMedPts
    idxMed = idxMed(idxMedPtsRight); %only used those idxs which belong to the med. side
    idxLat = (1:1:length(idxLatPtsRight))';
    idxLat = idxLat(idxLatPtsRight);
end

for i = 1 : nFrames
    
    %% Get pressure values per *.vtp file and calculate med/lat components
    frame = i-1; % *.vtp files start with 0
    
    % Read vtk file with toolbox:
    currentFile = fullfile(folder,filesVtpPatContact(find(contains(filesVtpPatContact,strcat('_', num2str(frame), '.vtp')))));
    outTmp = vtkRead(char(currentFile));

    % Get vertices with vtk toolbox:
    vertices = double(outTmp.points);
    
    % Get faces with vtk toolbox:
    faces = double(outTmp.cells);

    % Get tf pressure with vtk toolbox
    pressure = double(outTmp.cellData.casting_triangle_pressure_pf_contact);
    idxIsPressureTot = pressure > 0;
    
    %% Pressure medial/lateral
    % Get all faces for pat cart *.stl file
    x_faces = stl_pat_cart.ConnectivityList(:,1);
    y_faces = stl_pat_cart.ConnectivityList(:,2);
    z_faces = stl_pat_cart.ConnectivityList(:,3);
    
    % Only get faces which hold pressure values but for med & lat!
    facesPress = [];
    facesPress(:,1) = x_faces(idxIsPressureTot);
    facesPress(:,2) = y_faces(idxIsPressureTot);
    facesPress(:,3) = z_faces(idxIsPressureTot);
    
    % Only pressure
    onlyPressure = pressure(idxIsPressureTot);


    % Get medial faces & pressure
    medFaces = [];
    medPressure = [];
    idx = 1;
    for k = 1 : size(facesPress,1)
        if ismember(facesPress(k, :), idxMed)
            medFaces(idx, :) = facesPress(k,:);
            medPressure(idx,:) = onlyPressure(k);
            idx = idx + 1;
        end
    end

    
    % Get lateral faces & pressure
    latFaces = [];
    latPressure = [];
    idx = 1;
    for k = 1 : size(facesPress,1)
        if ismember(facesPress(k, :), idxLat)
            latFaces(idx, :) = facesPress(k,:);
            latPressure(idx,:) = onlyPressure(k);
            idx = idx + 1;
        end
    end
 

    % Compute the total area: the vtp files have a bug, the varnames have '.' in
    % it, so that vtk toolbox can't read i.e. the triangle size ....
    [~, areaTri] = heronsFormula(stl_pat_cart.Points, faces); %checked, similar to vtp files

    % Compute the areas of the triangles using heron's formula where
    % pressure is active
    if isempty(medFaces)
        areaTriMed = 0;
    else
        [~,areaTriMed] = heronsFormula(stl_pat_cart.Points, medFaces);
    end
    
    if isempty(latFaces)
        areaTriLat = 0;
    else
        [~,areaTriLat] = heronsFormula(stl_pat_cart.Points, latFaces);
    end
    
    % Get triangle normal vectors % This also looks good.
    normalsTot = compute_triangle_normal(vertices, faces);
    normalsMed = compute_triangle_normal(vertices, medFaces);
    normalsLat = compute_triangle_normal(vertices, latFaces);
   
    
    %% Compute variables    
    areaCart(i) = sum(areaTri);
    
    % Total
    if ~isempty(pressure(idxIsPressureTot))
    peakP_Pat(i) = max(pressure(idxIsPressureTot));
    avgP(i) = mean(pressure(idxIsPressureTot));
    contactF_tot(i) = sum(pressure.*areaTri.*normalsTot(:,1));
    contactAreaTot(i) = sum(areaTri(idxIsPressureTot));
    end
    
    % Medial
    if ~isempty(medPressure)
    peakPmed_Pat(i) = max(medPressure);
    avgPmed(i) = mean(medPressure);
    contactF_med(i) = sum(medPressure.*areaTriMed.*normalsMed(:,1));
    contactAreaMed(i) = sum(areaTriMed);
    end
    
    % Lateral
    if ~isempty(latPressure)
    peakPlat_Pat(i) = max(latPressure);
    avgPlat(i) = mean(latPressure);
    contactF_lat(i) = sum(latPressure.*areaTriLat.*normalsLat(:,1));
    contactAreaLat(i) = sum(areaTriLat);
    end
    
    % Frames
    frames2Save(i) = frame;
    
end
%% Normalize values to BW

% Pressure
peakP_BW = peakP_Pat/bodyWeight;
peakPmed_BW = peakPmed_Pat/bodyWeight;
peakPlat_BW = peakPlat_Pat/bodyWeight;

avgP_BW = avgP/bodyWeight;
avgPmed_BW = avgPmed/bodyWeight;
avgPlat_BW = avgPlat/bodyWeight;

% Contact force
contactF_tot_BW = contactF_tot/bodyWeight;
contactF_med_BW = contactF_med/bodyWeight;
contactF_lat_BW = contactF_lat/bodyWeight;

%% Put results in table for output
T.('PFpeakPressureTot') = peakP_Pat;
T.('PFpeakPressureMed') = peakPmed_Pat;
T.('PFpeakPressureLat') = peakPlat_Pat;
T.('PFpeakPressureTot_BW') = peakP_BW;
T.('PFpeakPressureMed_BW') = peakPmed_BW;
T.('PFpeakPressureLat_BW') = peakPlat_BW;


T.('PFavgPressureTot') = avgP;
T.('PFavgPressureMed') = avgPmed;
T.('PFavgPressureLat') = avgPlat;
T.('PFavgPressureTot_BW') = avgP_BW;
T.('PFavgPressureMed_BW') = avgPmed_BW;
T.('PFavgPressureLat_BW') = avgPlat_BW;

T.('PFcontactForceTot') = contactF_tot;
T.('PFcontactForceMed') = contactF_med;
T.('PFcontactForceLat') = contactF_lat;
T.('PFcontactForceTot_BW') = contactF_tot_BW;
T.('PFcontactForceMed_BW') = contactF_med_BW;
T.('PFcontactForceLat_BW') = contactF_lat_BW;

T.('PFareaCartilage') = areaCart;
T.('PFContactAreaTot') = contactAreaTot;
T.('PFContactAreaMed') = contactAreaMed;
T.('PFContactAreaLat') = contactAreaLat;

%% Finally fill 'nan' gaps here
% In case of single Nans interpolate them
MaxGapSize = 10; %hardcoded size which I tolerate
[TF, A] = fillmissing(T, 'linear', 'MaxGap', MaxGapSize);

if max(sum(ismissing(T))) < MaxGapSize && sum(sum(ismissing(T))) > 0
    disp(['>>>>> A total of ', num2str(sum(sum(A))) , 'NaNs filled with linear interpolation. Maximum number of NaNs per variable:', num2str(max(sum(A)))]);
end

% Rise warning when MaxGapSize is reached for one variable
if max(sum(ismissing(T))) >= MaxGapSize 
    warning([num2str(max(sum(ismissing(T)))), ' NaN values in one variable! NaNs were NOT interpolated but will results in a NaN-only variable because of the interpft function! Check that! See <vtp_get_pressure_vars>'])
end

% TF is new output
T = TF;

end