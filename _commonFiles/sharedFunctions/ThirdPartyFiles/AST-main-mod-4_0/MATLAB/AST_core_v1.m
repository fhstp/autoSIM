%% AST core
import org.opensim.modeling.*;

while true % iterate until the stop condition at the end of script
    if k>1 && flag==0 && flag2==0
        tolerance=0.64; % defaukte = 0.08; tolerance of scaling factor to determine the ScaleFactorRange
        if RMSErr(k-1)>=ManualScaleErr % if the RMSE is too highm, the tolerance is 0 and the tool will scale all segments using Mean Scale Factor
            tolerance=0;
        end
        if MeanScaleFactY > 1 % The ScaleFactorRange is not symmetrical with the Mean Scaling Factor
            ScaleFactorsRange=[MeanScaleFactY-tolerance/4, MeanScaleFactY+tolerance*3/4];%define the scaling factors range to consider the scaling reliable : it's translated toward higher values
        elseif MeanScaleFactY < 1
            ScaleFactorsRange=[MeanScaleFactY-tolerance*3/4, MeanScaleFactY+tolerance/4];%define the scaling factors range to consider the scaling reliable: it's translated toward lower values
        elseif MeanScaleFactY == 1
            ScaleFactorsRange=[MeanScaleFactY-tolerance/2, MeanScaleFactY+tolerance/2];%define the scaling factors range to consider the scaling reliable: it's translated toward lower values
        end
    end
    if k==1 % at  first iteration
        ScaleTool(fullfile(modelFolder,SetupFile)).run; % Run scale setup file at first iteration
    elseif k > 1 % for iterations > 1
        if RMSErr(k-1) > ManualScaleErr && flag==0 || flag2==1 % condition to scale with all segments = mean scaling factor.
            ind2=ind2+1;
            flag2=1;
            if ind2>=rep2 % flag2 returns to 0 just when the tool performs Rep2 times in a row the manual scaling for all the segments
                flag2=0;
                ind2=0;
            end
            ScaleTool(path_manualScale).run;
            disp('................................................................................................................');
            disp('Manual scaling for all bodies');
            disp('................................................................................................................');
        elseif ManualBodies>1 || flag==1 % condition to scale with just some segments = mean scaling factor
            ind=ind+1;
            flag=1;
            if ind>=rep % flag returns to 0 just when the tool performs Rep2 times in a row the manual scaling for the detected segments
                flag=0;
                ind=0;
            end
            ScaleTool(path_ScalerMix).run; % Run manual setup file
            for l=1:ManualBodies
                Str(l)=convertCharsToStrings(NamesBodyManual{l,1}.toString.toCharArray);
            end
            if exist('Str','var')
                disp('................................................................................................................');
                disp(append('Manual scale for body: ', strjoin(unique(Str)))); % displaying the segments been scaled manually
                disp('................................................................................................................');
            end
            clear NamesBodyManual
        else
            ScaleTool(fullfile(modelFolder, SetupFile)).run;% scale the model with the settings in the scale tool setup
        end
    end
    % IK
    ikTool=InverseKinematicsTool(path_ik_static); % call back IK tool for static trial
    ScaledModel= Model(fullfile(modelFolder,ScaledFileName)); % calling the Scaled model
    ikTool.setModel(ScaledModel); % set the Scaled model in the IK tool
    ikTool.setDebugLevel(0);
    ikTool.run; % run IK
    [MarkerLocation, HeadSTO]=load_sto(fullfile(modelFolder,'_ik_model_marker_locations.sto')); %loading the Model markers locations in output from IK
    %% obtaining the list of selected markers for the scaling (Markers True in IK Setup file)
    for i = 1:length(HeadSTO)-1 % do not consider time position in the HeadSTO vector
        SelectedMarker(i)=HeadSTO(i+1); % storing markers names
        SelectedMarker(i)=erase(SelectedMarker(i), '_tx'); %erasing extensions
        SelectedMarker(i)=erase(SelectedMarker(i), '_ty'); %erasing extensions
        SelectedMarker(i)=erase(SelectedMarker(i), '_tz'); %erasing extensions
    end
    SelectedMarkerList=unique(SelectedMarker,'stable'); %list of markers used
    %% determining the error between Experimental and Model markers using IK
    LengthMarkErr=min(length(MarkerLocation(:,1)),length(StatTRC(:,1))); % min length between TRC frames and Ik frames
    t1=time_range.get(0); % 
    t2=time_range.get(1);
    [~, TRC_start] = min((abs(StatTRC(:,2) - t1))); %searching in the time column % >>>> bhorsak (08/2024): made this line more generic.
    [~, TRC_end] = min((abs(StatTRC(:,2) - t2))); %searching in the time column   % >>>> bhorsak (08/2024): made this line more generic.
    for i=1:length(SelectedMarkerList) % for every marker
        p_TRC=contains(HeadTRC,SelectedMarkerList{i});%position of a certain marker data in .TRC file
        pos_TRC=(find(p_TRC==1)-2)*3; % first coordinate of .TRC marker
        pos_STO=(i-1)*3+2:(i*3)+1;% all 3 coordinates of .STO marker
        MarkErr(1:LengthMarkErr,3*(i-1)+1:3*i)=StatTRC(TRC_start:TRC_end,pos_TRC:pos_TRC+2)/trcScale-MarkerLocation(1:LengthMarkErr,pos_STO); % Error matrix between .TRC and ModelMarkerLocation  % >>>> bhorsak (08/2024): changed from "1000" to a "trcScale"  user can specify in the start up script..
        ModuleErr(i)=sqrt(sum(mean(MarkErr(:,3*(i-1)+1:3*i)).^2));% modulus of error for every marker
    end
    MeanMarkErr=mean(MarkErr); % mean of Markers errors for every direction
    [MaxErr,PosMaxErr]=max(abs(MeanMarkErr)); % finding max error
    RMSErr(k)=rms(ModuleErr);% finding RMS error
    %% determining what marker and in what direction to change
    ToChange=HeadSTO(PosMaxErr+1); % +1 beacuse the time is also considered  in HeadSTO (retrieving marker name from headSTO)
    if contains(ToChange,'_tx')
        dirToChange=1;
        MarkerToChange=erase(ToChange,'_tx');
    elseif contains( ToChange ,'_ty')
        dirToChange=2;
        MarkerToChange=erase(ToChange,'_ty');
    elseif contains( ToChange ,'_tz')
        dirToChange=3;
        MarkerToChange=erase(ToChange,'_tz');
    end
    direction(k)=dirToChange;
    %% storing errors and Marker name errors at every cycle
    err(k)=MaxErr;
    Merr(k)=MarkerToChange;
    %% control for scale factor: if scale factor is > or < of a certain threshold
    ManualBodies=0; % number of bodies scaled with Mean scaling factor
    if k~=1 % for every iteration > 1
        for n=0:Scaler.getModelScaler.getMeasurementSet.getSize-1 % for every measurement set
            for z=0:Scaler.getModelScaler.getMeasurementSet.get(n).getMarkerPairSet.getSize-1 %for every markerPair

                % Only get measurement set if set as true in scale setup file. % >>>> bhorsak (08/2024): added this line to make sure only measurements sets are used that were specified as true.
                if Scaler.getModelScaler.getMeasurementSet.get(n).get.Apply

                    FirstMarker=Scaler.getModelScaler.getMeasurementSet.get(n).getMarkerPairSet.get(z).getMarkerName(0); % first marker of the MarkerPair
                    SecondMarker=Scaler.getModelScaler.getMeasurementSet.get(n).getMarkerPairSet.get(z).getMarkerName(1);% second marker of the MarkerPair
                    % find position on both 2 markers in .TRC file
                    for p=1:length(HeadTRC)
                        if  strcmp(HeadTRC(p),FirstMarker)
                            PosA_exp=(p-2)*3;
                        elseif  strcmp(HeadTRC(p),SecondMarker)
                            PosB_exp=(p-2)*3;
                        end
                    end
                    %
                    state=model.initSystem();
                    % find position on both 2 markers in generic model file
                    LocA_mod = osimVec3ToArray(markerset.get(FirstMarker).getLocationInGround(state))'; % >>>> bhorsak (08/2024): to make compatible with API 4.0%.getAsMat; % get first model marker on the starting model
                    LocB_mod = osimVec3ToArray(markerset.get(SecondMarker).getLocationInGround(state))'; % >>>> bhorsak (08/2024): to make compatible with API 4.0 %.getAsMat;% get first model markers on the starting model
                    %
                    dist_mod(z+1)=pdist([LocA_mod';LocB_mod'],'euclidean'); %distance of marker pairs in  model file
                    dist_exp(z+1)=pdist([mean(StatTRC(:,PosA_exp:PosA_exp+2));mean(StatTRC(:,PosB_exp:PosB_exp+2))],'euclidean')/trcScale;%distance of marker pairs in .TRC file  % >>>> bhorsak (08/2024): changed from "1000" to a "trcScale"  user can specify in the start up script..

                    ScaleFact(n+1)=mean(dist_exp./dist_mod);%definition of scale factor
                    clear dist_mod dist_exp
                    if ScaleFact(n+1) < ScaleFactorsRange(1) || ScaleFact(n+1) > ScaleFactorsRange(2) % if scale factor is outside the range
                        if Scaler.getModelScaler.getMeasurementSet.get(n).getApply()==1 % considering only applied measurement sets
                            NamesBodyManual{ManualBodies+1,1}=Scaler.getModelScaler.getMeasurementSet.get(n).getBodyScaleSet.get(0).getName;%saving names of bodies affected
                            ManualBodies=ManualBodies+1;
                        end
                    end
                end
            end
        end
        clear n
    end
    %% creation of mixed ( measurement set and manual) scaling setup
    if ManualBodies > 0 % if detected
        ScalerMix=Scaler;
        % write inside array object
        ScalerMix.getModelScaler.setScalingOrder(Array2);% in the scaling order we want the array just created ("manualScale")
        for m=0:NumBodies-1
            scale.setSegmentName(model.getBodySet.get(m).getName); % set body name
            scale.setApply(0)
            ScalerMix.getModelScaler.getScaleSet.cloneAndAppend(scale);% appending the scale factor on the scale set
        end
        for m=0:NumBodies-1
            for l=1:ManualBodies
                scale=Scale(); % defining scale object
                scale.setScaleFactors(MeanScaleFact)% set the same scale factor for each body
                scale.setSegmentName(model.getBodySet.get(m).getName); % set body name
                if strcmp(model.getBodySet.get(m).getName, convertCharsToStrings(NamesBodyManual{l,1}.toString.toCharArray))
                    scale.setApply(1);% apply: true
                    ScalerMix.getModelScaler.getScaleSet.cloneAndAppend(scale);% appending the scale factor on the scale set
                end
            end
        end
        path_ScalerMix=fullfile(modelFolder,'ScalerMix.xml');
        ScalerMix.print(path_ScalerMix);
    end
    %% increment to change the marker coordinate affected by the max error
    if flag==1 % in case of manual scaling just for some bodies
        step(k)=MaxErr;
    elseif flag2==1 % in case of manual scaling for all bodies
        step(k) = MaxErr/2;
    else% in case of normal scaling
        step(k)= MaxErr/6;
    end
    state=model.initSystem();%take the state of the model
    GenericMarker=markerset.get(Merr(k));%considering the marker to chaange
    SockName=GenericMarker.getParentFrame;%retrieve Socket name (local)
    GroundFrame=model.getJointSet.get(0).getParentFrame;%define the ground coordinate system
    GenericMarker.changeFramePreserveLocation(state,GroundFrame); %change from local to ground coodinate systems
    currentCoordAbs = osimVec3ToArray(GenericMarker.get_location)';%.getAsMat;%from Vec3 to vector % >>>> bhorsak (08/2024): to make compatible with API 4.0
    %% control to sign and step size: if err rises, go back of 2 step
    if  k~=1
        if  err(k) >= err(k-1) && strcmp(Merr(k),Merr(k-1)) && direction(k)==direction(k-1) %&&  if the Max error arises at next cycle to the same marker and to the same direction it means that the right sign is the opposite one
            s=-s;
            disp('sign changed');
            step(k)=2*step(k-1);% doubling this step because the step at iteration k-1 was wrong
        end
    end
    %% setting new coordinate
    switch dirToChange % add increment according to the direction to change
        case 1
            NewCoord=currentCoordAbs + [s*step(k);0;0];
        case 2
            NewCoord=currentCoordAbs + [0;s*step(k);0];
        case 3
            NewCoord=currentCoordAbs + [0;0;s*step(k)];
    end
    NewLocationGenericMarker=osimVec3FromArray(NewCoord'); % >>>> bhorsak (08/2024): to make compatible with API 4.0% new coordinate in Vec3
    GenericMarker.set_location(NewLocationGenericMarker);%Upload new coordinates
    GenericMarker.changeFramePreserveLocation(state,SockName);%change from ground to local coordinate system
    model.updModel; %updt model
    markerset.connectToModel(model) % maybe this line is not useful
    model.print(fullfile(modelFolder,modelFile));%Save the model with new markerset
    disp('================================================================================================================');
    disp(['cycle #', num2str(k),' RMS Error:', num2str(RMSErr(k)), ' Max Error:', num2str(err(k)), ' Coord:', ToChange{1,1}, ' increment:' num2str(s*step(k))]);
    disp('================================================================================================================');
    %% end condition of while loop
    if RMSErr(k) < EndErr || k>Km && RMSErr(k) > RMSErr(k-1) && flag==0 && flag2==0 % end codition for RMSE and iterations
        break
    end
    k=k+1;
    if flag==0
        clear Str
    end
    save(fullfile(modelFolder, 'err.mat'), 'err')
    save(fullfile(modelFolder, 'RMSErr.mat'), 'RMSErr')
end
