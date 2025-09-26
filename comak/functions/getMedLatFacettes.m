function [idxMedPts, idxLatPts] = getMedLatFacettes(path_patella,matrix_accuracy,side)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get medial and lateral facette of patella %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to evaluate patella cartilage surface and split it into
% medial and lateral facetteThis function computes a set of pressure 
% variables from the *.vtp files.
%
%
% Example: 
% patella = stlread('C:\COMAK\Version080822\osimjam-master\setupFiles\Models\Geometry\lenhart2015-R-patella-cartilage.stl');
% side = 'l';
% matrix_accuracy = 0.0001;
%  
%
% [idxMedPts, idxLatPts] = getMedLatFacettes(patella,matrix_accuracy,side);
%
% Input: 
%         patella =          path to the patella cartilage stl file
%         matrix_accuracy =  paramter to set the resolution of the surface
%                            matrice in m 
%         side =             'l', or 'r' to indicate if it is a left or 
%                            right trial
%
% Output: Two logical vectors showing which point of the patella cartilage 
%         stl belongs to the lateral or medial facette.
%
% Written by:           Bernhard Guggenberger - bernhard.guggenberger2@fh-joanneum.at
%
% Last changed:         08/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

patella = stlread(path_patella);    

    %% Pressure medial/lateral
    % Get all faces for pat cart *.stl file
    
    %Get Minimum and Maximum values of the cartilage STL file
    %proximal/distal
    y_max = max(patella.Points(:,2));
    y_min = max(patella.Points(:,2)*(-1))*(-1);
    %medial/lateral
    z_max = max(patella.Points(:,3));
    z_min = max(patella.Points(:,3)*(-1))*(-1);
    
    %Calculate the size of the patella cartilage and then create a matrix for
    %the cartilage profile using the set matrix accuracy
    y_range = (y_max-y_min)/matrix_accuracy;
    z_range = (z_max-z_min)/matrix_accuracy;
    
    cartilage_profile = zeros(ceil(y_range)+1,ceil(z_range)+1);
    
    %Calculate the surface profile of the cartilage by interpolation of the
    %nearest 3 points
    for i = 1:length(cartilage_profile(1,:))
        
        z_distance = z_min+(i-1)*matrix_accuracy;
        
        for k = 1:length(cartilage_profile(:,1))
            
            y_distance = y_min+(k-1)*matrix_accuracy;
            
            distancesToPoint = sqrt((y_distance - patella.Points(:,2)).^2 + (z_distance - patella.Points(:,3)).^2);
            [~,index] = min(distancesToPoint);
    
            %Find second nearest point
            patellaPointsCopy = [patella.Points(:,1) patella.Points(:,2) patella.Points(:,3)];
            patellaPointsCopy(index,:) = [NaN NaN NaN];
    
            distancesToPoint2 = sqrt((y_distance - patellaPointsCopy(:,2)).^2 + (z_distance - patellaPointsCopy(:,3)).^2);
            [~,index2] = min(distancesToPoint2);
    
            %Find third nearest point
            patellaPointsCopy2 = [patellaPointsCopy(:,1) patellaPointsCopy(:,2) patellaPointsCopy(:,3)];
            patellaPointsCopy2(index2,:) = [NaN NaN NaN];
    
            distancesToPoint3 = sqrt((y_distance - patellaPointsCopy2(:,2)).^2 + (z_distance - patellaPointsCopy2(:,3)).^2);
            [~,index3] = min(distancesToPoint3);
            
            %Interpolate the three points
            profileInterpolation = ((patella.Points(index,1)+patella.Points(index2,1)+patella.Points(index3,1))/3)*(-1);
            cartilage_profile(k,i) = profileInterpolation;
    
        end
    end
    
    %Get the ridge, consisting of the line of highest points
    ridge = zeros(length(cartilage_profile(:,1)),1);
    ridge_index = zeros(length(cartilage_profile(:,1)),1);
    
    for i = 1:length(cartilage_profile(:,1))
        [ridge(i),ridge_index(i)] = max(cartilage_profile(i,:));
        
    end

    %Smooth the ridge
    smoothRidgeIndex = smoothdata(ridge_index);
    
    %Creating a dividing line to check if points ly on medial or lateral side
    dividingLine = zeros(length(ridge),3);
    
    for i = 1:length(ridge)
        dividingLine(i,1) = ridge(i)*(-1); 
        dividingLine(i,2) = (y_min + (i-1)*matrix_accuracy);
        dividingLine(i,3) = (z_min + (smoothRidgeIndex(i) * matrix_accuracy));
    
    end
    
    %Check which point of the patella cartilage stl belongs to which facette
    idxMedPts = nan(length(patella.Points(:,1)),1);
    idxLatPts = nan(length(patella.Points(:,1)),1);
    
    for i = 1:length(patella.Points(:,1))
        
        distancesToDividingLine = sqrt((patella.Points(i,2) - dividingLine(:,2)).^2 + (patella.Points(i,3) - dividingLine(:,3)).^2);
        [~,index] = min(distancesToDividingLine);
    
        if strcmp(side, 'r')
            if patella.Points(i,3) < dividingLine(index,3)
                idxMedPts(i) = 1;
                idxLatPts(i) = 0;
            else
                idxMedPts(i) = 0;
                idxLatPts(i) = 1;
            end
        
        else
            if patella.Points(i,3) < dividingLine(index,3)
                idxMedPts(i) = 0;
                idxLatPts(i) = 1;
            else
                idxMedPts(i) = 1;
                idxLatPts(i) = 0;
            end
        end
    end

    idxMedPts = logical(idxMedPts);
    idxLatPts = logical(idxLatPts);
end