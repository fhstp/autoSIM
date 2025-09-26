% Function to calculate the depth of multi-variable functional data.
% [D, i_D, Outliers, LRT, F] = MultiFunDepth(Array, OThreshold)
%
% INPUT
% Array is a nStride x nEpoch x nVariable matrix
% For the typical gait profile and 7 strides, Array is a 
% [7 Strides x 51 or 101 Epochs x 9 Variables]
% OThreshold is the zscore threshold to decide that a stride is outlying,
% default value is 3
%
% OUTPUT
% D is the sorted depth vector averaged over epoch and Variable
% i_D is the stride number of Array that correspond to D, i_D(1) is the
% number of the deepest stride (most representative)
% Outliers is a logical array
% LRT is the zscore calculated to detect outliers for each stride
% F is the matrix of variables' depth

% The function calls a uni-variable functional depth function

% Implementation of the method presented in:
% A simple method to choose the most representative stride and detect outliers
% Gait & Posture 2014 DOI: 10.1016/j.gaitpost.2014.12.004
% by Morgan Sangeux and Julia Polak
% Copyright Morgan Sangeux 2016 - The Royal Children's Hospital

function [D, i_D, Outliers, LRT, F] = MultiFunDepth(Array, OThreshold)
    
    [nStride,nEpoch,nVariable] = size(Array);
    
    if nargin < 2
        OThreshold = 3;
    end
    
    % Check input
    if nEpoch ~= 51 && nEpoch ~= 101
        error('Unexpected number of points per curve. Only 51 or 101 points allowed!');
    end    
    
    % Prepare individual dimension depth
    F = zeros(nStride,nVariable);
    
    % Calculate depth for each variable
    for d = 1:nVariable
        uniArray = Array(:,:,d);
        [UD, i_UD] = UniFunDepth(uniArray);
        F(i_UD,d) = UD;
    end
    
    % The overall depth is the average depth over all variables 
    D = mean(F,2);
    
    % Calculate dispersion of the depth for each stride and determine if
    % one is an outlier
    MedD = median(D);
	MAD = median(abs(D-MedD));
	LRT = 0.6745*(D-MedD)/MAD;
    
    Outliers = LRT > OThreshold;
    
    % Depth values are sorted to provide the number of the deepest stride
    % in i_D(1)
    [D,i_D] = sort(D);	
end

% Function to calculate the depth of a one variable functional data.
% [D, i_D, F, time] = UniFunDepth(Array,option)
% INPUT
% Array is a nStride x nEpoch matrix
% option is an option structure with the fields:
% - method: FM[default]
% OUTPUT
% D is the sorted depth vector averaged over epoch
% i_D is the line number of Array that correspond to D
% F is the unsorted matrix of ecdf
% dF is the distance to the median for each epoch and each curve
% time is the time it took to calculate the depth matrix

function [D, i_D, time] = UniFunDepth(Array)
    
tic
    % Check inputs
    SA = size(Array);
    if numel(SA) > 2
        error('UniFunDepth is meant for samples in one dimension only');
    end
    if SA(1) >  SA(2)
        warning('First dimension greater than second');
    end
    
    F = zeros(SA);
    % Calculate depth for each epoch
    for t = 1:SA(2)
        med = median(Array(:,t));
        F(:,t) = abs(Array(:,t) - med);
    end
    
    % Integrate across epochs
    D = (sum(F,2)-(F(:,1)+F(:,end))/2)/(SA(2)-1);
    
    % Sort according to depth
    [D,i_D] = sort(D);
    
time = toc;    
end