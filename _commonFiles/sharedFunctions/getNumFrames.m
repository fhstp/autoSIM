function nFrames = getNumFrames(dat)

for i = 1 : length(dat)
    [dataTmp, ~, ~] = read_opensim_mot(dat{i});
    nFramesTmp(i) = size(dataTmp,1); %height(dataTmp);
end
nFrames = max(nFramesTmp);

end