function [avgCurve, avgPlusSDCurve, avgMinusSDCurve] = getMeanSDCurves(dat)
% This file simply calculates a mean and mean +/- std waveform. 
% Note: columns are trials, rows frames.

avgCurve = mean(dat,2);
avgPlusSDCurve = avgCurve + std(dat,0,2);
avgMinusSDCurve = avgCurve - std(dat,0,2);


