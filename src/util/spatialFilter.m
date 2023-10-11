function [timecourse] = spatialFilter(subjData, channels, weights)
% SPATIALFILTER applies given weights (an array of double-values) to the
% given channels (a cell-array of channel-names, same length)
timecourse = zeros(size(subjData.signal.(channels{1}).data));
for i = 1:length(channels)
    timecourse = timecourse + weights(i)*subjData.signal.(channels{i}).data;
end
end