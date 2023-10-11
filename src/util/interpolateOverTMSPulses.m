function [interpolated, removed] = interpolateOverTMSPulses(data, window, markers)
%INTERPOLATEOVERTMSPULSES linearly interpolate over the TMS-pulses, defined
%by the window (relative, in indices) and the markers (from NeurOne data
%structure)
%  data: a row-array

triggerMetaIdcs = find(strcmp(markers.type, 'Out'));
triggerIdcs = markers.index(triggerMetaIdcs);

removed = triggerIdcs + window;
xq = 1:size(data, 1);
x = setdiff(xq, removed(:));
interpolated = interp1(x, data(x)', xq)';
end

