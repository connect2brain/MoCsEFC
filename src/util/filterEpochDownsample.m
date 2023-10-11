function Y = filterEpochDownsample(data, filterChannels, filterWeights, interpolationWindow, window, downSampleStep)
% FILTEREPOCHDOWNSAMPLE applies given spatial filter to given NeurOne data
% structure, then interpolates over TMS artifacts (as specififed by 
% <interpolationWindow>, filters out 50Hz noise, epochs as specified by 
% <window>, detrends each epoch, low-pass filters to under 250Hz, and 
% downsamples by factor <downsampleStep>
%
% Returns processed datamatrix
interpolationWindow = interpolationWindow(1):interpolationWindow(end);
Y = spatialFilter(data, filterChannels, filterWeights);
Y = interpolateOverTMSPulses(Y, interpolationWindow, data.markers);

samplingrate = data.properties.samplingRate;

D = designfilt('bandstopfir', 'FilterOrder', 5000, ...
'CutoffFrequency1', 49, 'CutoffFrequency2', 51, ...
'SampleRate', samplingrate);
Y = filtfilt(D, Y);

stimuliFired = data.markers.index(strcmpi(data.markers.type, 'Out'));
windows = (stimuliFired + window)';
Y = detrend(Y(windows), 1); %real-time does only demean (0-order)
D = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 250, 'SampleRate', samplingrate);
Y = filtfilt(D, Y);

Y = Y(1:downSampleStep:end, :);
end
