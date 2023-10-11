function [bandpower_l, bandpower_r, signals] = filterEpochBandpower(signals, data, windowCfg)
samplingrate = data.properties.samplingRate;
D = designfilt('bandstopfir', 'FilterOrder', 5000, ...
    'CutoffFrequency1', 49, 'CutoffFrequency2', 51, ...
    'SampleRate', samplingrate);
signals.M1_l_filt = filtfilt(D, signals.M1_l);
signals.M1_r_filt = filtfilt(D, signals.M1_r);

window = round(windowCfg.raw(1)*samplingrate):round(windowCfg.raw(2)*samplingrate);

stimuliFired = data.markers.index(strcmpi(data.markers.type, 'Out'));
windows = (stimuliFired + window)';
X = detrend(signals.M1_r_filt(windows), 1);
Y = detrend(signals.M1_l_filt(windows), 1);


bandpower_r = bandpower(X, samplingrate, [8 13.5]);
bandpower_l = bandpower(Y, samplingrate, [8 13.5]);
end