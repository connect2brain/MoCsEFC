function [analyticSignalPerTrial, powerPerTrial, phase] = analyticSignalPowerAndPhasePerTrial(signal, samplingrate, triggerIndices, windowCfg, downsampleStep)
% Filter out 50 Hz noise:
% D = designfilt('bandstopfir', 'FilterOrder', 2500, ...
%     'CutoffFrequency1', 49, 'CutoffFrequency2', 51, ...
%     'SampleRate', samplingrate);
% signal = filtfilt(D, signal);

powerwindow = round(windowCfg.power(1) * samplingrate):round(windowCfg.power(2)*samplingrate);
windows = (triggerIndices + powerwindow)';
X = detrend(signal(windows), 1);
powerPerTrial = bandpower(X, samplingrate, [8 13.5]);

edge = 65;
iterationOverhang = 50; % has to be at least 1
iterations = edge + iterationOverhang;
ARorder = 20; % 10 also seems to work

% Here: Extend the window by the length of the edge that Phastimat will cut off to the left. On the
% right, phastimate will fill in this edge with predicted samples.
plvWindow = (round(windowCfg.analytic(1)*samplingrate)-(edge*downsampleStep)):round(windowCfg.analytic(2)*samplingrate);
windows = (triggerIndices + plvWindow)';
X = detrend(signal(windows), 1);
% Note: Lower filter order due to small window sizes (not identical to real-time!)
DownsampleFilter = designfilt('lowpassfir', 'FilterOrder', 300, 'CutoffFrequency', samplingrate/(2*downsampleStep), ...
    'SampleRate', samplingrate);
X = filtfilt(DownsampleFilter, X);
X = X(1:downsampleStep:end, :);

downsampledFs = samplingrate/downsampleStep;

D = designfilt('bandpassfir', 'FilterOrder', 150, ...
            'CutoffFrequency1', 8, 'CutoffFrequency2', 13.5, ...
            'SampleRate', downsampledFs);
[~, ~, analyticSignalPerTrial] = phastimate(X, D, edge, ARorder, round(diff(windowCfg.analytic)*downsampledFs) + iterationOverhang, 0, iterations);
% time x trials
analyticSignalPerTrial = analyticSignalPerTrial(1:(end-iterationOverhang),:);

phaseWindow = round(windowCfg.phase(1))*samplingrate:round(windowCfg.phase(2)*samplingrate);
windows = (triggerIndices + phaseWindow)';
X = detrend(signal(windows), 1);
% Note: Lower filter order due to small window sizes (not identical to real-time!)
DownsampleFilter = designfilt('lowpassfir', 'FilterOrder', 120, 'CutoffFrequency', samplingrate/(2*downsampleStep), ...
    'SampleRate', samplingrate);
X = filtfilt(DownsampleFilter, X);
X = X(1:downsampleStep:end, :);

phase = phastimate(X, D, 65, 25, 128, abs(round(windowCfg.phase(2) * downsampledFs)));
end
