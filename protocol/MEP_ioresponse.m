function [meps, fig] = MEP_ioresponse(bd, stimulator,sampleIntensities, nSamplesEach)
%MEP_IORESPONSE Summary of this function goes here
% bd: boss-device instance
% stimulator: MAGIC-toolbox-stimulator instance
% sampleIntensities: which intensities to sample at; use linspace(min, max, nSamplePoints)
% nSamplesEach: how many samples to get at each sampleIntensity
%
% Requires M:\2020-02 SCREEN3 Experiment\code to be in path. (for
% createTimeAxis and extractMEP)

sc = [];

sc.isi = 2.5;

sc.samplingFrequency = 5000; % 5000 sampling frequency (Hz)
sc.emgTimeWindow = [-100 100]; % [-100 100] time window of the EMG trial (ms)

sc.emgTimeAxis = createTimeAxis(sc.emgTimeWindow,sc.samplingFrequency);
sc.mepTimeWindow = [18 55]; % [18 55] time window for MEP detection (ms)
sc.baselineTimeWindow = [-100 0];

nSamplePoints = length(sampleIntensities);
meps = nan(nSamplesEach, nSamplePoints, 2);
nSamplesCollected = zeros(1,nSamplePoints);

trialIntensityIdx = repmat(1:nSamplePoints, nSamplesEach, 1);
trialIntensityIdx = trialIntensityIdx(:);
nTrials = length(trialIntensityIdx);
trialIntensityIdx = trialIntensityIdx(randperm(nTrials));


fig = figure;
axAPBEMG = subplot(2,2,1);
axFDIEMG = subplot(2,2,3);
axAPBMEP = subplot(2,2,2);
axFDIMEP = subplot(2,2,4);

title(axAPBEMG, 'APB EMG');
title(axFDIEMG, 'FDI EMG');
title(axAPBMEP, 'MEPs APB');
title(axFDIMEP, 'MEPs FDI');

hold(axAPBMEP, 'on');
axAPBMEP.XLim = [min(sampleIntensities) - 5, max(sampleIntensities) + 5];
axAPBMEP.XTick = sort(sampleIntensities);

hold(axFDIMEP, 'on');
axFDIMEP.XLim = [min(sampleIntensities) - 5, max(sampleIntensities) + 5];
axFDIMEP.XTick = sort(sampleIntensities);



tic;
stimulator.arm;
for iTrial = 1:nTrials
    whichInt = trialIntensityIdx(iTrial);
    intensityNow = sampleIntensities(whichInt);
    stimulator.setAmplitude(intensityNow, true);
    pause(0.6) % from calibration_RMT.m
    marker = 100 + trialIntensityIdx(iTrial); % label by index of intensity
    bd.configure_time_port_marker([0 1 marker]);
    while toc < sc.isi
    end
    
    bd.manualTrigger();
    tic;
    
    
    % Get EMG data
    emgAPB = bd.mep(1); % channel 73
    emgFDI = bd.mep(2);
    
    disp(size(emgAPB))
    
    plot(axAPBEMG, emgAPB);
    plot(axFDIEMG, emgFDI);
    
    emgAPB = emgAPB - mean(emgAPB(sc.baselineTimeWindow(1) <= sc.emgTimeAxis & sc.emgTimeAxis <= sc.baselineTimeWindow(2)));
    mepAPB = extractMEP(emgAPB, sc.emgTimeAxis, sc.mepTimeWindow, sc.baselineTimeWindow);
    
    emgFDI = emgFDI - mean(emgFDI(sc.baselineTimeWindow(1) <= sc.emgTimeAxis & sc.emgTimeAxis <= sc.baselineTimeWindow(2)));
    mepFDI = extractMEP(emgFDI, sc.emgTimeAxis, sc.mepTimeWindow, sc.baselineTimeWindow);
    
    plot(axAPBMEP, intensityNow, mepAPB.amplitude, 'kx')
    plot(axFDIMEP, intensityNow, mepFDI.amplitude, 'kx')
    drawnow
    
    nSamplesCollected(whichInt) = nSamplesCollected(whichInt) + 1;
    
    meps(nSamplesCollected(whichInt), whichInt, 1) = mepAPB.amplitude;
    meps(nSamplesCollected(whichInt), whichInt, 2) = mepFDI.amplitude;
    
    fprintf('  MEP-ioresponse:  trial %d/%d  intensity = %d  ->  APB = %d   FDI = %d \n', ...
        iTrial, nTrials, intensityNow, mepAPB.amplitude, mepFDI.amplitude)
end

end

