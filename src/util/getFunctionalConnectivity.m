function [FC] = getFunctionalConnectivity(signals, data, windowCfg, downSampleStep, doPhastimate)
% signals must include fields M1_l_filt/M1_r_filt (e.g. from
% filterEpochBandPower)
%

samplingrate = data.properties.samplingRate;
window = round((windowCfg.plv(1)-windowCfg.plv_overhang)*samplingrate):round((windowCfg.plv(2)+windowCfg.plv_overhang)*samplingrate);

fprintf('  Obtaining FC and phase from [%f , %f]ms\n', window(1)/(0.001*samplingrate), window(end)/(0.001*samplingrate))

stimuliFired = data.markers.index(strcmpi(data.markers.type, 'Out'));
windows = (stimuliFired + window)';
% detrending windows is very important!
C3 = detrend(signals.M1_l_filt(windows), 1);
C4 = detrend(signals.M1_r_filt(windows), 1);

% Low-pass below NEW Nyquist frequency before downsampling to avoid
% aliasing!
D = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', samplingrate/(2*downSampleStep), 'SampleRate', samplingrate);
C3 = filtfilt(D, C3);
C3 = C3(1:downSampleStep:end, :);
C4 = filtfilt(D, C4);
C4 = C4(1:downSampleStep:end, :);
samplingrate = samplingrate/downSampleStep;

if doPhastimate    
    D = designfilt('bandpassfir', 'FilterOrder', 100, 'CutoffFrequency1', 9, 'CutoffFrequency2', 13, 'SampleRate', samplingrate, 'DesignMethod', 'window');
    FC = [];
    % data=C3, D=D, edge=65, ord=25, hilbertwindow=128,
    % [offset_correction=4]
    % in MoCsEFC, the offset of the pre-stim EEG window is 5ms (seen from
    % windowCfg.plv(2)), which at 500Hz samplingrate corresponds to 2.5
    % samples. Thus, an offset-correction of 2 or 3 is appropriate.
    % results___2023-03-08T150505 is with 4 samples offset correction.
    offset_correction = floor(abs(windowCfg.plv(2)) * samplingrate);
    FC.phaseC3 = phastimate(C3, D, 65, 25, 128, offset_correction);
    FC.phaseC4 = phastimate(C4, D, 65, 25, 128, offset_correction);

    nTrial = length(FC.phaseC3);
    FC.cPLV = nan(1, nTrial);
    FC.PLV  = nan(1, nTrial);
    FC.iPLV = nan(1, nTrial);
    FC.powerCorr = nan(1, nTrial);
else
    %Dhp = designfilt('highpassfir', 'FilterOrder', 200, ...
    %    'CutoffFrequency', 6, 'SampleRate', samplingrate);

    D = designfilt('bandpassfir', 'FilterOrder', 200, ...
        'CutoffFrequency1', 8, 'CutoffFrequency2', 13.5, ...
        'SampleRate', samplingrate);

    %hpfilteredC3 = filtfilt(Dhp, C3);
    %hpfilteredC4 = filtfilt(Dhp, C4);

    filteredC3 = filtfilt(D, C3); %hpfilteredC3);
    filteredC4 = filtfilt(D, C4); %hpfilteredC4);

    hilbertC3 = hilbert(filteredC3);
    hilbertC4 = hilbert(filteredC4);

    overhangLength = round(windowCfg.plv_overhang * samplingrate);
    plvWindow = overhangLength:(size(hilbertC3,1)-overhangLength);
    phasesC3 = angle(hilbertC3(plvWindow,:));
    phasesC4 = angle(hilbertC4(plvWindow,:));
    phaseDifferences = phasesC3 - phasesC4;

    FC = [];
    FC.phaseC3 = phasesC3(end,:);
    FC.phaseC4 = phasesC4(end,:);
    FC.cPLV = mean(exp(1i .* phaseDifferences));
    FC.PLV  = abs(FC.cPLV);
    FC.iPLV = abs(imag(FC.cPLV));

    normedC3 = abs(hilbertC3(plvWindow,:));
    normedC3 = (normedC3 - mean(normedC3, 1)) ./ std(normedC3, 0, 1);
    normedC4 = abs(hilbertC4(plvWindow,:));
    normedC4 = (normedC4 - mean(normedC4, 1)) ./ std(normedC4, 0, 1);
    FC.powerCorr = mean(normedC3 .* normedC4, 1);
end
end

