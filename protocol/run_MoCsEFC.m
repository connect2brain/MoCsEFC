% MoCsEFC-experiment script.
% Run section by section and pay attention to comments and info printed to
% command window!
%
% Experiment consist of:
% - Cap preparation (+ EEG-marker definition, EMG-preparation)
% - Hotspot search
% - resting state EEG
% - Stimulation Block 1 (~25min)
% - Short break (subject should move a bit (stand up), then reposition
%   coil)
% - Stimulation Block 2 (~25min)
%

%% Imports
fprintf('Importing\n')
addpath(genpath('../../libraries/bossdevice-api-matlab-master/release 2017b'));
addpath(genpath('../MoCsEFC/src'));
addpath(genpath('../MAGIC-toolbox'));
addpath(genpath('M:/2020-02 SCREEN3 Experiment/code'));


%% Documentation: Always fill this in!
record = [];
record.rmt = 50; % these are default values
record.SIfactor = 1.1; % previously: 1.2 (starting on April 12th: 1.1)
record.singlePulseSI = round(record.rmt * record.SIfactor);
record.minISI = 2;
record.timeout = 6;



% F I L L  I N :
subject_number = 7;
record.session = 3; 
record.subject = sprintf('MoCsEFC_%0.3d', subject_number);





fprintf('%s\n\n MoCsEFC-experiment\n session %d on %s\n\n', repmat('_', 1, 80), record.session, record.subject)

OUT_ = ['X:\2022-01 MoCsEFC\participants\' record.subject '\session_' num2str(record.session)];
success = mkdir(OUT_);

fprintf('> [%d] created or found %s\n\n', success, OUT_)


%% Start up Boss-device
bd=bossdevice;
% Note: Potentially need to restart Bossdevice before this works

%% Configure BOSS-device
% (based strongly on SIMPAC_experiment_V2):
nChannels = 64;
bd.eeg_channels = nChannels;

c3hjorth = zeros(nChannels, 1);
c4hjorth = zeros(nChannels, 1);

c3hjorth([5 21 23 25 27]) = [1 -0.25 -0.25 -0.25 -0.25];
c4hjorth([6 22 24 26 28]) = [1 -0.25 -0.25 -0.25 -0.25];



setparam(bd.tg, 'QLY', 'eeg_artifact_threshold', [100 100]);
setparam(bd.tg, 'QLY', 'eye_artifact_threshold', 1e6);
max_instability = 1e6; 
setparam(bd.tg, 'QLY', 'inst_freq_max_instability', max_instability*ones(1,6))
bd.alpha.amplitude_min(1) = 0; bd.alpha.amplitude_max(1) = 1e6; 
bd.alpha.amplitude_min(2) = 0; bd.alpha.amplitude_max(2) = 1e6;
bd.spatial_filter_weights = [c3hjorth c4hjorth];

bd.configure_time_port_marker([0 1 100]);

bd.theta.ignore;
bd.beta.ignore;

bd.min_inter_trig_interval = 0.7;
bd.sample_and_hold_period = 0.1; % what is this?

% TODO: set frequency filter weights to extract subject-specific mu-rhythm
samplingFrequencyAlpha = 500; %Hz (for alpha, that is the samplingrate)
nyquistF = samplingFrequencyAlpha/2;
% Here can use individual mu-rhythm from Screen3
frequencyBandEdges            = [0 6 8 15 16 nyquistF] ./ nyquistF;
frequencyBandTargetAmplitudes = [0 0 1  1  0 0];
errorWeighting = [1 1 1]; % all bands equally important
bd.alpha.bpf_fir_coeffs = firls(70, frequencyBandEdges, ...
    frequencyBandTargetAmplitudes, errorWeighting);


%% Connect to Stimulators
% Taken from SIMPAC_experiment_V2

% Remote control the stimulator to set the correct burst ISI automatically
if (~exist('stimulator', 'var') || isempty('stimulator'))
    try
        stimulator = [];
        stimulator = magventure('COM3');
        stimulator.connect();
    catch
        fprintf(2, 'Problem with the MAGIC toolbox, please restart Matlab.\n');
    end
end

% Copied from SIMPAC_experiment_V2 (David E Baur, Christoph Zrenner)
% Test stimulator remote control:
fprintf('Testing and preparing stimulator remote control ...\n')
[errorOrSuccess, deviceResponse] = stimulator.setAmplitude(record.rmt, true);
assert(errorOrSuccess == 0, 'Error setting stimulus intensity')
assert(deviceResponse.amplitudePercentage_A == record.rmt, 'Unexpected stimulator response')
pause(0.5);
fprintf('\nSetting stimulator to single pulse condition....\n')
stimulator.setWaveform('Biphasic', 'Normal', 2, 1, 1, false);
pause(1);
[errorOrSuccess, deviceResponse] = stimulator.getStatus;
assert(errorOrSuccess == 0, 'Error getting status from stimulator')
assert(strcmp(deviceResponse.Waveform, 'Biphasic'), 'Stimulator response ''Biphasic'' expected')
% stimulator.arm; % TODO: why was this here?
[errorOrSuccess, deviceResponse] = stimulator.setAmplitude(record.singlePulseSI, true);
assert(errorOrSuccess == 0, 'Error setting stimulus intensity')
assert(deviceResponse.amplitudePercentage_A == record.singlePulseSI, 'Unexpected stimulator response')
fprintf('\nDone setting stimulator to single pulse condition.\n')
fprintf('\nDone testing stimulator remote control.\n')

%%
fprintf('%s\n Now sending a test-pulse -- listen for it!\n', repmat(' _', 1, 25))
pause(0.5)
% For demo-purposes
stimulator.arm;
stimulator.setAmplitude(record.singlePulseSI, true);
stimulator.fire;
pause(1)



%% Open Loop Stimulation - Hotspot search
%stimulators.right.arm;
%stimulators.right.setAmplitude(50);
stimulator.arm;
stimulator.setAmplitude(55);

% Hotspot search with triggered pulses, observe EMG
% Abort using Ctrl+C
fprintf('\nHotspot search. Press Ctrl+C to abort\n')
while(true), bd.sendPulse(1), fprintf('.'), pause(2), end



%%
record.rmt = 49; % <- fill in
record.singlePulseSI = round(record.rmt * record.SIfactor);

assert(record.singlePulseSI <= 100, 'Stimulus intensity has to be at most 100% stimulator output')

%% MEP-input-output curve
bd.scope_emg.stop();
pause(0.3);
bd.scope_emg.NumSamples = 1000;
bd.scope_emg.NumPrePostSamples = -500;
intensities = round(linspace(record.rmt, 1.4*record.rmt, 5));
[meps, fig] = MEP_ioresponse(bd, stimulator, intensities, 5);
fprintf('Done collecting MEP-io\n\n')
save([OUT_ filesep sprintf('mep_io_%s.mat', record.subject)], 'intensities', 'meps')

set(fig, 'Renderer','Painters') %export vectorized
set(fig, 'PaperUnits', 'centimeter', 'PaperSize', [16 16]) % set size
set(fig, 'PaperPositionMode', 'manual', 'PaperUnits', 'normalized', 'PaperPosition',[0 0 1 1]); % fill page
set(fig, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard
print(fig, [OUT_ filesep 'mep-io'], '-dpdf', '-r0')


%% Documentation and logging
record.starttime = datestr(now);
writetable(struct2table(record), [OUT_ filesep 'record.csv']);%, ...
    % this is apparently newer that 2017 'WriteMode', 'Append', 'WriteVariableNames', false);

fileConditions = fopen([OUT_ filesep sprintf('conditions_%s.txt', record.subject)], 'w');
fileCriteria   = fopen([OUT_ filesep sprintf('criteria_%s.txt', record.subject)], 'w');
filePLVs       = fopen([OUT_ filesep sprintf('plvs_%s.txt', record.subject)], 'w');
fileEvents     = fopen([OUT_ filesep sprintf('events_%s.txt', record.subject)], 'w');
fileJointSignal= fopen([OUT_ filesep sprintf('signal_%s.txt', record.subject)], 'w');


%% Set up experiment-data-structure: Fill in number of trials and size of PLV buffer

nInitialPLVs = 1000;
nTrials = 450; % note: this is PER condition (i.e. total: 2*nTrials)

breakAt = nTrials+1;


experiment = Experiment(nTrials, nTrials, nInitialPLVs);

plvPeriod = 0.5; %s
samplingFrequencyAlpha = 500; % recall: alpha has 0.5kHz, i.e. sample every 2ms
plvWindowLengthInSamples = round(plvPeriod*samplingFrequencyAlpha);
plvWindowIndices = (1:plvWindowLengthInSamples)-1;
% 2022-04-28: SHIFT THIS TO THE END OF THE BUFFER!
plvWindowIndices = plvWindowIndices + (samplingFrequencyAlpha - plvWindowLengthInSamples);
fprintf('  PLV-indices are: %s\n', sprintf('%d ', plvWindowIndices));

restingPLVPause = 0.5;

bufferIds = [];
bufferIds.all = getsignalid(bd.tg, 'alpha_window_buffer') + int32([plvWindowIndices 500+plvWindowIndices]);
bufferIds.s1 = 1:plvWindowLengthInSamples;
bufferIds.s2 = (plvWindowLengthInSamples+1):(2*plvWindowLengthInSamples);

fprintf('\n%s\n Setting up experiment:\n   # trials = %d\n   # PLVs in buffer = %d\n\n This will take at least %d (%d+%d) min\n', repmat(' _', 1, 25), 2*nTrials, nInitialPLVs, floor((nInitialPLVs*restingPLVPause + 2*nTrials*record.minISI)/60), floor(nInitialPLVs*restingPLVPause/60), floor(2*nTrials*record.minISI/60))

%% Check example oscillation
jointSignal = getsignal(bd.tg, bufferIds.all);
signal1 = jointSignal(bufferIds.s1);
signal2 = jointSignal(bufferIds.s2);

checkCosines = figure;
subplot(2,1,1)
plot(cos(signal1)); hold on; plot(movmean(cos(signal1), 20));
title('Check approximate sinusoid, and signals differ')
subplot(2,1,2)
plot(cos(signal2)); hold on; plot(movmean(cos(signal2), 20));


fprintf('\n Showing example-oscillations. Check them, potentially rerun.\n\n If satisfied, [\bstart NeurOne recording]\b (resting state) and continue.\n\n')

%% > Collect PLVs from resting-state EEG and log final result
close all;

fprintf(' Collecting PLV-distribution from resting state ...')

experiment.logEvent(fileEvents, sprintf('%s, Begin resting state PLV-acquisition', datetime));

nRounds = nInitialPLVs;
times = NaN(1, nRounds);
tic
for iSample = 1:nRounds
    times(iSample) = toc; % the first one is nonsense;
    tic
    jointSignal = getsignal(bd.tg, bufferIds.all);
    plv = abs(mean(exp(1i .* (jointSignal(bufferIds.s1) - jointSignal(bufferIds.s2)))));
    experiment.storePLV(plv);
    
    pause(restingPLVPause);
    % The loop body takes ~7.5ms instead of the desired 2ms!
    % This also means that the samples are discontinuous!
end
times(1) = toc; % replace first one with meaningful

experiment.logEvent(fileEvents, sprintf('%s, End resting state PLV-acquisition', datetime));
fprintf(' done.\n Now [\bstop NeurOne recording]\b!\n')

fig = figure;
h = histogram(experiment.getAllPLVs());
experiment.log(fileConditions, fileCriteria, filePLVs);

set(fig, 'Renderer','Painters') %export vectorized
set(fig, 'PaperUnits', 'centimeter', 'PaperSize', [24 8]) % set size
set(fig, 'PaperPositionMode', 'manual', 'PaperUnits', 'normalized', 'PaperPosition',[0 0 1 1]); % fill page
set(fig, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard
print(fig, [OUT_ filesep 'rsPLV'], '-dpdf', '-r0')



fprintf('\n Instruct participant for main session,\n then [\bstart NeurOne recording]\b!\n')

%% > Realtime stimulation: Block 1
conditionCodes = [1, 2]; 
timeoutCode = 100;

waitTimes = NaN(1,2*nTrials);

stimulator.arm;
stimulator.setAmplitude(record.singlePulseSI, true);
bd.prepareForFastManualFiring();

experiment.logEvent(fileEvents, sprintf('%s, Begin stimulation block 1', datetime));

tWait = tic;
tic
while experiment.CurrentTrial < breakAt && ~experiment.isDone()
    toc % To track how long processing of PLV takes   
    tic
    
    jointSignal = getsignal(bd.tg, bufferIds.all);
    plv = abs(mean(exp(1i .* (jointSignal(bufferIds.s1) - jointSignal(bufferIds.s2)))));
    
    if experiment.fire(plv)
        condition = experiment.Conditions(experiment.CurrentTrial);
        fprintf('<< Fire pulse (trial %d/%d)\n', experiment.CurrentTrial, 2*nTrials)
        
        bd.configure_time_port_marker([0, 1, conditionCodes(condition+1)]);
        bd.fire();
        fprintf('\t Pulse should have been delivered, took approx. %d.\n', toc)
        waitTimes(experiment.CurrentTrial) = toc(tWait);
        experiment.logEvent(fileEvents, sprintf('%s, trial %d: fire condition %d', datetime, experiment.CurrentTrial, condition));
        experiment.logEvent(fileJointSignal, sprintf('%d,', jointSignal));
        experiment.next();
        %
        pause(record.minISI);
        tWait = tic;
    elseif toc(tWait) > record.timeout
        fprintf('<< Fire TIMEOUT pulse\n')
        bd.configure_time_port_marker([0, 1, timeoutCode]);
        bd.fire();
        experiment.logEvent(fileEvents, sprintf('%s, timeout pulse', datetime));
        pause(record.minISI);
        tWait = tic;
    else
        experiment.logEvent(fileEvents, 'waiting');
    end
    experiment.storePLV(plv);
    experiment.log([], fileCriteria, filePLVs);
    
    pause(0.1); % to reduce overlap between windows (which skews the plv-distribution)
end
experiment.logEvent(fileEvents, sprintf('%s, Break begin', datetime));
fprintf('\nAbout 10 min Break.\nStop NeurOne recording.\nInstruct Subject to move a little bit, then place coil again.\n[\bStart NeurOne recording]\b!, then run next block\n\n')

%% > Realtime stimulation: Block 2
stimulator.arm;
stimulator.setAmplitude(record.singlePulseSI, true);
bd.prepareForFastManualFiring();

experiment.logEvent(fileEvents, sprintf('%s, Break end', datetime));
tWait = tic;
tic
while ~experiment.isDone()
    toc % To track how long processing of PLV takes   
    tic
    
    jointSignal = getsignal(bd.tg, bufferIds.all);
    plv = abs(mean(exp(1i .* (jointSignal(bufferIds.s1) - jointSignal(bufferIds.s2)))));
    
    if experiment.fire(plv)
        condition = experiment.Conditions(experiment.CurrentTrial);
        fprintf('<< Fire pulse (trial %d/%d)\n', experiment.CurrentTrial, 2*nTrials)
        
        bd.configure_time_port_marker([0, 1, conditionCodes(condition+1)]);
        bd.fire();
        fprintf('\t Pulse should have been delivered.\n')
        waitTimes(experiment.CurrentTrial) = toc(tWait);
        experiment.logEvent(fileEvents, sprintf('%s, trial %d: fire condition %d', datetime, experiment.CurrentTrial, condition));
        experiment.logEvent(fileJointSignal, sprintf('%d,', jointSignal));
        experiment.next();
        pause(record.minISI);
        tWait = tic;
    elseif toc(tWait) > record.timeout
        fprintf('<< Fire TIMEOUT pulse\n')
        bd.configure_time_port_marker([0, 1, timeoutCode]);
        bd.fire();
        experiment.logEvent(fileEvents, sprintf('%s, timeout pulse', datetime));
        pause(record.minISI);
        tWait = tic;
    else
        experiment.logEvent(fileEvents, 'waiting');
    end
    experiment.storePLV(plv);
    experiment.log([], fileCriteria, filePLVs);
    
    pause(0.1); % to reduce overlap between windows (which skews the plv-distribution)
end



timeComparison = [];
timeComparison.acquire = times;
timeComparison.wait = waitTimes;

experiment.logEvent(fileEvents, sprintf('%s, Finished', datetime));
save([OUT_ filesep 'times'], 'timeComparison')

fprintf('\n. Done\n  Don''t forget clean-up!\n\n')

%% Wrap-up

fclose(filePLVs);
fclose(fileCriteria);
fclose(fileConditions);
fclose(fileEvents);
fclose(fileJointSignal);

stimulator.disconnect;
fprintf('\n -  -  -  -  -  -  -  -  -  -  - \n Performed Clean-up\n\n')
