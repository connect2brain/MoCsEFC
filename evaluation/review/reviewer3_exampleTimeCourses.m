% In this script, retrieve some C3/C4-Hjorth-time courses from some of the participants to show to
% reviewer 3 in the rebuttal letter

%% Imports
addpath(genpath('../../libraries/phastimate_mod'));
addpath(genpath('../../libraries/neurone_tools_for_matlab_1.1.3.11_mod'));
addpath(genpath('../../src'));

rng('shuffle')

%%
dlDir = fullfile(getenv('USERPROFILE'), 'Downloads');
OUT_ = dlDir;

IN_ = 'B:/Experimental Data/2022-01 MoCsEFC/participants';
T = readtable([IN_ '\Sessions.xlsx'],'Format','auto');
T = T(logical(T.Valid_),:);

subjects = unique(T.subject);
nSubjects = length(subjects);


%% Pick a random subject:
% Take resting state EEG
subject = subjects(randsample(1:nSubjects, 1))

%%
validRows = ismember(T.subject, subject) & logical(T.Valid_);

% Ignore EEG-blocks for this analysis.
eegRows = T(validRows & strcmpi(T.condition, 'eeg'),:);
subjName = eegRows.subject{:};
fprintf('\nParticipant: %s\n', subjName)

row = eegRows(end,:);
session_file = row.neuroneFile{:};
neurone_idx  = row.neuroneIndex;
subjData = module_read_neurone(sprintf('%s/%s/session_%d/%s', IN_, row.subject{:}, row.session, row.neuroneFile{:}), sessionPhaseNumber=row.neuroneIndex);
samplingrate = subjData.properties.samplingRate;

C3H = spatialFilter(subjData, {'C3', 'FC1', 'FC5', 'CP1', 'CP5'}, [1 -1/4 -1/4 -1/4 -1/4]);
C4H = spatialFilter(subjData, {'C4', 'FC2', 'FC6', 'CP2', 'CP6'}, [1 -1/4 -1/4 -1/4 -1/4]);

timeAxis = (1:length(C3H)) / samplingrate;
%% If needed: Filter out 50Hz noise
C3H = spatialFilter(subjData, {'C3', 'FC1', 'FC5', 'CP1', 'CP5'}, [1 -1/4 -1/4 -1/4 -1/4]);
D = designfilt('bandstopfir', 'FilterOrder', 5000, ...
            'CutoffFrequency1', 49, 'CutoffFrequency2', 51, ...
            'SampleRate', samplingrate);
C3H = filtfilt(D, C3H);

%%

% Pick a random time within the time-axis:
windowStart = rand() * (timeAxis(end) - 1);
windowEnd = windowStart + 1; % in seconds

fig = figure('Position', [50 50 900 150]);
subplot(1,2,1)
plot(timeAxis, C3H, 'Color', '#0000BB')
xlim([windowStart windowEnd])
title('C3-Hjorth signal')
subplot(1,2,2)
plot(timeAxis, C4H, 'Color', '#EE0055')
xlim([windowStart windowEnd])
title('C4-Hjorth signal')

exportgraphics(fig, sprintf('%s/%s-EEG-example.pdf', OUT_, subjName), ...
    'BackgroundColor', 'none', 'ContentType', 'vector')