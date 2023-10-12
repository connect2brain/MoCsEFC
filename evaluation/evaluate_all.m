%% Imports
addpath(genpath('../libraries/phastimate'));
addpath('../libraries/neurone_tools_for_matlab_1.1.3.11_mod');
addpath(genpath('../src'));

%% Configuration and Setup:
LOG_MEP = false;
ANTIMEDIAN_MEP = false;
offset = 5/1000; % in seconds, this is the offset of the windows from the TMS-pulse: 5ms

interpolate = false; % Previous attempt at analysis, questionable results, thus turned off.
interpolationWindow = [-40, 800];
windowCfg = [];
% For phastimate alone:
windowCfg.raw = [-0.5 0] - offset; % for bandpower
windowCfg.plv = [-0.6 0] - offset;
windowCfg.plv_overhang = 0;

% Write files to Downloads folder
dlDir = fullfile(getenv('USERPROFILE'), 'Downloads');
OUT_ = dlDir;

IN_ = 'B:/Experimental Data/2022-01 MoCsEFC/participants';
T = readtable([IN_ '\Sessions.xlsx'],'Format','auto');
T = T(logical(T.Valid_),:);

subjects = unique(T.subject);
nSubjects = length(subjects);
fprintf(' Evaluating all n=%d complete participants\n\n', nSubjects)

ResponseAPB  = [];
ResponseFDI  = [];
Response     = [];
Subject      = [];
Session      = [];
Condition    = [];
MuPowerL     = [];
MuPowerR     = [];
PosthocPLV   = [];
PosthociPLV  = [];
PosthocRealPart = [];
PosthocImagPart = [];
PosthocPowCor   = [];
PosthocPhaseC3  = [];
PosthocPhaseC4  = [];
Waittime     = [];
PreRangeAPB  = [];
PreRangeFDI  = [];


%% Iterate over all subjects
% Collect all data into the above vectors
for iSubject = 1:nSubjects
    validRows = ismember(T.subject, subjects(iSubject)) & logical(T.Valid_);
    if sum(validRows) < 6
        continue
    end

    subjectRows = T(validRows,:);
    subjName = subjectRows.subject{:};
    fprintf('\nParticipant: %s\n', subjName)

    %% Determine which rows of the table contain which information
    indices = [];
    sessionOffset = min(subjectRows.session) - 1;
    for ses = 1:2
        indices(ses).eeg = find(strcmpi(subjectRows.condition, 'eeg') & subjectRows.session == ses+sessionOffset);
        indices(ses).stim_1 = find(strcmpi(subjectRows.condition, 'stim 1') & subjectRows.session == ses+sessionOffset);
        indices(ses).stim_2 = find(strcmpi(subjectRows.condition, 'stim 2') & subjectRows.session == ses+sessionOffset);
    end

    %% Read all the data
    IN_1 = [IN_ filesep subjName filesep 'session_' num2str(subjectRows(indices(1).eeg,:).session)];
    load([IN_1 filesep 'times'])
    wait_12 = timeComparison.wait;

    if interpolate
        [MEP_1, signals_1, bandpower_1_l, bandpower_1_r, FC_1] = readAndEvaluate(...
            IN_1, subjectRows(indices(1).stim_1,:), windowCfg, ...
            'InterpolationWindow', interpolationWindow);

        [MEP_2, signals_2, bandpower_2_l, bandpower_2_r, FC_2] = readAndEvaluate(...
            IN_1, subjectRows(indices(1).stim_2,:), windowCfg, ...
            'InterpolationWindow', interpolationWindow);
    else
        [MEP_1, signals_1, bandpower_1_l, bandpower_1_r, FC_1] = readAndEvaluate(...
            IN_1, subjectRows(indices(1).stim_1,:), windowCfg);

        [MEP_2, signals_2, bandpower_2_l, bandpower_2_r, FC_2] = readAndEvaluate(...
            IN_1, subjectRows(indices(1).stim_2,:), windowCfg);
    end

    [MEP_12, FC_12] = processSession(MEP_1, MEP_2, LOG_MEP, ANTIMEDIAN_MEP, FC_1, FC_2);
    % input to PCA should be: nTrials x nChannels
    [~, projected_12] = pca([MEP_12.APB' MEP_12.FDI']);

    bandpower_12_l = [bandpower_1_l bandpower_2_l];
    bandpower_12_r = [bandpower_1_r bandpower_2_r];



    IN_2 = [IN_ filesep subjName filesep 'session_' num2str(subjectRows(indices(2).eeg,:).session)];
    load([IN_2 filesep 'times'])
    wait_34 = timeComparison.wait;

    if interpolate
        [MEP_3, signals_3, bandpower_3_l, bandpower_3_r, FC_3] = readAndEvaluate(...
            IN_2, subjectRows(indices(2).stim_1,:), windowCfg, ...
            'InterpolationWindow', interpolationWindow);

        [MEP_4, signals_4, bandpower_4_l, bandpower_4_r, FC_4] = readAndEvaluate(...
            IN_2, subjectRows(indices(2).stim_2,:), windowCfg, ...
            'InterpolationWindow', interpolationWindow);
    else
        [MEP_3, signals_3, bandpower_3_l, bandpower_3_r, FC_3] = readAndEvaluate(...
            IN_2, subjectRows(indices(2).stim_1,:), windowCfg);

        [MEP_4, signals_4, bandpower_4_l, bandpower_4_r, FC_4] = readAndEvaluate(...
            IN_2, subjectRows(indices(2).stim_2,:), windowCfg);
    end


    % Subject 12 requires particular treatment due to a NeurOne-crash (device failure), that makes it necessary to
    % stitch together two measurement blocks (run ad hoc in response to the crash to get the
    % required number of trials):
    if strcmpi(subjName, 'MoCsEFC_012')
        fprintf('  Patching %s\n', subjName)
        IN_3 = [IN_ filesep subjName filesep 'session_' num2str(subjectRows(9,:).session)];
        [subjData_6] = readNeuroneFromTable(IN_3, subjectRows(9,:));
        if interpolate
            [MEP_6, signals_6] = evaluateBlock(subjData_6, 'InterpolationWindow', interpolationWindow);
        else
            [MEP_6, signals_6] = evaluateBlock(subjData_6);
        end
        [bandpower_6_l, bandpower_6_r, signals_6] = filterEpochBandpower(signals_6, subjData_6, windowCfg);
        FC_6 = getFunctionalConnectivity(signals_6, subjData_6, windowCfg, 10, ~interpolate);

        bandpower_4_l = [bandpower_4_l bandpower_6_l];
        bandpower_4_r = [bandpower_4_r bandpower_6_r];

        MEP_4.high = [MEP_4.high; MEP_6.high];
        MEP_4.low  = [MEP_4.low;  MEP_6.low];
        MEP_4.APB  = [MEP_4.APB   MEP_6.APB];
        MEP_4.FDI  = [MEP_4.FDI   MEP_6.FDI];
        MEP_4.waited = [MEP_4.waited; MEP_6.waited];
        MEP_4.pre.APB = [MEP_4.pre.APB MEP_6.pre.APB];
        MEP_4.pre.FDI = [MEP_4.pre.FDI MEP_6.pre.FDI];

        signals_4.M1_l = [signals_4.M1_l; signals_6.M1_l];
        signals_4.M1_r = [signals_4.M1_r; signals_6.M1_r];
        signals_4.APB = [signals_4.APB; signals_6.APB];
        signals_4.FDI = [signals_4.FDI; signals_6.FDI];
        signals_4.APB_windows = [signals_4.APB_windows signals_6.APB_windows];
        signals_4.FDI_windows = [signals_4.FDI_windows signals_6.FDI_windows];
        names = fieldnames(FC_6);
        nFieldnames = size(names, 1);
        for i = 1:nFieldnames
            f = names{i};
            FC_4.(f) = [FC_4.(f) FC_6.(f)];
        end
    end



    [MEP_34, FC_34] = processSession(MEP_3, MEP_4, LOG_MEP, ANTIMEDIAN_MEP, FC_3, FC_4);
    bandpower_34_l = [bandpower_3_l bandpower_4_l];
    bandpower_34_r = [bandpower_3_r bandpower_4_r];



    if strcmpi(subjName, 'MoCsEFC_012')
        fprintf('  Patching MoCsEFC_012 (removing trials)\n')
        % NeurOne crashed around trial 670, so i replace the adjacent
        % trials with the additionally recorded "session 3". (and some
        % others to get to 450:450)
        selectedTrials = setdiff(1:size(MEP_34.APB, 2), [630:760 763 764 767 628 623]);

        MEP_34.high = MEP_34.high(selectedTrials);
        MEP_34.low  = MEP_34.low(selectedTrials);
        MEP_34.APB  = MEP_34.APB(selectedTrials);
        MEP_34.FDI  = MEP_34.FDI(selectedTrials);
        MEP_34.waited = MEP_34.waited(selectedTrials);
        MEP_34.pre.APB = MEP_34.pre.APB(selectedTrials);
        MEP_34.pre.FDI = MEP_34.pre.FDI(selectedTrials);
        bandpower_34_r = bandpower_34_r(selectedTrials);
        bandpower_34_l = bandpower_34_l(selectedTrials);
        names = fieldnames(FC_1);
        nFieldnames = size(names, 1);
        for i = 1:nFieldnames
            f = names{i};
            FC_34.(f) = FC_34.(f)(selectedTrials);
        end
    end

    [~, projected_34] = pca([MEP_34.APB' MEP_34.FDI']);

    %% Format the data into the arrays, but
    % exclude timeout trials:
    valid_12 = MEP_12.high | MEP_12.low;
    valid_34 = MEP_34.high | MEP_34.low;
    ResponseFDI = [ResponseFDI; ...
        MEP_12.FDI(valid_12)';
        MEP_34.FDI(valid_34)'];
    ResponseAPB = [ResponseAPB; ...
        MEP_12.APB(valid_12)';
        MEP_34.APB(valid_34)'];
    Response = [Response;
        projected_12(valid_12);
        projected_34(valid_34)];

    Waittime = [Waittime; MEP_12.waited(valid_12); MEP_34.waited(valid_34)];
    PreRangeAPB = [PreRangeAPB; MEP_12.pre.APB(valid_12)'; MEP_34.pre.APB(valid_34)'];
    PreRangeFDI = [PreRangeFDI; MEP_12.pre.FDI(valid_12)'; MEP_34.pre.FDI(valid_34)'];

    MuPowerR = [MuPowerR, bandpower_12_r(valid_12), bandpower_34_r(valid_34)];
    MuPowerL = [MuPowerL, bandpower_12_l(valid_12), bandpower_34_l(valid_34)];

    PosthocPLV  = [PosthocPLV  FC_12.PLV(valid_12)  FC_34.PLV(valid_34)];
    PosthociPLV = [PosthociPLV FC_12.iPLV(valid_12) FC_34.iPLV(valid_34)];
    PosthocRealPart = [PosthocRealPart real(FC_12.cPLV(valid_12)) real(FC_34.cPLV(valid_34))];
    PosthocImagPart = [PosthocImagPart imag(FC_12.cPLV(valid_12)) imag(FC_34.cPLV(valid_34))];
    PosthocPowCor = [PosthocPowCor FC_12.powerCorr(valid_12) FC_34.powerCorr(valid_34)];
    PosthocPhaseC3 = [PosthocPhaseC3 FC_12.phaseC3(valid_12) FC_34.phaseC3(valid_34)];
    PosthocPhaseC4 = [PosthocPhaseC4 FC_12.phaseC4(valid_12) FC_34.phaseC4(valid_34)];

    nTrials_s1 = sum(MEP_12.low) +  sum(MEP_12.high);
    nTrials_s2 = sum(MEP_34.low) +  sum(MEP_34.high);

    fprintf('   #low  = %d  (%d + %d)\n   #high = %d  (%d + %d)\n\n', ...
        sum(MEP_12.low) + sum(MEP_34.low), sum(MEP_12.low), sum(MEP_34.low), ...
        sum(MEP_12.high) + sum(MEP_34.high), sum(MEP_12.high), sum(MEP_34.high))

    Subject = [Subject; ones(nTrials_s1 + nTrials_s2, 1) * iSubject];
    Session = [Session; ones(nTrials_s1, 1); 2*ones(nTrials_s2, 1)];

    Condition = [Condition;
        MEP_12.high(valid_12);
        MEP_34.high(valid_34)];
end

%% Format the collected arrays
Subject = nominal(Subject);
Condition = nominal(Condition);
Session = nominal(Session);
MuPowerR = MuPowerR';
MuPowerL = MuPowerL';
PosthocPLV = PosthocPLV';
PosthociPLV = PosthociPLV';
PosthocPowCor = PosthocPowCor';
PosthocRealPart = PosthocRealPart';
PosthocImagPart = PosthocImagPart';
PosthocPhaseC3 = PosthocPhaseC3';
PosthocPhaseC4 = PosthocPhaseC4';

%% And write them into a joint table
fprintf('\n\n%s', repmat('_', 1, 60))
fprintf('\n Test for FDI:\n')
dataset = table(Subject, Session, Condition, MuPowerL, MuPowerR, ...
    PosthocPLV, PosthociPLV, PosthocPowCor, PosthocRealPart, PosthocImagPart, PosthocPhaseC3, PosthocPhaseC4, Waittime, PreRangeAPB, PreRangeFDI, ResponseFDI);

writetable(dataset, sprintf('%s/results___%s.txt', OUT_, datestr(datetime, 'yyyy-mm-ddTHHMMSS')))
fprintf('\n%s\n\n', repmat('_', 1, 60))





%% Functions
function [MEP, FC] = processSession(MEP_1, MEP_2, doLog, doAntimedian, FC_1, FC_2)
% Helper function to concatenate the data from the two sessions, and optionally process the
% MEP-amplitudes within the session:
%   MEP_1 and FC_1 are the MEP and FC structures of session 1, MEP_2 and FC_2 of session 2
%   doLog: Whether to log-transform the MEP-amplitudes
%   doAntimedian: Whether to subtract the moving window median (window = 50) from the
%   (log-transformed) MEP-amplitudes. Not recommended, do this in R.
MEP = [];
MEP.high = [MEP_1.high; MEP_2.high];
MEP.low = [MEP_1.low; MEP_2.low];
MEP.raw.APB  = [MEP_1.APB MEP_2.APB];
MEP.raw.FDI  = [MEP_1.FDI MEP_2.FDI];
MEP.waited = [MEP_1.waited; MEP_2.waited];
MEP.pre = [];
MEP.pre.APB = [MEP_1.pre.APB MEP_2.pre.APB];
MEP.pre.FDI = [MEP_1.pre.FDI MEP_2.pre.FDI];

if doLog
    MEP.APB = log(MEP.raw.APB);
    MEP.FDI = log(MEP.raw.FDI);
else
    MEP.APB = MEP.raw.APB;
    MEP.FDI = MEP.raw.FDI;
end

if doAntimedian
    MEP.APB = MEP.APB - movmedian(MEP.APB, 50);
    MEP.FDI = MEP.FDI - movmedian(MEP.FDI, 50);
end

FC = [];
names = fieldnames(FC_1);
nFieldnames = size(names, 1);
for i = 1:nFieldnames
    f = names{i};
    FC.(f) = [FC_1.(f) FC_2.(f)];
end
end






