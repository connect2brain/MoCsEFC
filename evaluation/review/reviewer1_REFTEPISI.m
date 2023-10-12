%% Imports
addpath(genpath('../../libraries/phastimate_mod'));
addpath(genpath('../../libraries/neurone_tools_for_matlab_1.1.3.11_mod'));
addpath(genpath('../../src'));

%%
dlDir = fullfile(getenv('USERPROFILE'), 'Downloads');
OUT_ = dlDir;

IN_ = 'B:\Experimental Data\2018-06 REFTEP\NeurOne Raw Data';
T = readtable('REFTEP_sessions.xlsx','Format','auto');
T = T(strcmpi(T.Condition, 'tep'),:);

%%
downSampleStep = 10;
windowCfg = [];
windowCfg.power    = [-0.5 0] - 0.005;
windowCfg.phase    = [-1 0] - 0.005;
windowCfg.analytic = [-1 0] - 0.005;
% Downsample-step is 10: result Fs=500Hz
width = 0.5;
windowCfg.indicesForPlv = (round((diff(windowCfg.analytic) - width)*(5000/downSampleStep))+1):round(diff(windowCfg.analytic)*(5000/downSampleStep));



%%
channels = {'C3', 'FC1', 'FC5', 'CP1', 'CP5', 'C4', 'FC2', 'FC6', 'CP2', 'CP6', 'FDIr', 'APBr'};
nRows = size(T,1);
ResultTable = [];

subjectName = '';
iBlock = 1;

for iRow = 1:nRows
    row = T(iRow,:);
    if ~strcmpi(row.Subject{:}, subjectName)
        iBlock = 1;
    end
    subjectName = row.Subject{:};

    fprintf('\n >> %s [%d]\n', subjectName, iBlock)

    session_file = row.NeurOneFolder{:};
    neurone_idx  = row.NeurOneIndex;
    subjData = module_read_neurone(sprintf('%s/%s/%s', IN_, row.Subject{:}, session_file), sessionPhaseNumber=neurone_idx);
    samplingrate = subjData.properties.samplingRate;
    triggerMetaIdcs = find(strcmp(subjData.markers.type, 'Out'));
    triggerIndices = subjData.markers.index(triggerMetaIdcs);
    triggerTimes = subjData.markers.time(triggerMetaIdcs);
    ISI = [Inf; diff(triggerTimes)];
    nTrials = length(triggerIndices);
    
    % MEP:
    window = round(0.02 * samplingrate):round(0.04 * samplingrate);
    preInnervationWindow = round(-0.505 * samplingrate):round(-0.005 * samplingrate);
    
    meps = range(subjData.signal.FDIr.data(triggerIndices + window)');
    PreRangeAPB = range(detrend(subjData.signal.APBr.data(triggerIndices + preInnervationWindow)', 1))';
    PreRangeFDI = range(detrend(subjData.signal.FDIr.data(triggerIndices + preInnervationWindow)', 1))';
    
    signal_C3H = spatialFilter(subjData, {'C3', 'FC1', 'FC5', 'CP1', 'CP5'}, [1 -1/4 -1/4 -1/4 -1/4]);
    [analytic_C3H, PowerC3H, PhaseC3H] = analyticSignalPowerAndPhasePerTrial(signal_C3H, samplingrate, triggerIndices, windowCfg, downSampleStep);
    PowerC3H = PowerC3H';
    PhaseC3H = PhaseC3H';
    
    signal_C4H = spatialFilter(subjData, {'C4', 'FC2', 'FC6', 'CP2', 'CP6'}, [1 -1/4 -1/4 -1/4 -1/4]);
    [analytic_C4H, PowerC4H, PhaseC4H] = analyticSignalPowerAndPhasePerTrial(signal_C4H, samplingrate, triggerIndices, windowCfg, downSampleStep);
    PowerC4H = PowerC4H';
    PhaseC4H = PhaseC4H';
    
    % time x trial
    cPLV = mean(exp(1i.*(angle(analytic_C3H(windowCfg.indicesForPlv,:)) - angle(analytic_C4H(windowCfg.indicesForPlv,:)))), 1)';
    PLV  = abs(cPLV);
    iPLV = abs(imag(cPLV));
    Subject = repmat(row.Subject, nTrials, 1);
    Block = repmat(iBlock, nTrials, 1);
    ResponseFDI = meps(:);
    BlockTable = table(Subject, Block, ISI, PowerC3H, PhaseC3H, PowerC4H, PhaseC4H, cPLV, PLV, iPLV, PreRangeFDI, PreRangeAPB, ResponseFDI);
    ResultTable = [ResultTable; BlockTable];
    writetable(ResultTable, sprintf('%s/REFTEP++_ISI_TEMP.txt', OUT_))
    
    iBlock = iBlock + 1;
end

filename = sprintf('%s/REFTEP++_ISI__%s.txt', OUT_, datestr(datetime, 'yyyy-mm-ddTHHMMSS'));
writetable(ResultTable, filename)
fprintf(' << Final ResultTable written to %s\n\n', filename)
