function [MEP, signals] = evaluateBlock(subjData, varargin)
% EVALUATEBLOCK retrieves (raw) MEP-amplitudes and lM1/rM1-signals from raW
%   neurone-datastructure "subjData".
%
%  Optional argument: 
%     'InterpolationWindow': a 2-element row vector specifying the relative
%     indices of the window around the TMS-pulse to be removed and
%     interpolated (e.g. 'InterpolationWindow', [-5, 10])
%
% RETURNS:
%   MEP: a structure with fields
%     MEP.window : specifies window used, as indices
%     MEP.low : boolean mask -- each trial is marked for whether it is a
%     low-condition trial
%     MEP.high: boolean mask -- each trial is marked for whether it is a
%     high-condition trial. NOTE: this is NOT generally ~MEP.low, as there
%     can be time-out-trials, which are neither high nor low condition.
%     MEP.APB: APB-MEP-amplitudes in all trials (including timeout)
%     MEP.FDI: FDI-MEP-amplitudes in all trials (including timeout)
%     MEP.waited: time waited since the last pulse (set to 10s for the
%     first pulse)
%   
%   Example: MEP.APB(MEP.low) to retrieve the MEP-amplitudes in APB muscle
%   in the low-condition trials
%
%   signals: a structure containing the full, raw left M1 (M1_l), right M1 
%   (M1_r), APB and FDI time-courses (not epoched), and the MEP-windows 
%   (epoched EMG: FDI_windows, APB_windows)

conditionCodes = {'1', '2'};

samplingrate = subjData.properties.samplingRate;
signals = [];
signals.M1_l = spatialFilter(subjData, {'C3', 'FC1', 'FC5', 'CP1', 'CP5'}, ...
    [1 -0.25 -0.25 -0.25 -0.25]);
signals.M1_r = spatialFilter(subjData, {'C4', 'FC2', 'FC6', 'CP2', 'CP6'}, ...
    [1 -0.25 -0.25 -0.25 -0.25]);

signals.APB = subjData.signal.APBr.data;
signals.FDI = subjData.signal.FDIr.data;


triggerMetaIdcs = find(strcmp(subjData.markers.type, 'Out'));
triggerIdcs = subjData.markers.index(triggerMetaIdcs);

MEP = [];
MEP.window = round(0.02 * samplingrate):round(0.04 * samplingrate);
MEP.preInnervationWindow = round(-0.505 * samplingrate):round(-0.005 * samplingrate);
MEP.pre = [];

MEP.low  = strcmp(subjData.markers.type(triggerMetaIdcs-1), conditionCodes{1});
MEP.high = strcmp(subjData.markers.type(triggerMetaIdcs-1), conditionCodes{2});
APB_windows    = signals.APB(triggerIdcs + MEP.window)';
APB_prewindows = signals.APB(triggerIdcs + MEP.preInnervationWindow)';
signals.APB_windows = APB_windows;
MEP.APB = max(APB_windows) - min(APB_windows);
MEP.pre.APB = range(detrend(APB_prewindows, 1));

FDI_windows = signals.FDI(triggerIdcs + MEP.window)';
FDI_prewindows = signals.FDI(triggerIdcs + MEP.preInnervationWindow)';
signals.FDI_windows = FDI_windows;
MEP.FDI = range(FDI_windows);
MEP.pre.FDI = range(detrend(FDI_prewindows, 1));

MEP.waited = [10; diff(subjData.markers.time(triggerMetaIdcs))];


if nargin > 0
    assert(rem(nargin-1, 2) == 0, sprintf('Optional arguments must be given as pairs of name and value (received %d arguments)', nargin));
    interpolationWindow = false;
    for i = 1:2:(nargin-1)
        switch varargin{i}
            case 'InterpolationWindow'
                w = varargin{i+1};
                interpolationWindow = w(1):w(2);
        end
    end
    if ~islogical(interpolationWindow)
        remove = triggerIdcs + interpolationWindow;
        fprintf('Interpolating\n')
        xq = 1:size(signals.M1_l, 1);
        x = setdiff(xq, remove(:));
        signals.raw.M1_l = signals.M1_l;
        signals.M1_l = interpolateOverTMSPulses(signals.M1_l, ...
            interpolationWindow, subjData.markers);

        signals.raw.M1_r = signals.M1_r;
        signals.M1_r = interpolateOverTMSPulses(signals.M1_r, ...
            interpolationWindow, subjData.markers);

        signals.interpolated = remove;
    end
end


end
