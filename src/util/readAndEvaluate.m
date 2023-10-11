function [MEP, signals, bandpower_l, bandpower_r, FC] = readAndEvaluate(dataPath, rows, windowCfg, varargin)
%READANDEVALUATE Summary of this function goes here
%  rows may be multiple rows, if so, their results will be concatenated
%
%  Optional argument:
%     'InterpolationWindow': a 2-element row vector specifying the relative
%     indices of the window around the TMS-pulse to be removed and
%     interpolated (e.g. 'InterpolationWindow', [-5, 10])
%     If no interpolationWindow is given, there will not be any
%     interpolation; AND PHASTIMATE WILL BE USED for phase-estimation

interpolate = false;
nFixedArgs = 3;
if nargin > 0
    assert(rem(nargin-nFixedArgs, 2) == 0, sprintf('Optional arguments must be given as pairs of name and value (received %d arguments)', nargin));
    for i = 1:2:(nargin-nFixedArgs)
        switch varargin{i}
            case 'InterpolationWindow'
                interpolate = true;
                interpolationWindow = varargin{i+1};
        end
    end
end

MEP = [];
MEP.high = [];
MEP.low = [];
MEP.APB = [];
MEP.FDI = [];
MEP.waited = [];
MEP.pre = [];
MEP.pre.APB = [];
MEP.pre.FDI = [];

signals = [];
signals.M1_l = [];
signals.M1_r = [];
signals.M1_l_filt = [];
signals.M1_r_filt = [];
signals.APB = [];
signals.FDI = [];
signals.APB_windows = [];
signals.FDI_windows = [];

bandpower_r = [];
bandpower_l = [];

FC = [];
FC.PLV = [];
FC.iPLV = [];
FC.cPLV = [];
FC.powerCorr = [];
FC.phaseC3 = [];
FC.phaseC4 = [];

nRows = size(rows, 1);
for iRow = 1:nRows
    subjData = readNeuroneFromTable(dataPath, rows(iRow,:));
    if interpolate
        [MEP_, signals_] = evaluateBlock(subjData, 'InterpolationWindow', interpolationWindow);
    else
        [MEP_, signals_] = evaluateBlock(subjData);
    end

    [bandpower_l_, bandpower_r_, signals_] = filterEpochBandpower(signals_, subjData, windowCfg);

    % PLV-computation
    doPhastimate = ~interpolate;
    FC_ = getFunctionalConnectivity(signals_, subjData, windowCfg, 10, doPhastimate);
    
    
    bandpower_r = [bandpower_r bandpower_r_];
    bandpower_l = [bandpower_l bandpower_l_];

    MEP.high   = [MEP.high;   MEP_.high];
    MEP.low    = [MEP.low;    MEP_.low];
    MEP.APB    = [MEP.APB     MEP_.APB];
    MEP.FDI    = [MEP.FDI     MEP_.FDI];
    MEP.waited = [MEP.waited; MEP_.waited];
    MEP.pre.APB = [MEP.pre.APB MEP_.pre.APB];
    MEP.pre.FDI = [MEP.pre.FDI MEP_.pre.FDI];

    signals.M1_l = [signals.M1_l; signals_.M1_l];
    signals.M1_r = [signals.M1_r; signals_.M1_r];
    signals.M1_l_filt = [signals.M1_l_filt; signals_.M1_l_filt];
    signals.M1_r_filt = [signals.M1_r_filt; signals_.M1_r_filt];
    signals.APB = [signals.APB; signals_.APB];
    signals.FDI = [signals.FDI; signals_.FDI];
    signals.APB_windows = [signals.APB_windows signals_.APB_windows];
    signals.FDI_windows = [signals.FDI_windows signals_.FDI_windows];

    FC.PLV  = [FC.PLV FC_.PLV];
    FC.iPLV = [FC.iPLV FC_.iPLV];
    FC.cPLV = [FC.cPLV FC_.cPLV];
    FC.powerCorr = [FC.powerCorr FC_.powerCorr];
    FC.phaseC3 = [FC.phaseC3 FC_.phaseC3];
    FC.phaseC4 = [FC.phaseC4 FC_.phaseC4];
    %disp(FC_.powerCorr)
end
end

