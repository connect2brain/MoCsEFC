% With this script, test the MoCsEFC-signal-extraction-pipeline
% 
% MoCsEFC uses two simple Hjorth/Surface-Laplacian-Montages (C3-, 
% C4-Hjorth), the signal from which is spectrally filtered to 8-15Hz, and
% functional connectivity is then estimated by the PLV (from the Hilbert
% transform of the 500Hz-signal).


%% Imports
fprintf('Importing\n')
addpath(genpath('P:\2020-12 VETTERPHD Project\libraries\mambolab'));
dlDir = fullfile(getenv('USERPROFILE'), 'Downloads');


%% Configure:
cfg = [];
cfg.fsample = 500; % 500Hz as in Boss-device-"alpha"

cfg.NumParcel = 2; % left and right primary motor cortex
cfg.NumBCS = 2; % ? 2 states? -- high and low connectivity
cfg.NumNodesBCS = [2 2];
cfg.FreqBandInter = [8 13];

% Define state time series:
cfg.SeqBCS    = [1 2   1 2 1 2   1 2   1];
cfg.LenBCSloc = [3 1.5 4 2 3 2.5 2 4.2 0.8];
cfg.PhaseDelayInter = {pi/2 + [-pi/4 pi/4]; [-pi pi]};

cfg.LengthTS = sum(cfg.LenBCSloc);

% This is what actually defines the states:
cfg.C = {};
% in state 1: no connection
cfg.C{1}      = sparse(cfg.NumParcel, cfg.NumParcel);
% in state 2: parcels 1 and 2 (lM1, rM1) are connected
cfg.C{2}      = sparse(cfg.NumParcel, cfg.NumParcel);       
cfg.C{2}(1,2) = 1;



% Generate data
data = simulatedata(cfg);

plot(data.trial{1}')


%% Translate into sensor space by forward model
%% Load headmodel
load('P:/2020-12 VETTERPHD Project/ny-head/sa_nyhead.mat')

headmodel = [];
headmodel.label     = sa.clab_electrodes;
% V_fem: channels x dipoles x directions
% V_fem_normal: channels x dipoles
headmodel.leadfield = sa.cortex75K.V_fem_normal(:, sa.cortex10K.in_from_cortex75K, :);
headmodel.smesh.e   = sa.cortex10K.tri;
headmodel.smesh.p   = sa.cortex75K.vc(sa.cortex10K.in_from_cortex75K,:);

%% Pick source locations
zPos = 70;
% Slightly asymmetric sources:
targetPos = [-40,  -5, zPos; 
              35, -25, zPos];
% Symmetric sources: -35 -20 70
%targetPos = [-35 -20 zPos;
%              35 -20 zPos];

nTargets = size(targetPos, 1);
targetDistances = NaN(nTargets, size(headmodel.smesh.p, 1));
targetMasks = false(nTargets, size(headmodel.smesh.p, 1));
maxDistance = 10;
for iTarget = 1:nTargets
    targetDistances(iTarget,:) = sqrt(sum((headmodel.smesh.p - targetPos(iTarget, :)) .^2, 2));
    targetMasks(iTarget,:) = targetDistances(iTarget,:) <= maxDistance;
end

fig = figure('Position', [100, 100, 400, 400]);
selection = sqrt(sum(headmodel.smesh.p .^2, 2));
h = trisurf(headmodel.smesh.e,headmodel.smesh.p(:,1), ...
    headmodel.smesh.p(:,2),headmodel.smesh.p(:,3), ...
    selection, 'EdgeAlpha',0);
axis vis3d
xlabel('x'); ylabel('y'), zlabel('z');
hold on;
scatter3(targetPos(2,1), targetPos(2,2), targetPos(2,3), 'bx');
scatter3(targetPos(1,1), targetPos(1,2), targetPos(1,3), 'rx');

srcColors = [0 0 0.8; 
    1 0 0];

for iTarget = 1:nTargets
    scatter3(headmodel.smesh.p(targetMasks(iTarget,:), 1), ...
        headmodel.smesh.p(targetMasks(iTarget,:), 2), ...
        headmodel.smesh.p(targetMasks(iTarget,:), 3), ...
        1, srcColors(iTarget,:))
end

set(gca, 'Color', 'none')
view([0 90])
exportgraphics(fig, sprintf('%s/MoCsEFC-anatomy.pdf', dlDir), ...
    'BackgroundColor', 'none', 'ContentType', 'vector')

%% Apply forward model
sourceTimeCourses = [];
for iTarget = 1:nTargets
    dataForTarget = nan(size(data.trial{1}, 2), sum(targetMasks(iTarget, :)));
    nSources = sum(targetMasks(iTarget, :));
    for iSource = 1:nSources
        dataForTarget(:, iSource) = data.trial{1}(iTarget, :)';
    end
    sourceTimeCourses = [sourceTimeCourses dataForTarget];
end

sensorTimeCourses = (headmodel.leadfield(:, any(targetMasks, 1)) * sourceTimeCourses')';

%% TODO: apply C3/C4-Hjorth
c3hjorth = zeros(length(headmodel.label),1);
c3vsRef1   = zeros(length(headmodel.label),1);
c4hjorth = zeros(length(headmodel.label),1);
c4vsRef2 = zeros(length(headmodel.label),1);

c3hjorth(strcmpi(headmodel.label, 'C3')) = 1;
c3hjorth(ismember(headmodel.label, {'CP1', 'CP5', 'FC1', 'FC5'})) = -1/4;
c4hjorth(strcmpi(headmodel.label, 'C4')) = 1;
c4hjorth(ismember(headmodel.label, {'CP2', 'CP6', 'FC2', 'FC6'})) = -1/4;

c3vsRef1(strcmpi(headmodel.label, 'C3')) = 1;
c3vsRef1(strcmpi(headmodel.label, 'Cz')) = -1;
c4vsRef2(strcmpi(headmodel.label, 'C4')) = 1;
c4vsRef2(strcmpi(headmodel.label, 'FT9')) = -1/2;
c4vsRef2(strcmpi(headmodel.label, 'FT10')) = -1/2;

hjorthFilters = [c3hjorth c4hjorth];
differentReference = [c3vsRef1 c4vsRef2];



% Try out these two (by commenting/uncommenting right now)
timeCourses = sensorTimeCourses * hjorthFilters;
%timeCourses = sensorTimeCourses * differentReference;

% Observe: Even for different symmetric reference (e.g. C3 vs Cz, C4 vs FT9
% and FT10), the single-electrode-montages do not work! (as Marzetti
% suspected)
% But, at least for symmetric (and slightly asymmetric) sources, Hjorths work

%% Compute Hilbert transform
h = hilbert(timeCourses);
phaseDiffs = angle(h(:,1)) - angle(h(:,2));
complexRepresentation = exp(1i .* phaseDiffs);


windowLengthInS = 0.5;
windowLengthInSamples = windowLengthInS * cfg.fsample;
window = 1:round(windowLengthInSamples);

protoPLVs = movmean(complexRepresentation, windowLengthInSamples); 
PLVs = abs(protoPLVs);
iPLVs = abs(imag(protoPLVs));



%% Visualize
times = data.time{1};

stateColors = [0 0 0.3; 1 0 0];
fig = figure('Position', [50, 50, 1100, 500]);
subplot(2,3,[1 2])
hold on

label = nan(size(data.time{1}));
start = 0;
for iBCS = 1:length(cfg.SeqBCS)
    timeMask = start <= times & times <= start+cfg.LenBCSloc(iBCS);
    label(timeMask) = cfg.SeqBCS(iBCS);

    plot(times(timeMask), PLVs(timeMask), 'Color', stateColors(cfg.SeqBCS(iBCS),:));
    plot(times(timeMask), iPLVs(timeMask), ':', 'Color', stateColors(cfg.SeqBCS(iBCS),:));

    start = start+cfg.LenBCSloc(iBCS);
end

every100ms = 1:round(0.1*cfg.fsample):length(label);
PLVevery100ms   = PLVs(every100ms);
iPLVevery100ms  = iPLVs(every100ms);
classEvery100ms = label(every100ms);
timesEvery100ms = times(every100ms);

plot(timesEvery100ms(classEvery100ms == 1), PLVevery100ms(classEvery100ms == 1), 'o', 'Color', stateColors(1,:))
plot(timesEvery100ms(classEvery100ms == 2), PLVevery100ms(classEvery100ms == 2), 'o', 'Color', stateColors(2,:))
ylabel('PLV')
xlabel('time [s]')
set(gca, 'Color', 'none')

subplot(2,3,[4 5])
plot(data.time{1}, label)
ylim([0.5 2.5])
xlabel('time [s]')
yticks([1 2])
yticklabels({'off', 'on'})
ylabel('True Connectivity State')
set(gca, 'Color', 'none')

subplot(2,3,3)
hold on
kernelShape = 'triangle';
bandwidth = 0.05;
ps = linspace(0,1,100);
[f, yi] = ksdensity(PLVevery100ms(classEvery100ms == 1), ps, 'Bandwidth',bandwidth, 'kernel', kernelShape);
p_magPLVlow = plot(f, yi, 'Color', stateColors(1,:));
[f, yi] = ksdensity(PLVevery100ms(classEvery100ms == 2), ps, 'Bandwidth',bandwidth, 'kernel', kernelShape);
p_magPLVhigh = plot(f, yi, 'Color', stateColors(2,:));

[f, yi] = ksdensity(iPLVevery100ms(classEvery100ms == 1), ps, 'Bandwidth',bandwidth, 'kernel', kernelShape);
p_imPLVlow = plot(f, yi, ':', 'Color', stateColors(1,:));
[f, yi] = ksdensity(iPLVevery100ms(classEvery100ms == 2), ps, 'Bandwidth',bandwidth, 'kernel', kernelShape);
p_imPLVhigh = plot(f, yi, ':', 'Color', stateColors(2,:));

set(gca, 'Color', 'none')
legend([p_magPLVhigh, p_magPLVlow, p_imPLVhigh, p_imPLVlow], {'PLV for high', 'PLV for low', 'ImPLV for high', 'ImPLV for low'}, 'Location', 'best')

KS_PLV = kolmogorovSmirnov(PLVevery100ms(classEvery100ms == 1), PLVevery100ms(classEvery100ms == 2), bandwidth, kernelShape);
KS_iPLV = kolmogorovSmirnov(iPLVevery100ms(classEvery100ms == 1), iPLVevery100ms(classEvery100ms == 2), bandwidth, kernelShape);

fprintf('Kolmogorov-Smirnov:\n KS(PLV)  = %d \n KS(iPLV) = %d\n\n', KS_PLV, KS_iPLV)

exportgraphics(fig, sprintf('%s/MoCsEFC-sim.pdf', dlDir), ...
    'BackgroundColor', 'none', 'ContentType', 'vector')

function [dist] = kolmogorovSmirnov(x, y, bandwidth, kernel)
% KOLMOGOROVSMIRNOV Estimates Kolmogorov-Smirnov distance between
% distributions of data x and y (vectors of realizations)
ps = linspace(0,1,100);
fx = ksdensity(x, ps, 'Bandwidth',bandwidth, 'kernel', kernel, 'Function', 'cdf');
fy = ksdensity(y, ps, 'Bandwidth',bandwidth, 'kernel', kernel, 'Function', 'cdf');

dist = max(abs(fx - fy));
end