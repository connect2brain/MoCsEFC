%% Imports
addpath(genpath('../libraries/phastimate'));
addpath('../libraries/neurone_tools_for_matlab_1.1.3.11_mod');
addpath(genpath('../src'));


%%
lowColor = [20, 0, 255] / 255;
highColor = [1, 0, 0];
gold =  [255, 179, 0] / 255;
alpha = 0.5;

dlDir = fullfile(getenv('USERPROFILE'), 'Downloads');
OUT_ = dlDir;
IN = 'B:\Experimental Data\2022-01 MoCsEFC\participants';

% For the method figure, we use example-data for illustration:
subjName = 'MoCsEFC_014';
session = 1;



WRITE_NAME = sprintf('%s_%d', subjName, session);


T = readtable([IN '\Sessions.xlsx']);

IN_ = [IN filesep subjName filesep 'session_' num2str(session)];

T_ = readtable([IN_ filesep 'signal_' subjName '.txt'],'ReadVariableNames', false);
T_ = T_(:,1:500); % cut off empty cell array column due to formatting (1,2,)
loggedsignal = table2array(T_);
clear T_


rawEvents = readtable([IN_ filesep 'events_' subjName '.txt'],'ReadVariableNames', false);
if strcmpi(subjName, 'MoCsEFC_012') && session == 2
    duplicateBreakEnds = find(startsWith(rawEvents.Var2, ' Break end'));
    rawEvents = rawEvents(setdiff(1:size(rawEvents, 1), duplicateBreakEnds(2)+1), :);
end

events = rawEvents(startsWith(rawEvents.Var2, ' trial') | startsWith(rawEvents.Var2, ' timeout') | strcmpi(rawEvents.Var1, 'waiting'),:);
triggermask = startsWith(events.Var2, ' trial'); % here only trials, not timeout!

if strcmpi(subjName, 'MoCsEFC_005') && session == 3
    triggermask = triggermask(setdiff(1:length(triggermask), 4997), :);
end

condition = [];
condition.high = endsWith(events.Var2(triggermask), '1');
condition.low = endsWith(events.Var2(triggermask), '0');

plvs = table2array(readtable([IN_ filesep 'plvs_' subjName '.txt'],'ReadVariableNames', false));

%%
iTrial = 2;
phC3 = loggedsignal(iTrial, 1:250);
phC4 = loggedsignal(iTrial, 251:end);

fig = figure('Renderer', 'painters', 'Position', [100 100 230 200]);
subplot(3,1,1)
plot(phC3, 'k.')
%ylabel('\phi_1')
yticks([-pi 0 pi])
yticklabels({'-180°', '0°', '180°'})
ylim([-pi pi])
xlim([1 250])
ax = gca;
set(ax, 'Color', 'none', 'XColor', 'none', 'TickDir', 'out')

subplot(3,1,2)
plot(phC4, 'k.')
%ylabel('\phi_2')
yticks([-pi 0 pi])
yticklabels({'-180°', '0°', '180°'})
ylim([-pi pi])
xlim([1 250])
ax = gca;
set(ax, 'Color', 'none', 'XColor', 'none', 'TickDir', 'out')

subplot(3,1,3)
plot(phC3-phC4, '.', 'Color', gold)
%ylabel('\Delta')
yticks([-pi 0 pi])
yticklabels({'-180°', '0°', '180°'})
ylim([-pi pi])
xlim([1 250])
ax = gca;
set(ax, 'Color', 'none', 'XColor', 'none', 'TickDir', 'out')

exportgraphics(fig, sprintf('%s/method-fig-phases-%s.pdf', OUT_, WRITE_NAME), ...
    'BackgroundColor', 'none', 'ContentType', 'vector')

%%
cPhasors = exp(1i.*(phC3-phC4));
cPLV = mean(cPhasors);


fig = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
plvLineWidth = 1;
polarplot((0.9 + 0.2*rand(1,250)) .* cPhasors, '.', 'Color', gold)
hold on

plvColor = 'k';
polarplot([cPLV 0], 'Color', plvColor, 'LineWidth', plvLineWidth)
rticklabels([])
thetaticklabels({'0°', '', '', '90°', '', '', '180°', '', '', '-90°', '', ''})

endBar = 1i * cPLV;
endBar = 0.08*(endBar ./ abs(endBar));
polarplot([endBar -endBar], 'Color', plvColor, 'LineWidth', plvLineWidth)
polarplot([endBar -endBar]+cPLV, 'Color', plvColor, 'LineWidth', plvLineWidth)
ax = gca;
set(ax, 'Color', 'none')
exportgraphics(fig, sprintf('%s/method-fig-plv-%s.pdf', OUT_, WRITE_NAME), ...
    'BackgroundColor', 'none', 'ContentType', 'vector')



%% Colored distribution
criteriaDistributions = plvs(1:end-1,:);
criteriaDistributions = criteriaDistributions(triggermask,:);
plvsForCurrentTrial = criteriaDistributions(iTrial,:);



Q1 = quantile(plvsForCurrentTrial, 0.25);
Q3 = quantile(plvsForCurrentTrial, 0.75);
xi = sort([linspace(0,1,100) Q1 Q3]);
lowMask = xi <= Q1;
highMask = xi >= Q3;
midMask = ~(lowMask | highMask);


df = ksdensity(plvsForCurrentTrial, xi, 'Bandwidth', 0.08);
fig = figure('Renderer', 'painters', 'Position', [100 100 400 150]);
hold on
area(xi(lowMask), df(lowMask), ...
    'FaceColor', lowColor, 'EdgeColor', lowColor, 'FaceAlpha', alpha)
area(xi(highMask), df(highMask), ...
    'FaceColor', highColor, 'EdgeColor', highColor, 'FaceAlpha', alpha)
plot(xi, df, 'k--')

ax = gca;
set(ax, 'Color', 'none', 'yColor', 'none', 'TickDir', 'out')

xticks([0 Q1 Q3 1])
xticklabels({'0', 'Q1', 'Q3', '1'})
xlabel('stPLV')

exportgraphics(fig, sprintf('%s/method-fig-criteria-%s.pdf', OUT_, WRITE_NAME), ...
    'BackgroundColor', 'none', 'ContentType', 'vector')