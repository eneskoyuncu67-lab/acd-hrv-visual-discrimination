%% ============================
%   FIGURE: ACD per subject and condition average
%   Single plot — all conditions overlaid
%   X: beat positions -3, -2, -1 (discrete)
%   Y: mean ΔIBI (ms)
%   Thin lines = subjects; thick lines = condition average
%   Dashed vertical line after -1 = stimulus onset
% ============================

saving_figures = 1;
figuredir      = '.';
FS             = 100;

targetKeys = {'cond_0_6667','cond_0_7692','cond_1_0000','cond_1_4286','cond_2_0000'};
gridLabels = {'66.67%','76.92%','100%','142.86%','200%'};
nConds     = numel(targetKeys);

% Color per condition
condColors = [
    0.85  0.20  0.20;   % 66.67%  — red
    0.95  0.60  0.10;   % 76.92%  — orange
    0.20  0.65  0.25;   % 100%    — green
    0.10  0.45  0.85;   % 142.86% — blue
    0.50  0.15  0.70;   % 200%    — purple
];

beatPositions = [-3, -2, -1];
xStim         = -0.5; % dashed line position just after -1

%% Compute per-subject mean ΔIBI at each beat position

% Storage: condSubMeans{cIdx} = [nSubs_in_cond x 3] matrix
condSubMeans = cell(nConds, 1);
condSubNames = cell(nConds, 1);

for cIdx = 1:nConds
    cName = targetKeys{cIdx};
    if ~isfield(organizedData, cName), continue; end
    condData = organizedData.(cName);

    nSubsHere = length(condData.subjects);
    subMat    = nan(nSubsHere, 3);

    for sIdx = 1:nSubsHere
        sname = condData.subjects{sIdx};
        if iscell(sname), sname = sname{1}; end

        try
            rpeak_row  = squeeze(condData.RPeak(sIdx, :));
            timing_row = squeeze(condData.TrialTiming(sIdx, :));

            peak_positions = find(rpeak_row == 1);
            stim_positions = find(timing_row > 0);
            ibi_ms         = diff(peak_positions) * (1000 / FS);

            delta_by_beat = nan(length(stim_positions), 3);

            for tIdx = 1:length(stim_positions)
                stim_samp      = stim_positions(tIdx);
                pre_stim_peaks = find(peak_positions < stim_samp);

                if numel(pre_stim_peaks) < 4, continue; end

                k_stim  = pre_stim_peaks(end);
                k_stim1 = pre_stim_peaks(end-1);
                k_stim2 = pre_stim_peaks(end-2);

                if k_stim < 2 || k_stim1 < 2 || k_stim2 < 2, continue; end
                if k_stim > length(ibi_ms) + 1,               continue; end

                ibi_at_stim    = ibi_ms(k_stim  - 1);
                ibi_at_stim1   = ibi_ms(k_stim1 - 1);
                ibi_at_stim2   = ibi_ms(k_stim2 - 1);
                ibi_prev_stim  = ibi_ms(k_stim  - 2);
                ibi_prev_stim1 = ibi_ms(k_stim1 - 2);
                ibi_prev_stim2 = ibi_ms(k_stim2 - 2);

                all_ibis = [ibi_at_stim, ibi_at_stim1, ibi_at_stim2, ...
                            ibi_prev_stim, ibi_prev_stim1, ibi_prev_stim2];
                if any(all_ibis < 444) || any(all_ibis > 1333), continue; end

                delta_by_beat(tIdx, 1) = ibi_at_stim2 - ibi_prev_stim2; % beat -3
                delta_by_beat(tIdx, 2) = ibi_at_stim1 - ibi_prev_stim1; % beat -2
                delta_by_beat(tIdx, 3) = ibi_at_stim  - ibi_prev_stim;  % beat -1
            end

            subMat(sIdx, :) = nanmean(delta_by_beat, 1);

        catch
            continue;
        end
    end

    condSubMeans{cIdx} = subMat;
    condSubNames{cIdx} = condData.subjects;
end

%% Plot
fig = figure('Color','w','Units','inches','Position',[1 1 7 5]);
ax  = axes('Units','normalized','Position',[0.12 0.13 0.72 0.78]);
hold on; grid on;

for cIdx = 1:nConds
    subMat = condSubMeans{cIdx};
    if isempty(subMat), continue; end
    col = condColors(cIdx, :);

    % Thin lines — one per subject
    for sIdx = 1:size(subMat, 1)
        row = subMat(sIdx, :);
        if all(isnan(row)), continue; end
        plot(beatPositions, row, '-', ...
            'Color', [col, 0.35], ...   % same color, low alpha
            'LineWidth', 0.9, ...
            'HandleVisibility', 'off');
    end

    % Thick line — condition average
    condMean = nanmean(subMat, 1);
    plot(beatPositions, condMean, '-o', ...
        'Color',           col, ...
        'LineWidth',       2.8, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'w', ...
        'MarkerSize',      8, ...
        'DisplayName',     gridLabels{cIdx});

    % Linear slope of the condition average across the 3 beat positions
    p          = polyfit(beatPositions, condMean, 1);
    condSlope  = p(1);

    % Annotate next to the beat -1 point, staggered vertically per condition
    yOffset = (cIdx - 3) * 2.5;
    text(-1 + 0.06, condMean(3) + yOffset, ...
        sprintf('\\beta = %.2f', condSlope), ...
        'Color', col, 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'Clipping', 'on');
end

% Stimulus onset dashed line
xline(xStim, '--k', 'LineWidth', 1.6, 'HandleVisibility','off');
text(xStim + 0.04, ax.YLim(1), 'Stim', ...
    'FontSize', 9, 'Color', [0.2 0.2 0.2], 'VerticalAlignment','bottom');

% Zero reference
yline(0, ':', 'Color',[0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility','off');

% Axes formatting
xticks(beatPositions);
xticklabels({'-3','-2','-1'});
xlim([-3.4, 0.0]);
xlabel('Beat position relative to stimulus onset', 'FontSize',12, 'FontWeight','bold');
ylabel('\DeltaIBI (ms)', 'FontSize',12, 'FontWeight','bold');
title('Anticipatory cardiac deceleration by condition', 'FontSize',13, 'FontWeight','bold');

set(ax, 'Box','off', 'TickDir','out', 'LineWidth', 1.3);
ax.XAxis.LineWidth = 1.3;
ax.YAxis.LineWidth = 1.3;

lgd = legend('Location','eastoutside', 'FontSize',10, 'Box','off');
lgd.Title.String = 'Condition';

if saving_figures
    print(fig, '-dpng', '-r400', ...
        sprintf('%s/figure_ACD_allConds_singlePlot.png', figuredir));
    fprintf('Saved figure_ACD_allConds_singlePlot.png\n');
end
