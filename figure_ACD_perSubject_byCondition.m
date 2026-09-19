%% ============================
%   FIGURE: ACD per subject, all conditions overlaid
%   X-axis: beat position relative to stimulus (-3, -2, -1)
%   Y-axis: mean ΔIBI (ms) per beat position, per subject
%   Each subject = one line; each condition = one subplot
%   Dashed vertical line at stimulus onset (position 0)
% ============================

saving_figures = 1;
figuredir      = '.';

FS = 100; % Hz

targetKeys = {'cond_0_6667','cond_0_7692','cond_1_0000','cond_1_4286','cond_2_0000'};
gridLabels = {'66.67%','76.92%','100%','142.86%','200%'};

% All unique subjects across all conditions
allSubjects = {'Sub01','Sub03','Sub04','Sub05','Sub06','Sub08'};
nSubs       = numel(allSubjects);

% Color per subject — fixed across conditions
subColors = [
    0.12  0.47  0.71;   % Sub01 — blue
    0.20  0.63  0.17;   % Sub03 — green
    0.89  0.10  0.11;   % Sub04 — red
    1.00  0.50  0.05;   % Sub05 — orange
    0.42  0.24  0.60;   % Sub06 — purple
    0.65  0.34  0.16;   % Sub08 — brown
];

beatPositions = [-3, -2, -1]; % relative to stimulus onset

%% Layout
fig = figure('Color','w','Units','inches','Position',[1 1 14 4]);

nConds      = numel(targetKeys);
left_margin = 0.06;
right_margin = 0.02;
top_margin   = 0.12;
bot_margin   = 0.18;
gap          = 0.025;
usable_w     = 1 - left_margin - right_margin;
sub_w        = (usable_w - (nConds-1)*gap) / nConds;
sub_x        = left_margin + (0:nConds-1) .* (sub_w + gap);
panel_h      = 1 - top_margin - bot_margin;

%% Loop over conditions
for cIdx = 1:nConds
    cName = targetKeys{cIdx};
    if ~isfield(organizedData, cName), continue; end

    condData = organizedData.(cName);
    ax = axes('Units','normalized', ...
              'Position',[sub_x(cIdx), bot_margin, sub_w, panel_h]);
    hold on; grid on;

    % --- collect per-subject mean ΔIBI at each beat position ---
    plotted_any = false;

    for sIdx = 1:length(condData.subjects)
        sname = condData.subjects{sIdx};
        if iscell(sname), sname = sname{1}; end
        sname = char(sname);

        % Match to global subject list for color
        colorIdx = find(strcmp(allSubjects, sname));
        if isempty(colorIdx), colorIdx = 1; end

        try
            rpeak_row  = squeeze(condData.RPeak(sIdx, :));
            timing_row = squeeze(condData.TrialTiming(sIdx, :));

            peak_positions = find(rpeak_row == 1);
            stim_positions = find(timing_row > 0);

            ibi_ms = diff(peak_positions) * (1000 / FS);

            % Accumulate ΔIBI per beat position across trials
            delta_by_beat = nan(length(stim_positions), 3);

            for tIdx = 1:length(stim_positions)
                stim_samp     = stim_positions(tIdx);
                pre_stim_peaks = find(peak_positions < stim_samp);

                if numel(pre_stim_peaks) < 4, continue; end

                k_stim  = pre_stim_peaks(end);
                k_stim1 = pre_stim_peaks(end-1);
                k_stim2 = pre_stim_peaks(end-2);

                if k_stim < 2 || k_stim1 < 2 || k_stim2 < 2, continue; end
                if k_stim > length(ibi_ms)+1, continue; end

                ibi_at_stim  = ibi_ms(k_stim  - 1);
                ibi_at_stim1 = ibi_ms(k_stim1 - 1);
                ibi_at_stim2 = ibi_ms(k_stim2 - 1);
                ibi_prev_stim  = ibi_ms(k_stim  - 2);
                ibi_prev_stim1 = ibi_ms(k_stim1 - 2);
                ibi_prev_stim2 = ibi_ms(k_stim2 - 2);

                all_ibis = [ibi_at_stim, ibi_at_stim1, ibi_at_stim2, ...
                            ibi_prev_stim, ibi_prev_stim1, ibi_prev_stim2];
                if any(all_ibis < 444) || any(all_ibis > 1333), continue; end

                % Beat -3, -2, -1 = k_stim2, k_stim1, k_stim
                delta_by_beat(tIdx, 1) = ibi_at_stim2 - ibi_prev_stim2; % beat -3
                delta_by_beat(tIdx, 2) = ibi_at_stim1 - ibi_prev_stim1; % beat -2
                delta_by_beat(tIdx, 3) = ibi_at_stim  - ibi_prev_stim;  % beat -1
            end

            % Mean across valid trials
            sub_mean = nanmean(delta_by_beat, 1);

            if all(isnan(sub_mean)), continue; end

            plot(beatPositions, sub_mean, '-o', ...
                'Color', subColors(colorIdx,:), ...
                'LineWidth', 1.8, ...
                'MarkerFaceColor', subColors(colorIdx,:), ...
                'MarkerSize', 6, ...
                'DisplayName', sname);

            plotted_any = true;

        catch
            continue;
        end
    end

    % Stimulus onset line
    xline(0, '--k', 'LineWidth', 1.5, 'HandleVisibility','off');

    % Zero reference
    yline(0, ':', 'Color',[0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility','off');

    title(gridLabels{cIdx}, 'FontSize', 12, 'FontWeight','bold');
    xlabel('Beat position relative to stimulus', 'FontSize', 10);
    if cIdx == 1
        ylabel('Mean \DeltaIBI (ms)', 'FontSize', 11, 'FontWeight','bold');
    end

    xticks([-3, -2, -1, 0]);
    xticklabels({'-3','-2','-1','Stim'});
    xlim([-3.4, 0.4]);

    set(ax, 'Box','off', 'TickDir','out', 'LineWidth', 1.2);
    ax.XAxis.LineWidth = 1.2;
    ax.YAxis.LineWidth = 1.2;

    % Legend only on last panel
    if cIdx == nConds && plotted_any
        legend('Location','northeast','FontSize',8,'Box','off');
    end
end

%% Supertitle
annotation('textbox',[0 0.92 1 0.08], ...
    'String','Mean \DeltaIBI per beat position by subject and condition', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',13,'FontWeight','bold');

if saving_figures
    print(fig, '-dpng', '-r400', ...
        sprintf('%s/figure_ACD_perSubject_byCondition.png', figuredir));
    fprintf('Saved figure.\n');
end
