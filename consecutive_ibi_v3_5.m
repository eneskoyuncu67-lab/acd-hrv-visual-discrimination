%% ============================
%   SETTINGS: Saving
% ============================
saving_figures = 1;
figuredir      = '.';

%% Load organizedData
if ~exist('organizedData','var')
    if exist('organizedData_by_condition','var')
        organizedData = organizedData_by_condition;
    elseif isfile('organizedData_by_condition.mat')
        S = load('organizedData_by_condition.mat');
        if isfield(S,'organizedData')
            organizedData = S.organizedData;
        elseif isfield(S,'organizedData_by_condition')
            organizedData = S.organizedData_by_condition;
        else
            error('organizedData not found in .mat file.');
        end
    else
        error('organizedData not found.');
    end
end
if ~exist('organizedData_by_condition','var')
    organizedData_by_condition = organizedData;
end

%% ============================
%   1. Data Processing
%
%   ΔIBI definition (heart-rate deceleration indexed per Sroufe, 1971):
%   For a given R-peak k, ΔIBI_k = IBI_after - IBI_before
%   where IBI_after  = time from peak k   to peak k+1  = ibi_ms(k)
%         IBI_before = time from peak k-1 to peak k    = ibi_ms(k-1)
%   Positive ΔIBI = interval lengthening = deceleration.
%
%   Outlier exclusion:
%   z ± 3 applied per beat position per subject before trial averaging.
%   Two-pass: (1) collect raw ΔIBI per beat, (2) z-score per column,
%   (3) NaN outliers, (4) average surviving values per trial.
% ============================

targetKeys = {'cond_0_6667', 'cond_0_7692', 'cond_1_0000', 'cond_1_4286', 'cond_2_0000'};
gridLabels = {'66.67%', '76.92%', '100%', '142.86%', '200%'};

existingFields     = fieldnames(organizedData);
validIdx           = ismember(targetKeys, existingFields);
gridConditionNames = targetKeys(validIdx);
gridLabels         = gridLabels(validIdx);

FS = 100;

all_trial_data_cell = {};

for cIdx = 1:length(gridConditionNames)
    cName    = gridConditionNames{cIdx};
    cLabel   = gridLabels{cIdx};
    condData = organizedData.(cName);

    for sIdx = 1:length(condData.subjects)
        sname = condData.subjects{sIdx};
        if iscell(sname), sname = sname{1}; end

        try
            rt_vec  = condData.RT{sIdx};
            acc_vec = condData.Acc{sIdx};
            if iscell(rt_vec),  rt_vec  = cell2mat(rt_vec);  end
            if iscell(acc_vec), acc_vec = cell2mat(acc_vec); end

            rpeak_row  = squeeze(condData.RPeak(sIdx, :));
            timing_row = squeeze(condData.TrialTiming(sIdx, :));

            stim_positions = find(timing_row > 0);
            peak_positions = find(rpeak_row == 1);

            % ibi_ms(k) = interval from peak k to peak k+1
            % so ibi_ms has length N-1 for N peaks
            % IBI_after  peak k = ibi_ms(k)
            % IBI_before peak k = ibi_ms(k-1)
            % ΔIBI at peak k   = ibi_ms(k) - ibi_ms(k-1)
            ibi_ms = diff(peak_positions) * (1000 / FS);

            nTrials       = length(stim_positions);
            % raw_deltas: [nTrials x 3], columns = beat -3, -2, -1
            raw_deltas    = nan(nTrials, 3);
            trial_meta    = nan(nTrials, 3); % [tIdx, rt, acc]

            % --- PASS 1: collect raw ΔIBI per trial per beat ---
            for tIdx = 1:nTrials
                stim_samp      = stim_positions(tIdx);
                pre_stim_peaks = find(peak_positions < stim_samp);

                % Need peaks k-1, k, k+1 for each of the 3 beats.
                % For beat at index k: need ibi_ms(k) and ibi_ms(k-1),
                % and ibi_ms(k) requires peak k+1 to exist.
                % For the 3 pre-stim beats (k_stim2, k_stim1, k_stim):
                %   k_stim  needs ibi_ms(k_stim)   -> peak k_stim+1 must exist
                %   k_stim1 needs ibi_ms(k_stim1)  -> peak k_stim1+1 must exist
                %   k_stim2 needs ibi_ms(k_stim2)  -> peak k_stim2+1 must exist
                % Also need ibi_ms(k-1) for each -> k >= 2
                % Minimum: 4 pre-stim peaks + 1 post-peak peak available

                if numel(pre_stim_peaks) < 4, continue; end

                k_stim  = pre_stim_peaks(end);
                k_stim1 = pre_stim_peaks(end-1);
                k_stim2 = pre_stim_peaks(end-2);

                % Need k-1 >= 1 (so k >= 2) and k+1 <= length(ibi_ms)
                % ibi_ms(k) = diff at index k, valid for k = 1..N-1
                % where N = numel(peak_positions)
                if k_stim  < 2 || k_stim1 < 2 || k_stim2 < 2, continue; end
                if k_stim  > length(ibi_ms), continue; end % need ibi_ms(k_stim)
                if k_stim1 > length(ibi_ms), continue; end
                if k_stim2 > length(ibi_ms), continue; end

                % Raw IBIs for plausibility check
                % IBI_before and IBI_after for each of the 3 beats
                ibi_before_stim  = ibi_ms(k_stim  - 1);
                ibi_after_stim   = ibi_ms(k_stim);
                ibi_before_stim1 = ibi_ms(k_stim1 - 1);
                ibi_after_stim1  = ibi_ms(k_stim1);
                ibi_before_stim2 = ibi_ms(k_stim2 - 1);
                ibi_after_stim2  = ibi_ms(k_stim2);

                all_ibis = [ibi_before_stim,  ibi_after_stim, ...
                            ibi_before_stim1, ibi_after_stim1, ...
                            ibi_before_stim2, ibi_after_stim2];
                if any(all_ibis < 444) || any(all_ibis > 1333), continue; end

                % ΔIBI = IBI_after - IBI_before
                raw_deltas(tIdx, 1) = ibi_after_stim2  - ibi_before_stim2;  % beat -3
                raw_deltas(tIdx, 2) = ibi_after_stim1  - ibi_before_stim1;  % beat -2
                raw_deltas(tIdx, 3) = ibi_after_stim   - ibi_before_stim;   % beat -1

                if tIdx <= numel(rt_vec) && tIdx <= numel(acc_vec)
                    trial_meta(tIdx, :) = [tIdx, rt_vec(tIdx), acc_vec(tIdx)];
                end
            end

            % --- PASS 1b: absolute ΔIBI threshold based on subject IBI SD ---
            % Any ΔIBI exceeding ±3 * SD(IBI) for this subject/condition
            % is flagged as implausible before z-scoring.
            % IBI SD computed from full valid IBI series for this block.
            valid_ibis = ibi_ms(ibi_ms >= 444 & ibi_ms <= 1333);
            if numel(valid_ibis) >= 2
                ibi_sd_thresh = 3 * std(valid_ibis);
                raw_deltas(abs(raw_deltas) > ibi_sd_thresh) = NaN;
            end

            % --- PASS 2: z ± 3 per beat column, per subject ---
            for b = 1:3
                col   = raw_deltas(:, b);
                valid = ~isnan(col);
                if sum(valid) < 2, continue; end
                mu_b  = mean(col(valid));
                sd_b  = std(col(valid));
                if sd_b == 0, continue; end
                z_b   = (col - mu_b) / sd_b;
                col(abs(z_b) > 3) = NaN;
                raw_deltas(:, b) = col;
            end

            % --- PASS 3: trial-level ACD = mean of surviving beat deltas ---
            for tIdx = 1:nTrials
                row = raw_deltas(tIdx, :);
                if all(isnan(row)), continue; end
                if any(isnan(trial_meta(tIdx, :))), continue; end

                acd_trial = nanmean(row);
                rt_t      = trial_meta(tIdx, 2);
                acc_t     = trial_meta(tIdx, 3);

                all_trial_data_cell = [all_trial_data_cell; ...
                    {cLabel, char(sname), tIdx, acd_trial, rt_t, acc_t}];
            end

        catch
            continue;
        end
    end
end

% Build table
all_trial_pooled = cell2table(all_trial_data_cell, ...
    'VariableNames', {'Condition','SubjectID','TrialIndex','ACD','RT','Acc'});
all_trial_pooled.Condition = categorical(all_trial_pooled.Condition, gridLabels);
all_trial_pooled.Acc       = categorical(all_trial_pooled.Acc, [0,1], {'Miss','Hit'});
all_trial_pooled           = rmmissing(all_trial_pooled);

condNumMap = containers.Map({'66.67%','76.92%','100%','142.86%','200%'}, ...
                             [66.67, 76.92, 100, 142.86, 200]);
all_trial_pooled.ConditionVal = cellfun(@(x) condNumMap(x), ...
    cellstr(all_trial_pooled.Condition));

all_trial_pooled_full = all_trial_pooled;

fprintf('Total trials extracted: %d\n', height(all_trial_pooled_full));
fprintf('ACD stats: mean=%.2f ms, SD=%.2f ms, min=%.2f, max=%.2f\n', ...
    mean(all_trial_pooled_full.ACD), std(all_trial_pooled_full.ACD), ...
    min(all_trial_pooled_full.ACD), max(all_trial_pooled_full.ACD));

%% Downsample for Panel C (Hit/Miss balance)
balanced_data = table();
rng('default');
for c = 1:numel(gridLabels)
    sub  = all_trial_pooled(all_trial_pooled.Condition == gridLabels{c}, :);
    mIdx = find(sub.Acc == 'Miss');
    hIdx = find(sub.Acc == 'Hit');
    if length(hIdx) > length(mIdx) && ~isempty(mIdx)
        hIdx = hIdx(randperm(length(hIdx), length(mIdx)));
    end
    balanced_data = [balanced_data; sub(sort([mIdx; hIdx]), :)];
end
all_trial_pooled_downsampled = balanced_data;


%% ============================
%   FIGURE 1: PANELS A + B
% ============================
figAB = figure('Color','w','Units','inches','Position',[1 1 6.5 5]);

panelA_ylim = [-150 150];
nConds      = numel(gridLabels);
condColors  = [0.7 0.85 0.9; 0.5 0.7 0.85; 0.3 0.55 0.8; 0.1 0.35 0.7; 0.05 0.15 0.5];

left_margin  = 0.10;
right_margin = 0.03;
row_top_bot  = 0.57;
row_top_h    = 0.35;
row_bot_bot  = 0.10;
row_bot_h    = 0.35;

usable_w = 1 - left_margin - right_margin;
gap      = 0.025;
sub_w    = (usable_w - (nConds-1)*gap) / nConds;
sub_x    = left_margin + (0:nConds-1) .* (sub_w + gap);

%% PANEL A — RT vs ACD scatter + LME fit
fprintf('\n=== ANALYSIS A: ACD ~ RT + (1|SubjectID), per condition ===\n');

for i = 1:nConds
    cLabel  = gridLabels{i};
    subData = all_trial_pooled_full(all_trial_pooled_full.Condition == cLabel, :);

    ax = axes('Units','normalized', ...
              'Position',[sub_x(i), row_top_bot, sub_w, row_top_h]);
    hold on; grid on;

    if ~isempty(subData)
        lme  = fitlme(subData, 'ACD ~ RT + (1|SubjectID)', 'FitMethod','ML');
        beta = lme.Coefficients.Estimate(2);
        r2   = lme.Rsquared.Ordinary;
        pval = lme.Coefficients.pValue(2);

        fprintf('\n--- Condition: %s (n=%d) ---\n', cLabel, height(subData));
        disp(lme.Coefficients);
        anova_A = anova(lme, 'DFMethod','Satterthwaite');
        fprintf('ANOVA (Type 3):\n'); disp(anova_A);
        fprintf('R2=%.4f | AIC=%.1f | BIC=%.1f | LogLik=%.1f\n', ...
            r2, lme.ModelCriterion.AIC, lme.ModelCriterion.BIC, lme.LogLikelihood);

        scatter(subData.RT, subData.ACD, 10, [0.5 0.5 0.5], ...
            'filled', 'MarkerFaceAlpha', 0.25);
        xl = linspace(min(subData.RT), max(subData.RT), 30);
        yl = lme.Coefficients.Estimate(1) + beta * xl;
        plot(xl, yl, 'r', 'LineWidth', 2);

        if pval < 0.001, pStr = 'p < 0.001';
        else,            pStr = sprintf('p = %.3f', pval); end
        statStr = sprintf('R^2 = %.3f\n\\beta = %.3f\n%s', r2, beta, pStr);
        text(0.05, 0.95, statStr, 'Units','normalized', 'FontSize',9, ...
            'VerticalAlignment','top', 'BackgroundColor','none', 'EdgeColor','none');
    end

    ylim(panelA_ylim);
    title(cLabel);
    xlabel('RT (ms)');
    if i == 1, ylabel('ACD (ms)'); end
    set(ax, 'Box','off', 'TickDir','out', 'LineWidth', 1.5);
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.LineWidth = 1.5;
    if i == 1
        annotation('textbox', ...
            [ax.Position(1)-0.04, ax.Position(2)+ax.Position(4)+0.02, 0.04, 0.04], ...
            'String','A','LineStyle','none','FontSize',12,'FontWeight','bold');
    end
end

%% PANEL B — Violin: ACD by Condition
axB = axes('Units','normalized', ...
           'Position',[left_margin, row_bot_bot, usable_w, row_bot_h]);
hold on; grid on;

fprintf('\n=== ANALYSIS B: ACD ~ RT + ConditionVal + (1|SubjectID) ===\n');

lme_B = fitlme(all_trial_pooled_full, 'ACD ~ RT + ConditionVal + (1|SubjectID)', 'FitMethod','ML');
disp(lme_B);
anovaResults_B = anova(lme_B, 'DFMethod', 'Satterthwaite');
fprintf('\nANOVA (Type 3):\n'); disp(anovaResults_B);
fprintf('R2=%.4f | AIC=%.1f | BIC=%.1f | LogLik=%.1f\n', ...
    lme_B.Rsquared.Ordinary, lme_B.ModelCriterion.AIC, ...
    lme_B.ModelCriterion.BIC, lme_B.LogLikelihood);

beta_RT_B   = lme_B.Coefficients.Estimate(2);
pval_RT_B   = lme_B.Coefficients.pValue(2);
r2_B        = lme_B.Rsquared.Ordinary;
pval_Cond_B = anovaResults_B.pValue(3);

for i = 1:nConds
    vals = all_trial_pooled_full.ACD(all_trial_pooled_full.Condition == gridLabels{i});
    [f,xi] = ksdensity(vals);
    f = f / max(f) * 0.35;
    fill(i + [f -fliplr(f)], [xi fliplr(xi)], condColors(i,:), ...
        'FaceAlpha',0.7,'EdgeColor','k','LineWidth',1);
    scatter(i + (rand(size(vals))-0.5)*0.1, vals, 8, 'k', ...
        'filled', 'MarkerFaceAlpha',0.2);
end

xticks(1:nConds); xticklabels(gridLabels);
ylabel('ACD (ms)');
xlabel('Experimental Condition');
yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2);

if pval_RT_B < 0.001,   pStr_RT   = 'p_{RT} < 0.001';
else,                    pStr_RT   = sprintf('p_{RT} = %.3f', pval_RT_B); end
if pval_Cond_B < 0.001, pStr_Cond = 'p_{Cond} < 0.001';
else,                    pStr_Cond = sprintf('p_{Cond} = %.3f', pval_Cond_B); end

statStr_B = sprintf('R^2 = %.3f\n\\beta_{RT} = %.3f\n%s\n%s', ...
    r2_B, beta_RT_B, pStr_RT, pStr_Cond);
text(0.02, 0.98, statStr_B, 'Units','normalized', 'FontSize',11, ...
    'VerticalAlignment','top', 'BackgroundColor','none', 'EdgeColor','none');

set(axB, 'Box','off', 'TickDir','out', 'LineWidth', 1.5);
axB.XAxis.LineWidth = 1.5;
axB.YAxis.LineWidth = 1.5;
annotation('textbox', ...
    [axB.Position(1)-0.04, axB.Position(2)+axB.Position(4)+0.02, 0.04, 0.04], ...
    'String','B','LineStyle','none','FontSize',12,'FontWeight','bold');

if saving_figures == 1
    print(figAB, '-dpng', '-r400', sprintf('%s/figure_AB_ACD_RT_Condition.png', figuredir));
end


%% ============================
%   FIGURE 2: PANEL C (Accuracy x Condition)
% ============================
figC = figure('Color','w','Units','inches','Position',[1 1 6.5 5]);
axC  = axes('Units','normalized','Position',[0.10 0.12 0.82 0.80]);
hold on; grid on;

fprintf('\n=== ANALYSIS C: ACD ~ Acc * Condition + (1|SubjectID) ===\n');
fprintf('    (downsampled for hit/miss balance — accuracy analysis only)\n');

lme_C = fitlme(all_trial_pooled_downsampled, 'ACD ~ Acc * Condition + (1|SubjectID)', 'FitMethod','ML');
disp(lme_C);
anovaResults_C = anova(lme_C, 'DFMethod', 'Satterthwaite');
fprintf('\nANOVA (Type 3):\n'); disp(anovaResults_C);
fprintf('R2=%.4f | AIC=%.1f | BIC=%.1f | LogLik=%.1f\n', ...
    lme_C.Rsquared.Ordinary, lme_C.ModelCriterion.AIC, ...
    lme_C.ModelCriterion.BIC, lme_C.LogLikelihood);

r2_C       = lme_C.Rsquared.Ordinary;
acc_row    = strcmp(lme_C.Coefficients.Name, 'Acc_Hit');
beta_Acc_C = lme_C.Coefficients.Estimate(acc_row);
pval_Acc_C = lme_C.Coefficients.pValue(acc_row);
pval_Interaction_C = anovaResults_C.pValue(4);

colors     = [0.9 0.4 0.4; 0.2 0.6 0.2];
edgeColors = [0.5 0.1 0.1; 0.1 0.3 0.1];
accLevels  = {'Miss','Hit'};

for i = 1:nConds
    x_center = i;
    for aIdx = 1:2
        data = all_trial_pooled_downsampled.ACD( ...
            all_trial_pooled_downsampled.Condition == gridLabels{i} & ...
            all_trial_pooled_downsampled.Acc == accLevels{aIdx});

        if ~isempty(data)
            [f,xi] = ksdensity(data);
            f = f / max(f) * 0.38;
            side = (aIdx*2-3);
            fill(x_center + side*f, xi, colors(aIdx,:), ...
                'FaceAlpha',0.45,'EdgeColor',edgeColors(aIdx,:),'LineWidth',1.5);
            scatter(x_center + side*0.15 + (rand(size(data))-0.5)*0.15, ...
                data, 10, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha',0.3);
        end
    end
end

xticks(1:nConds); xticklabels(gridLabels);
ylabel('ACD (ms)', 'FontSize',12, 'FontWeight','bold');
xlabel('Experimental Condition', 'FontSize',12, 'FontWeight','bold');
yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2);

if pval_Acc_C < 0.001,         pStr_Acc = 'p_{Acc} < 0.001';
else,                           pStr_Acc = sprintf('p_{Acc} = %.3f', pval_Acc_C); end
if pval_Interaction_C < 0.001, pStr_Int = 'p_{Int} < 0.001';
else,                           pStr_Int = sprintf('p_{Int} = %.3f', pval_Interaction_C); end

statStr_C = sprintf('R^2 = %.3f\n\\beta_{Acc} = %.3f\n%s\n%s', ...
    r2_C, beta_Acc_C, pStr_Acc, pStr_Int);
text(0.02, 0.98, statStr_C, 'Units','normalized', 'FontSize',11, ...
    'VerticalAlignment','top', 'BackgroundColor','none', 'EdgeColor','none');

set(axC, 'Box','off', 'TickDir','out', 'LineWidth', 1.5);

annotation('rectangle',[0.845 0.66 0.025 0.045], ...
    'FaceColor',[0.2 0.6 0.2],'EdgeColor',[0.1 0.3 0.1],'LineWidth',1.5);
annotation('textbox',[0.875 0.655 0.10 0.055], ...
    'String','Hit','EdgeColor','none','FontSize',12,'FontWeight','bold', ...
    'VerticalAlignment','middle');
annotation('rectangle',[0.845 0.59 0.025 0.045], ...
    'FaceColor',[0.9 0.4 0.4],'EdgeColor',[0.5 0.1 0.1],'LineWidth',1.5);
annotation('textbox',[0.875 0.585 0.10 0.055], ...
    'String','Miss','EdgeColor','none','FontSize',12,'FontWeight','bold', ...
    'VerticalAlignment','middle');

if saving_figures == 1
    print(figC, '-dpng', '-r400', sprintf('%s/figure_C_HitMiss_ACD.png', figuredir));
end


%% ============================
%   FIGURE 3: Block-level HRV (Time + Frequency Domain)
%
%   One value per subject per condition, computed over the full block.
%   Time domain: mean IBI, mean HR, SDNN, RMSSD, SDSD, pNN50.
%   Frequency domain: LF (0.04-0.15 Hz), HF (0.15-0.40 Hz), LF/HF,
%   total power — via cubic spline interpolation at 4 Hz + Welch PSD.
%   z ± 3 applied per participant across their condition values
%   for each metric before group summary.
%   PSD surfaces stored separately for Figure 3B (3D plots).
% ============================

hrv_keys      = {'cond_0_6667','cond_0_7692','cond_1_0000','cond_1_4286','cond_2_0000'};
hrv_labels    = {'66.67%','76.92%','100%','142.86%','200%'};
hrv_numeric   = [66.67, 76.92, 100, 142.86, 200];
nHRVConds     = numel(hrv_keys);

% Welch parameters for frequency domain
fs_interp  = 4;        % Hz — interpolation rate for IBI series
welch_win  = 60 * fs_interp;   % 60-second Welch window (samples at 4 Hz)
welch_ovlp = round(welch_win * 0.5);
nfft       = 2^nextpow2(welch_win * 2);
freq_res   = fs_interp / nfft;
freq_axis  = (0 : nfft/2) * freq_res;

lf_band = [0.04, 0.15];
hf_band = [0.15, 0.40];

% Storage
hrv_cell = {};   % {condLabel, subjectID, meanIBI, meanHR, SDNN, RMSSD, SDSD, pNN50, LF, HF, LFHF, TotalPow}

% PSD surfaces: cell array per condition, each entry = [nFreq x nWindows_across_subs]
psd_surfaces = cell(nHRVConds, 1);

fprintf('\n=== BLOCK-LEVEL HRV COMPUTATION ===\n');

for cIdx = 1:nHRVConds
    cName = hrv_keys{cIdx};
    if ~isfield(organizedData, cName), continue; end
    condData = organizedData.(cName);
    cLabel   = hrv_labels{cIdx};

    psd_all_windows = []; % accumulate Welch windows across subjects

    for sIdx = 1:length(condData.subjects)
        sname = condData.subjects{sIdx};
        if iscell(sname), sname = sname{1}; end
        sname = char(sname);

        try
            rpeak_row      = squeeze(condData.RPeak(sIdx, :));
            peak_positions = find(rpeak_row == 1);

            if numel(peak_positions) < 10, continue; end

            % Full block IBI series in ms
            rr_ms = diff(peak_positions) * (1000 / FS);

            % z ± 3 on raw IBI per participant before anything else
            rr_mean = mean(rr_ms);
            rr_sd   = std(rr_ms);
            rr_ms   = rr_ms(abs(rr_ms - rr_mean) <= 3 * rr_sd);

            % Physiological plausibility (45-135 bpm)
            rr_ms = rr_ms(rr_ms >= 444 & rr_ms <= 1333);

            if numel(rr_ms) < 10, continue; end

            % ---- TIME DOMAIN ----
            mean_ibi = mean(rr_ms);
            mean_hr  = 60000 / mean_ibi;
            sdnn     = std(rr_ms);
            rmssd    = sqrt(mean(diff(rr_ms).^2));
            sdsd     = std(diff(rr_ms));
            nn50     = sum(abs(diff(rr_ms)) > 50);
            pnn50    = 100 * nn50 / (numel(rr_ms) - 1);

            % ---- FREQUENCY DOMAIN ----
            % Cumulative time axis for IBI series (in seconds)
            % Each IBI ends at the cumulative sum of preceding intervals
            t_ibi = cumsum(rr_ms) / 1000; % seconds
            t_ibi = t_ibi - t_ibi(1);     % start at 0

            % Cubic spline interpolation at fs_interp Hz
            t_uniform = (0 : 1/fs_interp : t_ibi(end))';
            if numel(t_uniform) < welch_win + 1, continue; end

            rr_interp = interp1(t_ibi, rr_ms, t_uniform, 'spline');

            % Welch PSD
            [pxx, f_welch] = pwelch(rr_interp - mean(rr_interp), ...
                hanning(welch_win), welch_ovlp, nfft, fs_interp);

            % Store PSD windows for 3D surface
            % pwelch returns one PSD averaged over windows; for the 3D
            % surface we need per-window PSDs — use spectrogram instead
            [~, ~, ~, psd_matrix] = spectrogram(rr_interp - mean(rr_interp), ...
                hanning(welch_win), welch_ovlp, nfft, fs_interp);
            % psd_matrix: [nfft/2+1 x nWindows], power in ms^2/Hz
            psd_matrix = abs(psd_matrix).^2 / (fs_interp * welch_win);
            psd_all_windows = [psd_all_windows, psd_matrix];

            % Band powers (trapezoidal integration)
            lf_idx  = f_welch >= lf_band(1) & f_welch <= lf_band(2);
            hf_idx  = f_welch >= hf_band(1) & f_welch <= hf_band(2);
            lf_pow  = trapz(f_welch(lf_idx), pxx(lf_idx));
            hf_pow  = trapz(f_welch(hf_idx), pxx(hf_idx));
            tot_pow = trapz(f_welch, pxx);
            lf_hf   = lf_pow / max(hf_pow, eps);

            hrv_cell = [hrv_cell; {cLabel, sname, cIdx, ...
                mean_ibi, mean_hr, sdnn, rmssd, sdsd, pnn50, ...
                lf_pow, hf_pow, lf_hf, tot_pow}];

        catch ME
            fprintf('  Skipped %s cond %s: %s\n', sname, cLabel, ME.message);
            continue;
        end
    end

    % Trim PSD surface to match freq_axis length
    nFreqRows = numel(freq_axis);
    if ~isempty(psd_all_windows)
        psd_surfaces{cIdx} = psd_all_windows(1:min(nFreqRows,end), :);
    end
end

% Build HRV table
hrv_table = cell2table(hrv_cell, 'VariableNames', ...
    {'Condition','SubjectID','CondIdx', ...
     'MeanIBI','MeanHR','SDNN','RMSSD','SDSD','pNN50', ...
     'LF','HF','LFHF','TotalPow'});
hrv_table.Condition = categorical(hrv_table.Condition, hrv_labels);

% z ± 3 per participant across conditions for each metric
hrv_metrics = {'MeanIBI','MeanHR','SDNN','RMSSD','SDSD','pNN50','LF','HF','LFHF','TotalPow'};
subjects_hrv = unique(hrv_table.SubjectID);

for m = 1:numel(hrv_metrics)
    metric = hrv_metrics{m};
    for s = 1:numel(subjects_hrv)
        rows = strcmp(hrv_table.SubjectID, subjects_hrv{s});
        vals = hrv_table.(metric)(rows);
        if numel(vals) < 2, continue; end
        mu_s = mean(vals);
        sd_s = std(vals);
        if sd_s == 0, continue; end
        z_s  = (vals - mu_s) / sd_s;
        vals(abs(z_s) > 3) = NaN;
        hrv_table.(metric)(rows) = vals;
    end
end

fprintf('Block-level HRV table: %d rows\n', height(hrv_table));
disp(grpstats(hrv_table, 'Condition', {'mean','std'}, 'DataVars', hrv_metrics));

%% HRV LME MODELS: metric ~ ConditionVal + (1|SubjectID)
fprintf('\n=== HRV LME MODELS: metric ~ ConditionVal + (1|SubjectID) ===\n');

condValMap = [66.67, 76.92, 100, 142.86, 200];

for m = 1:numel(hrv_metrics)
    metric = hrv_metrics{m};

    tbl = hrv_table(:, {'SubjectID','CondIdx', metric});
    tbl = tbl(~isnan(tbl.(metric)), :);
    tbl.ConditionVal = condValMap(tbl.CondIdx)';

    if height(tbl) < 4, continue; end

    try
        lme_h  = fitlme(tbl, [metric ' ~ ConditionVal + (1|SubjectID)'], ...
                        'FitMethod','ML');
        anov_h = anova(lme_h, 'DFMethod','Satterthwaite');

        fprintf('\n--- %s ---\n', metric);
        disp(lme_h.Coefficients);
        fprintf('ANOVA (Type 3):\n'); disp(anov_h);
        fprintf('R2=%.4f | AIC=%.1f | BIC=%.1f | LogLik=%.1f\n', ...
            lme_h.Rsquared.Ordinary, lme_h.ModelCriterion.AIC, ...
            lme_h.ModelCriterion.BIC, lme_h.LogLikelihood);
    catch ME
        fprintf('--- %s --- ERROR: %s\n', metric, ME.message);
    end
end

%% FIGURE 3A — Time-domain HRV bar chart per condition
fig3A = figure('Color','w','Units','inches','Position',[1 1 12 8]);

td_metrics  = {'MeanIBI','MeanHR','SDNN','RMSSD','SDSD','pNN50'};
td_labels   = {'Mean IBI (ms)','Mean HR (bpm)','SDNN (ms)','RMSSD (ms)','SDSD (ms)','pNN50 (%)'};
nTD         = numel(td_metrics);

for m = 1:nTD
    ax = subplot(2, 3, m);
    hold on; grid on;

    for cIdx = 1:nHRVConds
        rows = hrv_table.Condition == hrv_labels{cIdx};
        vals = hrv_table.(td_metrics{m})(rows);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        mu_c  = mean(vals);
        sem_c = std(vals) / sqrt(numel(vals));

        bar(cIdx, mu_c, 0.6, 'FaceColor', condColors(cIdx,:), ...
            'EdgeColor', condColors(cIdx,:)*0.7, 'FaceAlpha', 0.8);
        errorbar(cIdx, mu_c, sem_c, 'k', 'LineWidth', 1.5, 'CapSize', 6);

        % Individual subject points
        scatter(cIdx + (rand(size(vals))-0.5)*0.15, vals, 30, 'k', ...
            'filled', 'MarkerFaceAlpha', 0.5);
    end

    xticks(1:nHRVConds); xticklabels(hrv_labels);
    ylabel(td_labels{m}, 'FontSize', 10, 'FontWeight','bold');
    set(ax, 'Box','off', 'TickDir','out', 'LineWidth', 1.2);
    title(td_labels{m}, 'FontSize', 10);
end

sgtitle('Block-level Time-Domain HRV by Condition', 'FontSize', 13, 'FontWeight','bold');

if saving_figures
    print(fig3A, '-dpng', '-r400', sprintf('%s/figure_HRV_TimeDomain.png', figuredir));
end

%% FIGURE 3B — Frequency-domain HRV bar chart
fig3B = figure('Color','w','Units','inches','Position',[1 1 10 4]);

fd_metrics = {'LF','HF','LFHF','TotalPow'};
fd_labels  = {'LF Power (ms^2)','HF Power (ms^2)','LF/HF Ratio','Total Power (ms^2)'};
nFD        = numel(fd_metrics);

for m = 1:nFD
    ax = subplot(1, 4, m);
    hold on; grid on;

    for cIdx = 1:nHRVConds
        rows = hrv_table.Condition == hrv_labels{cIdx};
        vals = hrv_table.(fd_metrics{m})(rows);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        mu_c  = mean(vals);
        sem_c = std(vals) / sqrt(numel(vals));

        bar(cIdx, mu_c, 0.6, 'FaceColor', condColors(cIdx,:), ...
            'EdgeColor', condColors(cIdx,:)*0.7, 'FaceAlpha', 0.8);
        errorbar(cIdx, mu_c, sem_c, 'k', 'LineWidth', 1.5, 'CapSize', 6);
        scatter(cIdx + (rand(size(vals))-0.5)*0.15, vals, 30, 'k', ...
            'filled', 'MarkerFaceAlpha', 0.5);
    end

    xticks(1:nHRVConds); xticklabels(hrv_labels);
    xtickangle(30);
    ylabel(fd_labels{m}, 'FontSize', 9, 'FontWeight','bold');
    set(ax, 'Box','off', 'TickDir','out', 'LineWidth', 1.2);
    title(fd_labels{m}, 'FontSize', 9);
end

sgtitle('Block-level Frequency-Domain HRV by Condition', 'FontSize', 13, 'FontWeight','bold');

if saving_figures
    print(fig3B, '-dpng', '-r400', sprintf('%s/figure_HRV_FreqDomain.png', figuredir));
end

%% FIGURE 3C — 3D PSD surfaces per condition (subplots)
% X = frequency (0-0.5 Hz), Y = window index (time x subject),
% Z = log10 power. Each condition is one subplot.

fig3C = figure('Color','w','Units','inches','Position',[2 1 20 5]);

freq_plot_idx = freq_axis <= 0.5;
freq_plot     = freq_axis(freq_plot_idx);
nFreqPts      = sum(freq_plot_idx);

for cIdx = 1:nHRVConds
    ax = subplot(1, nHRVConds, cIdx);

    psd_surf = psd_surfaces{cIdx};
    if isempty(psd_surf)
        title(hrv_labels{cIdx}); continue;
    end

    % Trim to 0-0.5 Hz before any further processing
    psd_surf = psd_surf(1:min(nFreqPts, size(psd_surf,1)), :);
    freq_use = freq_plot(1:size(psd_surf,1));

    % Log10 transform
    psd_log = log10(psd_surf + eps);

    nWin  = size(psd_log, 2);
    nFreq = size(psd_log, 1);

    % meshgrid: X = frequency (nFreq points), Y = window index (nWin points)
    % surf(X,Y,Z) expects X,Y,Z all [nWin x nFreq]
    [XX, YY] = meshgrid(freq_use, 1:nWin);
    ZZ       = psd_log';   % transpose: [nWin x nFreq]

    surf(XX, YY, ZZ, 'EdgeColor','none', 'FaceAlpha', 0.95);

    hold on;
    % LF/HF band boundary lines at top of surface
    z_top = max(ZZ(:));
    plot3([lf_band(1) lf_band(1)], [1 nWin], [z_top z_top], ...
        '--w', 'LineWidth', 1.5);
    plot3([lf_band(2) lf_band(2)], [1 nWin], [z_top z_top], ...
        '--w', 'LineWidth', 1.5);
    plot3([hf_band(2) hf_band(2)], [1 nWin], [z_top z_top], ...
        '--w', 'LineWidth', 1.5);

    % Band labels
    text(mean(lf_band),  nWin*1.05, z_top, 'LF',  'Color','w', ...
        'FontSize',8, 'FontWeight','bold', 'HorizontalAlignment','center');
    text(mean(hf_band),  nWin*1.05, z_top, 'HF',  'Color','w', ...
        'FontSize',8, 'FontWeight','bold', 'HorizontalAlignment','center');

    colormap(ax, parula);
    shading interp;
    view([225, 30]);   % azimuth 225 = looking from front-left, elevation 30
    xlim([0 0.5]);
    ylim([1 nWin]);

    xlabel('Freq (Hz)',    'FontSize', 8, 'FontWeight','bold');
    ylabel('Window',       'FontSize', 8, 'FontWeight','bold');
    zlabel('log_{10} Pwr', 'FontSize', 8, 'FontWeight','bold');
    title(hrv_labels{cIdx}, 'FontSize', 11, 'FontWeight','bold');

    % Tighten tick marks
    xticks([0 0.04 0.15 0.40 0.5]);
    xticklabels({'0','0.04','0.15','0.40','0.5'});
    xtickangle(ax, 30);

    set(ax, 'LineWidth', 1.0, 'FontSize', 7);
    grid on; box off;
end

sgtitle('IBI Power Spectral Density by Condition (Welch, pooled windows)', ...
    'FontSize', 13, 'FontWeight','bold');

if saving_figures
    print(fig3C, '-dpng', '-r400', sprintf('%s/figure_HRV_PSD_3D.png', figuredir));
end


%% ============================
%   FIGURE 4: Histograms (A) + ACD Violin (B)
% ============================
gridLabels = categories(all_trial_pooled_full.Condition);
nConds     = numel(gridLabels);

figHist = figure('Color','w','Units','inches','Position',[1 1 6.5 5]);

hist_sub_w   = (usable_w - (nConds-1)*gap) / nConds;
hist_sub_x   = left_margin + (0:nConds-1) .* (hist_sub_w + gap);
hist_row_bot = 0.57;
hist_row_h   = 0.35;

for i = 1:nConds
    ax = axes('Units','normalized', ...
              'Position',[hist_sub_x(i), hist_row_bot, hist_sub_w, hist_row_h]);
    hold on; grid on;

    vals = all_trial_pooled_full.ACD(all_trial_pooled_full.Condition == gridLabels{i});
    histogram(vals, 40, 'FaceColor',condColors(i,:), ...
        'EdgeColor',condColors(i,:)*0.6, 'FaceAlpha',0.8, 'LineWidth',0.8);

    mu  = mean(vals);
    sig = std(vals);
    xl  = linspace(min(vals), max(vals), 200);
    yl  = numel(vals) * (xl(2)-xl(1)) * normpdf(xl, mu, sig);
    plot(xl, yl, 'k-', 'LineWidth',1.8);
    xline(mu, '--r', 'LineWidth',1.5);

    text(0.97, 0.97, sprintf('M = %.1f\nSD = %.1f\nn = %d', mu, sig, numel(vals)), ...
        'Units','normalized', 'HorizontalAlignment','right', ...
        'VerticalAlignment','top', 'FontSize',9);

    title(gridLabels{i}, 'FontSize',12);
    xlabel('ACD (ms)', 'FontSize',10);
    if i == 1, ylabel('Count', 'FontSize',10); end
    set(ax, 'Box','off', 'TickDir','out', 'LineWidth', 1.5);

    if i == 1
        annotation('textbox', ...
            [ax.Position(1)-0.04, ax.Position(2)+ax.Position(4)+0.02, 0.04, 0.04], ...
            'String','A','LineStyle','none','FontSize',12,'FontWeight','bold');
    end
end

lme_cond   = fitlme(all_trial_pooled_full, 'ACD ~ ConditionVal + (1|SubjectID)', 'FitMethod','ML');
anova_cond = anova(lme_cond, 'DFMethod', 'Satterthwaite');
r2_cond    = lme_cond.Rsquared.Ordinary;
pval_cond  = anova_cond.pValue(2);
beta_cond  = lme_cond.Coefficients.Estimate(2);

fprintf('\n=== ANALYSIS: ACD ~ ConditionVal + (1|SubjectID) ===\n');
disp(lme_cond.Coefficients);
fprintf('ANOVA (Type 3):\n'); disp(anova_cond);
fprintf('R2=%.4f | AIC=%.1f | BIC=%.1f | LogLik=%.1f\n', ...
    r2_cond, lme_cond.ModelCriterion.AIC, ...
    lme_cond.ModelCriterion.BIC, lme_cond.LogLikelihood);

axB4 = axes('Units','normalized', ...
            'Position',[left_margin, row_bot_bot, usable_w, row_bot_h]);
hold on; grid on;

for i = 1:nConds
    vals = all_trial_pooled_full.ACD(all_trial_pooled_full.Condition == gridLabels{i});
    if ~isempty(vals)
        [f, xi] = ksdensity(vals, 'BoundaryCorrection','reflection');
        f = f / max(f) * 0.35;
        fill(i + [f, -fliplr(f)], [xi, fliplr(xi)], condColors(i,:), ...
            'FaceAlpha',0.7,'EdgeColor','k','LineWidth',1);
        scatter(i + (rand(size(vals))-0.5)*0.1, vals, 8, 'k', ...
            'filled', 'MarkerFaceAlpha',0.2);
        scatter(i, mean(vals), 60, 'w', 'filled', ...
            'MarkerEdgeColor',[0.8 0 0],'LineWidth',2);
    end
end

yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2);

if pval_cond < 0.001, pStr = 'p_{Cond} < 0.001';
else,                  pStr = sprintf('p_{Cond} = %.3f', pval_cond); end

text(0.02, 0.98, sprintf('R^2 = %.3f\n\\beta = %.3f\n%s', r2_cond, beta_cond, pStr), ...
    'Units','normalized', 'FontSize',11, 'VerticalAlignment','top');

xticks(1:nConds); xticklabels(gridLabels);
xlabel('Experimental Condition', 'FontSize',12, 'FontWeight','bold');
ylabel('ACD (ms)', 'FontSize',12, 'FontWeight','bold');
xlim([0.4, nConds + 0.6]);
set(axB4, 'Box','off', 'TickDir','out', 'LineWidth', 1.5);

annotation('textbox', ...
    [axB4.Position(1)-0.04, axB4.Position(2)+axB4.Position(4)+0.02, 0.04, 0.04], ...
    'String','B','LineStyle','none','FontSize',12,'FontWeight','bold');

if saving_figures == 1
    print(figHist, '-dpng', '-r400', sprintf('%s/figure_Hist_ACDViolin.png', figuredir));
end


%% ============================
%   FIGURE 5: Per-subject ACD trajectory
%   Same ΔIBI convention and z±3 per beat per subject as main loop
% ============================

fig5condColors = [
    0.85  0.20  0.20;
    0.95  0.60  0.10;
    0.20  0.65  0.25;
    0.10  0.45  0.85;
    0.50  0.15  0.70;
];

beatPositions = [-3, -2, -1];
xStim         = -0.5;

fig5Keys   = {'cond_0_6667','cond_0_7692','cond_1_0000','cond_1_4286','cond_2_0000'};
fig5Labels = {'66.67%','76.92%','100%','142.86%','200%'};
nFig5Conds = numel(fig5Keys);

condSubMeans = cell(nFig5Conds, 1);

for cIdx = 1:nFig5Conds
    cName = fig5Keys{cIdx};
    if ~isfield(organizedData, cName), continue; end
    condData  = organizedData.(cName);
    nSubsHere = length(condData.subjects);
    subMat    = nan(nSubsHere, 3);

    for sIdx = 1:nSubsHere
        try
            rpeak_row  = squeeze(condData.RPeak(sIdx, :));
            timing_row = squeeze(condData.TrialTiming(sIdx, :));

            peak_positions = find(rpeak_row == 1);
            stim_positions = find(timing_row > 0);
            ibi_ms         = diff(peak_positions) * (1000 / FS);

            nTrials       = length(stim_positions);
            delta_by_beat = nan(nTrials, 3);

            for tIdx = 1:nTrials
                stim_samp      = stim_positions(tIdx);
                pre_stim_peaks = find(peak_positions < stim_samp);

                if numel(pre_stim_peaks) < 4, continue; end

                k_stim  = pre_stim_peaks(end);
                k_stim1 = pre_stim_peaks(end-1);
                k_stim2 = pre_stim_peaks(end-2);

                if k_stim < 2 || k_stim1 < 2 || k_stim2 < 2, continue; end
                if k_stim  > length(ibi_ms), continue; end
                if k_stim1 > length(ibi_ms), continue; end
                if k_stim2 > length(ibi_ms), continue; end

                ibi_before_stim  = ibi_ms(k_stim  - 1);
                ibi_after_stim   = ibi_ms(k_stim);
                ibi_before_stim1 = ibi_ms(k_stim1 - 1);
                ibi_after_stim1  = ibi_ms(k_stim1);
                ibi_before_stim2 = ibi_ms(k_stim2 - 1);
                ibi_after_stim2  = ibi_ms(k_stim2);

                all_ibis = [ibi_before_stim,  ibi_after_stim, ...
                            ibi_before_stim1, ibi_after_stim1, ...
                            ibi_before_stim2, ibi_after_stim2];
                if any(all_ibis < 444) || any(all_ibis > 1333), continue; end

                delta_by_beat(tIdx, 1) = ibi_after_stim2  - ibi_before_stim2;
                delta_by_beat(tIdx, 2) = ibi_after_stim1  - ibi_before_stim1;
                delta_by_beat(tIdx, 3) = ibi_after_stim   - ibi_before_stim;
            end

            % Absolute ΔIBI threshold: ±3 * SD(IBI) for this subject/condition
            valid_ibis_f5 = ibi_ms(ibi_ms >= 444 & ibi_ms <= 1333);
            if numel(valid_ibis_f5) >= 2
                ibi_sd_thresh_f5 = 3 * std(valid_ibis_f5);
                delta_by_beat(abs(delta_by_beat) > ibi_sd_thresh_f5) = NaN;
            end

            % z ± 3 per beat column per subject
            for b = 1:3
                col   = delta_by_beat(:, b);
                valid = ~isnan(col);
                if sum(valid) < 2, continue; end
                mu_b  = mean(col(valid));
                sd_b  = std(col(valid));
                if sd_b == 0, continue; end
                z_b   = (col - mu_b) / sd_b;
                col(abs(z_b) > 3) = NaN;
                delta_by_beat(:, b) = col;
            end

            subMat(sIdx, :) = nanmean(delta_by_beat, 1);

        catch
            continue;
        end
    end

    condSubMeans{cIdx} = subMat;
end

fig5 = figure('Color','w','Units','inches','Position',[1 1 7 5]);
ax5  = axes('Units','normalized','Position',[0.12 0.13 0.72 0.78]);
hold on; grid on;

for cIdx = 1:nFig5Conds
    subMat = condSubMeans{cIdx};
    if isempty(subMat), continue; end
    col = fig5condColors(cIdx, :);

    for sIdx = 1:size(subMat, 1)
        row = subMat(sIdx, :);
        if all(isnan(row)), continue; end
        plot(beatPositions, row, '-', ...
            'Color', [col, 0.30], 'LineWidth', 0.9, 'HandleVisibility','off');
    end

    condMean = nanmean(subMat, 1);
    plot(beatPositions, condMean, '-o', ...
        'Color', col, 'LineWidth', 2.8, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 8, 'DisplayName', fig5Labels{cIdx});

    p         = polyfit(beatPositions, condMean, 1);
    condSlope = p(1);
    yOffset   = (cIdx - 3) * 2.5;
    text(-1 + 0.06, condMean(3) + yOffset, ...
        sprintf('\\beta = %.2f', condSlope), ...
        'Color', col, 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'Clipping', 'on');
end

xline(xStim, '--k', 'LineWidth', 1.6, 'HandleVisibility','off');
text(xStim + 0.04, min(ylim), 'Stim', ...
    'FontSize', 9, 'Color', [0.2 0.2 0.2], 'VerticalAlignment','bottom');
yline(0, ':', 'Color',[0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility','off');

xticks(beatPositions);
xticklabels({'-3','-2','-1'});
xlim([-3.4, 0.0]);
xlabel('Beat position relative to stimulus onset', 'FontSize',12, 'FontWeight','bold');
ylabel('\DeltaIBI (ms)', 'FontSize',12, 'FontWeight','bold');
title('Anticipatory cardiac deceleration by condition', 'FontSize',13, 'FontWeight','bold');

set(ax5, 'Box','off', 'TickDir','out', 'LineWidth', 1.3);
ax5.XAxis.LineWidth = 1.3;
ax5.YAxis.LineWidth = 1.3;

lgd = legend('Location','eastoutside', 'FontSize',10, 'Box','off');
lgd.Title.String = 'Condition';

if saving_figures
    print(fig5, '-dpng', '-r400', ...
        sprintf('%s/figure_ACD_perSubject_byCondition.png', figuredir));
    fprintf('Saved figure_ACD_perSubject_byCondition.png\n');
end

fprintf('\nAnalysis Complete.\n');
