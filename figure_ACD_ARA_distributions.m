%% ============================
%   FIGURE_ACD_ARA_DISTRIBUTIONS.M  --  STANDALONE
%
%   Histogram distributions of ACD and ARA per condition.
%   ACD = HR(stim-3) - HR(stim)   [bpm]
%   ARA = HR(stim+3) - HR(stim+1) [bpm]
%
%   Reuses the identical trial-extraction and 7-criterion exclusion
%   pipeline as acd_ara_analysis.m (Steps 1-3, Criterion 7), so the
%   trials shown here are exactly the ones the LME/correlation models
%   are fit to. See README.md for the ACD/ARA definitions and exclusion
%   criteria.
%
%   Requires: organizedData in workspace (or organizedData_by_condition.mat)
% ============================

%% ---- Load organizedData ----
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
        error('organizedData not found in workspace or on disk.');
    end
end

%% ---- Settings ----
FS           = 100;
stim_dur_smp = round(16.7 / (1000/FS));  % ~2 samples

targetKeys  = {'cond_0_6667','cond_0_7692','cond_1_0000','cond_1_4286','cond_2_0000'};
gridLabels  = {'66.67%','76.92%','100%','142.86%','200%'};
condNumVals = [66.67, 76.92, 100, 142.86, 200];

existingFields     = fieldnames(organizedData);
validIdx           = ismember(targetKeys, existingFields);
gridConditionNames = targetKeys(validIdx);
gridLabels         = gridLabels(validIdx);
condNumVals        = condNumVals(validIdx);
nConds             = numel(gridLabels);

condColors = [0.7 0.85 0.9; 0.5 0.7 0.85; 0.3 0.55 0.8; 0.1 0.35 0.7; 0.05 0.15 0.5];
condNumMap = containers.Map(gridLabels, condNumVals);

%% ============================
%   STEP 1: COLLECT PAIRED ACD + ARA PER TRIAL
%   Only trials with valid values for BOTH measures are retained.
%   (Identical logic to acd_ara_analysis.m.)
% ============================
fprintf('\n=== DATA PROCESSING ===\n');
fprintf('Computing paired ACD and ARA per trial.\n');
fprintf('ACD = HR(stim-3) - HR(stim)  [bpm]\n');
fprintf('ARA = HR(stim+3) - HR(stim+1)  [bpm]\n\n');

all_data_cell = {};

for cIdx = 1:nConds
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
            ibi_ms         = diff(peak_positions) * (1000 / FS);

            nTrials        = length(stim_positions);

            raw_d_acd  = nan(nTrials, 3);
            ibi_sroufe = nan(nTrials, 2);

            raw_d_ara  = nan(nTrials, 3);
            ibi_ara    = nan(nTrials, 2);

            trial_meta = nan(nTrials, 3);

            for tIdx = 1:nTrials
                stim_samp = stim_positions(tIdx);
                stim_off  = stim_samp + stim_dur_smp;

                pre_peaks = find(peak_positions < stim_samp);
                if numel(pre_peaks) >= 4
                    k0  = pre_peaks(end);
                    k1p = pre_peaks(end-1);
                    k2p = pre_peaks(end-2);

                    if k0 >= 2 && k1p >= 2 && k2p >= 2 && ...
                       k0 <= length(ibi_ms) && k1p <= length(ibi_ms) && k2p <= length(ibi_ms)

                        ib_s   = ibi_ms(k0);    ib_sb  = ibi_ms(k0-1);
                        ib_s1  = ibi_ms(k1p);   ib_s1b = ibi_ms(k1p-1);
                        ib_s2  = ibi_ms(k2p);   ib_s2b = ibi_ms(k2p-1);

                        all_pre = [ib_sb,ib_s,ib_s1b,ib_s1,ib_s2b,ib_s2];
                        if ~(any(all_pre < 444) || any(all_pre > 1333))
                            raw_d_acd(tIdx,1) = ib_s2  - ib_s2b;
                            raw_d_acd(tIdx,2) = ib_s1  - ib_s1b;
                            raw_d_acd(tIdx,3) = ib_s    - ib_sb;
                            ibi_sroufe(tIdx,1) = ib_s2;
                            ibi_sroufe(tIdx,2) = ib_s;
                        end
                    end
                end

                post_peaks = find(peak_positions > stim_off);
                if numel(post_peaks) >= 3
                    kp1 = post_peaks(1);
                    kp2 = post_peaks(2);
                    kp3 = post_peaks(3);

                    if kp1 >= 2 && kp2 >= 2 && kp3 >= 2 && ...
                       kp1 <= length(ibi_ms) && kp2 <= length(ibi_ms) && kp3 <= length(ibi_ms)

                        ib1=ibi_ms(kp1); ib0=ibi_ms(kp1-1);
                        ib2=ibi_ms(kp2); ib_b2=ibi_ms(kp2-1);
                        ib3=ibi_ms(kp3); ib_b3=ibi_ms(kp3-1);

                        all_post = [ib0,ib1,ib_b2,ib2,ib_b3,ib3];
                        if ~(any(all_post < 444) || any(all_post > 1333))
                            raw_d_ara(tIdx,1) = ib1 - ib0;
                            raw_d_ara(tIdx,2) = ib2 - ib_b2;
                            raw_d_ara(tIdx,3) = ib3 - ib_b3;
                            ibi_ara(tIdx,1)   = ib1;
                            ibi_ara(tIdx,2)   = ib3;
                        end
                    end
                end

                if tIdx <= numel(rt_vec) && tIdx <= numel(acc_vec)
                    trial_meta(tIdx,:) = [tIdx, rt_vec(tIdx), acc_vec(tIdx)];
                end
            end

            valid_ibis = ibi_ms(ibi_ms >= 444 & ibi_ms <= 1333);
            if numel(valid_ibis) >= 2
                thresh = 3 * std(valid_ibis);
                raw_d_acd(abs(raw_d_acd) > thresh) = NaN;
                raw_d_ara(abs(raw_d_ara) > thresh) = NaN;
            end

            for b = 1:3
                col = raw_d_acd(:,b); valid=~isnan(col);
                if sum(valid)>=2
                    mu_b=mean(col(valid)); sd_b=std(col(valid));
                    if sd_b>0, col(abs((col-mu_b)/sd_b)>3)=NaN; raw_d_acd(:,b)=col; end
                end
                col = raw_d_ara(:,b); valid=~isnan(col);
                if sum(valid)>=2
                    mu_b=mean(col(valid)); sd_b=std(col(valid));
                    if sd_b>0, col(abs((col-mu_b)/sd_b)>3)=NaN; raw_d_ara(:,b)=col; end
                end
            end

            for tIdx = 1:nTrials
                if all(isnan(raw_d_acd(tIdx,:))), continue; end
                if all(isnan(raw_d_ara(tIdx,:))), continue; end
                if any(isnan(trial_meta(tIdx,:))), continue; end

                is3 = ibi_sroufe(tIdx,1); is0 = ibi_sroufe(tIdx,2);
                ip1 = ibi_ara(tIdx,1);    ip3 = ibi_ara(tIdx,2);
                if isnan(is3)||isnan(is0)||isnan(ip1)||isnan(ip3), continue; end

                acd_val = (60000/is3) - (60000/is0);
                ara_val = (60000/ip3) - (60000/ip1);

                all_data_cell = [all_data_cell; ...
                    {cLabel, char(sname), tIdx, acd_val, ara_val, ...
                     trial_meta(tIdx,2), trial_meta(tIdx,3)}];
            end

        catch
            continue;
        end
    end
end

tbl = cell2table(all_data_cell, ...
    'VariableNames',{'Condition','SubjectID','TrialIndex','ACD','ARA','RT','Acc'});
tbl.Condition    = categorical(tbl.Condition, gridLabels);
tbl.Acc          = categorical(tbl.Acc, [0,1], {'Miss','Hit'});
tbl              = rmmissing(tbl);
tbl.ConditionVal = cellfun(@(x) condNumMap(x), cellstr(tbl.Condition));

fprintf('After criteria 1-6 (paired ACD+ARA): %d trials retained\n', height(tbl));

for i = 1:nConds
    rows = tbl.Condition == gridLabels{i};
    vals = tbl.ACD(rows);
    mu_c=mean(vals,'omitnan'); sd_c=std(vals,'omitnan');
    if sd_c==0, continue; end
    vals(abs((vals-mu_c)/sd_c)>3) = NaN;
    tbl.ACD(rows) = vals;
end
tbl = tbl(~isnan(tbl.ACD),:);

for i = 1:nConds
    rows = tbl.Condition == gridLabels{i};
    vals = tbl.ARA(rows);
    mu_c=mean(vals,'omitnan'); sd_c=std(vals,'omitnan');
    if sd_c==0, continue; end
    vals(abs((vals-mu_c)/sd_c)>3) = NaN;
    tbl.ARA(rows) = vals;
end
tbl = tbl(~isnan(tbl.ARA),:);

fprintf('After criterion 7 (z+-3 on ACD and ARA): %d trials retained\n\n', height(tbl));

%% ============================
%   FIGURE: ACD and ARA distributions per condition
% ============================
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 13.3 5.5]);

for i = 1:nConds
    rows = tbl.Condition == gridLabels{i};

    ax1 = subplot(2, nConds, i);
    vals = tbl.ACD(rows); vals = vals(~isnan(vals));
    histogram(vals, 'FaceColor', condColors(i,:), 'FaceAlpha', 0.8, ...
        'EdgeColor', 'k', 'Normalization', 'pdf');
    xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
    title(gridLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
    if i == 1, ylabel('ACD density', 'FontSize', 9, 'FontWeight', 'bold'); end
    text(0.97, 0.95, sprintf('M=%.2f\nSD=%.2f\nn=%d', mean(vals), std(vals), numel(vals)), ...
        'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top');
    set(ax1, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);

    ax2 = subplot(2, nConds, nConds + i);
    vals = tbl.ARA(rows); vals = vals(~isnan(vals));
    histogram(vals, 'FaceColor', condColors(i,:), 'FaceAlpha', 0.8, ...
        'EdgeColor', 'k', 'Normalization', 'pdf');
    xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
    xlabel('bpm', 'FontSize', 8);
    if i == 1, ylabel('ARA density', 'FontSize', 9, 'FontWeight', 'bold'); end
    text(0.97, 0.95, sprintf('M=%.2f\nSD=%.2f\nn=%d', mean(vals), std(vals), numel(vals)), ...
        'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top');
    set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);
end

annotation('textbox', [0.02, 0.96, 0.9, 0.04], ...
    'String', 'ACD (top) and ARA (bottom) distributions by condition', ...
    'EdgeColor', 'none', 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(fig, './figure_ACD_ARA_distributions.png', 'Resolution', 300);
fprintf('Saved figure_ACD_ARA_distributions.png\n');
