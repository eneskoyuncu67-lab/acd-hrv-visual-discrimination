%% ============================
%   ACD EXCLUSION BREAKDOWN
%   Reruns the exact processing logic from consecutive_ibi_v3.m
%   with counters at every exclusion step, per condition.
%   Requires: organizedData in workspace.
% ============================

if ~exist('organizedData','var')
    error('Run consecutive_ibi_v3.m first to load organizedData.');
end

targetKeys = {'cond_0_6667','cond_0_7692','cond_1_0000','cond_1_4286','cond_2_0000'};
gridLabels = {'66.67%','76.92%','100%','142.86%','200%'};

existingFields     = fieldnames(organizedData);
validIdx           = ismember(targetKeys, existingFields);
gridConditionNames = targetKeys(validIdx);
gridLabels         = gridLabels(validIdx);

FS = 100;

% Per-condition counters
nConds = numel(gridLabels);
cnt_total          = zeros(1,nConds); % all stimulus onsets
cnt_few_peaks      = zeros(1,nConds); % < 4 pre-stim peaks
cnt_ibi_plaus      = zeros(1,nConds); % IBI outside 444-1333 ms
cnt_abs_thresh     = zeros(1,nConds); % beats removed by abs ΔIBI threshold
cnt_z3_beats       = zeros(1,nConds); % beats removed by z±3 per beat column
cnt_all_nan        = zeros(1,nConds); % trials lost because all beats NaN
cnt_missing_meta   = zeros(1,nConds); % trials lost because RT/Acc missing
cnt_retained       = zeros(1,nConds); % final retained trials

for cIdx = 1:nConds
    cName    = gridConditionNames{cIdx};
    if ~isfield(organizedData, cName), continue; end
    condData = organizedData.(cName);

    for sIdx = 1:length(condData.subjects)
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

            nTrials    = length(stim_positions);
            raw_deltas = nan(nTrials, 3);
            trial_meta = nan(nTrials, 3);

            cnt_total(cIdx) = cnt_total(cIdx) + nTrials;

            %% PASS 1: per-trial exclusions
            for tIdx = 1:nTrials
                stim_samp      = stim_positions(tIdx);
                pre_stim_peaks = find(peak_positions < stim_samp);

                % Criterion 1: fewer than 4 pre-stimulus peaks
                if numel(pre_stim_peaks) < 4
                    cnt_few_peaks(cIdx) = cnt_few_peaks(cIdx) + 1;
                    continue;
                end

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

                all_ibis = [ibi_before_stim, ibi_after_stim, ...
                            ibi_before_stim1, ibi_after_stim1, ...
                            ibi_before_stim2, ibi_after_stim2];

                % Criterion 2: IBI outside physiological range
                if any(all_ibis < 444) || any(all_ibis > 1333)
                    cnt_ibi_plaus(cIdx) = cnt_ibi_plaus(cIdx) + 1;
                    continue;
                end

                raw_deltas(tIdx, 1) = ibi_after_stim2  - ibi_before_stim2;
                raw_deltas(tIdx, 2) = ibi_after_stim1  - ibi_before_stim1;
                raw_deltas(tIdx, 3) = ibi_after_stim   - ibi_before_stim;

                if tIdx <= numel(rt_vec) && tIdx <= numel(acc_vec)
                    trial_meta(tIdx, :) = [tIdx, rt_vec(tIdx), acc_vec(tIdx)];
                end
            end

            %% PASS 1b: absolute ΔIBI threshold
            valid_ibis = ibi_ms(ibi_ms >= 444 & ibi_ms <= 1333);
            if numel(valid_ibis) >= 2
                ibi_sd_thresh = 3 * std(valid_ibis);
                before_abs    = sum(~isnan(raw_deltas(:)));
                raw_deltas(abs(raw_deltas) > ibi_sd_thresh) = NaN;
                after_abs     = sum(~isnan(raw_deltas(:)));
                cnt_abs_thresh(cIdx) = cnt_abs_thresh(cIdx) + (before_abs - after_abs);
            end

            %% PASS 2: z±3 per beat column per subject
            for b = 1:3
                col   = raw_deltas(:, b);
                valid = ~isnan(col);
                if sum(valid) < 2, continue; end
                mu_b = mean(col(valid));
                sd_b = std(col(valid));
                if sd_b == 0, continue; end
                z_b  = (col - mu_b) / sd_b;
                n_before = sum(valid);
                col(abs(z_b) > 3) = NaN;
                n_after  = sum(~isnan(col));
                cnt_z3_beats(cIdx) = cnt_z3_beats(cIdx) + (n_before - n_after);
                raw_deltas(:, b) = col;
            end

            %% PASS 3: trial-level ACD
            for tIdx = 1:nTrials
                row = raw_deltas(tIdx, :);

                % Criterion 3: all beats NaN
                if all(isnan(row))
                    cnt_all_nan(cIdx) = cnt_all_nan(cIdx) + 1;
                    continue;
                end

                % Criterion 4: missing RT or accuracy
                if any(isnan(trial_meta(tIdx, :)))
                    cnt_missing_meta(cIdx) = cnt_missing_meta(cIdx) + 1;
                    continue;
                end

                cnt_retained(cIdx) = cnt_retained(cIdx) + 1;
            end

        catch
            continue;
        end
    end
end

%% Print results
fprintf('\n');
fprintf('=====================================================================================================\n');
fprintf('  ACD EXCLUSION BREAKDOWN PER CONDITION\n');
fprintf('=====================================================================================================\n');
fprintf('%-12s %8s %10s %10s %12s %12s %10s %10s %9s\n', ...
    'Condition','Total','<4 peaks','IBI range', ...
    'Abs thresh*','z±3 beats*','All NaN','No meta','Retained');
fprintf('%s\n', repmat('-',1,101));

for i = 1:nConds
    fprintf('%-12s %8d %10d %10d %12d %12d %10d %10d %9d\n', ...
        gridLabels{i}, ...
        cnt_total(i), ...
        cnt_few_peaks(i), ...
        cnt_ibi_plaus(i), ...
        cnt_abs_thresh(i), ...
        cnt_z3_beats(i), ...
        cnt_all_nan(i), ...
        cnt_missing_meta(i), ...
        cnt_retained(i));
end

fprintf('%s\n', repmat('-',1,101));
fprintf('%-12s %8d %10d %10d %12d %12d %10d %10d %9d\n', ...
    'TOTAL', ...
    sum(cnt_total), sum(cnt_few_peaks), sum(cnt_ibi_plaus), ...
    sum(cnt_abs_thresh), sum(cnt_z3_beats), ...
    sum(cnt_all_nan), sum(cnt_missing_meta), sum(cnt_retained));
fprintf('=====================================================================================================\n');
fprintf('* Beat-level counts (not trial-level): one trial has 3 beats, so these can exceed trial counts.\n');
fprintf('\n');

%% Summary percentages
fprintf('  Retention summary:\n');
fprintf('  %-30s %6.1f%%\n','Excluded (<4 pre-stim peaks):', ...
    100*sum(cnt_few_peaks)/sum(cnt_total));
fprintf('  %-30s %6.1f%%\n','Excluded (IBI outside range):', ...
    100*sum(cnt_ibi_plaus)/sum(cnt_total));
fprintf('  %-30s %6.1f%%\n','Excluded (all beats NaN):', ...
    100*sum(cnt_all_nan)/sum(cnt_total));
fprintf('  %-30s %6.1f%%\n','Excluded (missing metadata):', ...
    100*sum(cnt_missing_meta)/sum(cnt_total));
fprintf('  %-30s %6.1f%%\n','Final retained:', ...
    100*sum(cnt_retained)/sum(cnt_total));
fprintf('\n');
