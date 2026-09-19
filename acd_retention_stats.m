%% ACD retention stats per condition
% Requires: all_trial_pooled_full in workspace
% Shows: trials before and after each exclusion step, per condition

if ~exist('all_trial_pooled_full','var')
    error('Run consecutive_ibi_v3.m first.');
end

gridLabels  = categories(all_trial_pooled_full.Condition);
nConds      = numel(gridLabels);

fprintf('\n');
fprintf('===========================================================================\n');
fprintf('  ACD RETENTION AFTER EXCLUSION CRITERIA\n');
fprintf('===========================================================================\n');
fprintf('%-12s %10s %10s %10s %10s %10s\n', ...
    'Condition','After main','After z±3','Removed','Retained%','Mean ACD');
fprintf('%s\n', repmat('-',1,62));

total_main = 0;
total_filt = 0;

for i = 1:nConds
    cond = gridLabels{i};

    % Step 1: after main loop (all_trial_pooled_full)
    rows_main = all_trial_pooled_full.Condition == cond;
    vals_main = all_trial_pooled_full.ACD(rows_main);
    n_main    = sum(~isnan(vals_main));

    % Step 2: apply z±3 per condition (same as presentation_figures.m)
    mu_c = mean(vals_main, 'omitnan');
    sd_c = std(vals_main,  'omitnan');
    if sd_c > 0
        z_c      = (vals_main - mu_c) / sd_c;
        retained = vals_main(abs(z_c) <= 3 & ~isnan(vals_main));
    else
        retained = vals_main(~isnan(vals_main));
    end

    n_filt   = numel(retained);
    n_removed = n_main - n_filt;
    pct      = 100 * n_filt / n_main;
    mu_acd   = mean(retained);

    fprintf('%-12s %10d %10d %10d %9.1f%% %10.2f\n', ...
        cond, n_main, n_filt, n_removed, pct, mu_acd);

    total_main = total_main + n_main;
    total_filt = total_filt + n_filt;
end

fprintf('%s\n', repmat('-',1,62));
fprintf('%-12s %10d %10d %10d %9.1f%%\n', ...
    'TOTAL', total_main, total_filt, ...
    total_main - total_filt, 100*total_filt/total_main);
fprintf('===========================================================================\n\n');
