%% ============================
%   ACD_ARA_ANALYSIS.M  —  STANDALONE
%
%   PART 1: LME model
%     ARA ~ ConditionVal + (1|SubjectID)
%     ARA = HR(stim+3) - HR(stim+1)  [bpm]
%
%   PART 2: ACD–ARA correlation per condition
%     For each condition separately:
%       Compute per-trial ACD and ARA for trials that pass ALL exclusion
%       criteria for BOTH measures.
%       Pool across subjects within condition.
%       Pearson correlation between ACD and ARA.
%       One correlation coefficient per condition.
%
%   All exclusion criteria applied identically to both ACD and ARA.
%   Only trials with valid ACD AND valid ARA are included in correlations.
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

            % --- ACD exclusion screening (pre-stimulus) ---
            raw_d_acd  = nan(nTrials, 3);
            ibi_sroufe = nan(nTrials, 2);

            % --- ARA exclusion screening (post-stimulus) ---
            raw_d_ara  = nan(nTrials, 3);
            ibi_ara    = nan(nTrials, 2);

            trial_meta = nan(nTrials, 3);

            %% PASS 1: per-trial collection for both ACD and ARA
            for tIdx = 1:nTrials
                stim_samp = stim_positions(tIdx);
                stim_off  = stim_samp + stim_dur_smp;

                %% ACD: pre-stimulus peaks
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
                            ibi_sroufe(tIdx,1) = ib_s2;   % IBI at beat -3
                            ibi_sroufe(tIdx,2) = ib_s;    % IBI at beat -1 (stim)
                        end
                    end
                end

                %% ARA: post-stimulus peaks
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
                            ibi_ara(tIdx,1)   = ib1;  % IBI at beat +1
                            ibi_ara(tIdx,2)   = ib3;  % IBI at beat +3
                        end
                    end
                end

                if tIdx <= numel(rt_vec) && tIdx <= numel(acc_vec)
                    trial_meta(tIdx,:) = [tIdx, rt_vec(tIdx), acc_vec(tIdx)];
                end
            end

            %% PASS 1b: Criterion 3 — abs ΔIBI threshold
            valid_ibis    = ibi_ms(ibi_ms >= 444 & ibi_ms <= 1333);
            if numel(valid_ibis) >= 2
                thresh = 3 * std(valid_ibis);
                raw_d_acd(abs(raw_d_acd) > thresh) = NaN;
                raw_d_ara(abs(raw_d_ara) > thresh) = NaN;
            end

            %% PASS 2: Criterion 4 — z±3 per beat column
            for b = 1:3
                for raw_d = {raw_d_acd, raw_d_ara}
                    col   = raw_d{1}(:,b);
                    valid = ~isnan(col);
                    if sum(valid) < 2, continue; end
                    mu_b = mean(col(valid)); sd_b = std(col(valid));
                    if sd_b == 0, continue; end
                    col(abs((col-mu_b)/sd_b) > 3) = NaN;
                    if isequal(raw_d{1}, raw_d_acd)
                        raw_d_acd(:,b) = col;
                    else
                        raw_d_ara(:,b) = col;
                    end
                end
            end
            % Re-apply separately (loop above is read-only on raw_d)
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

            %% PASS 3: Criteria 5 & 6 — collect valid paired trials
            for tIdx = 1:nTrials
                % Criteria 5: both ACD and ARA beats must not be all-NaN
                if all(isnan(raw_d_acd(tIdx,:))), continue; end
                if all(isnan(raw_d_ara(tIdx,:))), continue; end
                % Criterion 6: metadata must be present
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

%% Build paired table
tbl = cell2table(all_data_cell, ...
    'VariableNames',{'Condition','SubjectID','TrialIndex','ACD','ARA','RT','Acc'});
tbl.Condition    = categorical(tbl.Condition, gridLabels);
tbl.Acc          = categorical(tbl.Acc, [0,1], {'Miss','Hit'});
tbl              = rmmissing(tbl);
tbl.ConditionVal = cellfun(@(x) condNumMap(x), cellstr(tbl.Condition));

fprintf('After criteria 1-6 (paired ACD+ARA): %d trials retained\n', height(tbl));

%% Criterion 7: z±3 on ACD per condition
for i = 1:nConds
    rows = tbl.Condition == gridLabels{i};
    vals = tbl.ACD(rows);
    mu_c=mean(vals,'omitnan'); sd_c=std(vals,'omitnan');
    if sd_c==0, continue; end
    vals(abs((vals-mu_c)/sd_c)>3) = NaN;
    tbl.ACD(rows) = vals;
end
tbl = tbl(~isnan(tbl.ACD),:);

%% Criterion 7: z±3 on ARA per condition
for i = 1:nConds
    rows = tbl.Condition == gridLabels{i};
    vals = tbl.ARA(rows);
    mu_c=mean(vals,'omitnan'); sd_c=std(vals,'omitnan');
    if sd_c==0, continue; end
    vals(abs((vals-mu_c)/sd_c)>3) = NaN;
    tbl.ARA(rows) = vals;
end
tbl = tbl(~isnan(tbl.ARA),:);

fprintf('After criterion 7 (z±3 on ACD and ARA): %d trials retained\n\n', height(tbl));

%% RT units
rt_check = tbl.RT(~isnan(tbl.RT)&tbl.RT>0);
if median(rt_check) < 10, tbl.RT = tbl.RT * 1000; end

%% Downsample for accuracy model
balanced_data = table();
rng('default');
for c = 1:nConds
    sub  = tbl(tbl.Condition==gridLabels{c},:);
    mIdx = find(sub.Acc=='Miss'); hIdx = find(sub.Acc=='Hit');
    if length(hIdx)>length(mIdx) && ~isempty(mIdx)
        hIdx = hIdx(randperm(length(hIdx),length(mIdx)));
    end
    balanced_data = [balanced_data; sub(sort([mIdx;hIdx]),:)];
end
fprintf('N paired: %d | N downsampled: %d\n\n', height(tbl), height(balanced_data));

%% ============================
%   PART 1: LME MODELS — ARA
% ============================

%% MODEL 1: ARA ~ ConditionVal + (1|SubjectID)
fprintf('================================================================\n');
fprintf('  MODEL 1: ARA ~ ConditionVal + (1|SubjectID)\n');
fprintf('================================================================\n');
lme1  = fitlme(tbl,'ARA ~ ConditionVal + (1|SubjectID)','FitMethod','ML');
anov1 = anova(lme1,'DFMethod','Satterthwaite');
disp(lme1); fprintf('Type III ANOVA:\n'); disp(anov1);
fprintf('R2=%.4f | AIC=%.2f | BIC=%.2f | LogLik=%.4f | N=%d\n\n', ...
    lme1.Rsquared.Ordinary,lme1.ModelCriterion.AIC, ...
    lme1.ModelCriterion.BIC,lme1.LogLikelihood,lme1.NumObservations);

%% MODEL 2: ARA ~ RT + ConditionVal + (1|SubjectID)
fprintf('================================================================\n');
fprintf('  MODEL 2: ARA ~ RT + ConditionVal + (1|SubjectID)\n');
fprintf('================================================================\n');
lme2  = fitlme(tbl,'ARA ~ RT + ConditionVal + (1|SubjectID)','FitMethod','ML');
anov2 = anova(lme2,'DFMethod','Satterthwaite');
disp(lme2); fprintf('Type III ANOVA:\n'); disp(anov2);
fprintf('R2=%.4f | AIC=%.2f | BIC=%.2f | LogLik=%.4f | N=%d\n\n', ...
    lme2.Rsquared.Ordinary,lme2.ModelCriterion.AIC, ...
    lme2.ModelCriterion.BIC,lme2.LogLikelihood,lme2.NumObservations);

%% MODEL 3: ARA ~ Acc * ConditionVal + (1|SubjectID)  [downsampled]
fprintf('================================================================\n');
fprintf('  MODEL 3: ARA ~ Acc * ConditionVal + (1|SubjectID)\n');
fprintf('================================================================\n');
lme3  = fitlme(balanced_data,'ARA ~ Acc * ConditionVal + (1|SubjectID)','FitMethod','ML');
anov3 = anova(lme3,'DFMethod','Satterthwaite');
disp(lme3); fprintf('Type III ANOVA:\n'); disp(anov3);
fprintf('R2=%.4f | AIC=%.2f | BIC=%.2f | LogLik=%.4f | N=%d\n\n', ...
    lme3.Rsquared.Ordinary,lme3.ModelCriterion.AIC, ...
    lme3.ModelCriterion.BIC,lme3.LogLikelihood,lme3.NumObservations);

%% ============================
%   PART 2: ACD–ARA CORRELATION PER CONDITION
%   Pool all trials within each condition (across subjects).
%   Pearson r with 95% CI and p-value.
% ============================
fprintf('================================================================\n');
fprintf('  PART 2: ACD–ARA Pearson correlation per condition\n');
fprintf('  Trials pooled across subjects within each condition.\n');
fprintf('================================================================\n');
fprintf('%-10s %8s %10s %10s %10s %8s\n','Condition','N','r','p','95%CI_lo','95%CI_hi');
fprintf('%s\n',repmat('-',1,60));

corr_r   = zeros(1,nConds);
corr_p   = zeros(1,nConds);
corr_n   = zeros(1,nConds);
corr_ci  = zeros(nConds,2);

for i = 1:nConds
    rows = tbl.Condition == gridLabels{i};
    acd  = tbl.ACD(rows);
    ara  = tbl.ARA(rows);
    valid = ~isnan(acd) & ~isnan(ara);
    acd  = acd(valid); ara = ara(valid);
    n    = numel(acd);
    corr_n(i) = n;

    if n < 5
        fprintf('%-10s %8d  (insufficient data)\n', gridLabels{i}, n);
        corr_r(i) = NaN; corr_p(i) = NaN; corr_ci(i,:) = [NaN NaN];
        continue;
    end

    [r, p] = corr(acd, ara, 'Type', 'Pearson');
    corr_r(i) = r; corr_p(i) = p;

    % Fisher z 95% CI
    z     = atanh(r);
    se    = 1 / sqrt(n - 3);
    z_ci  = z + [-1 1] * 1.96 * se;
    ci    = tanh(z_ci);
    corr_ci(i,:) = ci;

    if p < 0.001, ps = '< .001'; else, ps = sprintf('= %.3f', p); end
    fprintf('%-10s %8d %10.4f %10s %10.4f %8.4f\n', ...
        gridLabels{i}, n, r, ps, ci(1), ci(2));
end
fprintf('%s\n\n',repmat('-',1,60));

%% ============================
%   FIGURES
% ============================

%% Figure 1 — ARA violin per condition (Model 1)
int_c  = lme1.Coefficients.Estimate(1);
beta_c = lme1.Coefficients.Estimate(2);
pval_c = lme1.Coefficients.pValue(2);
r2_c   = lme1.Rsquared.Ordinary;
if pval_c<0.001, pc_s='< .001'; else, pc_s=sprintf('= %.3f',pval_c); end

fig1  = figure('Color','w','Units','inches','Position',[0.5 0.5 13.3 5.5]);
ax1   = axes('Units','normalized','Position',[0.07, 0.13, 0.89, 0.78]);
hold on; grid on;
cond_mu1=zeros(1,nConds); cond_sd1=zeros(1,nConds); cond_n1=zeros(1,nConds);
for i = 1:nConds
    vals=tbl.ARA(tbl.Condition==gridLabels{i}); vals=vals(~isnan(vals));
    if numel(vals)<3, continue; end
    [f,xi]=ksdensity(vals,'BoundaryCorrection','reflection'); f=f/max(f)*0.35;
    fill(i+[f,-fliplr(f)],[xi,fliplr(xi)],condColors(i,:),'FaceAlpha',0.65,'EdgeColor','k','LineWidth',0.8,'HandleVisibility','off');
    scatter(i+(rand(size(vals))-0.5)*0.08,vals,5,'k','filled','MarkerFaceAlpha',0.12,'HandleVisibility','off');
    mu_i=mean(vals); sd_i=std(vals); n_i=numel(vals);
    cond_mu1(i)=mu_i; cond_sd1(i)=sd_i; cond_n1(i)=n_i;
    scatter(i,mu_i,90,'w','filled','MarkerEdgeColor',[0.15 0.15 0.15],'LineWidth',1.8,'HandleVisibility','off');
end
plot(1:nConds, int_c+beta_c*condNumVals,'--k','LineWidth',2.2,'HandleVisibility','off');
yline(0,'--','Color',[0.45 0.45 0.45],'LineWidth',1.0,'HandleVisibility','off');
yl=ylim; y_lab=yl(1)+0.03*diff(yl);
for i=1:nConds
    text(i,y_lab,sprintf('M=%.1f\nSD=%.1f\nn=%d',cond_mu1(i),cond_sd1(i),cond_n1(i)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7,'Color',[0.15 0.15 0.15]);
end
xticks(1:nConds); xticklabels(gridLabels); xlim([0.4 nConds+0.6]);
xlabel('Experimental Condition','FontSize',11,'FontWeight','bold');
ylabel('ARA (bpm)','FontSize',11,'FontWeight','bold');
title('ARA \sim Condition','FontSize',12,'FontWeight','bold');
text(0.02,0.97,sprintf('\\beta=%.3f,  p %s,  R^2=%.3f',beta_c,pc_s,r2_c),'Units','normalized','FontSize',10,'VerticalAlignment','top','BackgroundColor','w','EdgeColor',[0.7 0.7 0.7]);
set(ax1,'Box','off','TickDir','out','LineWidth',1.2,'FontSize',9);
exportgraphics(fig1,'./acd_ara_figure1_ARA_condition.png','Resolution',300);
fprintf('Saved acd_ara_figure1_ARA_condition.png\n');

%% Figure 2 — ACD vs ARA scatter per condition (5 panels)
fig2 = figure('Color','w','Units','inches','Position',[0.5 0.5 13.3 5.5]);
panel_b=0.13; panel_h=0.75;
left_m=0.07; right_m=0.02; uw=1-left_m-right_m;
gap_c=0.018; pw=(uw-(nConds-1)*gap_c)/nConds;
px=left_m+(0:nConds-1).*(pw+gap_c);

for i = 1:nConds
    ax = axes('Units','normalized','Position',[px(i),panel_b,pw,panel_h]);
    hold on; grid on;

    rows = tbl.Condition==gridLabels{i};
    acd  = tbl.ACD(rows); ara = tbl.ARA(rows);
    valid = ~isnan(acd)&~isnan(ara);
    acd  = acd(valid); ara = ara(valid);

    scatter(acd, ara, 6, condColors(i,:), 'filled', 'MarkerFaceAlpha', 0.25, ...
        'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');

    if numel(acd) >= 5
        p_fit = polyfit(acd, ara, 1);
        x_fit = linspace(min(acd), max(acd), 100);
        plot(x_fit, polyval(p_fit,x_fit), '-', 'Color', condColors(i,:)*0.6, ...
            'LineWidth', 1.8, 'HandleVisibility', 'off');

        r   = corr_r(i); p_v = corr_p(i);
        if isnan(r), rstr = 'r=NaN'; else
            if p_v<0.001, pstr='p<.001'; else, pstr=sprintf('p=%.3f',p_v); end
            rstr = sprintf('r=%.3f\n%s\nn=%d', r, pstr, corr_n(i));
        end
        text(0.05,0.97,rstr,'Units','normalized','FontSize',8,'VerticalAlignment','top', ...
            'BackgroundColor','w','EdgeColor',[0.7 0.7 0.7],'Color',condColors(i,:)*0.5);
    end

    xline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',0.7,'HandleVisibility','off');
    yline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',0.7,'HandleVisibility','off');
    title(gridLabels{i},'FontSize',10,'FontWeight','bold');
    if i==1
        xlabel('ACD (bpm)','FontSize',8.5); ylabel('ARA (bpm)','FontSize',8.5);
    else
        xlabel('ACD (bpm)','FontSize',8.5);
    end
    set(ax,'Box','off','TickDir','out','LineWidth',1.0,'FontSize',8);
end
annotation('textbox',[left_m+0.02,panel_b+panel_h+0.02,0.5,0.05], ...
    'String','ACD–ARA Correlation per Condition (trials pooled across subjects)', ...
    'EdgeColor','none','FontSize',11,'FontWeight','bold');
exportgraphics(fig2,'./acd_ara_figure2_correlation.png','Resolution',300);
fprintf('Saved acd_ara_figure2_correlation.png\n');

%% Figure 3 — r values across conditions bar chart
fig3 = figure('Color','w','Units','inches','Position',[0.5 0.5 7 5]);
ax3  = axes('Units','normalized','Position',[0.13,0.13,0.83,0.72]);
hold on; grid on;

valid_conds = ~isnan(corr_r);
b_h = bar(find(valid_conds), corr_r(valid_conds), 0.55, 'FaceColor','flat');
for i = 1:sum(valid_conds)
    ci = find(valid_conds); b_h.CData(i,:) = condColors(ci(i),:);
end
errorbar(find(valid_conds), corr_r(valid_conds), ...
    corr_r(valid_conds)-corr_ci(valid_conds,1), ...
    corr_ci(valid_conds,2)-corr_r(valid_conds), ...
    'k.','LineWidth',1.2,'CapSize',5,'HandleVisibility','off');
yline(0,'--','Color',[0.4 0.4 0.4],'LineWidth',1.0,'HandleVisibility','off');

for i = find(valid_conds)
    if corr_p(i) < 0.05
        yl = ylim;
        text(i, corr_r(i) + sign(corr_r(i))*0.02*diff(yl), '*', ...
            'HorizontalAlignment','center','FontSize',14,'Color',[0.1 0.1 0.1]);
    end
end

xticks(1:nConds); xticklabels(gridLabels); xlim([0.3 nConds+0.7]);
xlabel('Experimental Condition','FontSize',10,'FontWeight','bold');
ylabel('Pearson r  (ACD ~ ARA)','FontSize',10,'FontWeight','bold');
title('ACD–ARA Correlation per Condition','FontSize',11,'FontWeight','bold');
set(ax3,'Box','off','TickDir','out','LineWidth',1.2,'FontSize',9);
exportgraphics(fig3,'./acd_ara_figure3_r_barchart.png','Resolution',300);
fprintf('Saved acd_ara_figure3_r_barchart.png\n');

fprintf('\nAll done.\n');
