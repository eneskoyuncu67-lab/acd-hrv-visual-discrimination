%% ============================
%   ACD_ARA_MODERATION.M  —  STANDALONE
%
%   Tests whether accuracy and RT moderate the ACD->ARA slope.
%
%   MODEL A: ARA ~ ACD * Acc + ConditionVal + (1|SubjectID)
%     ACD x Acc coefficient: does the slope of ARA on ACD differ
%     between Hit and Miss trials, controlling for condition?
%
%   MODEL B: ARA ~ ACD * RT + ConditionVal + (1|SubjectID)
%     ACD x RT coefficient: does the slope of ARA on ACD change
%     with response latency, controlling for condition?
%
%   All 7 exclusion criteria applied to both ACD and ARA simultaneously.
%   Only paired trials (valid ACD AND ARA) are included.
%
%   FIGURE LAYOUT (mirrors ACD accuracy model):
%     Top row (Model A): scatter ACD vs ARA per condition,
%       Hit (green) and Miss (red) overlaid with separate fitted lines.
%     Bottom row (Model B): scatter ACD vs ARA per condition,
%       points coloured by RT (cool=fast, warm=slow), continuous fitted
%       lines at RT mean ± 1 SD.
%
%   Requires: organizedData in workspace (or organizedData_by_condition.mat)
% ============================

%% ---- Load organizedData ----
if ~exist('organizedData','var')
    if exist('organizedData_by_condition','var')
        organizedData = organizedData_by_condition;
    elseif isfile('organizedData_by_condition.mat')
        S = load('organizedData_by_condition.mat');
        if isfield(S,'organizedData'), organizedData = S.organizedData;
        elseif isfield(S,'organizedData_by_condition'), organizedData = S.organizedData_by_condition;
        else, error('organizedData not found.'); end
    else, error('organizedData not found.'); end
end

%% ---- Settings ----
FS           = 100;
stim_dur_smp = round(16.7 / (1000/FS));

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
hit_col    = [0.20 0.60 0.20];
miss_col   = [0.90 0.40 0.40];
hit_edg    = [0.10 0.30 0.10];
miss_edg   = [0.50 0.10 0.10];

%% ============================
%   STEP 1: DATA PROCESSING — PAIRED ACD + ARA
%   Same 7 exclusion criteria applied to both measures simultaneously.
% ============================
fprintf('\n=== DATA PROCESSING ===\n');
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

            nTrials    = length(stim_positions);
            raw_d_acd  = nan(nTrials, 3);
            raw_d_ara  = nan(nTrials, 3);
            ibi_sroufe = nan(nTrials, 2);
            ibi_ara    = nan(nTrials, 2);
            trial_meta = nan(nTrials, 3);

            %% PASS 1
            for tIdx = 1:nTrials
                stim_samp  = stim_positions(tIdx);
                stim_off   = stim_samp + stim_dur_smp;
                pre_peaks  = find(peak_positions < stim_samp);
                post_peaks = find(peak_positions > stim_off);

                % ACD
                if numel(pre_peaks) >= 4
                    k0=pre_peaks(end); k1p=pre_peaks(end-1); k2p=pre_peaks(end-2);
                    if k0>=2 && k1p>=2 && k2p>=2 && ...
                       k0<=length(ibi_ms) && k1p<=length(ibi_ms) && k2p<=length(ibi_ms)
                        ib_s=ibi_ms(k0); ib_sb=ibi_ms(k0-1);
                        ib_s1=ibi_ms(k1p); ib_s1b=ibi_ms(k1p-1);
                        ib_s2=ibi_ms(k2p); ib_s2b=ibi_ms(k2p-1);
                        all_pre=[ib_sb,ib_s,ib_s1b,ib_s1,ib_s2b,ib_s2];
                        if ~(any(all_pre<444)||any(all_pre>1333))
                            raw_d_acd(tIdx,1)=ib_s2-ib_s2b;
                            raw_d_acd(tIdx,2)=ib_s1-ib_s1b;
                            raw_d_acd(tIdx,3)=ib_s-ib_sb;
                            ibi_sroufe(tIdx,1)=ib_s2;
                            ibi_sroufe(tIdx,2)=ib_s;
                        end
                    end
                end

                % ARA
                if numel(post_peaks) >= 3
                    kp1=post_peaks(1); kp2=post_peaks(2); kp3=post_peaks(3);
                    if kp1>=2 && kp2>=2 && kp3>=2 && ...
                       kp1<=length(ibi_ms) && kp2<=length(ibi_ms) && kp3<=length(ibi_ms)
                        ib1=ibi_ms(kp1); ib0=ibi_ms(kp1-1);
                        ib2=ibi_ms(kp2); ib_b2=ibi_ms(kp2-1);
                        ib3=ibi_ms(kp3); ib_b3=ibi_ms(kp3-1);
                        all_post=[ib0,ib1,ib_b2,ib2,ib_b3,ib3];
                        if ~(any(all_post<444)||any(all_post>1333))
                            raw_d_ara(tIdx,1)=ib1-ib0;
                            raw_d_ara(tIdx,2)=ib2-ib_b2;
                            raw_d_ara(tIdx,3)=ib3-ib_b3;
                            ibi_ara(tIdx,1)=ib1;
                            ibi_ara(tIdx,2)=ib3;
                        end
                    end
                end

                if tIdx<=numel(rt_vec) && tIdx<=numel(acc_vec)
                    trial_meta(tIdx,:)=[tIdx, rt_vec(tIdx), acc_vec(tIdx)];
                end
            end

            %% PASS 1b: Criterion 3
            valid_ibis = ibi_ms(ibi_ms>=444 & ibi_ms<=1333);
            if numel(valid_ibis) >= 2
                thresh = 3*std(valid_ibis);
                raw_d_acd(abs(raw_d_acd)>thresh) = NaN;
                raw_d_ara(abs(raw_d_ara)>thresh) = NaN;
            end

            %% PASS 2: Criterion 4
            for b = 1:3
                col=raw_d_acd(:,b); valid=~isnan(col);
                if sum(valid)>=2
                    mu_b=mean(col(valid)); sd_b=std(col(valid));
                    if sd_b>0, col(abs((col-mu_b)/sd_b)>3)=NaN; raw_d_acd(:,b)=col; end
                end
                col=raw_d_ara(:,b); valid=~isnan(col);
                if sum(valid)>=2
                    mu_b=mean(col(valid)); sd_b=std(col(valid));
                    if sd_b>0, col(abs((col-mu_b)/sd_b)>3)=NaN; raw_d_ara(:,b)=col; end
                end
            end

            %% PASS 3: Criteria 5 & 6
            for tIdx = 1:nTrials
                if all(isnan(raw_d_acd(tIdx,:))), continue; end
                if all(isnan(raw_d_ara(tIdx,:))), continue; end
                if any(isnan(trial_meta(tIdx,:))), continue; end
                is3=ibi_sroufe(tIdx,1); is0=ibi_sroufe(tIdx,2);
                ip1=ibi_ara(tIdx,1);    ip3=ibi_ara(tIdx,2);
                if isnan(is3)||isnan(is0)||isnan(ip1)||isnan(ip3), continue; end
                acd_val = (60000/is3)-(60000/is0);
                ara_val = (60000/ip3)-(60000/ip1);
                all_data_cell = [all_data_cell; ...
                    {cLabel, char(sname), tIdx, acd_val, ara_val, ...
                     trial_meta(tIdx,2), trial_meta(tIdx,3)}];
            end

        catch, continue; end
    end
end

%% Build table
tbl = cell2table(all_data_cell, ...
    'VariableNames',{'Condition','SubjectID','TrialIndex','ACD','ARA','RT','Acc'});
tbl.Condition    = categorical(tbl.Condition, gridLabels);
tbl.Acc          = categorical(tbl.Acc, [0,1], {'Miss','Hit'});
tbl              = rmmissing(tbl);
tbl.ConditionVal = cellfun(@(x) condNumMap(x), cellstr(tbl.Condition));

%% Criterion 7: z±3 on ACD and ARA per condition
for i = 1:nConds
    rows = tbl.Condition==gridLabels{i};
    for fld = {'ACD','ARA'}
        v=tbl.(fld{1})(rows); mu=mean(v,'omitnan'); sd=std(v,'omitnan');
        if sd>0, v(abs((v-mu)/sd)>3)=NaN; tbl.(fld{1})(rows)=v; end
    end
end
tbl = tbl(~isnan(tbl.ACD) & ~isnan(tbl.ARA), :);
fprintf('Paired trials after all criteria: %d\n', height(tbl));

%% RT units
if median(tbl.RT(tbl.RT>0)) < 10, tbl.RT = tbl.RT*1000; end

%% Downsample for Model A (Hit/Miss balance per condition)
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
fprintf('N paired: %d | N downsampled (Model A): %d\n\n', height(tbl), height(balanced_data));

%% ============================
%   STEP 2: FIT MODELS
% ============================

%% MODEL A: ARA ~ ACD * Acc + ConditionVal + (1|SubjectID)  [downsampled]
fprintf('================================================================\n');
fprintf('  MODEL A: ARA ~ ACD * Acc + ConditionVal + (1|SubjectID)\n');
fprintf('  Does accuracy moderate the ACD->ARA slope?\n');
fprintf('  ACD x Acc: slope difference Hit vs Miss\n');
fprintf('================================================================\n');
lme_A  = fitlme(balanced_data, ...
    'ARA ~ ACD * Acc + ConditionVal + (1|SubjectID)', 'FitMethod','ML');
anov_A = anova(lme_A,'DFMethod','Satterthwaite');
disp(lme_A); fprintf('Type III ANOVA:\n'); disp(anov_A);
fprintf('R2=%.4f | AIC=%.2f | BIC=%.2f | LogLik=%.4f | N=%d\n\n', ...
    lme_A.Rsquared.Ordinary, lme_A.ModelCriterion.AIC, ...
    lme_A.ModelCriterion.BIC, lme_A.LogLikelihood, lme_A.NumObservations);

%% MODEL B: ARA ~ ACD * RT + ConditionVal + (1|SubjectID)
fprintf('================================================================\n');
fprintf('  MODEL B: ARA ~ ACD * RT + ConditionVal + (1|SubjectID)\n');
fprintf('  Does RT moderate the ACD->ARA slope?\n');
fprintf('  ACD x RT: slope change per unit RT\n');
fprintf('================================================================\n');
lme_B  = fitlme(tbl, ...
    'ARA ~ ACD * RT + ConditionVal + (1|SubjectID)', 'FitMethod','ML');
anov_B = anova(lme_B,'DFMethod','Satterthwaite');
disp(lme_B); fprintf('Type III ANOVA:\n'); disp(anov_B);
fprintf('R2=%.4f | AIC=%.2f | BIC=%.2f | LogLik=%.4f | N=%d\n\n', ...
    lme_B.Rsquared.Ordinary, lme_B.ModelCriterion.AIC, ...
    lme_B.ModelCriterion.BIC, lme_B.LogLikelihood, lme_B.NumObservations);

%% Extract Model A coefficients
int_A     = lme_A.Coefficients.Estimate(strcmp(lme_A.Coefficients.Name,'(Intercept)'));
b_ACD_A   = lme_A.Coefficients.Estimate(strcmp(lme_A.Coefficients.Name,'ACD'));
b_Acc_A   = lme_A.Coefficients.Estimate(strcmp(lme_A.Coefficients.Name,'Acc_Hit'));
b_Cond_A  = lme_A.Coefficients.Estimate(strcmp(lme_A.Coefficients.Name,'ConditionVal'));
b_Int_A   = lme_A.Coefficients.Estimate(strcmp(lme_A.Coefficients.Name,'ACD:Acc_Hit'));
p_ACD_A   = lme_A.Coefficients.pValue(strcmp(lme_A.Coefficients.Name,'ACD'));
p_Acc_A   = lme_A.Coefficients.pValue(strcmp(lme_A.Coefficients.Name,'Acc_Hit'));
p_Int_A   = lme_A.Coefficients.pValue(strcmp(lme_A.Coefficients.Name,'ACD:Acc_Hit'));
r2_A      = lme_A.Rsquared.Ordinary;

F_Int_A   = anov_A.FStat(strcmp(anov_A.Term,'ACD:Acc'));
df2_Int_A = anov_A.DF2(strcmp(anov_A.Term,'ACD:Acc'));

if p_ACD_A<0.001, sp_ACD_A='< .001'; else, sp_ACD_A=sprintf('= %.3f',p_ACD_A); end
if p_Acc_A<0.001, sp_Acc_A='< .001'; else, sp_Acc_A=sprintf('= %.3f',p_Acc_A); end
if p_Int_A<0.001, sp_Int_A='< .001'; else, sp_Int_A=sprintf('= %.3f',p_Int_A); end

%% Extract Model B coefficients
int_B      = lme_B.Coefficients.Estimate(strcmp(lme_B.Coefficients.Name,'(Intercept)'));
b_ACD_B    = lme_B.Coefficients.Estimate(strcmp(lme_B.Coefficients.Name,'ACD'));
b_RT_B     = lme_B.Coefficients.Estimate(strcmp(lme_B.Coefficients.Name,'RT'));
b_Cond_B   = lme_B.Coefficients.Estimate(strcmp(lme_B.Coefficients.Name,'ConditionVal'));
b_Int_B    = lme_B.Coefficients.Estimate(strcmp(lme_B.Coefficients.Name,'ACD:RT'));
p_ACD_B    = lme_B.Coefficients.pValue(strcmp(lme_B.Coefficients.Name,'ACD'));
p_RT_B     = lme_B.Coefficients.pValue(strcmp(lme_B.Coefficients.Name,'RT'));
p_Int_B    = lme_B.Coefficients.pValue(strcmp(lme_B.Coefficients.Name,'ACD:RT'));
r2_B       = lme_B.Rsquared.Ordinary;

F_Int_B    = anov_B.FStat(strcmp(anov_B.Term,'ACD:RT'));
df2_Int_B  = anov_B.DF2(strcmp(anov_B.Term,'ACD:RT'));

if p_ACD_B<0.001, sp_ACD_B='< .001'; else, sp_ACD_B=sprintf('= %.3f',p_ACD_B); end
if p_RT_B<0.001,  sp_RT_B='< .001';  else, sp_RT_B=sprintf('= %.3f',p_RT_B);   end
if p_Int_B<0.001, sp_Int_B='< .001'; else, sp_Int_B=sprintf('= %.3f',p_Int_B); end

mean_RT  = mean(tbl.RT);
sd_RT    = std(tbl.RT);

%% ============================
%   STEP 3: FIGURE
%   Top row:    Model A — scatter per condition, Hit vs Miss
%   Bottom row: Model B — scatter per condition, coloured by RT
%   Layout mirrors figure3 from lme_all_models.m
% ============================

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 13.3 9.0]);

panel_h  = 0.36;
gap_row  = 0.10;
row_b_b  = 0.07;                          % bottom of Model B row
row_a_b  = row_b_b + panel_h + gap_row;   % bottom of Model A row

left_m = 0.07; right_m = 0.02;
uw     = 1 - left_m - right_m;
gap_c  = 0.018;
pw     = (uw - (nConds-1)*gap_c) / nConds;
px     = left_m + (0:nConds-1).*(pw+gap_c);

%% --- ROW A: Model A — Hit vs Miss scatter per condition ---
annotation('textbox',[left_m-0.01, row_a_b+panel_h+0.005, 0.04, 0.04], ...
    'String','A','LineStyle','none','FontSize',13,'FontWeight','bold');
annotation('textbox',[left_m+0.02, row_a_b+panel_h+0.005, 0.55, 0.04], ...
    'String','ARA \sim ACD \times Accuracy + Condition  (downsampled)', ...
    'EdgeColor','none','FontSize',10,'FontWeight','bold');

% Condition-specific ACD range for fitted lines
acd_fit_x = linspace(min(tbl.ACD), max(tbl.ACD), 100);

% Fitted slopes for Miss (Acc_Hit=0) and Hit (Acc_Hit=1) at mean ConditionVal
mean_cond = mean(tbl.ConditionVal);
slope_miss = b_ACD_A;                  % slope for Miss (reference)
slope_hit  = b_ACD_A + b_Int_A;       % slope for Hit
int_miss   = int_A + b_Cond_A*mean_cond;
int_hit    = int_A + b_Acc_A + b_Cond_A*mean_cond;

for i = 1:nConds
    ax = axes('Units','normalized','Position',[px(i), row_a_b, pw, panel_h]);
    hold on; grid on;

    for aIdx = 1:2
        aLabel = {'Miss','Hit'}; acol = {miss_col, hit_col}; aedg = {miss_edg, hit_edg};
        rows = balanced_data.Condition==gridLabels{i} & balanced_data.Acc==aLabel{aIdx};
        acd_v = balanced_data.ACD(rows); ara_v = balanced_data.ARA(rows);
        valid = ~isnan(acd_v)&~isnan(ara_v);
        acd_v = acd_v(valid); ara_v = ara_v(valid);
        if numel(acd_v) < 3, continue; end

        scatter(acd_v+(rand(size(acd_v))-0.5)*0.2, ara_v, 5, acol{aIdx}, 'filled', ...
            'MarkerFaceAlpha',0.18,'MarkerEdgeColor','none','HandleVisibility','off');
    end

    % Fitted lines at condition-specific intercept adjustment
    cond_adj = b_Cond_A * (condNumVals(i) - mean_cond);
    y_miss = (int_miss + cond_adj) + slope_miss * acd_fit_x;
    y_hit  = (int_hit  + cond_adj) + slope_hit  * acd_fit_x;
    plot(acd_fit_x, y_miss, '--', 'Color', miss_col, 'LineWidth', 2.0, 'HandleVisibility','off');
    plot(acd_fit_x, y_hit,  '--', 'Color', hit_col,  'LineWidth', 2.0, 'HandleVisibility','off');

    % Mean markers per accuracy level
    for aIdx = 1:2
        aLabel = {'Miss','Hit'}; acol = {miss_col, hit_col};
        rows = balanced_data.Condition==gridLabels{i} & balanced_data.Acc==aLabel{aIdx};
        acd_v=balanced_data.ACD(rows); ara_v=balanced_data.ARA(rows);
        valid=~isnan(acd_v)&~isnan(ara_v); acd_v=acd_v(valid); ara_v=ara_v(valid);
        if numel(acd_v)<3, continue; end
        scatter(mean(acd_v), mean(ara_v), 75, acol{aIdx}, 'filled', ...
            'MarkerEdgeColor','k','LineWidth',1.2,'HandleVisibility','off');
    end

    xline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',0.7,'HandleVisibility','off');
    yline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',0.7,'HandleVisibility','off');
    title(gridLabels{i},'FontSize',10,'FontWeight','bold');
    xlabel('ACD (bpm)','FontSize',8);
    if i==1
        ylabel('ARA (bpm)','FontSize',9,'FontWeight','bold');
        scatter(nan,nan,30,hit_col,'filled','MarkerEdgeColor','k','DisplayName','Hit');
        scatter(nan,nan,30,miss_col,'filled','MarkerEdgeColor','k','DisplayName','Miss');
        legend('Location','best','FontSize',7,'Box','off');
    end
    set(ax,'Box','off','TickDir','out','LineWidth',1.0,'FontSize',8);

    % Annotation box — only on first panel to avoid clutter
    if i == 1
        text(0.03,0.97, sprintf('\\beta_{ACD}=%.3f, p %s\n\\beta_{Acc}=%.3f, p %s\n\\beta_{ACD\\timesAcc}=%.3f, p %s\nR^2=%.3f', ...
            b_ACD_A,sp_ACD_A, b_Acc_A,sp_Acc_A, b_Int_A,sp_Int_A, r2_A), ...
            'Units','normalized','FontSize',6.5,'VerticalAlignment','top', ...
            'BackgroundColor','w','EdgeColor',[0.6 0.6 0.6]);
    end
end

%% --- ROW B: Model B — RT-coloured scatter per condition ---
annotation('textbox',[left_m-0.01, row_b_b+panel_h+0.005, 0.04, 0.04], ...
    'String','B','LineStyle','none','FontSize',13,'FontWeight','bold');
annotation('textbox',[left_m+0.02, row_b_b+panel_h+0.005, 0.55, 0.04], ...
    'String','ARA \sim ACD \times RT + Condition', ...
    'EdgeColor','none','FontSize',10,'FontWeight','bold');

% RT-based fitted lines: at mean RT ± 1 SD
rt_lo = mean_RT - sd_RT;
rt_hi = mean_RT + sd_RT;

for i = 1:nConds
    ax = axes('Units','normalized','Position',[px(i), row_b_b, pw, panel_h]);
    hold on; grid on;

    rows = tbl.Condition==gridLabels{i};
    acd_v=tbl.ACD(rows); ara_v=tbl.ARA(rows); rt_v=tbl.RT(rows);
    valid=~isnan(acd_v)&~isnan(ara_v)&~isnan(rt_v);
    acd_v=acd_v(valid); ara_v=ara_v(valid); rt_v=rt_v(valid);
    if numel(acd_v) < 3, continue; end

    % Colour points by RT (cool = fast, warm = slow)
    rt_norm = (rt_v - min(rt_v)) / max(range(rt_v), eps);
    cmap    = cool(256);
    cidx    = max(1, min(256, round(rt_norm*255)+1));
    pt_cols = cmap(cidx,:);
    scatter(acd_v, ara_v, 5, pt_cols, 'filled', ...
        'MarkerFaceAlpha',0.25,'HandleVisibility','off');

    % Fitted lines at RT mean-1SD (cool/blue) and mean+1SD (warm/red)
    cond_adj_B = b_Cond_B*(condNumVals(i) - mean(tbl.ConditionVal));
    y_lo = (int_B + b_RT_B*rt_lo + cond_adj_B) + (b_ACD_B + b_Int_B*rt_lo)*acd_fit_x;
    y_hi = (int_B + b_RT_B*rt_hi + cond_adj_B) + (b_ACD_B + b_Int_B*rt_hi)*acd_fit_x;
    plot(acd_fit_x, y_lo, '--', 'Color', [0.20 0.55 0.85], 'LineWidth', 2.0, 'HandleVisibility','off');
    plot(acd_fit_x, y_hi, '--', 'Color', [0.85 0.30 0.20], 'LineWidth', 2.0, 'HandleVisibility','off');

    xline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',0.7,'HandleVisibility','off');
    yline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',0.7,'HandleVisibility','off');
    title(gridLabels{i},'FontSize',10,'FontWeight','bold');
    xlabel('ACD (bpm)','FontSize',8);
    if i==1
        ylabel('ARA (bpm)','FontSize',9,'FontWeight','bold');
        plot(nan,nan,'--','Color',[0.20 0.55 0.85],'LineWidth',2,'DisplayName','RT mean-1SD (fast)');
        plot(nan,nan,'--','Color',[0.85 0.30 0.20],'LineWidth',2,'DisplayName','RT mean+1SD (slow)');
        legend('Location','best','FontSize',7,'Box','off');
    end
    set(ax,'Box','off','TickDir','out','LineWidth',1.0,'FontSize',8);

    if i == 1
        text(0.03,0.97, sprintf('\\beta_{ACD}=%.3f, p %s\n\\beta_{RT}=%.3f, p %s\n\\beta_{ACD\\timesRT}=%.3f, p %s\nR^2=%.3f', ...
            b_ACD_B,sp_ACD_B, b_RT_B,sp_RT_B, b_Int_B,sp_Int_B, r2_B), ...
            'Units','normalized','FontSize',6.5,'VerticalAlignment','top', ...
            'BackgroundColor','w','EdgeColor',[0.6 0.6 0.6]);
    end
end

exportgraphics(fig,'./figure_ACD_ARA_moderation.png','Resolution',300);
fprintf('Saved figure_ACD_ARA_moderation.png\n\nAll done.\n');
