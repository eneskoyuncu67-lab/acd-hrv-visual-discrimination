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

%% 1. Data Processing
targetKeys = {'cond_0_6667','cond_0_7692','cond_1_0000','cond_1_4286','cond_2_0000'};
gridLabels = {'66.67%','76.92%','100%','142.86%','200%'};

existingFields     = fieldnames(organizedData);
validIdx           = ismember(targetKeys, existingFields);
gridConditionNames = targetKeys(validIdx);
gridLabels         = gridLabels(validIdx);

FS = 100;

all_trial_data_cell = {};
all_delta_cell      = {};

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

            rpeak_row  = squeeze(condData.RPeak(sIdx,:));
            timing_row = squeeze(condData.TrialTiming(sIdx,:));

            stim_positions = find(timing_row > 0);
            peak_positions = find(rpeak_row == 1);
            ibi_ms         = diff(peak_positions) * (1000/FS);

            for tIdx = 1:length(stim_positions)
                stim_samp      = stim_positions(tIdx);
                pre_stim_peaks = find(peak_positions < stim_samp);

                if numel(pre_stim_peaks) < 4, continue; end

                k_stim  = pre_stim_peaks(end);
                k_stim1 = pre_stim_peaks(end-1);
                k_stim2 = pre_stim_peaks(end-2);

                if k_stim < 2 || k_stim1 < 2 || k_stim2 < 2, continue; end
                if k_stim > length(ibi_ms)+1, continue; end

                ibi_at_stim    = ibi_ms(k_stim  - 1);
                ibi_at_stim1   = ibi_ms(k_stim1 - 1);
                ibi_at_stim2   = ibi_ms(k_stim2 - 1);
                ibi_prev_stim  = ibi_ms(k_stim  - 2);
                ibi_prev_stim1 = ibi_ms(k_stim1 - 2);
                ibi_prev_stim2 = ibi_ms(k_stim2 - 2);

                all_ibis = [ibi_at_stim, ibi_at_stim1, ibi_at_stim2, ...
                            ibi_prev_stim, ibi_prev_stim1, ibi_prev_stim2];
                if any(all_ibis < 444) || any(all_ibis > 1333), continue; end

                d2 = ibi_at_stim2 - ibi_prev_stim2;
                d1 = ibi_at_stim1 - ibi_prev_stim1;
                d0 = ibi_at_stim  - ibi_prev_stim;
                acd_trial = mean([d0, d1, d2]);

                if tIdx <= numel(rt_vec) && tIdx <= numel(acc_vec)
                    all_trial_data_cell = [all_trial_data_cell; ...
                        {cLabel, char(sname), tIdx, acd_trial, rt_vec(tIdx), acc_vec(tIdx)}];
                    all_delta_cell = [all_delta_cell; ...
                        {cLabel, char(sname), tIdx, d2, d1, d0, acc_vec(tIdx)}];
                end
            end
        catch
            continue;
        end
    end
end

%% Build tables
all_trial_pooled = cell2table(all_trial_data_cell, ...
    'VariableNames',{'Condition','SubjectID','TrialIndex','ACD','RT','Acc'});
all_trial_pooled.Condition = categorical(all_trial_pooled.Condition, gridLabels);
all_trial_pooled.Acc       = categorical(all_trial_pooled.Acc, [0,1], {'Miss','Hit'});
all_trial_pooled           = rmmissing(all_trial_pooled);

condNumMap = containers.Map({'66.67%','76.92%','100%','142.86%','200%'}, ...
                             [66.67, 76.92, 100, 142.86, 200]);
all_trial_pooled.ConditionVal = cellfun(@(x) condNumMap(x), ...
    cellstr(all_trial_pooled.Condition));
all_trial_pooled_full = all_trial_pooled;

delta_table = cell2table(all_delta_cell, ...
    'VariableNames',{'Condition','SubjectID','TrialIndex','D_stim2','D_stim1','D_stim','Acc'});
delta_table.Condition = categorical(delta_table.Condition, gridLabels);

fprintf('Total trials extracted: %d\n', height(all_trial_pooled_full));
fprintf('ACD stats: mean=%.2f ms, SD=%.2f ms, min=%.2f, max=%.2f\n', ...
    mean(all_trial_pooled_full.ACD), std(all_trial_pooled_full.ACD), ...
    min(all_trial_pooled_full.ACD), max(all_trial_pooled_full.ACD));

%% Downsample for Panel C
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

%% Layout constants
nConds       = numel(gridLabels);
condColors   = [0.7 0.85 0.9; 0.5 0.7 0.85; 0.3 0.55 0.8; 0.1 0.35 0.7; 0.05 0.15 0.5];
condColors_e = [0.8 0.2 0.2; 0.9 0.6 0.1; 0.2 0.6 0.2; 0.1 0.4 0.8; 0.5 0.1 0.7];
left_margin  = 0.10;
right_margin = 0.03;
usable_w     = 1 - left_margin - right_margin;
gap          = 0.025;
sub_w        = (usable_w - (nConds-1)*gap) / nConds;
sub_x        = left_margin + (0:nConds-1) .* (sub_w + gap);
row_top_bot  = 0.57;
row_top_h    = 0.35;
row_bot_bot  = 0.10;
row_bot_h    = 0.35;
panelA_ylim  = [-150 150];

%% ============================
%   FIT ALL MODELS
% ============================
fprintf('\n=== ANALYSIS A: RT -> DIBI ===\n');
lme_A = cell(nConds,1);
for i = 1:nConds
    subData  = all_trial_pooled_full(all_trial_pooled_full.Condition == gridLabels{i}, :);
    lme_A{i} = fitlme(subData, 'ACD ~ RT + (1|SubjectID)');
    fprintf('Condition %-10s | beta=%7.3f | R2=%.3f | p=%.4f | n=%d\n', ...
        gridLabels{i}, lme_A{i}.Coefficients.Estimate(2), ...
        lme_A{i}.Rsquared.Ordinary, lme_A{i}.Coefficients.pValue(2), height(subData));
end

fprintf('\n=== ANALYSIS B: DIBI ~ RT + ConditionVal ===\n');
lme_B = fitlme(all_trial_pooled_full, 'ACD ~ RT + ConditionVal + (1|SubjectID)');
disp(lme_B);
anovaResults_B = anova(lme_B, 'DFMethod','Satterthwaite');
fprintf('\nANOVA (Type 3):\n'); disp(anovaResults_B);

fprintf('\n=== ANALYSIS C (Downsampled): Accuracy x Condition ===\n');
lme_C = fitlme(all_trial_pooled_downsampled, 'ACD ~ Acc * Condition + (1|SubjectID)');
disp(lme_C);
anovaResults_C = anova(lme_C, 'DFMethod','Satterthwaite');
fprintf('\nANOVA (Type 3):\n'); disp(anovaResults_C);

%% Extract model values
beta_RT_B   = lme_B.Coefficients.Estimate(2);
pval_RT_B   = lme_B.Coefficients.pValue(2);
r2_B        = lme_B.Rsquared.Ordinary;
pval_Cond_B = anovaResults_B.pValue(3);
f_rt_B      = anovaResults_B.FStat(2);
f_cv_B      = anovaResults_B.FStat(3);
b_cv_B      = lme_B.Coefficients.Estimate(3);
se_cv_B     = lme_B.Coefficients.SE(3);
t_cv_B      = lme_B.Coefficients.tStat(3);

acc_row    = strcmp(lme_C.Coefficients.Name, 'Acc_Hit');
beta_Acc_C = lme_C.Coefficients.Estimate(acc_row);
se_Acc_C   = lme_C.Coefficients.SE(acc_row);
t_Acc_C    = lme_C.Coefficients.tStat(acc_row);
pval_Acc_C = lme_C.Coefficients.pValue(acc_row);
r2_C       = lme_C.Rsquared.Ordinary;
pval_Interaction_C = anovaResults_C.pValue(4);

%% p-value strings
if pval_RT_B         < 0.001, ps_RT_B   = '< .001'; else, ps_RT_B   = sprintf('%.3f',pval_RT_B);   end
if pval_Cond_B       < 0.001, ps_Cond_B = '< .001'; else, ps_Cond_B = sprintf('%.3f',pval_Cond_B); end
if pval_Acc_C        < 0.001, ps_Acc_C  = '< .001'; else, ps_Acc_C  = sprintf('%.3f',pval_Acc_C);  end
if pval_Interaction_C< 0.001, ps_Int_C  = '< .001'; else, ps_Int_C  = sprintf('%.3f',pval_Interaction_C); end

f_cond_C = anovaResults_C.FStat(strcmp(anovaResults_C.Term,'Condition'));
p_cond_C = anovaResults_C.pValue(strcmp(anovaResults_C.Term,'Condition'));
f_int_C  = anovaResults_C.FStat(strcmp(anovaResults_C.Term,'Condition:Acc'));
p_int_C  = anovaResults_C.pValue(strcmp(anovaResults_C.Term,'Condition:Acc'));
if p_cond_C < 0.001, ps_cond_C = '< .001'; else, ps_cond_C = sprintf('%.3f',p_cond_C); end
if p_int_C  < 0.001, ps_int_C  = '< .001'; else, ps_int_C  = sprintf('%.3f',p_int_C);  end

%% ============================
%   FIGURE 1: PANELS A + B
% ============================
figAB = figure('Color','w','Units','inches','Position',[1 1 6.5 5]);

for i = 1:nConds
    ax = axes('Units','normalized', ...
              'Position',[sub_x(i), row_top_bot, sub_w, row_top_h]);
    hold on; grid on;

    subData = all_trial_pooled_full(all_trial_pooled_full.Condition == gridLabels{i}, :);
    beta    = lme_A{i}.Coefficients.Estimate(2);
    r2      = lme_A{i}.Rsquared.Ordinary;
    pval    = lme_A{i}.Coefficients.pValue(2);
    if pval < 0.001, ps = '< .001'; else, ps = sprintf('%.3f',pval); end

    scatter(subData.RT, subData.ACD, 10, [0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.25);
    xl = linspace(min(subData.RT), max(subData.RT), 30);
    yl = lme_A{i}.Coefficients.Estimate(1) + beta * xl;
    plot(xl, yl, 'r', 'LineWidth',2);

    text(0.05, 0.95, sprintf('R^2=%.3f\n\\beta=%.3f\np=%s',r2,beta,ps), ...
        'Units','normalized','FontSize',9,'VerticalAlignment','top');

    ylim(panelA_ylim);
    title(gridLabels{i});
    xlabel('RT (s)');
    if i==1, ylabel('\DeltaIBI (ms)'); end
    set(ax,'Box','off','TickDir','out','LineWidth',1.5);
    ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;
    if i==1
        annotation('textbox',[ax.Position(1)-0.04,ax.Position(2)+ax.Position(4)+0.02,0.04,0.04], ...
            'String','A','LineStyle','none','FontSize',12,'FontWeight','bold');
    end
end

axB = axes('Units','normalized','Position',[left_margin,row_bot_bot,usable_w,row_bot_h]);
hold on; grid on;

for i = 1:nConds
    vals = all_trial_pooled_full.ACD(all_trial_pooled_full.Condition == gridLabels{i});
    [f,xi] = ksdensity(vals);
    f = f/max(f)*0.35;
    fill(i+[f -fliplr(f)],[xi fliplr(xi)],condColors(i,:),'FaceAlpha',0.7,'EdgeColor','k','LineWidth',1);
    scatter(i+(rand(size(vals))-0.5)*0.1,vals,8,'k','filled','MarkerFaceAlpha',0.2);
end

xticks(1:nConds); xticklabels(gridLabels);
ylabel('\DeltaIBI (ms)');
xlabel('Experimental Condition');
yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2);
text(0.02,0.98,sprintf('R^2=%.3f\n\\beta_{RT}=%.3f\np_{RT}=%s\np_{Cond}=%s', ...
    r2_B,beta_RT_B,ps_RT_B,ps_Cond_B), ...
    'Units','normalized','FontSize',11,'VerticalAlignment','top');
set(axB,'Box','off','TickDir','out','LineWidth',1.5);
axB.XAxis.LineWidth=1.5; axB.YAxis.LineWidth=1.5;
annotation('textbox',[axB.Position(1)-0.04,axB.Position(2)+axB.Position(4)+0.02,0.04,0.04], ...
    'String','B','LineStyle','none','FontSize',12,'FontWeight','bold');

figure(figAB);
if saving_figures==1
    print(figAB,'-dpng','-r300',sprintf('%s/figure_AB_DIBI_RT_Condition.png',figuredir));
    fprintf('Saved figure AB\n');
end

%% ============================
%   FIGURE 2: RMSSD
% ============================
targetKeys_r  = {'cond_0_6667','cond_0_7692','cond_1_0000','cond_1_4286','cond_2_0000'};
numericValues = [66.67, 76.92, 100, 142.86, 200];
window_samples = 10 * FS;
step_samples   = 2  * FS;
all_win_metrics = {};

for cIdx = 1:length(targetKeys_r)
    cName = targetKeys_r{cIdx};
    if ~isfield(organizedData,cName), continue; end
    condData = organizedData.(cName);
    for sIdx = 1:length(condData.subjects)
        subjectRPeak = squeeze(condData.RPeak(sIdx,:));
        for start_s = 1:step_samples:(length(subjectRPeak)-window_samples)
            end_s      = start_s + window_samples;
            rp_samples = find(subjectRPeak(start_s:end_s)==1);
            if numel(rp_samples) >= 6
                rr_ms = diff(rp_samples)*(1000/FS);
                rr_ms = rr_ms(rr_ms>=444 & rr_ms<=1333);
                if isempty(rr_ms), continue; end
                v_rmssd = sqrt(mean(diff(rr_ms).^2));
                if v_rmssd>5 && v_rmssd<100
                    all_win_metrics=[all_win_metrics; ...
                        {numericValues(cIdx),cIdx,char(condData.subjects{sIdx}),v_rmssd}];
                end
            end
        end
    end
end

metrics_table = cell2table(all_win_metrics, ...
    'VariableNames',{'ConditionVal','ConditionIdx','SubjectID','RMSSD'});
lme_hrv   = fitlme(metrics_table,'RMSSD ~ ConditionVal + (1|SubjectID)');
b0        = lme_hrv.Coefficients.Estimate(1);
b1        = lme_hrv.Coefficients.Estimate(2);
r2_hrv    = lme_hrv.Rsquared.Ordinary;
anova_hrv = anova(lme_hrv,'DFMethod','Satterthwaite');
pval_hrv  = anova_hrv.pValue(2);
if pval_hrv<0.001, ps_hrv='< .001'; else, ps_hrv=sprintf('%.3f',pval_hrv); end

figHRV = figure('Color','w','Units','inches','Position',[1 1 6.5 3]);
ax = axes('Units','normalized','Position',[0.12,0.18,0.83,0.68]);
hold on; grid on;

for i = 1:5
    vals = metrics_table.RMSSD(metrics_table.ConditionIdx==i);
    if ~isempty(vals)
        [f,xi] = ksdensity(vals,'BoundaryCorrection','reflection');
        f = f/max(f)*0.35;
        fill(i+[f,-fliplr(f)],[xi,fliplr(xi)],condColors(i,:),'FaceAlpha',0.7,'EdgeColor','k','LineWidth',1);
        scatter(i+(rand(size(vals))-0.5)*0.1,vals,5,'k','filled','MarkerFaceAlpha',0.1);
    end
end

yFit = b0+b1*numericValues;
scatter(1:5,yFit,80,'w','filled','MarkerEdgeColor',[0.8 0 0],'LineWidth',2);
text(0.98,0.95,sprintf('R^2=%.3f\n\\beta=%.4f\np=%s',r2_hrv,b1,ps_hrv), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',11);
ylabel('10s Window RMSSD (ms)','FontWeight','bold','FontSize',12);
xlabel('Experimental Condition','FontSize',12,'FontWeight','bold');
xticks(1:5); xticklabels(gridLabels); xlim([0.4 5.6]);
set(ax,'Box','off','TickDir','out','LineWidth',1.5);

figure(figHRV);
if saving_figures==1
    print(figHRV,'-dpng','-r300',sprintf('%s/figure_HRV_RMSSD.png',figuredir));
    fprintf('Saved figure HRV\n');
end

%% ============================
%   FIGURE 3: Histograms + Violin
% ============================
figHist = figure('Color','w','Units','inches','Position',[1 1 6.5 5]);

for i = 1:nConds
    ax = axes('Units','normalized', ...
              'Position',[sub_x(i),row_top_bot,sub_w,row_top_h]);
    hold on; grid on;
    vals = all_trial_pooled_full.ACD(all_trial_pooled_full.Condition==gridLabels{i});
    histogram(vals,40,'FaceColor',condColors(i,:),'EdgeColor',condColors(i,:)*0.6,'FaceAlpha',0.8,'LineWidth',0.8);
    mu=mean(vals); sig=std(vals);
    xl=linspace(min(vals),max(vals),200);
    yl=numel(vals)*(xl(2)-xl(1))*normpdf(xl,mu,sig);
    plot(xl,yl,'k-','LineWidth',1.8);
    xline(mu,'--r','LineWidth',1.5);
    text(0.97,0.97,sprintf('M=%.1f\nSD=%.1f\nn=%d',mu,sig,numel(vals)), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',9);
    title(gridLabels{i},'FontSize',12);
    xlabel('\DeltaIBI (ms)','FontSize',10);
    if i==1, ylabel('Count','FontSize',10); end
    set(ax,'Box','off','TickDir','out','LineWidth',1.5);
    if i==1
        annotation('textbox',[ax.Position(1)-0.04,ax.Position(2)+ax.Position(4)+0.02,0.04,0.04], ...
            'String','A','LineStyle','none','FontSize',12,'FontWeight','bold');
    end
end

lme_cond   = fitlme(all_trial_pooled_full,'ACD ~ ConditionVal + (1|SubjectID)');
anova_cond = anova(lme_cond,'DFMethod','Satterthwaite');
r2_cond    = lme_cond.Rsquared.Ordinary;
pval_cond  = anova_cond.pValue(2);
beta_cond  = lme_cond.Coefficients.Estimate(2);
if pval_cond<0.001, ps_cond='< .001'; else, ps_cond=sprintf('%.3f',pval_cond); end

axB4 = axes('Units','normalized','Position',[left_margin,row_bot_bot,usable_w,row_bot_h]);
hold on; grid on;

for i = 1:nConds
    vals = all_trial_pooled_full.ACD(all_trial_pooled_full.Condition==gridLabels{i});
    if ~isempty(vals)
        [f,xi] = ksdensity(vals,'BoundaryCorrection','reflection');
        f=f/max(f)*0.35;
        fill(i+[f,-fliplr(f)],[xi,fliplr(xi)],condColors(i,:),'FaceAlpha',0.7,'EdgeColor','k','LineWidth',1);
        scatter(i+(rand(size(vals))-0.5)*0.1,vals,8,'k','filled','MarkerFaceAlpha',0.2);
        scatter(i,mean(vals),60,'w','filled','MarkerEdgeColor',[0.8 0 0],'LineWidth',2);
    end
end

yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2);
text(0.02,0.98,sprintf('R^2=%.3f\n\\beta=%.3f\np=%s',r2_cond,beta_cond,ps_cond), ...
    'Units','normalized','FontSize',11,'VerticalAlignment','top');
xticks(1:nConds); xticklabels(gridLabels);
xlabel('Experimental Condition','FontSize',12,'FontWeight','bold');
ylabel('\DeltaIBI (ms)','FontSize',12,'FontWeight','bold');
xlim([0.4,nConds+0.6]);
set(axB4,'Box','off','TickDir','out','LineWidth',1.5);
annotation('textbox',[axB4.Position(1)-0.04,axB4.Position(2)+axB4.Position(4)+0.02,0.04,0.04], ...
    'String','B','LineStyle','none','FontSize',12,'FontWeight','bold');

figure(figHist);
if saving_figures==1
    print(figHist,'-dpng','-r300',sprintf('%s/figure_Hist_DIBIViolin.png',figuredir));
    fprintf('Saved figure Hist\n');
end

%% ============================
%   FIGURE B: Violin + Model B Table (presentation)
% ============================
figB_full = figure('Color','w','Units','inches','Position',[1 1 13 7]);

axB_new = axes('Units','normalized','Position',[0.08 0.40 0.88 0.53]);
hold on; grid on;

for i = 1:nConds
    vals = all_trial_pooled_full.ACD(all_trial_pooled_full.Condition==gridLabels{i});
    [f,xi] = ksdensity(vals);
    f=f/max(f)*0.35;
    fill(i+[f -fliplr(f)],[xi fliplr(xi)],condColors(i,:),'FaceAlpha',0.7,'EdgeColor','k','LineWidth',1);
    scatter(i+(rand(size(vals))-0.5)*0.1,vals,8,'k','filled','MarkerFaceAlpha',0.2);
    scatter(i,mean(vals),70,'w','filled','MarkerEdgeColor',[0.8 0 0],'LineWidth',2);
end

yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2);
xticks(1:nConds); xticklabels(gridLabels);
ylabel('\DeltaIBI (ms)','FontSize',13,'FontWeight','bold');
xlabel('Experimental Condition','FontSize',13,'FontWeight','bold');
title('Anticipatory Cardiac Deceleration by Condition','FontSize',13,'FontWeight','bold');
set(axB_new,'Box','off','TickDir','out','LineWidth',1.5,'FontSize',11);

axT_B = axes('Units','normalized','Position',[0 0.02 1 0.35]);
axis off; hold on;

b_rt  = lme_B.Coefficients.Estimate(2);
se_rt = lme_B.Coefficients.SE(2);
t_rt  = lme_B.Coefficients.tStat(2);
p_rt  = lme_B.Coefficients.pValue(2);
b_cv  = lme_B.Coefficients.Estimate(3);
se_cv = lme_B.Coefficients.SE(3);
t_cv  = lme_B.Coefficients.tStat(3);
p_cv  = lme_B.Coefficients.pValue(3);
if p_rt<0.001, ps_rt='< .001'; else, ps_rt=sprintf('%.3f',p_rt); end
if p_cv<0.001, ps_cv='< .001'; else, ps_cv=sprintf('%.3f',p_cv); end

B_headers = {'Predictor','\beta','SE','t','F','p','R^2'};
B_col_w   = [0.28 0.10 0.10 0.10 0.10 0.12 0.10];
B_col_x   = [0,cumsum(B_col_w(1:end-1))];
B_rows    = {
    {'RT',                 sprintf('%.3f',b_rt), sprintf('%.3f',se_rt), sprintf('%.3f',t_rt), sprintf('%.3f',f_rt_B), ps_rt, sprintf('%.3f',r2_B)};
    {'Condition (linear)', sprintf('%.3f',b_cv), sprintf('%.3f',se_cv), sprintf('%.3f',t_cv), sprintf('%.3f',f_cv_B), ps_cv, '—'}
};
B_row_bg = {[0.93 0.96 1.00],[0.86 0.92 0.98]};
rh=0.24; hh=0.28; ty=0.88;

text(0.5,0.98,'Model B: \DeltaIBI ~ RT + Condition + (1|SubjectID)', ...
    'Units','normalized','HorizontalAlignment','center','FontSize',11,'FontWeight','bold','Color',[0.1 0.1 0.1]);

rectangle('Position',[0.02 ty-hh 0.96 hh],'FaceColor',[0.15 0.25 0.50],'EdgeColor','none');
for h=1:numel(B_headers)
    text(0.02+B_col_x(h)+B_col_w(h)/2,ty-hh/2,B_headers{h}, ...
        'Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontSize',10,'FontWeight','bold','Color','w');
end

for r=1:numel(B_rows)
    ry=ty-hh-r*rh; rd=B_rows{r};
    rectangle('Position',[0.02 ry 0.96 rh],'FaceColor',B_row_bg{r},'EdgeColor',[0.82 0.82 0.82]);
    for h=1:numel(B_headers)
        if h==1, tx=0.02+B_col_x(h)+0.008; ha='left';
        else,    tx=0.02+B_col_x(h)+B_col_w(h)/2; ha='center'; end
        fw='normal';
        if h==6 && strcmp(rd{h},'< .001'), fw='bold'; end
        text(tx,ry+rh/2,rd{h},'Units','normalized','HorizontalAlignment',ha, ...
            'VerticalAlignment','middle','FontSize',10,'FontWeight',fw,'Color',[0.1 0.1 0.1]);
    end
end
rectangle('Position',[0.02 ty-hh-numel(B_rows)*rh 0.96 hh+numel(B_rows)*rh], ...
    'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5,'FaceColor','none');

figure(figB_full);
if saving_figures==1
    print(figB_full,'-dpng','-r300',sprintf('%s/figure_B_violin_table.png',figuredir));
    fprintf('Saved figure B with table\n');
end

%% ============================
%   FIGURE C: Hit/Miss Violin + Model C Table (presentation)
% ============================
figC_full = figure('Color','w','Units','inches','Position',[1 1 13 7]);

axC_new = axes('Units','normalized','Position',[0.08 0.40 0.80 0.53]);
hold on; grid on;

colors_c     = [0.9 0.4 0.4; 0.2 0.6 0.2];
edgeColors_c = [0.5 0.1 0.1; 0.1 0.3 0.1];
accLevels    = {'Miss','Hit'};

for i=1:nConds
    for aIdx=1:2
        data = all_trial_pooled_downsampled.ACD( ...
            all_trial_pooled_downsampled.Condition==gridLabels{i} & ...
            all_trial_pooled_downsampled.Acc==accLevels{aIdx});
        if ~isempty(data)
            [f,xi]=ksdensity(data);
            f=f/max(f)*0.38; side=(aIdx*2-3);
            fill(i+side*f,xi,colors_c(aIdx,:),'FaceAlpha',0.45,'EdgeColor',edgeColors_c(aIdx,:),'LineWidth',1.5);
            scatter(i+side*0.15+(rand(size(data))-0.5)*0.15,data,10,[0.2 0.2 0.2],'filled','MarkerFaceAlpha',0.3);
        end
    end
end

yline(0,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2);
xticks(1:nConds); xticklabels(gridLabels);
ylabel('\DeltaIBI (ms)','FontSize',13,'FontWeight','bold');
xlabel('Experimental Condition','FontSize',13,'FontWeight','bold');
title('\DeltaIBI by Accuracy and Condition','FontSize',13,'FontWeight','bold');
set(axC_new,'Box','off','TickDir','out','LineWidth',1.5,'FontSize',11);

annotation('rectangle',[0.845 0.87 0.02 0.04],'FaceColor',[0.2 0.6 0.2],'EdgeColor',[0.1 0.3 0.1],'LineWidth',1.5);
annotation('textbox',[0.868 0.865 0.08 0.05],'String','Hit','EdgeColor','none','FontSize',11,'FontWeight','bold','VerticalAlignment','middle');
annotation('rectangle',[0.845 0.82 0.02 0.04],'FaceColor',[0.9 0.4 0.4],'EdgeColor',[0.5 0.1 0.1],'LineWidth',1.5);
annotation('textbox',[0.868 0.815 0.08 0.05],'String','Miss','EdgeColor','none','FontSize',11,'FontWeight','bold','VerticalAlignment','middle');

axT_C = axes('Units','normalized','Position',[0 0.02 1 0.35]);
axis off; hold on;

C_headers = {'Term','\beta / F','SE','t','p','R^2'};
C_col_w   = [0.32 0.13 0.12 0.12 0.14 0.10];
C_col_x   = [0,cumsum(C_col_w(1:end-1))];
C_rows    = {
    {'Accuracy (Hit vs Miss)', sprintf('\\beta=%.3f',beta_Acc_C), sprintf('%.3f',se_Acc_C), sprintf('%.3f',t_Acc_C), ps_Acc_C,  sprintf('%.3f',r2_C)};
    {'Condition (ANOVA)',       sprintf('F=%.3f',f_cond_C),        '—',                       '—',                    ps_cond_C, '—'};
    {'Acc \times Condition',    sprintf('F=%.3f',f_int_C),         '—',                       '—',                    ps_int_C,  '—'}
};
C_row_bg = {[0.93 0.96 1.00],[0.86 0.92 0.98],[0.93 0.96 1.00]};
rh=0.20; hh=0.26; ty=0.88;

text(0.5,0.98,'Model C: \DeltaIBI ~ Accuracy \times Condition + (1|SubjectID)', ...
    'Units','normalized','HorizontalAlignment','center','FontSize',11,'FontWeight','bold','Color',[0.1 0.1 0.1]);

rectangle('Position',[0.02 ty-hh 0.96 hh],'FaceColor',[0.15 0.25 0.50],'EdgeColor','none');
for h=1:numel(C_headers)
    text(0.02+C_col_x(h)+C_col_w(h)/2,ty-hh/2,C_headers{h}, ...
        'Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontSize',10,'FontWeight','bold','Color','w');
end

for r=1:numel(C_rows)
    ry=ty-hh-r*rh; rd=C_rows{r};
    rectangle('Position',[0.02 ry 0.96 rh],'FaceColor',C_row_bg{r},'EdgeColor',[0.82 0.82 0.82]);
    for h=1:numel(C_headers)
        if h==1, tx=0.02+C_col_x(h)+0.008; ha='left';
        else,    tx=0.02+C_col_x(h)+C_col_w(h)/2; ha='center'; end
        fw='normal';
        if h==5 && strcmp(rd{h},'< .001'), fw='bold'; end
        text(tx,ry+rh/2,rd{h},'Units','normalized','HorizontalAlignment',ha, ...
            'VerticalAlignment','middle','FontSize',10,'FontWeight',fw,'Color',[0.1 0.1 0.1]);
    end
end
rectangle('Position',[0.02 ty-hh-numel(C_rows)*rh 0.96 hh+numel(C_rows)*rh], ...
    'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5,'FaceColor','none');

figure(figC_full);
if saving_figures==1
    print(figC_full,'-dpng','-r300',sprintf('%s/figure_C_violin_table.png',figuredir));
    fprintf('Saved figure C with table\n');
end

%% ============================
%   FIGURE E: DIBI across pre-stimulus R-peaks
% ============================
subjList_e = unique(delta_table.SubjectID);
nSubj_e    = numel(subjList_e);
nPeaks     = 3;
peak_labels = {'R:stim-2','R:stim-1','R:stim'};

subj_means = cell(nConds,1);
for c=1:nConds
    subj_means{c} = NaN(nSubj_e,nPeaks);
    for s=1:nSubj_e
        mask = strcmp(delta_table.SubjectID,subjList_e{s}) & ...
               delta_table.Condition==gridLabels{c};
        if sum(mask)>0
            subj_means{c}(s,1) = mean(delta_table.D_stim2(mask));
            subj_means{c}(s,2) = mean(delta_table.D_stim1(mask));
            subj_means{c}(s,3) = mean(delta_table.D_stim(mask));
        end
    end
end

figE = figure('Color','w','Units','inches','Position',[1 1 10 6]);
axE  = axes('Units','normalized','Position',[0.10 0.13 0.72 0.78]);
hold on; grid on;

x_pos = [1 2 3];

for c=1:nConds
    sm  = subj_means{c};
    col = condColors_e(c,:);
    for s=1:nSubj_e
        if ~any(isnan(sm(s,:)))
            plot(x_pos,sm(s,:),'-','Color',[col 0.25],'LineWidth',0.8);
        end
    end
    gm  = nanmean(sm,1);
    gse = nanstd(sm,0,1)./sqrt(sum(~isnan(sm),1));
    patch([x_pos fliplr(x_pos)],[gm+gse fliplr(gm-gse)],col,'FaceAlpha',0.15,'EdgeColor','none');
    plot(x_pos,gm,'-o','Color',col,'LineWidth',2.5,'MarkerSize',8, ...
        'MarkerFaceColor',col,'DisplayName',gridLabels{c});
end

yline(0,'--','Color',[0.4 0.4 0.4],'LineWidth',1.2);
xline(3.08,'--','Color',[0.15 0.15 0.15],'LineWidth',1.8);
text(3.11,0,'Stimulus Onset','FontSize',10,'Color',[0.15 0.15 0.15],'VerticalAlignment','middle');

xlim([0.7 3.5]);
xticks([1 2 3]); xticklabels(peak_labels);
xlabel('Pre-Stimulus R-Peak','FontSize',13,'FontWeight','bold');
ylabel('\DeltaIBI (ms)','FontSize',13,'FontWeight','bold');
title('Anticipatory Cardiac Deceleration across Pre-Stimulus Heartbeats','FontSize',12,'FontWeight','bold');
legend(gridLabels,'Location','northwest','FontSize',10,'Box','off');
set(axE,'Box','off','TickDir','out','LineWidth',1.5,'FontSize',11);

figure(figE);
if saving_figures==1
    print(figE,'-dpng','-r300',sprintf('%s/figure_E_ACD_slopes.png',figuredir));
    fprintf('Saved figure E\n');
end

fprintf('\nAll figures complete.\n');
