%% ============================
%   ACD_ARA_RT_MODEL.M  —  STANDALONE
%
%   MODEL: ARA ~ ACD * RT + ConditionVal + (1|SubjectID)
%
%   Tests whether RT moderates the ACD->ARA slope.
%   ACD x RT: does the slope of ARA on ACD change with response latency,
%   controlling for condition?
%
%   FIGURE: mirrors figure3 Model B from lme_all_models.m
%     Split violin per condition: Fast RT (left, blue) / Slow RT (right, orange)
%     Split defined by median RT within each condition.
%     Dashed LME fitted lines at RT mean-1SD (fast) and mean+1SD (slow).
%     M/SD labels at bottom.
%
%   All 7 exclusion criteria applied to both ACD and ARA simultaneously.
%   Only paired trials (valid ACD AND ARA) are included.
%   No downsampling — RT is continuous.
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
fast_col   = [0.20 0.45 0.75];   % blue  — fast RT
slow_col   = [0.90 0.55 0.10];   % orange — slow RT
fast_edg   = [0.10 0.25 0.50];
slow_edg   = [0.55 0.30 0.05];

%% ============================
%   STEP 1: DATA PROCESSING — PAIRED ACD + ARA
%   Same 7 exclusion criteria applied to both simultaneously.
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

%% Split RT into Fast / Slow by median per condition (for visualisation only)
tbl.RTGroup = categorical(repmat({'Slow'}, height(tbl), 1));
for i = 1:nConds
    rows   = tbl.Condition==gridLabels{i};
    med_rt = median(tbl.RT(rows), 'omitnan');
    tbl.RTGroup(rows & tbl.RT <= med_rt) = categorical({'Fast'});
end

fprintf('N paired total: %d\n\n', height(tbl));

%% ============================
%   STEP 2: FIT MODEL
%   ARA ~ ACD * RT + ConditionVal + (1|SubjectID)
%   Full table — no downsampling, RT is continuous
% ============================
fprintf('================================================================\n');
fprintf('  MODEL: ARA ~ ACD * RT + ConditionVal + (1|SubjectID)\n');
fprintf('  ACD x RT: does RT moderate the ACD->ARA slope?\n');
fprintf('================================================================\n');

lme  = fitlme(tbl, ...
    'ARA ~ ACD * RT + ConditionVal + (1|SubjectID)', 'FitMethod','ML');
anov = anova(lme, 'DFMethod','Satterthwaite');
disp(lme); fprintf('Type III ANOVA:\n'); disp(anov);
fprintf('R2=%.4f | AIC=%.2f | BIC=%.2f | LogLik=%.4f | N=%d\n\n', ...
    lme.Rsquared.Ordinary, lme.ModelCriterion.AIC, ...
    lme.ModelCriterion.BIC, lme.LogLikelihood, lme.NumObservations);

%% Extract coefficients
int_v   = lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'(Intercept)'));
b_ACD   = lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'ACD'));
b_RT    = lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'RT'));
b_Cond  = lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'ConditionVal'));
b_Int   = lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'ACD:RT'));
p_ACD   = lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'ACD'));
p_RT    = lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'RT'));
p_Int   = lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'ACD:RT'));
r2_v    = lme.Rsquared.Ordinary;

if p_ACD<0.001, sp_ACD='< .001'; else, sp_ACD=sprintf('= %.3f',p_ACD); end
if p_RT<0.001,  sp_RT='< .001';  else, sp_RT=sprintf('= %.3f',p_RT);   end
if p_Int<0.001, sp_Int='< .001'; else, sp_Int=sprintf('= %.3f',p_Int); end

mean_RT = mean(tbl.RT); sd_RT = std(tbl.RT);
rt_fast = mean_RT - sd_RT;   % fast RT probe value
rt_slow = mean_RT + sd_RT;   % slow RT probe value

% Slope of ARA on ACD at fast and slow RT
slope_fast = b_ACD + b_Int * rt_fast;
slope_slow = b_ACD + b_Int * rt_slow;

%% ============================
%   STEP 3: FIGURE — mirrors figure3 Model B from lme_all_models.m
%   Split violin per condition: Fast RT (left, blue) / Slow RT (right, orange)
%   Split defined by within-condition median RT
%   Dashed LME fitted lines at mean RT ± 1 SD
%   M/SD labels at bottom
% ============================
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 13.3 5.5]);
plot_b = 0.14; plot_h = 0.78; left_m = 0.07; right_m = 0.10;
uw     = 1 - left_m - right_m;
ax     = axes('Units','normalized','Position',[left_m, plot_b, uw, plot_h]);
hold on; grid on;

rtLevs  = {'Fast','Slow'};
rcols   = {fast_col; slow_col};
redgs   = {fast_edg; slow_edg};

cond_mu_fast = zeros(1,nConds); cond_sd_fast = zeros(1,nConds);
cond_mu_slow = zeros(1,nConds); cond_sd_slow = zeros(1,nConds);

for i = 1:nConds
    for rIdx = 1:2  % 1=Fast (left, side=-1), 2=Slow (right, side=+1)
        data = tbl.ARA( ...
            tbl.Condition==gridLabels{i} & tbl.RTGroup==rtLevs{rIdx});
        if numel(data) < 3, continue; end

        [f,xi] = ksdensity(data,'BoundaryCorrection','reflection');
        f    = f/max(f)*0.35;
        side = rIdx*2 - 3;  % Fast=-1, Slow=+1

        fill(i+side*f, xi, rcols{rIdx}, ...
            'FaceAlpha',0.45,'EdgeColor',redgs{rIdx},'LineWidth',1.1, ...
            'HandleVisibility','off');
        scatter(i+side*0.13+(rand(size(data))-0.5)*0.09, data, 5, ...
            [0.2 0.2 0.2],'filled','MarkerFaceAlpha',0.15,'HandleVisibility','off');

        mu_i = mean(data); sd_i = std(data);
        if rIdx==1, cond_mu_fast(i)=mu_i; cond_sd_fast(i)=sd_i;
        else,        cond_mu_slow(i)=mu_i; cond_sd_slow(i)=sd_i; end

        scatter(i+side*0.13, mu_i, 75, rcols{rIdx}, 'filled', ...
            'MarkerEdgeColor','k','LineWidth',1.2,'HandleVisibility','off');
    end
end

%% LME fitted lines: predicted ARA at mean ACD within RT group x condition cell
%  Mirrors ACD accuracy model: evaluate fitted value using cell mean ACD
fitted_fast = zeros(1,nConds);
fitted_slow = zeros(1,nConds);
for i = 1:nConds
    acd_fast = tbl.ACD(tbl.Condition==gridLabels{i} & tbl.RTGroup=='Fast');
    acd_slow = tbl.ACD(tbl.Condition==gridLabels{i} & tbl.RTGroup=='Slow');
    mu_acd_fast = mean(acd_fast(~isnan(acd_fast)));
    mu_acd_slow = mean(acd_slow(~isnan(acd_slow)));
    if isempty(mu_acd_fast), mu_acd_fast=0; end
    if isempty(mu_acd_slow), mu_acd_slow=0; end

    fitted_fast(i) = int_v + b_RT*rt_fast + b_Cond*condNumVals(i) + slope_fast*mu_acd_fast;
    fitted_slow(i) = int_v + b_RT*rt_slow + b_Cond*condNumVals(i) + slope_slow*mu_acd_slow;
end

plot(1:nConds, fitted_fast, '--','Color',fast_col,'LineWidth',2.1,'HandleVisibility','off');
plot(1:nConds, fitted_slow, '--','Color',slow_col,'LineWidth',2.1,'HandleVisibility','off');

%% Legend proxies
scatter(nan,nan,60,fast_col,'filled','MarkerEdgeColor','k','DisplayName','Fast RT (median split)');
scatter(nan,nan,60,slow_col,'filled','MarkerEdgeColor','k','DisplayName','Slow RT (median split)');

yline(0,'--','Color',[0.45 0.45 0.45],'LineWidth',1.0,'HandleVisibility','off');

%% M/SD labels at bottom
yl = ylim; ylim([yl(1)-0.22*diff(yl), yl(2)]); yl = ylim;
y_b  = yl(1) + 0.01*diff(yl);
dy_b = 0.11*diff(yl);
for i = 1:nConds
    text(i, y_b+dy_b, sprintf('F: %.1f (%.1f)', cond_mu_fast(i), cond_sd_fast(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',8,'Color',fast_col,'FontWeight','bold');
    text(i, y_b, sprintf('S: %.1f (%.1f)', cond_mu_slow(i), cond_sd_slow(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',8,'Color',slow_col,'FontWeight','bold');
end

%% Axes formatting
xticks(1:nConds); xticklabels(gridLabels); xlim([0.4 nConds+0.6]);
xlabel('Experimental Condition','FontSize',10,'FontWeight','bold');
ylabel('ARA (bpm)','FontSize',11,'FontWeight','bold');
title('ARA \sim ACD \times RT + Condition','FontSize',11,'FontWeight','bold');

text(0.03, 0.97, sprintf('\\beta_{ACD}=%.3f, p %s\n\\beta_{RT}=%.3f, p %s\n\\beta_{ACD\\timesRT}=%.3f, p %s\nR^2=%.3f', ...
    b_ACD,sp_ACD, b_RT,sp_RT, b_Int,sp_Int, r2_v), ...
    'Units','normalized','FontSize',9,'VerticalAlignment','top', ...
    'BackgroundColor','w','EdgeColor',[0.6 0.6 0.6]);

legend('Location','eastoutside','FontSize',9,'Box','off');
set(ax,'Box','off','TickDir','out','LineWidth',1.2,'FontSize',9);

exportgraphics(fig,'./figure_ARA_ACD_RT_model.png','Resolution',300);
fprintf('Saved figure_ARA_ACD_RT_model.png\n\nAll done.\n');
