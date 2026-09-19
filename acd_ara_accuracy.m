%% ============================
%   ACD_ARA_ACCURACY.M  —  STANDALONE
%
%   Computes Pearson r between ACD and ARA, split by accuracy (Hit / Miss)
%   within each condition. Trials pooled across subjects per accuracy x
%   condition cell.
%
%   ACD = HR(stim-3) - HR(stim)   [bpm]
%   ARA = HR(stim+3) - HR(stim+1) [bpm]
%
%   All 7 exclusion criteria applied to both measures simultaneously.
%   Only paired trials (valid ACD AND ARA) are included.
%
%   OUTPUT:
%     Console: r, p, 95% CI, N per accuracy x condition cell
%     Figure A: scatter panels per condition, Hit vs Miss overlaid
%     Figure B: r bar chart per condition, Hit vs Miss side by side
%
%   Requires: organizedData in workspace
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

%% ============================
%   STEP 1: COLLECT PAIRED ACD + ARA
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

%% Criterion 7: z±3 on ACD and ARA separately per condition
for i = 1:nConds
    rows=tbl.Condition==gridLabels{i};
    for fld={'ACD','ARA'}
        v=tbl.(fld{1})(rows); mu=mean(v,'omitnan'); sd=std(v,'omitnan');
        if sd>0, v(abs((v-mu)/sd)>3)=NaN; tbl.(fld{1})(rows)=v; end
    end
end
tbl = tbl(~isnan(tbl.ACD) & ~isnan(tbl.ARA), :);
fprintf('Paired trials after all criteria: %d\n\n', height(tbl));

%% RT units
if median(tbl.RT(tbl.RT>0)) < 10, tbl.RT = tbl.RT*1000; end

%% ============================
%   STEP 2: CORRELATION PER ACCURACY x CONDITION CELL
% ============================
accLevels = {'Hit','Miss'};
nAcc      = 2;

% Storage: rows=conditions, cols=Hit/Miss
r_mat  = nan(nConds, nAcc);
p_mat  = nan(nConds, nAcc);
n_mat  = zeros(nConds, nAcc);
ci_mat = nan(nConds, nAcc, 2);  % dim3: lo/hi

fprintf('=== ACD–ARA CORRELATION BY ACCURACY x CONDITION ===\n');
fprintf('%-10s %-6s %8s %10s %10s %8s\n','Condition','Acc','N','r','p','95%%CI');
fprintf('%s\n',repmat('-',1,58));

for i = 1:nConds
    for aIdx = 1:nAcc
        rows  = tbl.Condition==gridLabels{i} & tbl.Acc==accLevels{aIdx};
        acd_v = tbl.ACD(rows);
        ara_v = tbl.ARA(rows);
        valid = ~isnan(acd_v) & ~isnan(ara_v);
        acd_v = acd_v(valid); ara_v = ara_v(valid);
        n     = numel(acd_v);
        n_mat(i,aIdx) = n;

        if n < 5
            fprintf('%-10s %-6s %8d  (insufficient)\n', gridLabels{i}, accLevels{aIdx}, n);
            continue;
        end

        [r,p] = corr(acd_v, ara_v, 'Type','Pearson');
        r_mat(i,aIdx) = r; p_mat(i,aIdx) = p;

        z    = atanh(r); se = 1/sqrt(n-3);
        ci   = tanh(z + [-1 1]*1.96*se);
        ci_mat(i,aIdx,:) = ci;

        if p<0.001, ps='< .001'; else, ps=sprintf('= %.3f',p); end
        fprintf('%-10s %-6s %8d %10.4f %10s  [%.3f, %.3f]\n', ...
            gridLabels{i}, accLevels{aIdx}, n, r, ps, ci(1), ci(2));
    end
end
fprintf('%s\n\n',repmat('-',1,58));

%% ============================
%   FIGURES
% ============================

%% Figure A — Scatter panels per condition, Hit vs Miss overlaid
fig1 = figure('Color','w','Units','inches','Position',[0.5 0.5 13.3 5.5]);
panel_b=0.14; panel_h=0.76;
left_m=0.07; right_m=0.02; uw=1-left_m-right_m;
gap_c=0.018; pw=(uw-(nConds-1)*gap_c)/nConds;
px=left_m+(0:nConds-1).*(pw+gap_c);

acols = {hit_col, miss_col};
alabels = {'Hit','Miss'};

for i = 1:nConds
    ax = axes('Units','normalized','Position',[px(i),panel_b,pw,panel_h]);
    hold on; grid on;

    for aIdx = 1:nAcc
        rows  = tbl.Condition==gridLabels{i} & tbl.Acc==accLevels{aIdx};
        acd_v = tbl.ACD(rows); ara_v = tbl.ARA(rows);
        valid = ~isnan(acd_v)&~isnan(ara_v);
        acd_v = acd_v(valid); ara_v = ara_v(valid);
        if numel(acd_v) < 3, continue; end

        scatter(acd_v, ara_v, 5, acols{aIdx}, 'filled', ...
            'MarkerFaceAlpha', 0.20, 'MarkerEdgeColor','none', 'HandleVisibility','off');

        if numel(acd_v) >= 5
            p_fit = polyfit(acd_v, ara_v, 1);
            x_fit = linspace(min(acd_v), max(acd_v), 100);
            plot(x_fit, polyval(p_fit,x_fit), '-', 'Color', acols{aIdx}, ...
                'LineWidth', 1.8, 'HandleVisibility','off');
        end
    end

    % r annotations for Hit and Miss
    yl = ylim; xl = xlim;
    for aIdx = 1:nAcc
        r = r_mat(i,aIdx); p = p_mat(i,aIdx);
        if isnan(r), continue; end
        if p<0.001, ps='p<.001'; else, ps=sprintf('p=%.3f',p); end
        rstr = sprintf('%s: r=%.2f (%s)', alabels{aIdx}(1), r, ps);
        ypos = yl(2) - (aIdx-1)*0.12*diff(yl);
        text(xl(1)+0.03*diff(xl), ypos, rstr, ...
            'FontSize', 7, 'Color', acols{aIdx}*0.75, 'FontWeight','bold', ...
            'VerticalAlignment','top');
    end

    xline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',0.7,'HandleVisibility','off');
    yline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',0.7,'HandleVisibility','off');
    title(gridLabels{i},'FontSize',10,'FontWeight','bold');
    xlabel('ACD (bpm)','FontSize',8);
    if i==1, ylabel('ARA (bpm)','FontSize',8.5); end
    set(ax,'Box','off','TickDir','out','LineWidth',1.0,'FontSize',8);
end

% Legend
scatter(nan,nan,30,hit_col,'filled','DisplayName','Hit');
scatter(nan,nan,30,miss_col,'filled','DisplayName','Miss');
legend('Location','eastoutside','FontSize',9,'Box','off');
annotation('textbox',[left_m+0.01,panel_b+panel_h+0.03,0.6,0.05], ...
    'String','ACD–ARA Correlation: Hit vs Miss per Condition', ...
    'EdgeColor','none','FontSize',11,'FontWeight','bold');
exportgraphics(fig1,'./acd_ara_acc_figure1_scatter.png','Resolution',300);
fprintf('Saved acd_ara_acc_figure1_scatter.png\n');

%% Figure B — r bar chart per condition, Hit vs Miss side by side
fig2 = figure('Color','w','Units','inches','Position',[0.5 0.5 13.3 5.5]);
ax2  = axes('Units','normalized','Position',[0.08, 0.14, 0.87, 0.76]);
hold on; grid on;

bar_w   = 0.32;
offsets = [-0.22, 0.22];  % Hit left, Miss right within each condition

for aIdx = 1:nAcc
    for i = 1:nConds
        r = r_mat(i,aIdx);
        if isnan(r), continue; end
        x = i + offsets(aIdx);
        b = bar(x, r, bar_w, 'FaceColor', acols{aIdx}, ...
            'EdgeColor', acols{aIdx}*0.6, 'LineWidth', 0.8, 'HandleVisibility','off');

        % CI error bar
        lo = ci_mat(i,aIdx,1); hi = ci_mat(i,aIdx,2);
        errorbar(x, r, r-lo, hi-r, 'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');

        % Significance marker
        if ~isnan(p_mat(i,aIdx)) && p_mat(i,aIdx) < 0.05
            yl = ylim;
            text(x, r+sign(r)*0.025*diff(yl), '*', ...
                'HorizontalAlignment','center','FontSize',13,'Color',[0.1 0.1 0.1]);
        end
    end
end

% Legend proxies
bar(nan, nan, 'FaceColor', hit_col,  'EdgeColor', hit_col*0.6,  'DisplayName','Hit');
bar(nan, nan, 'FaceColor', miss_col, 'EdgeColor', miss_col*0.6, 'DisplayName','Miss');

yline(0,'--','Color',[0.4 0.4 0.4],'LineWidth',1.0,'HandleVisibility','off');
xticks(1:nConds); xticklabels(gridLabels); xlim([0.4 nConds+0.6]);
xlabel('Experimental Condition','FontSize',11,'FontWeight','bold');
ylabel('Pearson r  (ACD ~ ARA)','FontSize',11,'FontWeight','bold');
title('ACD–ARA Correlation by Accuracy and Condition','FontSize',12,'FontWeight','bold');
legend('Location','best','FontSize',10,'Box','off');

% r and N labels inside/below bars
for aIdx = 1:nAcc
    for i = 1:nConds
        r = r_mat(i,aIdx); n = n_mat(i,aIdx);
        if isnan(r), continue; end
        x = i + offsets(aIdx);
        yl = ylim;
        text(x, yl(1)+0.01*diff(yl), sprintf('r=%.2f\nn=%d',r,n), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',6,'Color',[0.2 0.2 0.2]);
    end
end

set(ax2,'Box','off','TickDir','out','LineWidth',1.2,'FontSize',9);
exportgraphics(fig2,'./acd_ara_acc_figure2_r_barchart.png','Resolution',300);
fprintf('Saved acd_ara_acc_figure2_r_barchart.png\n');

fprintf('\nAll done.\n');
