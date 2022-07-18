clear all
close all

%---------------------------------------------------------------------------
% Pearon's corr. coef. (r)^2 = Coef. of determination, only for simple
% linear regression (only one predictor variable, or independent variable)
% https://stats.stackexchange.com/questions/83347/relationship-between-r2-and-correlation-coefficient
% https://towardsdatascience.com/r%C2%B2-or-r%C2%B2-when-to-use-what-4968eee68ed3
%---------------------------------------------------------------------------
% Drawn on 'static_analysis_reho.mat', which has 
% (i) ReHo, 
% (ii) IED rate mesured during long-term monitoring
% (iii) IED rate measured over entire scan, for each session of every
% subject
% This script makes scatterplot of (i) against (ii), and (i) against (iii)
% across all subjects
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter path where input struct. is stored at
% directname = ['/work/levan_lab/mtang/fmri_project/matrices', filesep 'static_analysis_reho'];   % direct. of output matrix
directname = 'C:\Users\siumichael.tang\Downloads\fmri_project\matrices\static_analysis_reho';
filename_input = 'static_analysis_reho.mat';   % filename of output file

% enter path where plots are saved at
% directname_op = '/work/levan_lab/mtang/fmri_project/plots/static_analysis_reho';
directname_op = 'C:\Users\siumichael.tang\Downloads\fmri_project\plots\static_analysis_reho';

% enter if user wants to write plots to file (1=yes, 0=no)
op_results = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% load structure
st = load(fullfile(directname, filename_input));
 
% set size of markers for scatterplot
sz = 12;

% get array of subjects in structure loaded
fd_name_sub = fieldnames(st.grand);

% f1 = figure;   % open figure window 
% f_grand = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
f_grand = figure('units','normalized','outerposition',[0 0 1 1]);

%-----------------------------------------------------------------------------
% (I): MAKE PLOT OF (Y) REHO AGAINST (X) IED RATE MEASUDRED DURING LONG-TERM 
% MONITORING FOR ALL SESSIONS

% initialize array for storing the required variables from all sessions
reho_all_sess = [];   % reho
ied_long_term_all_sess_log = [];   % ied_{long term} on log-scale

for sub_ind = 1:numel(fd_name_sub)   % for each subject
% for sub_ind = 1:1
    if ~isempty(st.grand.(fd_name_sub{sub_ind}))   % if value is not empty
        fd_run_cur = fieldnames(st.grand.(fd_name_sub{sub_ind}));   % get run fieldnames available
        for run_ind = 1:numel(fd_run_cur)   % for each run
            
            % access sub-structure of curr. sess.
            st_cur = st.grand.(fd_name_sub{sub_ind}).(fd_run_cur{run_ind});
            
            % get array of reho of every channel in curr. session from
            % sig_box_all, (rightmost col.)
            reho_cur = st_cur.sig_box_all(:, end);
            
            % concatenate vertically reho_cur to reho_all_sess 
            reho_all_sess = [reho_all_sess; cell2mat(reho_cur)];
            
            % get ied rate measured during long-term monitoring in curr.
            % sess., (rightmost col.)
            ied_long_term_cur = st_cur.ied_long_term(:, end);
             
            % concatenate vertically log10(ied_long_term_cur) to
            % ied_long_term_all_sess
            ied_long_term_all_sess_log = [ied_long_term_all_sess_log; log10(cell2mat(ied_long_term_cur))];
                   
        end
    end
end

% when applying log-scale on ied_{long term}, inf is obtained, filter out
% those entries, and obtain the corresponding entries in reho_all_sess
ind_req = ~isinf(ied_long_term_all_sess_log);   % get indices of non-inf entries in ied array
reho_all_sess_filt = reho_all_sess(ind_req);   % get corresponding entries in reho array
ied_long_term_all_sess_log_filt = ied_long_term_all_sess_log(ind_req);   % get non-inf entries in ied array

% after obtaining required arrays (non-inf entries) from all sessions, make plot
% f1 = figure;   % open figure window 
subplot(1, 2, 1);

%----------------------------------------
% make plot of (Y) ReHo against (X) IED_{LONG TERM} on log-scale
scatter(ied_long_term_all_sess_log_filt, reho_all_sess_filt, sz, 'filled', 'MarkerFaceColor', 'b'); hold on

%----------------------------------------
% add regression line to set of points:

% get coefficients of best-fit line
p_coef1 = polyfit(ied_long_term_all_sess_log_filt, reho_all_sess_filt, 1);

% get predicted y-values using coefficients above
y_pred1 = polyval(p_coef1, ied_long_term_all_sess_log_filt);

% add best-fit line in black to plot
plot(ied_long_term_all_sess_log_filt, y_pred1, 'k-', 'LineWidth', 1); hold on

%----------------------------------------
% add pearson coef (r) and its p-value to x-label:

% as best-fit line is a simple regression model, r^2 = R^2
[r_pear1, pval_pear1] = corr(ied_long_term_all_sess_log_filt, reho_all_sess_filt, 'Type', 'Pearson');

% format pearson r str, round value to 3 d.p.
r_pear_str1 = ['\(r = \) ', num2str(round(r_pear1, 3)), ', '];

% if pval < 0.001, then show 'p < 0.001' on x-label, otherwise, show value
% corrected to 1 sig. fig.
if pval_pear1 < 0.001
    pval_pear_str1 = '\(p < 0.001\)';
else
    pval_pear_str1 = ['\(p = \) ', num2str(round(pval_pear1, 1, 'significant'))];
end

% format str displaying stat. sig. btw. ReHo and IED_{long term}
stat_fig_str1 = [r_pear_str1, pval_pear_str1];
%----------------------------------------

box on;
% grid on;

% set x- and y-labels, 
xlabel(['IED\(_{\mathrm{long}}\)\(_{ }\)\(_{\mathrm{term}}\) on log scale [\(\mathrm{min}^{\mathrm{-1}}\)], ', ...
    stat_fig_str1], 'Interpreter','latex');
ylabel('ReHo','Interpreter','latex');

% set axis fontsize (default = 10)
ax1 = gca;
ax1.FontSize = 12;
 
% set interpreter of tick labels on both axes 
ax1.TickLabelInterpreter = 'latex';
            
% get y-axis min and max
ymin1 = ax1.YLim(1);
ymax1 = ax1.YLim(2);

% % get x-axis min and max
xmin1 = min(ied_long_term_all_sess_log_filt);
xmax1 = max(ied_long_term_all_sess_log_filt);

ax1.XTick = [floor(xmin1):1:ceil(xmax1)];

% set aspect ratio of plot
pbr = 1;
pbaspect([xmax1-xmin1, pbr*(xmax1-xmin1), 1]); % multiple y-axis by the factor

% place string on plot
text(0.93, 0.95, '(a)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 15);   
       
% END (I): MAKE PLOT OF (Y) REHO AGAINST (X) IED RATE MEASUDRED DURING LONG-TERM 
% MONITORING FOR ALL SESSIONS
%-----------------------------------------------------------------------------
% (II): MAKE PLOT OF (Y) REHO AGAINST (X) IED RATE MEASUDRED DURING SCAN FOR ALL SESSIONS

% initialize array for storing the required variables from all sessions
reho_ch_matched_all_sess = [];   % reho of channels with ied during scan registered
ied_during_scan_all_sess_log = [];   % ied_{during scan} on log-scale

for sub_ind = 1:numel(fd_name_sub)   % for each subject
% for sub_ind = 1:1
    if ~isempty(st.grand.(fd_name_sub{sub_ind}))   % if value is not empty
        fd_run_cur = fieldnames(st.grand.(fd_name_sub{sub_ind}));   % get run fieldnames available
        for run_ind = 1:numel(fd_run_cur)   % for each run
            
            % access sub-structure of curr. sess.
            st_cur = st.grand.(fd_name_sub{sub_ind}).(fd_run_cur{run_ind});
            
            % get array of reho of matched channel with ied_{during scan} recorded
            % in curr. session from ch_name_ied_reho (rightmost col.)
            reho_ch_matched_cur = st_cur.ch_name_matched_ied_reho(:, end);
            
            % concatenate vertically reho_ch_matched_cur to reho_ch_matched_all_sess 
            reho_ch_matched_all_sess = [reho_ch_matched_all_sess; cell2mat(reho_ch_matched_cur)];
            
            % get ied rate measured during scan in curr. sess. (2nd
            % rightmost col. of ch_name_matched_ied_reho)
            ied_during_scan_cur = st_cur.ch_name_matched_ied_reho(:, end-1);
             
            % concatenate vertically log(ied_during_scan_cur) to
            % ied_during_scan_all_sess_log
            ied_during_scan_all_sess_log = [ied_during_scan_all_sess_log; log10(cell2mat(ied_during_scan_cur))];
                   
        end
    end
end

% after obtaining required arrays (non-inf entries) from all sessions, make plot
% f2 = figure;   % open figure window 
subplot(1, 2, 2);

%----------------------------------------
% make plot of (Y) ReHo against (X) IED_{LONG TERM} on log-scale
scatter(ied_during_scan_all_sess_log, reho_ch_matched_all_sess, sz, 'filled', 'MarkerFaceColor', 'b'); hold on

%----------------------------------------
% add regression line to set of points:

% get coefficients of best-fit line
p_coef2 = polyfit(ied_during_scan_all_sess_log, reho_ch_matched_all_sess, 1);

% get predicted y-values using coefficients above
y_pred2 = polyval(p_coef2, ied_during_scan_all_sess_log);

% add best-fit line in black to plot
plot(ied_during_scan_all_sess_log, y_pred2, 'k-', 'LineWidth', 1); hold on

%----------------------------------------
% add pearson coef (r) and its p-value to x-label:

% as best-fit line is a simple regression model, r^2 = R^2
[r_pear2, pval_pear2] = corr(ied_during_scan_all_sess_log, reho_ch_matched_all_sess, 'Type', 'Pearson');

% format pearson r str, round value to 3 d.p.
r_pear_str2 = ['\(r = \) ', num2str(round(r_pear2, 3)), ', '];

% if pval < 0.001, then show 'p < 0.001' on x-label, otherwise, show value
% corrected to 1 sig. fig.
if pval_pear2 < 0.001
    pval_pear_str2 = '\(p < 0.001\)';
else
    pval_pear_str2 = ['\(p = \) ', num2str(round(pval_pear2, 1, 'significant'))];
end

% format str displaying stat. sig. btw. ReHo and IED_{during scan}
stat_fig_str2 = [r_pear_str2, pval_pear_str2];
%----------------------------------------

box on;
% grid on;

% set x- and y-labels, 
xlabel(['IED\(_{\mathrm{during}}\)\(_{ }\)\(_{\mathrm{scan}}\) on log scale [\(\mathrm{min}^{\mathrm{-1}}\)], ', ...
    stat_fig_str2], 'Interpreter','latex');
ylabel('ReHo','Interpreter','latex');

% set axis fontsize (default = 10)
ax2 = gca;
ax2.FontSize = 12;
 
% set interpreter of tick labels on both axes 
ax2.TickLabelInterpreter = 'latex';
            
% get y-axis min and max
ymin2 = ax2.YLim(1);
ymax2 = ax2.YLim(2);

% % get x-axis min and max
xmin2 = min(ied_during_scan_all_sess_log);
xmax2 = max(ied_during_scan_all_sess_log);

% format y-axis limits
ax2.YLim = [0, 1];

% set aspect ratio of plot
pbr = 1;
pbaspect([xmax2-xmin2, pbr*(xmax2-xmin2), 1]); % multiple y-axis by the factor

% place string on plot
text(0.93, 0.95, '(b)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 15);   

% END (II): MAKE PLOT OF (Y) REHO AGAINST (X) IED RATE MEASUDRED DURING SCAN FOR ALL SESSIONS
%-----------------------------------------------------------------------------

%----------------------------------------
% adjust positions of pannels:

% % hdiff = 0.001/10000;   % set horizontal diff. btw. panels
hdiff = (1/2.2)*(ax2.Position(1) - (ax1.Position(1) + ax1.Position(3)));

% set pannel positions
ax2.Position(1) = ax1.Position(1) + ax1.Position(3) + hdiff;
ax2.Position(1) = ax1.Position(1) + ax1.Position(3) + hdiff;

%----------------------------------------

% adjust y-axis limits on both figures, keep them the same
ymin_all = floor(min([ymin1, ymin2]));   % get global ymin
ymax_all = ceil(max([ymax1, ymax2]));   % get global ymax

% set y-axis limits, ticks, and labels on ax1
ax1.YLim = [ymin_all, ymax_all];   % limits of y-axis
ax2.YLim = [ymin_all, ymax_all];   % limits of y-axis

% adjust x-axis limits on both figures, keep them the same
xmin_all = floor(min([xmin1, xmin2]));   % get global xmin
xmax_all = ceil(max([xmax1, xmax2]));   % get global xmax

% format exponents on xlabels 
% though log scale is applied on x data, 
% set x-axis labels on ax1 as [0.01, 0.1, 1, 10, 100, ...] 
xlabel_expo = [xmin_all:1:xmax_all];   % format exponents

% set x-axis limits, ticks, and labels on ax1
ax1.XLim = [xmin_all, xmax_all];   % limits of x-axis
ax1.XTick = xmin_all:1:xmax_all;   % set xticks on ax1
ax1.XTickLabel = 10.^(xlabel_expo);   % format xtick labels

% set x-axis limits, ticks, and labels on ax2
ax2.XLim = [xmin_all, xmax_all];
ax2.XTick = xmin_all:1:xmax_all;   % set xticks on ax2
ax2.XTickLabel = 10.^(xlabel_expo);   % format xtick labels

% set aspect ratio on both ax1 and ax2
ax1.PlotBoxAspectRatio = [xmax_all-xmin_all, pbr*(xmax_all-xmin_all), 1];
ax2.PlotBoxAspectRatio = [xmax_all-xmin_all, pbr*(xmax_all-xmin_all), 1];

%----------------------------------------
% format filename for current plot
filename_op = ['static_analysis_reho'];
            
% save figure to file 
if op_results == 1
    saveas(f_grand, fullfile(directname_op, filename_op), 'epsc');
    saveas(f_grand, fullfile(directname_op, filename_op), 'jpeg');
end

