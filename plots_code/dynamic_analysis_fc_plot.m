close all
clear all

% This script requires 'dynamic_analysis_fc.mat' as input, and generate a
% plot of (Y) dynamic corr. of fc vs inst. ied rate against (X) coeff. of
% var. of inst. ied rate, for (left) ied onset-ied onset, (right) ied
% onset-ied propagation, channels pairs.
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

tic 

% enter path to directory where info. of all subjects are stored
subnum_dir = 'C:\Users\siumichael.tang\Downloads\fmri_project'; 
% subnum_dir = '/work/levan_lab/mtang/fmri_project/';

% format path where ouput struct. is stored at
fname_input_fc = [subnum_dir, filesep, 'matrices' filesep 'dynamic_analysis_fc']; 
filename_input_fc = 'dynamic_analysis_fc_fisher_z.mat';   % filename of fc file

% enter path where plots are saved at
% directname_op = '/work/levan_lab/mtang/fmri_project/plots/dynamic_analysis_fc_plot';
directname_op = 'C:\Users\siumichael.tang\Downloads\fmri_project\paper_plots';

% format output filename
filenameg_op = ['dynamic_analysis_fc_plot'];

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% obtain full path of input .mat file
file_path_fc = fullfile(fname_input_fc, filename_input_fc);

% load .mat file
fc_st = load(file_path_fc);

% use collect_ar, defined at the end of script, to combine array from 
% each subject together to obtain single layer array
fc_ar_all_onset_onset = collect_ar(fc_st.main, 'fc_ar_comb_onset_onset');
fc_ar_all_onset_propa = collect_ar(fc_st.main, 'fc_ar_comb_onset_propa');

%----------------------------------------------------
% make plot

% figure
f1 = figure('units','normalized','outerposition',[0 0 1 1]);

% specify circle size
sz = 15;

%-------------------------------------------------
% make plot for ied onset-ied onset channel pairs
sp1 = subplot(1, 2, 1);

coef_var_onset_onset = cell2mat(fc_ar_all_onset_onset(:, end-1));
fisher_z_fc_onset_onset = cell2mat(fc_ar_all_onset_onset(:, end));

% add scatterplot (X): coeff. of var., (Y): fisher z-score of fc vs ied
% rate for ied onset-ied onset channels
scatter(coef_var_onset_onset, fisher_z_fc_onset_onset, sz, 'filled', ...
    'MarkerFaceColor', [0 0 0.7]); hold on

% call get_regression_line to obtain predicted y data, pearson's corr.
% coeff., and its p-value based on (X) and (Y)
[y_pred_onset_onset, pear_rho_onset_onset, ...
    p_val_onset_onset] = get_regression_line(coef_var_onset_onset, ...
    fisher_z_fc_onset_onset);

% plot regression line
plot(coef_var_onset_onset, y_pred_onset_onset, '-k', 'LineWidth', 1.15); hold off

box on

% set x- and y-labels, 
xlabel_rho_1 = round(pear_rho_onset_onset, 3);   % corr. to 3 d. p.
xlabel_p_val_1 = round(p_val_onset_onset, 1, 'significant');   % corr. to 1 sig. fig.

if xlabel_p_val_1 >= 0.001
    xlabel(sprintf('\\((\\sigma/\\mu)_{\\mathrm{IED onset-IED onset}}\\), \\(r = %.3f, p = %.2f\\)', ...
        xlabel_rho_1, xlabel_p_val_1), 'Interpreter','latex');
else
    xlabel(sprintf('\\((\\sigma/\\mu)_{\\mathrm{IED onset-IED onset}}\\), \\(r = %.3f, p < 0.001\\)', ...
        xlabel_rho_1), 'Interpreter','latex');
end

ylabel('Dynamic correlation of FC','Interpreter','latex');

% format title:
fc_mean_fisher_z_onset_onset = mean(fisher_z_fc_onset_onset);   % get mean Fisher'z
[~, fc_mean_fisher_z_pval_onset_onset] = ttest(fisher_z_fc_onset_onset);   % get p-value of fisher z
fc_mean_fisher_pval_disp_onset_onset = round(fc_mean_fisher_z_pval_onset_onset, 1, 'significant');   % round numbers
t1 = title(sprintf('FC of IED onset-IED onset (\\(\\bar{z} = %.3f, p = %.1f\\))', fc_mean_fisher_z_onset_onset, ...
    fc_mean_fisher_pval_disp_onset_onset),'Interpreter','latex');

% get current gca
ax1 = gca;

% set line width of axes
ax1.LineWidth = 0.9;

% set x and y limits
xmin = 0;
xmax = 1.5;
xtick = 0.5;
    
ymin = -1.2;
ymax = 1;
ytick = 0.2;
    
% format axes
ax1.XLim = [xmin,xmax];   % x-limits
ax1.XTick = [xmin:xtick:xmax];
    
ax1.YLim = [ymin,ymax];   % y-limits
ax1.YTick = [ymin:ytick:ymax];
   
% set axis fontsize (default = 10)
ax1.FontSize = 10;   % axes font size
ax1.YLabel.FontSize = 16;   % y-label font size
ax1.XLabel.FontSize = 14;   % x-label font size
t1.FontSize = 14;   % title font size

% set properties of ticks
ax1.TickDir = 'out';   % direction of ticks
ax1.YMinorTick = 'on';   % show minor y-ticks
% ax1.TickLength = 2*ax1.TickLength;   % set length of ticks
ax1.TickLength = 2*[ax1.TickLength(1) 4*ax1.TickLength(1)];   % set length of ticks

% set interpreter of tick labels on both axes 
ax1.TickLabelInterpreter = 'latex';
    
% set aspect ratio of plot
pbr = 1;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

% place string on plot
text(0.92, 0.96, '(a)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 18);   

% END make plot for ied onset-ied onset channel pairs
%-------------------------------------------------

% make plot for ied onset-ied propa channel pairs
sp2 = subplot(1, 2, 2);

coef_var_onset_propa = cell2mat(fc_ar_all_onset_propa(:, end-1));
fisher_z_fc_onset_propa = cell2mat(fc_ar_all_onset_propa(:, end));

% add scatterplot (X): coeff. of var., (Y): fisher z-score of fc vs ied
% rate for ied onset-ied propa channels
scatter(coef_var_onset_propa, fisher_z_fc_onset_propa, sz, 'filled', ...
    'MarkerFaceColor', [0 0 0.7]); hold on

% call get_regression_line to obtain predicted y data, pearson's corr.
% coeff., and its p-value based on (X) and (Y)
[y_pred_onset_propa, pear_rho_onset_propa, ...
    p_val_onset_propa] = get_regression_line(coef_var_onset_propa, ...
    fisher_z_fc_onset_propa);

% plot regression line
plot(coef_var_onset_propa, y_pred_onset_propa, '-k', 'LineWidth', 1.15); hold off

box on

% set x- and y-labels, 
xlabel_rho_2 = round(pear_rho_onset_propa, 3);   % corr. to 3 d. p.
xlabel_p_val_2 = round(p_val_onset_propa, 1, 'significant');   % corr. to 1 sig. fig.

if xlabel_p_val_2 >= 0.001
    xlabel(sprintf('\\((\\sigma/\\mu)_{\\mathrm{IED onset-IED propa}}\\), \\(r = %.3f, p = %.2f\\)', ...
        xlabel_rho_2, xlabel_p_val_2), 'Interpreter','latex');
else
    xlabel(sprintf('\\((\\sigma/\\mu)_{\\mathrm{IED onset-IED propa}}\\), \\(r = %.3f, p < 0.001\\)', ...
        xlabel_rho_2), 'Interpreter','latex');
end

% format title:
fc_mean_fisher_z_onset_propa = mean(fisher_z_fc_onset_propa);   % get mean Fisher'z
[~, fc_mean_fisher_z_pval_onset_propa] = ttest(fisher_z_fc_onset_propa);   % get p-value of fisher z
fc_mean_fisher_pval_disp_onset_propa = round(fc_mean_fisher_z_pval_onset_propa, 1, 'significant');   % round numbers
t2 = title(sprintf('FC of IED onset-IED propa (\\(\\bar{z} = %.3f, p = %.3f\\))', fc_mean_fisher_z_onset_propa, ...
    fc_mean_fisher_pval_disp_onset_propa),'Interpreter','latex');

% get current gca
ax2 = gca;

% set line width of axes
ax2.LineWidth = 0.9;

% set x and y limits
xmin = 0;
xmax = 1.5;
xtick = 0.5;
    
ymin = -1.2;
ymax = 1;
ytick = 0.2;
    
% format axes
ax2.XLim = [xmin,xmax];   % x-limits
ax2.XTick = [xmin:xtick:xmax];
    
ax2.YLim = [ymin,ymax];   % y-limits
ax2.YTick = [ymin:ytick:ymax];
  
% set axis fontsize (default = 10)
ax2.FontSize = 10;   % axes font size
ax2.YLabel.FontSize = 16;   % y-label font size
ax2.XLabel.FontSize = 14;   % x-label font size
t2.FontSize = 14;   % title font size

% set properties of ticks
ax2.TickDir = 'out';   % direction of ticks
ax2.YMinorTick = 'on';   % show minor y-ticks
% ax2.TickLength = 2*ax2.TickLength;   % set length of ticks
ax2.TickLength = 2*[ax2.TickLength(1) 4*ax2.TickLength(1)];   % set length of ticks

% set interpreter of tick labels on both axes 
ax2.TickLabelInterpreter = 'latex';
    
% set aspect ratio of plot
pbr = 1;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

% place string on plot
text(0.92, 0.96, '(b)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 18);   

% END make plot for ied onset-ied propa channel pairs
%-------------------------------------------------
% adjust positions of subplots
hdiff = 0.1/2.5;
ax2.Position(1) = ax1.Position(1) + ax1.Position(3) + hdiff;
ax2.Position(2) = ax1.Position(2);

% end adjusting positions of scatterplots
%--------------------------------------------------

% output plots to path
if op_results == 1
saveas(f1, fullfile(directname_op, filenameg_op), 'epsc');
saveas(f1, fullfile(directname_op, filenameg_op), 'jpeg');
end


% END MAKING PLOTS 
%--------------------------------------------------

% This function intakes input struct (input_st) and name of array interested
% (ar_int_str), and generates output array which concatenates vertically array
% interested across subjects
function output_array = collect_ar(input_st, ar_int_str)

% initialize output
output_array = [];

% get fieldnames under input_st (all subjects)
fd_sub = fieldnames(input_st);

% vertically concatenate array interested across subjects
ar_int_all = [];   % initialize array, ar_int_all
for i = 1:numel(fd_sub)   % for each subject
    % obtain array interested for curr. subject
    ar_int_curr = input_st.(fd_sub{i}).(ar_int_str);   
    ar_int_all = [ar_int_all; ar_int_curr];   % vert. concat. arrays
end

% assign value to output array
output_array = ar_int_all;

end

%---------------------------------------

% This function requires x and y data as inputs, then generates the
% following outputs, 
% y_pred = predicted y-value based on x data in regression analysis
% pear_rho = pearson's corr. coeff. btw. x and y data
% p_val = p-value of pearson's corr. coeff.
function [y_pred, pear_rho, p_val] = get_regression_line(x_data, y_data)

% initialize output variables
y_pred = [];   % predicted y-values for plotting
pear_rho = [];   % pearson's corr. coeff. btw. x and y data
p_val = [];   % p-value of pearson's corr. coeff.

% get coefficient of straight line fitted to set of data points
% using polyfit with degree of 1:
% coef. found by polyfit are based on least-square principle
p_coef = polyfit(x_data, y_data, 1);

% get predicted y-values from coefficients obtained 
y_pred = polyval(p_coef, x_data);

% get pearson's corr. coeff. btw. x and y data, and its p-value
[pear_rho, p_val] = corr(x_data, y_data, ...
    'Type', 'Pearson');

end

%---------------------------------------