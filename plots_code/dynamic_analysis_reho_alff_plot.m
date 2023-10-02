close all
clear all

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

tic 

% enter path to directory where info. of all subjects are stored
subnum_dir = 'C:\Users\siumichael.tang\Downloads\fmri_project'; 
% subnum_dir = '/work/levan_lab/mtang/fmri_project/';

% format path where ouput struct. is stored at
fname_input_reho = [subnum_dir, filesep, 'matrices' filesep 'dynamic_analysis_reho']; 
filename_input_reho = 'dynamic_analysis_reho_fisher_z.mat';   % filename of reho file

fname_input_alff = [subnum_dir, filesep, 'matrices' filesep 'dynamic_analysis_alff']; 
filename_input_alff = 'dynamic_analysis_alff_fisher_z.mat';   % filename of alff file

% enter path where plots are saved at
% directname_op = '/work/levan_lab/mtang/fmri_project/plots/dynamic_analysis_reho_alff_plot';
directname_op = 'C:\Users\siumichael.tang\Downloads\fmri_project\paper_plots';

% format output filename
filenameg_op = ['dynamic_analysis_reho_alff_plot'];

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% obtain full path of input .mat file
file_path_reho = fullfile(fname_input_reho, filename_input_reho);
file_path_alff = fullfile(fname_input_alff, filename_input_alff);

% load .mat file
reho_st = load(file_path_reho);
alff_st = load(file_path_alff);

% use collect_ar, defined at the end of script, to combine array from 
% each subject together to obtain single layer array
ch_onset_ds_ied_reho_all = collect_ar(reho_st.main, 'ch_onset_ds_ied_reho_comb');
ch_propa_ds_ied_reho_all = collect_ar(reho_st.main, 'ch_propa_ds_ied_reho_comb');
ch_ext_ds_ied_reho_all = collect_ar(reho_st.main, 'ch_ext_ds_ied_reho_comb');

ch_onset_ds_ied_alff_all = collect_ar(alff_st.main, 'ch_onset_ds_ied_alff_comb');
ch_propa_ds_ied_alff_all = collect_ar(alff_st.main, 'ch_propa_ds_ied_alff_comb');
ch_ext_ds_ied_alff_all = collect_ar(alff_st.main, 'ch_ext_ds_ied_alff_comb');

% obtain fisher z-score and coeff. of variation (end-1 and end col.) from
% each subject and store them in single layer array, with entries of nan 
% and inf removed, for each channel category, i.e. onset, propa, and ext =
% onset+propa
% for reho, 
[ch_onset_coef_var_reho, ch_onset_fisher_z_reho] = ...
    remove_nan_inf(ch_onset_ds_ied_reho_all(:, end-1), ch_onset_ds_ied_reho_all(:, end));
[ch_propa_coef_var_reho, ch_propa_fisher_z_reho] = ...
    remove_nan_inf(ch_propa_ds_ied_reho_all(:, end-1), ch_propa_ds_ied_reho_all(:, end));
[ch_ext_coef_var_reho, ch_ext_fisher_z_reho] = ...
    remove_nan_inf(ch_ext_ds_ied_reho_all(:, end-1), ch_ext_ds_ied_reho_all(:, end));

% for alff, 
[ch_onset_coef_var_alff, ch_onset_fisher_z_alff] = ...
    remove_nan_inf(ch_onset_ds_ied_alff_all(:, end-1), ch_onset_ds_ied_alff_all(:, end));
[ch_propa_coef_var_alff, ch_propa_fisher_z_alff] = ...
    remove_nan_inf(ch_propa_ds_ied_alff_all(:, end-1), ch_propa_ds_ied_alff_all(:, end));
[ch_ext_coef_var_alff, ch_ext_fisher_z_alff] = ...
    remove_nan_inf(ch_ext_ds_ied_alff_all(:, end-1), ch_ext_ds_ied_alff_all(:, end));

% call make_plot defined below to obtain scatterplot for each category of channels, 
% provide op_results to signify if output should bed saved to path
% [f1, ax11, ax12] = make_plot(ch_onset_coef_var_reho, ch_onset_fisher_z_reho, ...
%     ch_onset_coef_var_alff, ch_onset_fisher_z_alff, 'onset', ...
%     op_results, directname_op, filenameg_op);   % onset channels
% [f2, ax21, ax22] = make_plot(ch_propa_coef_var_reho, ch_propa_fisher_z_reho, ...
%     ch_propa_coef_var_alff, ch_propa_fisher_z_alff, 'propa', ...
%     op_results, directname_op, filenameg_op);   % propa channels
[f3, ax31, ax32] = make_plot(ch_ext_coef_var_reho, ch_ext_fisher_z_reho, ...
    ch_ext_coef_var_alff, ch_ext_fisher_z_alff, 'ext', ...
    op_results, directname_op, filenameg_op);   % ext(onset+propa) channels


%--------------------------------------------
% make plot

% This function requires the following inputs to make scatterplot of (Y)
% dynamic correlation vs (X) coefficient of variation of reho (left) and 
% alff (right),
% coef_var_reho = coefficient of variation (cov) of ied rates from reho analysis
% fisher_z_reho = fisher-transformed dynamic corr. btw. reho and cov
% coef_var_reho = coefficient of variation (cov) of ied rates from alff analysis
% fisher_z_reho = fisher-transformed dynamic corr. btw. alff and cov
% ch_cat_str = str denoting type of channels interested (onset, propa, ext)
% op_results = 1 to denote save figures to path, 0 to denote otherwise
% directname_op = path to parent directory where outputs are saved
% filenameg_op = filename of output figure, to be appended with ch_cat_str,
% and generates the following outputs, 
% f1 = figure handle
% ax1 = gca for reho plot (left)
% ax2 = gca for alff plot (right)
function [f1, ax1, ax2] = make_plot(coef_var_reho, fisher_z_reho, coef_var_alff, ...
    fisher_z_alff, ch_cat_str, op_results, directname_op, filenameg_op)

f1 = figure('units','normalized','outerposition',[0 0 1 1]);

% specify circle size
sz = 15;

%----------------------------------------------
% subplot for reho

sp1 = subplot(1, 2, 1);

% add scatterplot (X): coeff. of var., (Y): fisher z-score of reho vs ied
% rate
scatter(coef_var_reho, fisher_z_reho, sz, 'filled', ...
    'MarkerFaceColor', [0 0 0.7]); hold on

% call get_regression_line to obtain predicted y data, pearson's corr.
% coeff., and its p-value based on (X) and (Y)
[y_pred_reho, pear_rho_reho, p_val_reho] = get_regression_line(coef_var_reho, fisher_z_reho);

% plot regression line
plot(coef_var_reho, y_pred_reho, '-k', 'LineWidth', 1.15); hold off

box on

% set x- and y-labels, 
xlabel_rho_1 = round(pear_rho_reho, 3);   % corr. to 3 d. p.
xlabel_p_val_1 = round(p_val_reho, 1, 'significant');   % corr. to 1 sig. fig.

if xlabel_p_val_1 >= 0.001
    xlabel(sprintf('\\((\\sigma/\\mu)_{\\mathrm{ReHo}}\\), \\(r = %.3f, p = %.3f\\)', ...
        xlabel_rho_1, xlabel_p_val_1), 'Interpreter','latex');
else
    xlabel(sprintf('\\((\\sigma/\\mu)_{\\mathrm{ReHo}}\\), \\(r = %.3f, p < 0.001\\)', ...
        xlabel_rho_1), 'Interpreter','latex');
end

ylabel('Dynamic correlation','Interpreter','latex');

% format title:
reho_mean_fisher_z = mean(fisher_z_reho);   % get mean Fisher'z
[~, reho_mean_fisher_z_pval] = ttest(fisher_z_reho);   % get p-value of fisher z
reho_mean_fisher_pval_disp = round(reho_mean_fisher_z_pval, 1, 'significant');   % round numbers
t1 = title(sprintf('ReHo (\\(\\bar{z} = %.3f, p = %.1f\\))', reho_mean_fisher_z, ...
    reho_mean_fisher_pval_disp),'Interpreter','latex');

% get current gca
ax1 = gca;

% set line width of axes
ax1.LineWidth = 0.9;

% set x and y limits
xmin = 0;
xmax = 1.5;
xtick = 0.5;
    
ymin = -1;
ymax = 1;
ytick = 0.2;
    
% format axes
ax1.XLim = [xmin,xmax];   % x-limits
ax1.XTick = [xmin:xtick:xmax];
    
ax1.YLim = [ymin,ymax];   % y-limits
ax1.YTick = [ymin:ytick:ymax];
   
% set axis fontsize (default = 10)
% ax1.FontSize = 12;
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
text(0.915, 0.95, '(a)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 18);   

% place string on plot denoting channel type (ied onset, propa, no ied,
% etc)
if strcmp(ch_cat_str, 'onset')
ch_cat_label = 'IED onset channels';
text(0.6, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14);  
end
if strcmp(ch_cat_str, 'propa')
ch_cat_label = 'IED propagation channels';
text(0.45, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14);
end
if strcmp(ch_cat_str, 'ext')
ch_cat_label = 'IED onset+propa channels';
text(0.45, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14);
end

% adjust title position after setting aspect ratio of plot
t1.Position(2) = t1.Position(2) + 0.04;

%----------------------------------------------
% subplot for alff

sp1 = subplot(1, 2, 2);

% add scatterplot (X): coeff. of var., (Y): fisher z-score of alff vs ied
% rate
scatter(coef_var_alff, fisher_z_alff, sz, 'filled', ...
    'MarkerFaceColor', [0 0 0.7]); hold on

% call get_regression_line to obtain predicted y data, pearson's corr.
% coeff., and its p-value based on (X) and (Y)
[y_pred_alff, pear_rho_alff, p_val_alff] = get_regression_line(coef_var_alff, fisher_z_alff);

% plot regression line
plot(coef_var_alff, y_pred_alff, '-k', 'LineWidth', 1.15); hold off

box on;
% grid on;

% set x- and y-labels, 
% append pearson's r and p-value to x-label, corr. to 3 decimal places
xlabel_rho_2 = round(pear_rho_alff, 3);   % corr. to 3 d. p.
xlabel_p_val_2 = round(p_val_alff, 1, 'significant');   % corr. to 1 sig. fig.

if xlabel_p_val_2 >= 0.001
    xlabel(sprintf('\\((\\sigma/\\mu)_{\\mathrm{ALFF}}\\), \\(r = %.3f, p = %.2f\\)', ...
        xlabel_rho_2, xlabel_p_val_2), 'Interpreter','latex');
else
    xlabel(sprintf('\\((\\sigma/\\mu)_{\\mathrm{ALFF}}\\), \\(r = %.3f, p < 0.001\\)', ...
        xlabel_rho_2), 'Interpreter','latex');
end

% ylabel('Dynamic correlation','Interpreter','latex');

% format title:
alff_mean_fisher_z = mean(fisher_z_alff);   % get mean Fisher'z
[~, alff_mean_fisher_z_pval] = ttest(fisher_z_alff);   % get p-value of fisher z
alff_mean_fisher_pval_disp = round(alff_mean_fisher_z_pval, 1, 'significant');   % round numbers
t2 = title(sprintf('ALFF (\\(\\bar{z} = %.3f, p = %.1f\\))', alff_mean_fisher_z, ...
    alff_mean_fisher_pval_disp),'Interpreter','latex');

% get current gca
ax2 = gca;

% set line width of axes
ax2.LineWidth = 0.9;

% set x and y limits
xmin = 0;
xmax = max(coef_var_alff);
xmax = 1.5;
xtick = 0.5;
    
ymin = -1;
ymax = 1;
ytick = 0.2;
    
% format axes
ax2.XLim = [xmin,xmax];   % x-limits
ax2.XTick = [xmin:xtick:xmax];
    
ax2.YLim = [ymin,ymax];   % y-limits
ax2.YTick = [ymin:ytick:ymax];
  
% set axis fontsize (default = 10)
ax2.FontSize = 10;
ax2.YLabel.FontSize = 16;
ax2.XLabel.FontSize = 14;
t2.FontSize = 14;   % title font size

% set properties of ticks
ax2.TickDir = 'out';   % direction of ticks
ax2.YMinorTick = 'on';   % show minor y-ticks
% ax2.TickLength = 2*ax1.TickLength;   % set length of ticks
ax2.TickLength = 2*[ax2.TickLength(1) 4*ax2.TickLength(1)];   % set length of ticks

% set interpreter of tick labels on both axes 
ax2.TickLabelInterpreter = 'latex';
    
% set aspect ratio of plot
pbr = 1;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

% place string on plot
text(0.915, 0.95, '(b)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 18);   

% place string on plot denoting channel type (ied onset, propa, no ied,
% etc)
if strcmp(ch_cat_str, 'onset')
ch_cat_label = 'IED onset channels';
text(0.6, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14);  
end
if strcmp(ch_cat_str, 'propa')
ch_cat_label = 'IED propagation channels';
text(0.45, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14);
end
if strcmp(ch_cat_str, 'ext')
ch_cat_label = 'IED onset+propa channels';
text(0.45, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14);
end

% adjust title position after setting aspect ratio of plot
t2.Position(2) = t2.Position(2) + 0.04;

% end subplot for alff
%----------------------------------------------

%--------------------------------------------------
% adjust positions of subplots
hdiff = 0.1/2.5;
ax2.Position(1) = ax1.Position(1) + ax1.Position(3) + hdiff;
ax2.Position(2) = ax1.Position(2);

% end adjusting positions of scatterplots
%--------------------------------------------------

% output plots to path
if op_results == 1
saveas(f1, fullfile(directname_op, [filenameg_op, '_', ch_cat_str]), 'epsc');
saveas(f1, fullfile(directname_op, [filenameg_op, '_', ch_cat_str]), 'jpeg');
end


end   % end function [f1, ax1] = make_plot(coef_var_reho, fisher_z_reho, coef_var_alff, fisher_z_alff)





% USER DEFINED FUNCTIONS:
%---------------------------------------
% This function intakes two arrays and output the same arrays with entries
% of NaN and Inf removed
% if an entry is removed in one input array, the corresponding entry in
% another array is also removed, so entries in input arrays are paired.
function [output_array_1, output_array_2] = remove_nan_inf(input_array_1, ...
    input_array_2)

input_array_1 = cell2mat(input_array_1);
input_array_2 = cell2mat(input_array_2);

ind_req_non_ran_non_inf = [];
for i = 1:size(input_array_1)
    if ~isnan(input_array_1(i)) && ~isnan(input_array_2(i)) ...
            && ~isinf(input_array_1(i)) && ~isinf(input_array_2(i))
        ind_req_non_ran_non_inf = [ind_req_non_ran_non_inf; i];
    end
end

% use indices obtained above to select required entries in arrays
output_array_1 = input_array_1(ind_req_non_ran_non_inf);
output_array_2 = input_array_2(ind_req_non_ran_non_inf);

end

%---------------------------------------

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

% This function intakes the following inputs, 
% pear_rho = pearon's corr. coeff.
% p_val = p-value of pearson's corr. ceoff.
% and generates the following outputs
% xlabel_command_str = command (str) to be used in sprintf for xlabel
% xlabel_rho = rounded (3 d. p.) pearson's corr. coeff. added to sprintf
% xlabel_p_val = p-value to be added to sprintf
function [xlabel_command_str, xlabel_rho, xlabel_p_val] = format_xlabel(pear_rho, p_val)

% round pearson's corr. coeff. to 3 decimal places
xlabel_rho = round(pear_rho, 3);
xlabel_p_val = round(p_val, 1, 'significant');

if p_val >= 0.001
    xlabel_command_str = '\\((\\sigma/\\mu)_{\\mathrm{ALFF}}\\), \\(r = %.3f, p = %.3f\\)';
else
    xlabel_command_str = '\\((\\sigma/\\mu)_{\\mathrm{ALFF}}\\), \\(r = %.3f, p < 0.001\\)';
end

end

%---------------------------------------


