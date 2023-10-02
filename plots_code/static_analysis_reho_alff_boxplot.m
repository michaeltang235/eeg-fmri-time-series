clear all
close all

%---------------------------------------------------------------------------
% Pearon's corr. coef. (r)^2 = Coef. of determination, only for simple
% linear regression (only one predictor variable, or independent variable)
% https://stats.stackexchange.com/questions/83347/relationship-between-r2-and-correlation-coefficient
% https://towardsdatascience.com/r%C2%B2-or-r%C2%B2-when-to-use-what-4968eee68ed3
%---------------------------------------------------------------------------
% Drawn on 'static_analysis_reho.mat' and 'static_analysis_alff.mat', which
% have
% (i) ReHo or ALFF, depending on filename 
% (ii) IED rate mesured during long-term monitoring
% (iii) IED rate measured over entire scan, for each session of every
% subject.
% The Spearman's correlation coefficient and p-value between 
% (i) and (ii), called it (rho_var_ied_lt), 
% (i) and (iii), called it (rho_var_ied_ds),
% where var denotes the fmri feature interested, was calculated for each session.
% This script first applies fisher transform on all rhos, then 
% makes boxplot of rho_var_ied_lt and rho_var_ied_ds across all subjects.
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter path where input struct. is stored at 
% directname = ['/work/levan_lab/mtang/fmri_project/matrices', filesep 'static_analysis_alff'];   % direct. of output matrix
directname = 'C:\Users\siumichael.tang\Downloads\fmri_project\matrices';
filename_input_reho = 'static_analysis_reho_1.mat';   % filename of input file of reho
filename_input_alff = 'static_analysis_alff_1.mat';   % filename of input file of alff

% enter path where plots are saved at
% directname_op = '/work/levan_lab/mtang/fmri_project/plots/static_analysis_alff';
directname_op = 'C:\Users\siumichael.tang\Downloads\fmri_project\paper_plots';

% format output filename
filenameg_op = ['static_analysis_reho_alff_boxplot'];

% enter if user wants to write plots to file (1=yes, 0=no)
op_results = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% load structure of reho and alff
st_reho = load(fullfile([directname, filesep, 'static_analysis_reho'], filename_input_reho));
st_alff = load(fullfile([directname, filesep, 'static_analysis_alff'], filename_input_alff));

% get array of subjects in each structure loaded, not all subjects have
% alff calculated due to length of time series
fd_name_sub_reho = fieldnames(st_reho.grand);
fd_name_sub_alff = fieldnames(st_alff.grand);

% (i) ReHo, (ii) ALFF, (iii), IED_{long term}, (iv) IED_{during scan}, were
% computed in each run for every subject, and their spearman's corr.
% coeff. was computed and classified into the following categories based on
% channel type. 
% For ReHo/ALFF - IED_{long term}, channel types are onset,
% propagation, no ied, and all = onset + propa + no ied
% For ReHo/ALFF - IED_{during scan}, channel types are onset,
% propagation, and ext = onset + propa.

% use process_var, defined at the end of script, to obtain a struct
% containing the following
% (1) fisher z-score of each rho in each session of every subject
% (2) mean fisher z-score
% (3) p-value of fisher z-score
% For instance,
% rho_reho_ied_lt_onset_st = strcut of rho btw. ReHo and IED_{long term} for ied onset
% channels

% for ReHo/ALFF - IED_{long term}, we have, 
rho_reho_ied_lt_onset_st = process_var(st_reho.grand, 'rho_reho_ied_lt_onset');
rho_reho_ied_lt_propa_st = process_var(st_reho.grand, 'rho_reho_ied_lt_propa');
rho_reho_ied_lt_noied_st = process_var(st_reho.grand, 'rho_reho_ied_lt_noied');
rho_reho_ied_lt_all_st = process_var(st_reho.grand, 'rho_reho_ied_lt_all');

rho_alff_ied_lt_onset_st = process_var(st_alff.grand, 'rho_alff_ied_lt_onset');
rho_alff_ied_lt_propa_st = process_var(st_alff.grand, 'rho_alff_ied_lt_propa');
rho_alff_ied_lt_noied_st = process_var(st_alff.grand, 'rho_alff_ied_lt_noied');
rho_alff_ied_lt_all_st = process_var(st_alff.grand, 'rho_alff_ied_lt_all');

% for ReHo/ALFF - IED_{during scan}, we have, 
rho_reho_ied_ds_onset_st = process_var(st_reho.grand, 'rho_reho_ied_ds_onset');
rho_reho_ied_ds_propa_st = process_var(st_reho.grand, 'rho_reho_ied_ds_propa');
rho_reho_ied_ds_ext_st = process_var(st_reho.grand, 'rho_reho_ied_ds_ext');

rho_alff_ied_ds_onset_st = process_var(st_alff.grand, 'rho_alff_ied_ds_onset');
rho_alff_ied_ds_propa_st = process_var(st_alff.grand, 'rho_alff_ied_ds_propa');
rho_alff_ied_ds_ext_st = process_var(st_alff.grand, 'rho_alff_ied_ds_ext');

% format channel category str required for make_plot function defined below
ch_cat_str_all = {'rho_reho_ied_lt_all', 'rho_alff_ied_lt_all', ...
    'rho_reho_ied_ds_ext', 'rho_alff_ied_ds_ext'};   % all channels

ch_cat_str_onset = {'rho_reho_ied_lt_onset', 'rho_alff_ied_lt_onset', ...
    'rho_reho_ied_ds_onset', 'rho_alff_ied_ds_onset'};   % ied onset channels

ch_cat_str_propa = {'rho_reho_ied_lt_propa', 'rho_alff_ied_lt_propa', ...
    'rho_reho_ied_ds_propa', 'rho_alff_ied_ds_propa'};   % ied propa channels

ch_cat_str_noied = {'rho_reho_ied_lt_noied', 'rho_alff_ied_lt_noied'};   % no ied channels

% call make_plot defined below to obtain boxplot for each category of channels, 
% provide op_results to signify if output should bed saved to path 
[f1, ax1, rho_reho_ied_lt_all_st, rho_reho_ied_ds_all_st, rho_alff_ied_lt_all_st, ...
    rho_alff_ied_ds_all_st] = make_plot(st_reho, st_alff, ch_cat_str_all, ...
    op_results, directname_op, filenameg_op);   % all channels

% [f2, ax2, rho_reho_ied_lt_onset_st, rho_reho_ied_ds_onset_st, rho_alff_ied_lt_onset_st, ...
%     rho_alff_ied_ds_onset_st] = make_plot(st_reho, st_alff, ch_cat_str_onset, ...
%     op_results, directname_op, filenameg_op);   % ied onset channels
% 
% [f3, ax3, rho_reho_ied_lt_propa_st, rho_reho_ied_ds_propa_st, rho_alff_ied_lt_propa_st, ...
%     rho_alff_ied_ds_propa_st] = make_plot(st_reho, st_alff, ch_cat_str_propa, ...
%     op_results, directname_op, filenameg_op);   % ied propa channels
% 
% [f4, ax4, rho_reho_ied_lt_noied_st, rho_reho_ied_ds_noied_st, rho_alff_ied_lt_noied_st, ...
%     rho_alff_ied_ds_noied_st] = make_plot(st_reho, st_alff, ch_cat_str_noied, ...
%     op_results, directname_op, filenameg_op);   % no ied channels


%---------------------------------------------------------------
%---------------------------------------------------------------
% FUNCTION [F1, AX1, RHO_REHO_IED_LT_CAT_ST, ...] = MAKE_PLOT(ST_REHO,
% ST_ALFF, ...)
% this function makes boxplot of static corr. btw. (1) reho - ied_{long term},
% (2) reho - ied_{during scan}, (3) alff - ied_{long term}, (1) alff - ied_{during scan},
% based on category of channels interested
% inputs required are,
% ST_REHO = struct obtained from static analysis of reho
% ST_ALFF = struct obtained from static analysis of alff
% CH_CAT_STR = cell array containing fieldnames interested inside each of the input
% struct (filednames specifying the category of channels required)
% OP_RESULTS = 1 when figures are saved to path, 0 when figures are not
% saved
% DIRECTNAME_OP = path of parent directory where output figures are saved
% FILENAME_OP = base filenames of output figures, additional str is
% appended depending on category of channels required, see below for
% details.
% outputs are, 
% F1 = figure handle
% AX1 = axes gca
% RHO_REHO_IED_LT_CAT_ST = struct listing variables required for boxplot,
% for static corr. btw. reho and ied_{long term}, see process_var for details, 
% RHO_REHO_IED_DS_CAT_ST = struct listing variables required for boxplot,
% for static corr. btw. reho and ied_{during scan}, see process_var for details, 
% RHO_ALFF_IED_LT_CAT_ST = struct listing variables required for boxplot,
% for static corr. btw. alff and ied_{long term}, see process_var for details, 
% RHO_ALFF_IED_DS_CAT_ST = struct listing variables required for boxplot,
% for static corr. btw. alff and ied_{during scan}, see process_var for details, 

function [f1, ax1, rho_reho_ied_lt_cat_st, rho_reho_ied_ds_cat_st, rho_alff_ied_lt_cat_st, ...
    rho_alff_ied_ds_cat_st] = make_plot(st_reho, st_alff, ch_cat_str, op_results, ...
    directname_op, filenameg_op)

% initialize output structs
rho_reho_ied_lt_cat_st = [];   % reho - ied_{long term}
rho_reho_ied_ds_cat_st = [];   % reho - ied_{during scan}
rho_alff_ied_lt_cat_st = [];   % alff - ied_{long term}
rho_alff_ied_ds_cat_st = [];   % alff - ied_{during scan}

% get struct containing all variables required for boxplot using
% process_var, see end of script for details
if numel(ch_cat_str) == 4
rho_reho_ied_lt_cat_st = process_var(st_reho.grand, ch_cat_str{1});   % reho - ied_{long term}
rho_alff_ied_lt_cat_st = process_var(st_alff.grand, ch_cat_str{2});   % reho - ied_{during scan}
rho_reho_ied_ds_cat_st = process_var(st_reho.grand, ch_cat_str{3});   % alff - ied_{long term}
rho_alff_ied_ds_cat_st = process_var(st_alff.grand, ch_cat_str{4});   % alff - ied_{during scan}

% from each struct, access variables required, fisher z-scores, mean fisher
% z score, and the corresponding p-value, for each column in boxplot (reho
% - ied_{long term}, reho - ied_{during scan}, etc.
rho_reho_ied_lt_cat_fisher_z = rho_reho_ied_lt_cat_st.fisher_z;
rho_reho_ied_ds_cat_fisher_z = rho_reho_ied_ds_cat_st.fisher_z;
rho_alff_ied_lt_cat_fisher_z = rho_alff_ied_lt_cat_st.fisher_z;
rho_alff_ied_ds_cat_fisher_z = rho_alff_ied_ds_cat_st.fisher_z;

rho_reho_ied_lt_cat_mean_fisher_z = rho_reho_ied_lt_cat_st.mean_fisher_z;
rho_reho_ied_ds_cat_mean_fisher_z = rho_reho_ied_ds_cat_st.mean_fisher_z;
rho_alff_ied_lt_cat_mean_fisher_z = rho_alff_ied_lt_cat_st.mean_fisher_z;
rho_alff_ied_ds_cat_mean_fisher_z = rho_alff_ied_ds_cat_st.mean_fisher_z;

rho_reho_ied_lt_cat_pval_fisher_z = rho_reho_ied_lt_cat_st.pval_fisher_z;
rho_reho_ied_ds_cat_pval_fisher_z = rho_reho_ied_ds_cat_st.pval_fisher_z;
rho_alff_ied_lt_cat_pval_fisher_z = rho_alff_ied_lt_cat_st.pval_fisher_z;
rho_alff_ied_ds_cat_pval_fisher_z = rho_alff_ied_ds_cat_st.pval_fisher_z;
end

% get struct containing all variables required for boxplot using
% process_var, see end of script for details, if ch_cat_str contains only 2
% entries, for the case that ied_{during scan} is missing 
% when considering channels with no ieds registered during fmri scan
if numel(ch_cat_str) == 2
rho_reho_ied_lt_cat_st = process_var(st_reho.grand, ch_cat_str{1});
rho_alff_ied_lt_cat_st = process_var(st_alff.grand, ch_cat_str{2});

% access variables from each struct
rho_reho_ied_lt_cat_fisher_z = rho_reho_ied_lt_cat_st.fisher_z;
rho_alff_ied_lt_cat_fisher_z = rho_alff_ied_lt_cat_st.fisher_z;

rho_reho_ied_lt_cat_mean_fisher_z = rho_reho_ied_lt_cat_st.mean_fisher_z;
rho_reho_ied_ds_cat_mean_fisher_z = [];
rho_alff_ied_lt_cat_mean_fisher_z = rho_alff_ied_lt_cat_st.mean_fisher_z;
rho_alff_ied_ds_cat_mean_fisher_z = [];

rho_reho_ied_lt_cat_pval_fisher_z = rho_reho_ied_lt_cat_st.pval_fisher_z;
rho_reho_ied_ds_cat_pval_fisher_z = [];
rho_alff_ied_lt_cat_pval_fisher_z = rho_alff_ied_lt_cat_st.pval_fisher_z;
rho_alff_ied_ds_cat_pval_fisher_z = [];
end

%---------------------------------------------------
% make boxplot for all channels
% figure
f1 = figure('units','normalized','outerposition',[0 0 1 1]);

% concatenate all var. vertically
if numel(ch_cat_str) == 4
    x_data = [rho_reho_ied_lt_cat_fisher_z; rho_reho_ied_ds_cat_fisher_z; ...
    rho_alff_ied_lt_cat_fisher_z; rho_alff_ied_ds_cat_fisher_z];
else
    x_data = [rho_reho_ied_lt_cat_fisher_z; rho_alff_ied_lt_cat_fisher_z];
end

% create grouping variable that assigns same value to rows corresponding to
% the same vector in x_data
if numel(ch_cat_str) == 4   % when there are 4 entries in ch_cat_str
g1 = repmat({'ReHo-IED_lt'}, length(rho_reho_ied_lt_cat_fisher_z), 1);   % reho - ied_{long term}
g2 = repmat({'ReHo-IED_ds'}, length(rho_reho_ied_ds_cat_fisher_z), 1);   % reho - ied_{during scan}
g3 = repmat({'ALFF-IED_lt'}, length(rho_alff_ied_lt_cat_fisher_z), 1);   % alff - ied_{long term}
g4 = repmat({'ALFF-IED_ds'}, length(rho_alff_ied_ds_cat_fisher_z), 1);   % alff - ied_{during scan}
else   % when there are not 4 entries (during scan ied rate is missing for the channel category)
g1 = repmat({'ReHo-IED_lt'}, length(rho_reho_ied_lt_cat_fisher_z), 1);   % reho - ied_{long term}
g2 = [];
g3 = repmat({'ALFF-IED_lt'}, length(rho_alff_ied_lt_cat_fisher_z), 1);   % alff - ied_{long term}
g4 = [];
end

g = [g1; g2; g3; g4];   % concatenate them together

% make boxplot
a = boxplot(x_data, g);

% format plot
set(a, 'LineWidth', 1.1);   % set linewidth of boxes
yline(0, '--', 'LineWidth', 1.5)   % place y=0 line 

box on 

% get current gca
ax1 = gca;

% set line width of axes
ax1.LineWidth = 0.9;

% set limits of y-axis
ax1.YLim = [-1.4 1.6];
ax1.YTick = [ax1.YLim(1):0.2:ax1.YLim(2)];

% y-label
ylabel('Static correlation','Interpreter','latex');

% set fontsize of axes
ax1.FontSize = 10;   % axes font size
ax1.YLabel.FontSize = 16;   % y-label font size
ax1.XAxis.FontSize = 12;   % x-axis font size (ticks)
% ax1.XLabel.FontSize = 13;   % x-label font size

% obtain axes limits
xmax = ax1.XLim(2);
xmin = ax1.XLim(1);

% set aspect ratio of plot
pbr = 1/1.2;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

% set properties of ticks
ax1.TickDir = 'out';   % direction of ticks
ax1.YMinorTick = 'on';   % show minor y-ticks
% ax1.TickLength = 2*ax1.TickLength;   % set length of ticks
ax1.TickLength = 2*[ax1.TickLength(1) 4*ax1.TickLength(1)];   % set length of ticks

ax1.TickLabelInterpreter = 'latex';

% use format_xlabel_str to obtain label for each column (each channel category)
% in boxplot, adding '*' to label if the corr. is significant (p-value ,=
% 0.05)
xlabel_str = format_xlabel_str(rho_reho_ied_lt_cat_pval_fisher_z, ...
    rho_reho_ied_ds_cat_pval_fisher_z, rho_alff_ied_lt_cat_pval_fisher_z, ...
    rho_alff_ied_ds_cat_pval_fisher_z);

% add xtick labels to plot
ax1.XTickLabel = xlabel_str;

% display messages on terminal showing mean fisher z and their p-values
display([rho_reho_ied_lt_cat_mean_fisher_z rho_reho_ied_ds_cat_mean_fisher_z, ...
    rho_alff_ied_lt_cat_mean_fisher_z, rho_alff_ied_ds_cat_mean_fisher_z])

display([rho_reho_ied_lt_cat_pval_fisher_z rho_reho_ied_ds_cat_pval_fisher_z, ...
    rho_alff_ied_lt_cat_pval_fisher_z, rho_alff_ied_ds_cat_pval_fisher_z])

% place string on plot denoting channel type (ied onset, propa, no ied,
% etc)
if contains(ch_cat_str{1}, 'all')
ch_cat_label = 'All channels';
text(0.8, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14); 
end
if contains(ch_cat_str{1}, 'onset')
ch_cat_label = 'IED onset channels';
text(0.72, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14); 
end
if contains(ch_cat_str{1}, 'propa')
ch_cat_label = 'IED propagation channels';
text(0.62, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14); 
end
if contains(ch_cat_str{1}, 'noied')
ch_cat_label = 'No IED channels';
text(0.75, 0.05, ch_cat_label, 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 14); 
end 

% format output figure name based on channel category
if contains(ch_cat_str{1}, 'all')
filenameg1 = [filenameg_op, '_all'];
end
if contains(ch_cat_str{1}, 'onset')
filenameg1 = [filenameg_op, '_onset'];
end
if contains(ch_cat_str{1}, 'propa')
filenameg1 = [filenameg_op, '_propa'];
end
if contains(ch_cat_str{1}, 'noied')
filenameg1 = [filenameg_op, '_noied']; 
end

% output plots to path
if op_results == 1
saveas(f1, fullfile(directname_op, filenameg1), 'epsc');
saveas(f1, fullfile(directname_op, filenameg1), 'jpeg');
end

end   % end function [f1, ax1] = make_plot(ch_cat_str)

%---------------------------------------
% This function intakes array of rho and output the same array with entries
% of NaN and Inf removed
function output_array = remove_nan_inf(input_array)

ind_req_non_nan = ~isnan(input_array);   % identify indices of entries with non-NaN value assigned
non_nan_array = input_array(ind_req_non_nan);

ind_req_non_inf = ~isinf(non_nan_array);
output_array = non_nan_array(ind_req_non_inf);   % filter array with indices found

end
%---------------------------------------

%---------------------------------------
% This function intakes array of rho and output the corresponding (i)
% fisher transformed rho, (ii) mean fisher transformed rho, and (iii)
% p-value of the fisher transformed rho
function [fisher_z_noninf, mean_fisher_z, pval_fisher_z] = fisher_transform(input_array)

% compute quantities interested,
fisher_z = atanh(input_array);   % fisher transformed rho

% remove entries of inf in fisher_z for cases with rho = 1 or -1
ind_non_inf = ~isinf(fisher_z);

% select non inf entries of fisher transformed z-score
fisher_z_noninf = fisher_z(ind_non_inf);

% get mean of fisher transformed rho
mean_fisher_z = mean(fisher_z(ind_non_inf)); 

% get p-value of fisher transformed rho
[~, pval_fisher_z] = ttest(fisher_z(ind_non_inf)); 

end
%---------------------------------------
%---------------------------------------

% This function requires the following inputs, 
% input_st = input structure
% var_int_str = name (str) of variable interested, stored inside input_st
% and generates the following output, 
% output_ar = output array, concatenation of var_int across all sessions
% and subjects
function output_ar = collect_var(input_st, var_int_str)

% initialize output array
output_ar = [];

% get array of subjects in input structure loaded, not all subjects have
% var_int calculated due to length of time series
fd_name_sub = fieldnames(input_st);

% access required array from input structure
for sub_ind = 1:numel(fd_name_sub)   % for each subject in input struct.
    if ~isempty(input_st.(fd_name_sub{sub_ind}))   % if value is not empty
        fd_run_cur = fieldnames(input_st.(fd_name_sub{sub_ind}));   % get run fieldnames available
        for run_ind = 1:numel(fd_run_cur)   % for each run
            
            % access sub-structure of curr. sess.
            st_cur = input_st.(fd_name_sub{sub_ind}).(fd_run_cur{run_ind});
            
            % get var_int in curr. session
            var_int_cur = st_cur.(var_int_str);
            
            % concatenate vertically each rho to arrays 
            output_ar = [output_ar; var_int_cur];                
        end
    end
end

end

%---------------------------------------
%---------------------------------------

% This function intakes a struct and variable interested and does the
% following processing steps:
% (1) vertically concatenate var. interested across all sessions and
% subjects into an array
% (2) remove entries of nan and inf in the array obtained in (1)
% (3) fisher transform all entries in arrays obtained in (2) to obtain 
% list of fisher z-scores, remove entries, compute mean fisher z-score and
% its p-value
% and generate an output struct containing the following, 
% (1) fisher z-score of each rho in each session of every subject
% (2) mean fisher z-score
% (3) p-value of fisher z-score
% (4) array of var. interested from all sessions and subjects
% (5) same as (4) but with entries of nan and inf removed
% note: this function requires functions defined above, collect_var,
% remove_nan_inf, and fisher_transform
% function output_st = process_var(input_st, var_int)
% input_st = input structure contaning variable interested
% var_int = name (str) of variable interested
% output_st = output structure containing processed variables, see above
% text for details
function output_st = process_var(input_st, var_int)

% initialize output structure
output_st = struct;

% call collect_var to gather var_int from all sessions and subjects
var_all_sess = collect_var(input_st, var_int);

% call remove_nan_inf to remove entries of nan and inf
var_all_sess_removed_nan_inf = remove_nan_inf(var_all_sess);

% call fisher_transform to obtain fisher z-score, mean fisher z-score, and
% its p-value for the set of var_int 
[fisher_z, mean_fisher_z, pval_fisher_z] = fisher_transform(var_all_sess_removed_nan_inf);

% assign variables computed above to output structrel
output_st.fisher_z = fisher_z;   % fisher z-score, with entries of inf removed
output_st.mean_fisher_z = mean_fisher_z;   % mean fisher z-score
output_st.pval_fisher_z = pval_fisher_z;   % its p-value

output_st.var_int_all_sess = var_all_sess;   % array containing var_int in all sess. and sub.
output_st.var_all_sess_removed_nan_inf = var_all_sess_removed_nan_inf;   % with entries of nan and inf removed

end

%---------------------------------------
%---------------------------------------
% This function requires array of p-value of mean fisher z score of each
% channel category and generates the corresponding xlabel str for plotting.
% '*' is added to the label if the p-value is <= threshold, see below for
% details.
function xlabel_str_output = format_xlabel_str(rho_reho_ied_lt_cat_pval_fisher_z, ...
    rho_reho_ied_ds_cat_pval_fisher_z, rho_alff_ied_lt_cat_pval_fisher_z, ...
    rho_alff_ied_ds_cat_pval_fisher_z)

% initialize xlabel str for plot
xlabel_str_output = {};

% format xlabel str for each column in boxplot
reho_ied_lt_xlabel_str = 'ReHo--IED\(_{\mathrm{long} \: \mathrm{term}}\)';
reho_ied_ds_xlabel_str = 'ReHo--IED\(_{\mathrm{during} \: \mathrm{scan}}\)';
alff_ied_lt_xlabel_str = 'ALFF--IED\(_{\mathrm{long} \: \mathrm{term}}\)';
alff_ied_ds_xlabel_str = 'ALFF--IED\(_{\mathrm{during} \: \mathrm{scan}}\)';

% set threshold of p-value
thresh = 0.05;

% check if p-value in each column is <= thresh, if so, add '*' to the label
if rho_reho_ied_lt_cat_pval_fisher_z <= thresh
    reho_ied_lt_xlabel_str = 'ReHo--IED\(_{\mathrm{long} \: \mathrm{term}}\)*';
end
if rho_reho_ied_ds_cat_pval_fisher_z <= thresh
    reho_ied_ds_xlabel_str = 'ReHo--IED\(_{\mathrm{during} \: \mathrm{scan}}\)*';
end
if rho_alff_ied_lt_cat_pval_fisher_z <= thresh
    alff_ied_lt_xlabel_str = 'ALFF--IED\(_{\mathrm{long} \: \mathrm{term}}\)*';
end
if rho_alff_ied_ds_cat_pval_fisher_z <= thresh
    alff_ied_ds_xlabel_str = 'ALFF--IED\(_{\mathrm{during} \: \mathrm{scan}}\)*';
end

% format output xlabel str based on input arrays of p-values of mean fisher
% z scores
if ~isempty(rho_reho_ied_ds_cat_pval_fisher_z) && ~isempty(rho_alff_ied_ds_cat_pval_fisher_z)
xlabel_str_output = {reho_ied_lt_xlabel_str, reho_ied_ds_xlabel_str, ...
    alff_ied_lt_xlabel_str, alff_ied_ds_xlabel_str};
else
xlabel_str_output = {reho_ied_lt_xlabel_str, alff_ied_lt_xlabel_str};
end

end


%---------------------------------------
%---------------------------------------