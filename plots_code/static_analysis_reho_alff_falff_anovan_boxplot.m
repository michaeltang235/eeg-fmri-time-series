clear all 
close all

% This script intakes 'static_analysis_reho.mat', 'static_analysis_alff.mat' 
% and 'static_analysis_falff.mat'
% and performs the following, 
% (1) n-way anova to identify if there are statistically diff. group means, here
% group denotes channel combination (onset-onset, onset-propa, etc., ...)
% (2) multiple comparison tests to identify which group is sig. diff. from
% one another
% (3) using results from mult. comp. tests, make box plot of 95% confidence
% interval of alff and falff of each group (mean +/- 2*standard error)
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

tic 

% enter path to directory where input data is stored
input_dir_reho = 'C:\Users\siumichael.tang\Downloads\fmri_project\matrices\static_analysis_reho';
input_dir_alff = 'C:\Users\siumichael.tang\Downloads\fmri_project\matrices\static_analysis_alff'; 
input_dir_falff = 'C:\Users\siumichael.tang\Downloads\fmri_project\matrices\static_analysis_falff'; 

% input_dir_falff = '/work/levan_lab/mtang/fmri_project/matrices/static_analysis_reho';
% input_dir_falff = '/work/levan_lab/mtang/fmri_project/matrices/static_analysis_falff';
% input_dir_alff = '/work/levan_lab/mtang/fmri_project/matrices/static_analysis_alff';

% enter filename of input data
filename_input_reho = 'static_analysis_reho_1.mat';
filename_input_alff = 'static_analysis_alff_1.mat';
filename_input_falff = 'static_analysis_falff_1.mat';

% enter path where plots are saved at
% directname_op = '/work/levan_lab/mtang/fmri_project/plots/static_analysis_alff_falff_boxplot';
directname_op = 'C:\Users\siumichael.tang\Downloads\fmri_project\plots\static_analysis_reho_alff_falff_boxplot';

% format output filename
filenameg_op = ['static_analysis_reho_alff_falff_anovan_boxplot'];

% enter if user wants to write plots to file (1=yes, 0=no)
op_results = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% get full path of input files
input_path_reho = fullfile(input_dir_reho, filename_input_reho);
input_path_alff = fullfile(input_dir_alff, filename_input_alff);
input_path_falff = fullfile(input_dir_falff, filename_input_falff);

% load input data using path
input_data_reho = load(input_path_reho);
input_data_alff = load(input_path_alff);
input_data_falff = load(input_path_falff);

% vertically concatenate array interested across sessions for reho, for
% each type of channels, onset, propa, no-ied, all = onset + propa + no ied
ch_onset_lt_ied_reho_all = collect_ar(input_data_reho, 'ch_onset_lt_ied_reho');
ch_propa_lt_ied_reho_all = collect_ar(input_data_reho, 'ch_propa_lt_ied_reho');
ch_noied_lt_ied_reho_all = collect_ar(input_data_reho, 'ch_noied_lt_ied_reho');
% ch_all_lt_ied_reho_all = collect_ar(input_data_reho, 'ch_all_lt_ied_reho');

% vertically concatenate array interested across sessions for alff, for
% each type of channels, onset, propa, no-ied, all = onset + propa + no ied
ch_onset_lt_ied_alff_all = collect_ar(input_data_alff, 'ch_onset_lt_ied_alff');
ch_propa_lt_ied_alff_all = collect_ar(input_data_alff, 'ch_propa_lt_ied_alff');
ch_noied_lt_ied_alff_all = collect_ar(input_data_alff, 'ch_noied_lt_ied_alff');
% ch_all_lt_ied_alff_all = collect_ar(input_data_alff, 'ch_all_lt_ied_alff');

% vertically concatenate array interested across sessions for falff, for
% each type of channels, onset, propa, no-ied, all = onset + propa + no ied
ch_onset_lt_ied_falff_all = collect_ar(input_data_falff, 'ch_onset_lt_ied_falff');
ch_propa_lt_ied_falff_all = collect_ar(input_data_falff, 'ch_propa_lt_ied_falff');
ch_noied_lt_ied_falff_all = collect_ar(input_data_falff, 'ch_noied_lt_ied_falff');
% ch_all_lt_ied_falff_all = collect_ar(input_data_falff, 'ch_all_lt_ied_falff');

% get arrays required for n-way anova for ReHo, ALFF, and fALFF
[op_reho, ch_type_all_reho, subject_all_reho, reho_all] = ...
    process_ar(ch_onset_lt_ied_reho_all, ch_propa_lt_ied_reho_all, ch_noied_lt_ied_reho_all);

[op_alff, ch_type_all_alff, subject_all_alff, alff_all] = ...
    process_ar(ch_onset_lt_ied_alff_all, ch_propa_lt_ied_alff_all, ch_noied_lt_ied_alff_all);

[op_falff, ch_type_all_falff, subject_all_falff, falff_all] = ...
    process_ar(ch_onset_lt_ied_falff_all, ch_propa_lt_ied_falff_all, ch_noied_lt_ied_falff_all);

% use anovan, (n-way) anova, to determine if there is sig. diff. btw.
% group means, channel type and subjects for ReHo, ALFF and fALFF
[p1_reho,tbl1_reho,stats1_reho, terms1_reho] = ...
    anovan(reho_all, {ch_type_all_reho subject_all_reho}, "Model", "linear", ...
    "random", 2, "Varnames", {'ch_type', 'subject'}, 'display', 'off');

[p1_alff,tbl1_alff,stats1_alff, terms1_alff] = ...
    anovan(alff_all, {ch_type_all_alff subject_all_alff}, "Model", "linear", ...
    "random", 2, "Varnames", {'ch_type', 'subject'}, 'display', 'off');

[p1_falff,tbl1_falff,stats1_falff, terms1_falff] = ...
    anovan(falff_all, {ch_type_all_falff subject_all_falff}, "Model", "linear", ...
    "random", 2, "Varnames", {'ch_type', 'subject'}, 'display', 'off');

% use multcompare to get mean and standard error of each group for
% pairwise comparison (which group means are sig. diff. from others), for
% ReHo, ALFF and fALFF
% see doc. for details
[c_reho,m_reho,h_reho,gnames_reho] = multcompare(stats1_reho, 'display', 'off');
[c_alff,m_alff,h_alff,gnames_alff] = multcompare(stats1_alff, 'display', 'off');
[c_falff,m_falff,h_falff,gnames_falff] = multcompare(stats1_falff, 'display', 'off');

% format of m is given below, 
% 1st col. = mean estimate of every group
% 2nd col. = standard error of the mean

%--------------------------

% open figure window
% f1 = figure;
% f1 = figure('units','normalized','outerposition',[0 0 1 1]);
f1 = figure('units','normalized','outerposition',[0 0 1 0.7]);

%----------------------------------
% make subplot, 1st one, ReHo
sp1 = subplot(1, 3, 1);

% obtain x data (mean ReHo of each group) from col. 1 of m_reho, generated by
% multcompare
x_data_reho = m_reho(:, 1);

% format name of every group, concatenate them vertically
g_reho = [{'IED onset'}; {'IED propa'}; {'no IED'}]; 

% make boxplot showing mean and 95% confidence interval (std error obtained 
% from col. 2 of m) of every group, 
a_reho = boxplot(x_data_reho, g_reho); hold on
e_reho = errorbar(x_data_reho, 2*m_reho(:, 2), 'LineStyle', 'none', 'Color', 'b'); hold on

% set length of cap at end of error bars
e_reho.CapSize = 10;

% add significant bridges on plot (manually)
% there are 2 sets of groups that are statistically different from each
% other

% enter positions of significant bridge for the 1st set of groups
xl1 = 1;   % x-pos, left 
xr1 = 3.05;   % x-pos, right
xh1 = 0.645;   % height of bridge
xd1 = 0.01/2;   % drop of bridge on each end

% draw sig. bridge using coordinates specified for 1st set of groups
line([0.5*(xl1+xr1), 0.5*(xl1+xr1)], [xh1+xd1/2, xh1+xd1/2], 'Color', 'k', 'Marker', '*', 'LineWidth', 1.1); hold on
line([xl1, xr1], [xh1, xh1], 'Color', 'k', 'LineWidth', 1.1); hold on
line([xl1, xl1], [xh1-xd1, xh1], 'Color', 'k', 'LineWidth', 1.1); hold on
line([xr1, xr1], [xh1-xd1, xh1], 'Color', 'k', 'LineWidth', 1.1); hold on

% enter positions of significant bridge for the 1st set of groups
xl2 = 2;   % x-pos, left 
xr2 = 2.95;   % x-pos, right
xh2 = 0.638;   % height of bridge
xd2 = 0.01/2;   % drop of bridge on each end

% draw sig. bridge using coordinates specified for 2nd set of groups
line([0.5*(xl2+xr2), 0.5*(xl2+xr2)], [xh2+xd2/2, xh2+xd2/2], 'Color', 'k', 'Marker', '*', 'LineWidth', 1.1); hold on
line([xl2, xr2], [xh2, xh2], 'Color', 'k', 'LineWidth', 1.1); hold on
line([xl2, xl2], [xh2-xd2, xh2], 'Color', 'k', 'LineWidth', 1.1); hold on
line([xr2, xr2], [xh2-xd2, xh2], 'Color', 'k', 'LineWidth', 1.1); hold on

hold off

% set linewidth of box plot
set(a_reho, 'LineWidth', 1.1);   % set linewidth of boxes
set(e_reho, 'LineWidth', 1.1);   % set linewidth of error bars

box on 

% get current gca
ax1 = gca;

% set line width of axes
ax1.LineWidth = 0.9;

% set limits of y-axis
ax1.YLim = [0.54 0.66];
ax1.YTick = [ax1.YLim(1):0.02:ax1.YLim(2)];

% set y-label
ylabel('ReHo','Interpreter','latex');

% set x-tick labels
ax1.XTickLabel = {'IED onset', 'IED propagation', 'no IED'};
ax1.XTickLabelRotation = 35;

% set fontsize of axes
% ax1.FontSize = 10;   % axes font size
ax1.XAxis.FontSize = 11.5;   % x-axis font size (ticks)
ax1.YAxis.FontSize = 9;   % y-axis font size (ticks)
ax1.YLabel.FontSize = 14;   % y-label font size
% ax1.XLabel.FontSize = 13;   % x-label font size

% obtain axes limits
xmax = ax1.XLim(2);
xmin = ax1.XLim(1);

% set aspect ratio of plot
pbr = 1;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

% set properties of ticks
ax1.TickDir = 'out';   % direction of ticks
ax1.YMinorTick = 'on';   % show minor y-ticks
% ax1.TickLength = 2*ax1.TickLength;   % set length of ticks
ax1.TickLength = 2*[ax1.TickLength(1) 4*ax1.TickLength(1)];   % set length of ticks

% set tick label interpreter
ax1.TickLabelInterpreter = 'latex';

% place string on plot
text(0.88, 0.935, '(a)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 17);  


%----------------------------------
% make subplot, 2nd one, ALFF
sp2 = subplot(1, 3, 2);

% obtain x data (mean alff of each group) from col. 1 of m_alff, generated by
% multcompare
x_data_alff = m_alff(:, 1);

% format name of every group, concatenate them vertically
g_alff = [{'IED onset'}; {'IED propa'}; {'no IED'}]; 

% make boxplot showing mean and 95% confidence interval (std error obtained 
% from col. 2 of m) of every group, 
a_alff = boxplot(x_data_alff, g_alff); hold on
e_alff = errorbar(x_data_alff, 2*m_alff(:, 2), 'LineStyle', 'none', 'Color', 'b'); hold on

% set length of cap at end of error bars
e_alff.CapSize = 10;

hold off

% set linewidth of box plot
set(a_alff, 'LineWidth', 1.1);   % set linewidth of boxes
set(e_alff, 'LineWidth', 1.1);   % set linewidth of error bars

box on 

% get current gca
ax2 = gca;

% set line width of axes
ax2.LineWidth = 0.9;

% set limits of y-axis
ax2.YLim = [0.7, 0.88];
ax2.YTick = [ax2.YLim(1):0.02:ax2.YLim(2)];

% set y-label
ylabel('ALFF','Interpreter','latex');

% set x-tick labels
ax2.XTickLabel = {'IED onset', 'IED propagation', 'no IED'};
ax2.XTickLabelRotation = 35;

% set fontsize of axes
% ax2.FontSize = 10;   % axes font size
ax2.XAxis.FontSize = 11.5;   % x-axis font size (ticks)
ax2.YAxis.FontSize = 9;   % y-axis font size (ticks)
ax2.YLabel.FontSize = 14;   % y-label font size
% ax2.XLabel.FontSize = 13;   % x-label font size

% obtain axes limits
xmax = ax2.XLim(2);
xmin = ax2.XLim(1);

% set aspect ratio of plot
pbr = 1;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

% set properties of ticks
ax2.TickDir = 'out';   % direction of ticks
ax2.YMinorTick = 'on';   % show minor y-ticks
% ax1.TickLength = 2*ax2.TickLength;   % set length of ticks
ax2.TickLength = 2*[ax2.TickLength(1) 4*ax2.TickLength(1)];   % set length of ticks

% set tick label interpreter
ax2.TickLabelInterpreter = 'latex';

% place string on plot
text(0.88, 0.935, '(b)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 17);  

%--------------------------
%----------------------------------
% make subplot, 2nd one, fALFF
sp2 = subplot(1, 3, 3);

% obtain x data (mean fALFF of each group) from col. 1 of m_falff, generated by
% multcompare
x_data_falff = m_falff(:, 1);

% format name of every group, concatenate them vertically
g_falff = [{'IED onset'}; {'IED propa'}; {'no IED'}]; 

% make boxplot showing mean and 95% confidence interval (std error obtained 
% from col. 2 of m) of every group, 
a_falff = boxplot(x_data_falff, g_falff); hold on
e_falff = errorbar(x_data_falff, 2*m_falff(:, 2), 'LineStyle', 'none', 'Color', 'b'); 

% add significant bridges on plot (manually)
% there are 2 sets of groups that are statistically different from each
% other

% enter positions of significant bridge for the 1st set of groups
xl1 = 1;   % x-pos, left 
xr1 = 3.05;   % x-pos, right
xh1 = 0.9675;   % height of bridge
xd1 = 0.01/2;   % drop of bridge on each end

% draw sig. bridge using coordinates specified for 1st set of groups
line([0.5*(xl1+xr1), 0.5*(xl1+xr1)], [xh1+xd1/2, xh1+xd1/2], 'Color', 'k', 'Marker', '*', 'LineWidth', 1.1); hold on
line([xl1, xr1], [xh1, xh1], 'Color', 'k', 'LineWidth', 1.1); hold on
line([xl1, xl1], [xh1-xd1, xh1], 'Color', 'k', 'LineWidth', 1.1); hold on
line([xr1, xr1], [xh1-xd1, xh1], 'Color', 'k', 'LineWidth', 1.1); hold on

% enter positions of significant bridge for the 1st set of groups
xl2 = 2;   % x-pos, left 
xr2 = 2.95;   % x-pos, right
xh2 = 0.9625;   % height of bridge
xd2 = 0.01/2;   % drop of bridge on each end

% draw sig. bridge using coordinates specified for 2nd set of groups
line([0.5*(xl2+xr2), 0.5*(xl2+xr2)], [xh2+xd2/2, xh2+xd2/2], 'Color', 'k', 'Marker', '*', 'LineWidth', 1.1); hold on
line([xl2, xr2], [xh2, xh2], 'Color', 'k', 'LineWidth', 1.1); hold on
line([xl2, xl2], [xh2-xd2, xh2], 'Color', 'k', 'LineWidth', 1.1); hold on
line([xr2, xr2], [xh2-xd2, xh2], 'Color', 'k', 'LineWidth', 1.1); hold on

% set length of cap at end of error bars
e_falff.CapSize = 10;

hold off

set(a_falff, 'LineWidth', 1.1);   % set linewidth of boxes
set(e_falff, 'LineWidth', 1.1);   % set linewidth of error bars

box on 

% get current gca
ax3 = gca;

% set line width of axes
ax3.LineWidth = 0.9;

% set limits of y-axis
ax3.YLim = [0.9, 0.98];
ax3.YTick = [ax3.YLim(1):0.02:ax3.YLim(2)];

% y-label
ylabel('fALFF','Interpreter','latex');
ax3.YLabel.FontSize = 14;

% set x-tick labels
ax3.XTickLabel = {'IED onset', 'IED propagation', 'no IED'};
ax3.XTickLabelRotation = 35;

% set fontsize of axes
% ax3.FontSize = 10;   % axes font size
ax3.XAxis.FontSize = 11.5;   % x-axis font size (ticks)
ax3.YAxis.FontSize = 9;   % y-axis font size (ticks)
ax3.YLabel.FontSize = 14;   % y-label font size
% ax3.XLabel.FontSize = 13;   % x-label font size

% obtain axes limits
xmax = ax3.XLim(2);
xmin = ax3.XLim(1);

% set aspect ratio of plot
pbr = 1;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

% set properties of ticks
ax3.TickDir = 'out';   % direction of ticks
ax3.YMinorTick = 'on';   % show minor y-ticks
% ax3.TickLength = 2*ax3.TickLength;   % set length of ticks
ax3.TickLength = 2*[ax3.TickLength(1) 4*ax3.TickLength(1)];   % set length of ticks

% set tick label interpreter
ax3.TickLabelInterpreter = 'latex';

% place string on plot
text(0.88, 0.935, '(c)', 'Interpreter', ...
    'latex','FontWeight','bold', 'Units', 'normalized', 'FontSize', 17);

%--------------------------

% % set spacing between subplots
hdiff = 0.1/1.4;   % horizontal spacing btw. subplots

ax2.Position(3) = ax1.Position(3);   % set width of 2nd subplot equal to 1st one
ax3.Position(3) = ax2.Position(3);   % set width of 3rd subplot equal to 2nd one

ax1.Position(1) = ax1.Position(1) - hdiff;   % adjust x-pos of 1st subplot
ax2.Position(1) = ax1.Position(1) + ax1.Position(3) + hdiff/1.3;
ax3.Position(1) = ax2.Position(1) + ax2.Position(3) + hdiff/1.3;


% output falff and alff plots to path
if op_results == 1
saveas(f1, fullfile(directname_op, filenameg_op), 'epsc');
saveas(f1, fullfile(directname_op, filenameg_op), 'jpeg');
end


% USER-DEFINED FUNCTION:
%---------------------------------------------------------------------------

%---------------------------------------------------------------
% This function requires the following inputs, 
% input_st = input structure
% ar_int_str = str of array to be concatenated vertically across sessions
% and generates 
% output_ar = output array (vertical concatenation of array interested across
% sessions)
% note: array being concatenated has the following format, 
% 1st col. = subject number
% 2nd col. = ied type 
% 3rd col. = channel name
% 4th col. = long term monitoring ied rate per minute
% 5th col. = time series feature (e.g. ALFF, fALFF)
function output_ar = collect_ar(input_st, ar_int_str)

output_ar = [];   % initialize output array

fd_name_sub = fieldnames(input_st.grand);   % obtain fieldnames under grand of input array 
for sub_ind = 1:numel(fd_name_sub)   % for every subject
    if ~isempty(input_st.grand.(fd_name_sub{sub_ind}))   % proceed if curr. sub. is not empty
        fd_name_run = fieldnames(input_st.grand.(fd_name_sub{sub_ind}));   % obtain run fieldnames
        for run_ind = 1:numel(fd_name_run)   % for every run
            % concatenate array vertically
            output_ar = [output_ar; ...
                input_st.grand.(fd_name_sub{sub_ind}).(fd_name_run{run_ind}).(ar_int_str)];   
        end
    end
end

end

%---------------------------------------------------------------

% This function does the following, 
% (1) vertically concatenate arrays of ied onset, propa, and no ied
% channels
% (2) prepares arrays required for n-way anova

% This function requires the following inputs, 
% onset_ar = array of ied onset channels
% propa_ar = array of ied propagation channels
% noied_ar = array of no-ied channels,
% and generates the following outputs
% output_ar = vertical concatenation of all input arrays
% ch_type_all = array specifying channel type in each row of output_ar
% subject_all = subject number of each row in output_ar
% feat_all = magnitude of feature (ALFF or fALFF) in each row in output_ar
function [output_ar, ch_type_all, subject_all, feat_all] = ...
    process_ar(onset_ar, propa_ar, noied_ar)

% get number of rows in each input array
num_row_onset = size(onset_ar, 1);   % ied onset channels
num_row_propa = size(propa_ar, 1);   % ied propa channels
num_row_noied = size(noied_ar, 1);   % no-ied channels

% vertically concatenate input arrays
output_ar = [onset_ar; propa_ar; noied_ar];

% assign channel type based on order at which entries in output_ar appears
ch_type_all = [repmat('onset', [num_row_onset, 1]); ...
    repmat('propa', [num_row_propa, 1]); ...
    repmat('noied', [num_row_noied, 1]);];

% get array of subject number
subject_all = output_ar(:, 1);

% get array of time series feature interested (ALFF or fALFF) from rightmost 
% column in output_ar
feat_all = cell2mat(output_ar(:, end));

end