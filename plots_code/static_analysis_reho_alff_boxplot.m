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
filename_input_reho = 'static_analysis_reho.mat';   % filename of input file of reho
filename_input_alff = 'static_analysis_alff.mat';   % filename of input file of alff

% enter path where plots are saved at
% directname_op = '/work/levan_lab/mtang/fmri_project/plots/static_analysis_alff';
directname_op = 'C:\Users\siumichael.tang\Downloads\fmri_project\plots';

% format output filename
filenameg1 = ['static_analysis_reho_alff_boxplot'];

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

% initialize arrays interested, spearman's corr. coef. (rho) btw. 
% (A): ReHo-IED_{long term}, (rho_reho_ied_lt),
% (B): ReHo-IED_{during scan}, (rho_reho_ied_ds),
% (C): ALFF-IED_{long term}, (rho_alff_ied_lt),
% (D): ALFF-IED_{during scan}, (rho_alff_ied_ds),

rho_reho_ied_lt_all_sess = [];   % (A) for all sessions
rho_reho_ied_ds_all_sess = [];   % (B) for all sessions
rho_alff_ied_lt_all_sess = [];   % (C) for all sessions
rho_alff_ied_ds_all_sess = [];   % (D) for all sessions

% access required arrays from structure, reho

for sub_ind = 1:numel(fd_name_sub_reho)   % for each subject in reho struct.
    if ~isempty(st_reho.grand.(fd_name_sub_reho{sub_ind}))   % if value is not empty
        fd_run_cur = fieldnames(st_reho.grand.(fd_name_sub_reho{sub_ind}));   % get run fieldnames available
        for run_ind = 1:numel(fd_run_cur)   % for each run
            
            % access sub-structure of curr. sess.
            st_reho_cur = st_reho.grand.(fd_name_sub_reho{sub_ind}).(fd_run_cur{run_ind});
            
            % get rho_reho_ied_lt and rho_reho_ied_ds in curr. session
            rho_reho_ied_lt_cur = st_reho_cur.rho_reho_ied_lt;
            rho_reho_ied_ds_cur = st_reho_cur.rho_reho_ied_ds;
            
            % concatenate vertically each rho to arrays 
            rho_reho_ied_lt_all_sess = [rho_reho_ied_lt_all_sess; rho_reho_ied_lt_cur];
            rho_reho_ied_ds_all_sess = [rho_reho_ied_ds_all_sess; rho_reho_ied_ds_cur];                   
        end
    end
end

% access required arrays from structure, alff

for sub_ind = 1:numel(fd_name_sub_alff)   % for each subject in alff struct.
    if ~isempty(st_alff.grand.(fd_name_sub_alff{sub_ind}))   % if value is not empty
        fd_run_cur = fieldnames(st_alff.grand.(fd_name_sub_alff{sub_ind}));   % get run fieldnames available
        for run_ind = 1:numel(fd_run_cur)   % for each run
            
            % access sub-structure of curr. sess.
            st_alff_cur = st_alff.grand.(fd_name_sub_alff{sub_ind}).(fd_run_cur{run_ind});
            
            % get rho_alff_ied_lt and rho_alff_ied_ds in curr. session
            rho_alff_ied_lt_cur = st_alff_cur.rho_alff_ied_lt;
            rho_alff_ied_ds_cur = st_alff_cur.rho_alff_ied_ds;
            
            % concatenate vertically each rho to arrays 
            rho_alff_ied_lt_all_sess = [rho_alff_ied_lt_all_sess; rho_alff_ied_lt_cur];
            rho_alff_ied_ds_all_sess = [rho_alff_ied_ds_all_sess; rho_alff_ied_ds_cur];              
        end
    end
end

% remove entries of NaN and Inf in each of the rho arrays
rho_reho_ied_lt_all_sess_nonnan = remove_nan_inf(rho_reho_ied_lt_all_sess);
rho_reho_ied_ds_all_sess_nonnan = remove_nan_inf(rho_reho_ied_ds_all_sess);
rho_alff_ied_lt_all_sess_nonnan = remove_nan_inf(rho_alff_ied_lt_all_sess);
rho_alff_ied_ds_all_sess_nonnan = remove_nan_inf(rho_alff_ied_ds_all_sess);

% for each of the rho array above, get
% (i)fisher transformed z,
% (ii) mean fisher z, and 
% (ii) p-value, function defined at the end of script
[rho_reho_ied_lt_all_sess_fisher_z, rho_reho_ied_lt_all_sess_fisher_z_mean, ...
    pval_rho_reho_ied_lt_all_sess_fisher_z] = fisher_transform(rho_reho_ied_lt_all_sess_nonnan);

[rho_reho_ied_ds_all_sess_fisher_z, rho_reho_ied_ds_all_sess_fisher_z_mean, ...
    pval_rho_reho_ied_ds_all_sess_fisher_z] = fisher_transform(rho_reho_ied_ds_all_sess_nonnan);

[rho_alff_ied_lt_all_sess_fisher_z, rho_alff_ied_lt_all_sess_fisher_z_mean, ...
    pval_rho_alff_ied_lt_all_sess_fisher_z] = fisher_transform(rho_alff_ied_lt_all_sess_nonnan);

[rho_alff_ied_ds_all_sess_fisher_z, rho_alff_ied_ds_all_sess_fisher_z_mean, ...
    pval_rho_alff_ied_ds_all_sess_fisher_z] = fisher_transform(rho_alff_ied_ds_all_sess_nonnan);

%---------------------------------------------------
% make boxplot

% figure
f1 = figure('units','normalized','outerposition',[0 0 1 1])

% concatenate all var. vertically
x_data = [rho_reho_ied_lt_all_sess_fisher_z; rho_reho_ied_ds_all_sess_fisher_z; ...
    rho_alff_ied_lt_all_sess_fisher_z; rho_alff_ied_ds_all_sess_fisher_z];

% create grouping variable that assigns same value to rows corresponding to
% the same vector in x_data
g1 = repmat({'ReHo-IED_lt'}, length(rho_reho_ied_lt_all_sess_fisher_z), 1);   % reho - ied_{long term}
g2 = repmat({'ReHo-IED_ds'}, length(rho_reho_ied_ds_all_sess_fisher_z), 1);   % reho - ied_{during scan}
g3 = repmat({'ALFF-IED_lt'}, length(rho_alff_ied_lt_all_sess_fisher_z), 1);   % alff - ied_{long term}
g4 = repmat({'ALFF-IED_ds'}, length(rho_alff_ied_ds_all_sess_fisher_z), 1);   % alff - ied_{during scan}
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
ax1.YLim = [-0.8 1.4];
ax1.YTick = [ax1.YLim(1):0.2:ax1.YLim(2)];

% set fontsize of axes
ax1.FontSize = 13;

% y-label
ylabel('Static correlation','Interpreter','latex');
% ax1.YLabel.FontSize = 20;

% obtain axes limits
xmax = ax1.XLim(2);
xmin = ax1.XLim(1);

% set aspect ratio of plot
pbr = 1/1.2;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

ax1.TickLabelInterpreter = 'latex';
% ax1.XTickLabel = {'ReHo-IED_lt', 'ReHo-IED_ds', 'ALFF-IED_lt', 'ALFF-IED_ds'};
% format x-tick label, '\:' is used to generate space in math mode in latex
ax1.XTickLabel = {'ReHo--IED\(_{\mathrm{long} \: \mathrm{term}}\)*', ...
    'ReHo--IED\(_{\mathrm{during} \: \mathrm{scan}}\)', ...
    'ALFF--IED\(_{\mathrm{long} \: \mathrm{term}}\)', ...
    'ALFF--IED\(_{\mathrm{during} \: \mathrm{scan}}\)*'};

% display messages on terminal showing mean fisher z and their p-values
display([rho_reho_ied_lt_all_sess_fisher_z_mean, rho_reho_ied_ds_all_sess_fisher_z_mean, ...
    rho_alff_ied_lt_all_sess_fisher_z_mean, rho_alff_ied_ds_all_sess_fisher_z_mean])

display([pval_rho_reho_ied_lt_all_sess_fisher_z, pval_rho_reho_ied_ds_all_sess_fisher_z, ...
    pval_rho_alff_ied_lt_all_sess_fisher_z, pval_rho_alff_ied_ds_all_sess_fisher_z])

% output plots to path
if op_results == 1
saveas(f1, fullfile(directname_op, filenameg1), 'epsc');
saveas(f1, fullfile(directname_op, filenameg1), 'jpeg');
end



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
function [fisher_z, mean_fisher_z, pval_fisher_z] = fisher_transform(input_array)

% compute quantities interested,
fisher_z = atanh(input_array);   % fisher transformed rho

% remove entries of inf in fisher_z for cases with rho = 1 or -1
ind_non_inf = ~isinf(fisher_z);

mean_fisher_z = mean(fisher_z(ind_non_inf));   % mean of fisher transformed rho

[~, pval_fisher_z] = ttest(fisher_z(ind_non_inf));   % p-value of fisher transformed rho

end
%---------------------------------------

