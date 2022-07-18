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
% The Spearman's correlation coefficient (rho) and p-value between 
% (i) and (ii), called it (rho_var_ied_lt), 
% (i) and (iii), called it (rho_var_ied_ds),
% where var denotes the fmri feature interested, was calculated for each session.
% This script intakes rho and p-value from every session and display them
% on a table (.csv)
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
directname_op_reho = 'C:\Users\siumichael.tang\Downloads\fmri_project\matrices\static_analysis_reho';
directname_op_alff = 'C:\Users\siumichael.tang\Downloads\fmri_project\matrices\static_analysis_alff';

% format output filename
filename_table_reho = ['static_analysis_reho_table.csv'];
filename_table_alff = ['static_analysis_alff_table.csv'];

% enter if user wants to write plots to file (1=yes, 0=no)
op_results = 1;

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

% initialize array with format given below
% col. 1 = subject number
% col. 2 = session number
% col. 3 = ReHo-IED_{long term}, (rho_reho_ied_lt) (A),
% col. 4 = p-value of col. 3
% col. 5 = ReHo-IED_{during scan}, (rho_reho_ied_ds) (B),
% col. 6 = p-value of col. 5
rho_reho_ied_all_sess = {};

% access required arrays from structure, reho
row_num = 1;

for sub_ind = 1:numel(fd_name_sub_reho)   % for each subject in reho struct.
    if ~isempty(st_reho.grand.(fd_name_sub_reho{sub_ind}))   % if value is not empty
        fd_run_cur = fieldnames(st_reho.grand.(fd_name_sub_reho{sub_ind}));   % get run fieldnames available
        for run_ind = 1:numel(fd_run_cur)   % for each run
            
            % access sub-structure of curr. sess.
            st_reho_cur = st_reho.grand.(fd_name_sub_reho{sub_ind}).(fd_run_cur{run_ind});
            
            % assign quantities to current row according to format given
            % above
            rho_reho_ied_all_sess{row_num, 1} = fd_name_sub_reho{sub_ind};   % sub. num.
            rho_reho_ied_all_sess{row_num, 2} = fd_run_cur{run_ind};   % sess. num.
            rho_reho_ied_all_sess{row_num, 3} = round(st_reho_cur.rho_reho_ied_lt, 3);   % (A)
            rho_reho_ied_all_sess{row_num, 4} = round(st_reho_cur.pval_reho_ied_lt, 3);   % p-value of (A)
            rho_reho_ied_all_sess{row_num, 5} = round(st_reho_cur.rho_reho_ied_ds, 3);   % (B)
            rho_reho_ied_all_sess{row_num, 6} = round(st_reho_cur.pval_reho_ied_ds, 3);   % p-value of (B)
            row_num = row_num + 1;   % increment row number by 1 for next sess. with non-empty value assigned               
        end
    end
end

% initialize array with format given below
% col. 1 = subject number
% col. 2 = session number
% col. 3 = ALFF-IED_{long term}, (rho_alff_ied_lt) (C),
% col. 4 = p-value of col. 3
% col. 5 = ALFF-IED_{during scan}, (rho_alff_ied_ds) (D),
% col. 6 = p-value of col. 5
rho_alff_ied_all_sess = {};

% access required arrays from structure, alff
row_num = 1;

for sub_ind = 1:numel(fd_name_sub_alff)   % for each subject in alff struct.
    if ~isempty(st_alff.grand.(fd_name_sub_alff{sub_ind}))   % if value is not empty
        fd_run_cur = fieldnames(st_alff.grand.(fd_name_sub_alff{sub_ind}));   % get run fieldnames available
        for run_ind = 1:numel(fd_run_cur)   % for each run
            
            % access sub-structure of curr. sess.
            st_alff_cur = st_alff.grand.(fd_name_sub_alff{sub_ind}).(fd_run_cur{run_ind});
            
            % assign quantities to current row according to format given
            % above
            rho_alff_ied_all_sess{row_num, 1} = fd_name_sub_alff{sub_ind};   % sub. num.
            rho_alff_ied_all_sess{row_num, 2} = fd_run_cur{run_ind};   % sess. num.
            rho_alff_ied_all_sess{row_num, 3} = round(st_alff_cur.rho_alff_ied_lt, 3);   % (A)
            rho_alff_ied_all_sess{row_num, 4} = round(st_alff_cur.pval_alff_ied_lt, 3);   % p-value of (A)
            rho_alff_ied_all_sess{row_num, 5} = round(st_alff_cur.rho_alff_ied_ds, 3);   % (B)
            rho_alff_ied_all_sess{row_num, 6} = round(st_alff_cur.pval_alff_ied_ds, 3);   % p-value of (B)
            row_num = row_num + 1;   % increment row number by 1 for next sess. with non-empty value assigned               
        end
    end
end

% convert cell arrays to tables
table_rho_reho_ied_all_sess = cell2table(rho_reho_ied_all_sess, "VariableNames", ...
    ["sub. num.", "sess. num.", "rho-reho-ied-lt", "p-val-reho-ied-lt", ...
    "rho-reho-ied-ds", "p-val-reho-ied-ds"]);

table_rho_alff_ied_all_sess = cell2table(rho_alff_ied_all_sess, "VariableNames", ...
    ["sub. num.", "sess. num.", "rho-alff-ied-lt", "p-val-alff-ied-lt", ...
    "rho-alff-ied-ds", "p-val-alff-ied-ds"]);

% format path where each table will be saved at
table_reho_path = fullfile(directname_op_reho, filename_table_reho);   % reho
table_alff_path = fullfile(directname_op_alff, filename_table_alff);   % alff

% output tables to their paths
if op_results == 1
    writetable(table_rho_reho_ied_all_sess, table_reho_path);   % reho
    writetable(table_rho_alff_ied_all_sess, table_alff_path);   % alff
end