clear all 
close all

% use this script to merge all tables generated from another script
% alff_spikes.m
tic 
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '14';

% enter path to directory where all input files are located
directname = ['/work/levan_lab/mtang/fmri_project/', 'sub', subnum, ...
    '/matrices/alff_spikes'];

% enter format of csv filename
filename_csv = ['sub', subnum, '*winsize*.csv'];

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% locate paths of func. images and their details 

% use get_path function to obtain full path of input files, output cell
% array containing paths of all matched files
[~, csv_file_path, ~] = get_path(directname, filename_csv); 

% identify session numbers embedded in paths
% which starts with run%d_winsize_%d.csv
sess_num = [];   % initialize array for storing sess. num.
for i = 1:numel(csv_file_path)   % for each entry in path cell array
    % use regexp to look for indices matching the patterns
    sind = regexp(csv_file_path{i}, 'run\d*');   % start index
    eind = regexp(csv_file_path{i}, '_winsize');   % end index
    % use indices to obtain current session number
    sess_num = [sess_num; str2num(csv_file_path{i}(sind+(eind-1-sind):eind-1))];
%     sess_num = [sess_num; str2num(csv_file_path{i}(sind+3:eind-1))];
end

% get unique session number of all paths
unique_sess_num = unique(sess_num);

% group paths by session number,
% initialize array, each cell correponds to paths within certain session
% number
path_grouped = {};   
for run = 1:numel(unique_sess_num)   % for each unique sess. num.
    % find cell indices with sess. num. equal to current unique sess. num,.
    ind_req = find(sess_num == unique_sess_num(run));   
    % select those entries and store them in current cell
    path_grouped{run} = csv_file_path(ind_req);   
end

% next, sort window size embedded in paths by their session number
winsize_num = cell(numel(path_grouped), 1);   % initialize array
for run = 1:numel(path_grouped)   % for each sess. num.
    for j = 1:numel(path_grouped{run})   % for each entry in current sess. num.
        % every path has pattern 'winsize_%d.csv', look for relevant
        % indices
        sind_ws = regexp(path_grouped{run}{j}, 'winsize\d*');   % start index, window size
        eind_ws = regexp(path_grouped{run}{j}, '.csv');   % end index, window size
        % use indices obtained to select window size for current entry
        winsize_num{run}{j} = str2num(path_grouped{run}{j}(sind_ws+8:eind_ws-1));
    end
end

% with paths grouped by session number and the corresponding window sizes
% identified, for each cell (sess. num.), sort window sizes in acending
% order, then get equivalent indices in unsorted array, use the indices to
% arrange paths 
path_sorted = cell(numel(path_grouped), 1);   % initialize array
for run = 1:numel(winsize_num)   % for each sess. num.
    % sort window sizes in current sess. in acending order
    [~, ind_ascend] = sort([winsize_num{run}{:}]);
    % use each index obtained to arrange paths 
    for i = 1:length(ind_ascend)
        path_sorted{run}{i, 1} = path_grouped{run}{ind_ascend(i)};
    end
end

% merge tables in each session together and store resultant product under a
% struct named t_struct
t_struct = struct;   % initialize structure

for run = 1:numel(path_sorted)   % for each session
    t_merged = [];   % initialize table for current session
    for item = 1:numel(path_sorted{run})   % for each path identified
        % get current table in current session
        cur_t = readtable(path_sorted{run}{item}, 'Delimiter', ',', ...
            'PreserveVariableNames', true);  
        % if current table is the first table in current session, append
        % all columns to t_merged, else, append 3rd and 4th col. to
        % t_merged (because first two columns are event type and channel
        % names, which are duplicate)
        if item == 1
            t_merged = [t_merged cur_t];
        else
            t_merged = [t_merged cur_t(:, 3:4)];
        end
    end
    % construct fieldname for current session
    fieldname = ['run', num2str(unique_sess_num(run))];
    % add table for current session to struct
    t_struct.(fieldname) = t_merged;
end

% create array for storing paths of merged tables (output of this script)
merged_table_path = cell(numel(path_sorted), 1);

% format filenames and full paths of tables
fname_op = directname;
for run = 1:numel(path_sorted)
    filename_table = ['sub', subnum, '_table_ns_mnalff_run', ...
        num2str(unique_sess_num(run)), '_merged', '.csv'];
    merged_table_path{run, 1} = fullfile(fname_op, filename_table); 
end

%---------------------------------------------------------------------------
% output table and structure created in current session
if op_results == 1
    fdname = fieldnames(t_struct);   % get all fieldnames
    for run = 1:numel(fdname)   % for each run, write table to path
        writetable(t_struct.(fdname{run}), merged_table_path{run});   % table of clin. det. ch.
    end
end
%---------------------------------------------------------------------------


