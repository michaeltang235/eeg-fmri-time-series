% This script extracts data from .mat file generated from the analysis of
% spearman's corr. coeff. btw. dynamic functional connectivity (dfc) and 
% number of spikes, and assemble an array with format given below for every subject,
% col. 1 = event type
% col. 2 = ref. channel name
% col. 3 = target channel name
% col. 4 = rho_all runs
% col. 5 = p-value
% col. 6 = Fisher's z all runs (computed in this script)
% col. 7 = mean number of spikes, all runs
% col. 8 = standard deviation of number of spikes, all runs

% the assembled array for each subject is then stored under a struct, with
% format given below
% struct.(subject number) = array assembled

clear all
close all

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter list of subject numbers interested
subnum_ar = {'14', '15', '18', '19', '20', '27', '30', '31', '32', '33', '34', ...
    '35', '41'};

% enter directory where ouput file is saved
op_dir = '/work/levan_lab/mtang/fmri_project';

% enter filename of output file
op_filename = 'fisher_z_dfc_all_subs.mat';

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% initialize ouput structure
opst = struct;

% loop through each sub. num. and extract relevant info.
for item = 1:numel(subnum_ar)
    
    % get current subject number
    subnum = subnum_ar{item};
    
    % format fieldname using current subject number
    fieldname = ['sub', subnum];

    % call function group_fisher_z to get array required, store var. under
    % current fieldname
    opst.(fieldname) = group_fisher_z(subnum);
end

% write output struct. to path
save(fullfile(op_dir, op_filename), 'opst', '-v7.3');

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

function array_req = group_fisher_z(subnum)

% obtain sub. num. from input
subnum = subnum;

% initialize output variable
array_req = [];

%-----------------------------------------------------
% PART (I): format paths of input files

% format path to directory where all input files are located
directname = ['/work/levan_lab/mtang/fmri_project/', 'sub', subnum];

% format path where struct. files are stored at
fname_dfc_struct = [directname, filesep, 'matrices' filesep 'dfc_spikes_fixed_winsize'];   % direct. of alff struct is located
filename_dfc_struct = 'dfc_spikes_fixed_windsize.mat';   % filename of alff struct. file

fname_spikes_stats_struct = [directname, filesep, 'matrices' filesep 'num_spikes_stats'];   % direct. of spikes stats matrix
filename_spikes_stats_struct = 'num_spikes_stats.mat';   % filename of spikes stats struct. file

% use get_path to get full paths of files interested
[~, dfc_struct_path, ~] = get_path(fname_dfc_struct, filename_dfc_struct);   % dfc
[~, spikes_stats_struct_path, ~] = get_path(fname_spikes_stats_struct, filename_spikes_stats_struct);   % spikes stats

% if any of the paths is empty, return function
if isempty(dfc_struct_path) || isempty(spikes_stats_struct_path)
    return
end

% END PART (I): format paths of input files
%-----------------------------------------------------

% PART (II): load var. interested from paths

% load sturct. from path to variable
st_dfc = load(dfc_struct_path{:});   % dfc struct.
st_spikes_stats = load(spikes_stats_struct_path{:});   % spikes stats struct.

% access fieldnames of struct. terms
fds_dfc = fieldnames(st_dfc.terms);   % dfc struct.
fds_spikes_stats = fieldnames(st_spikes_stats.terms);   % spikes stats struct.

% get non-empty fieldnames for dfc and spikes_stats struct.
fd_ind_dfc = [];   % initialize array for storing indices of non-empty fields
for sess_ind = 1:length(fds_dfc)   % for each field
    if ~isempty(st_dfc.terms.(fds_dfc{sess_ind}))   % check if field is empty
        fd_ind_dfc = [fd_ind_dfc sess_ind];   % if not, append index to array
    end
end

fd_ind_spikes_stats = [];   % initialize array for storing indices of non-empty fields
for sess_ind = 1:length(fds_spikes_stats)   % for each field
    if ~isempty(st_spikes_stats.terms.(fds_spikes_stats{sess_ind}))   % check if field is empty
        fd_ind_spikes_stats = [fd_ind_spikes_stats sess_ind];   % if not, append index to array
    end
end

% get common indices between the two arrays
common_ind = intersect(fd_ind_dfc, fd_ind_spikes_stats);

% access all non-empty fields using indices acquired
fdnames_dfc = fds_dfc(common_ind);   % dfc
fdnames_spikes_stats = fds_spikes_stats(common_ind);   % spikes stats

% variable interested: 
% spearman's rho and p-value btw. dfc and number of spikes (A),
% mean and standard deviation of number of spikes across all segments (B),
% in all runs combined
% note: all runs means combininig segments from all runs together and
% calculate A and B

% A is stored in cell array called rho_pval_comb_all, which is saved
% under EVERY field in struct., so access first field under struct., and
% obtain the cell array. Similar procedure is applied to B
rho_pval_comb_all = st_dfc.terms.(fdnames_dfc{1}).rho_pval_comb_all;   % spearman rho and p-value
spikes_stats_all = st_spikes_stats.terms.(fdnames_spikes_stats{1}).spikes_stats_all;   % stats of number of spikes 

% END PART (II): load var. interested from paths
%-----------------------------------------------------

% PART (III): assemble array required

% with A and B obtained, assemble the array, with format given below
% col. 1 = event type
% col. 2 = channel name
% col. 3 = rho_all runs
% col. 4 = p-value
% col. 5 = Fisher's z all runs
% col. 6 = mean number of spikes, all runs
% col. 7 = standard deviation of number of spikes, all runs

% assiagn first two columns of rho_pval_comb_all to array_req
array_req = rho_pval_comb_all(:, 1:2);   % assign event type and assoc. channel names

% get number of col. established
colct = size(array_req, 2);

% assign the 2 rightmost columns (rho and p-value from all runs combined) 
% to array req. 
array_req(:, colct+1:colct+2) = rho_pval_comb_all(:, end-1:end);

% get number of col. established
colct = size(array_req, 2);

% then, compute fisher's z for every spearman's rho obtained
for i = 1:size(array_req, 1)
    array_req{i, colct+1} = atanh(array_req{i, 3});
end

% get number of col. established
colct = size(array_req, 2);

% next, assign B (mean and standard dev.) to the next two col. of array
% req.
% note: B are stored in the 2 rightmost columns of spikes_stats_all
for i = 1:size(array_req, 1)   % for every row in array req.
    cur_ev = array_req{i, 1};   % get current event type
    for j = 1:size(spikes_stats_all, 1)   % for every row in spikes stats array
        % if event type in spikes stats equals to current event type
        % assign mean and standard dev. of all seg. from all runs combined
        if isequal(spikes_stats_all{j, 1}, cur_ev)
            array_req{i, colct+1} = spikes_stats_all{j, end-1};   % mean num. of spikes
            array_req{i, colct+2} = spikes_stats_all{j, end};   % standard dev. of num. of spikes
        end
    end
end

% print message to terminal
disp([subnum, ' done'])

% END PART (III): assemble array required
%-----------------------------------------------------


end   % end function array_req = group_fisher_z(input_array)

