clear all 
close all

% this script does the following:
% (i) combine instantaneous IED rate (A) and ALFF (B) computed in every time 
% segment of each channel across sessions ....
% (ii) compute Spearman's corr. coeff. btw. the combined (A) and (B)
% (iii) compute coeff. of variation of (A)
% (iv) compute Fisher transformed (z-score) of (iii)

% Note: not all channels have (A) registerd

% quantities calculated are saved in a struct with format given below, 
% grand --- subject number --- session number --- quantities
% e.g. to access reho calculated at session 7 for subject 14, 
% grand.sub14.run7.reho 

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

tic 

% enter path to directory where info. of all subjects are stored
% subnum_dir = 'C:\Users\siumichael.tang\Downloads\fmri_project'; 
subnum_dir = '/work/levan_lab/mtang/fmri_project';

% enter path where input struct. is stored at
fname_input = [subnum_dir, filesep, 'matrices' filesep 'dynamic_analysis_alff'];   % direct. of output matrix
filename_input = 'dynamic_analysis_alff.mat';   % filename of input file

% enter output filename
fname_output = fname_input;
filename_output = 'dynamic_analysis_alff_fisher_z.mat';

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% use get_path to get full paths of files interested
[~, alff_struct_path, ~] = get_path(fname_input, filename_input);   % reho

% load sturct. from path to variable
st_alff = load(alff_struct_path{:});   % alff struct.

% obtain fieldnames of subjects
fd_sub = fieldnames(st_alff.grand);

% initialize struct, main, for storing outputs
 main = struct;

% go through each subject 
for sub_ind = 1:numel(fd_sub)

    % call combine_ar, defined at end of script, to get arrays concatenated
    % across sessions for current subject, for func for details
    % onset = ied onset channels, propa = ied propagation channels
    % ext = ied onset + propa channels (those with ied registered during
    % fmri scan)
    ch_onset_ds_ied_alff_comb = combine_ar(st_alff.grand.(fd_sub{sub_ind}), 'ch_onset_ds_ied_alff');
    ch_propa_ds_ied_alff_comb = combine_ar(st_alff.grand.(fd_sub{sub_ind}), 'ch_propa_ds_ied_alff');
    ch_ext_ds_ied_alff_comb = combine_ar(st_alff.grand.(fd_sub{sub_ind}), 'ch_ext_ds_ied_alff');

    % store arrays to output struct under curr. subject fieldname
    main.(fd_sub{sub_ind}).ch_onset_ds_ied_alff_comb = ch_onset_ds_ied_alff_comb;
    main.(fd_sub{sub_ind}).ch_propa_ds_ied_alff_comb = ch_propa_ds_ied_alff_comb;
    main.(fd_sub{sub_ind}).ch_ext_ds_ied_alff_comb = ch_ext_ds_ied_alff_comb;

    % output structure created in current session
    if op_results == 1
        save(fullfile(fname_output , filename_output), 'main');
    end

end   % end for sub_ind = 1:numel(fd_sub)


% USER DEFINED FUNCTIONS BELOW:
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% FUNCTION OUTPUT_AR = COMBINE_AR(INPUT_AR, AR_INT_STR)
% this function requires the following inputs, 
% input_st = input struct containing data from all subjects,
% ar_int_st = name (str) of array interested to be concatenated across
% sessions,
% and generates the following output
% output_ar = output array with ied rate and feature (reho/alff) in each
% segment of time series combined across sessions for every pair of ied type and channel
function output_ar = combine_ar(input_st, ar_int_str)

% initialize output array
output_ar = [];

% check if input struct is empty, if so, return control to line calling
% this func
if isempty(input_st)
    return
end

% obtain fieldnames of runs in curr. subject
fd_run = fieldnames(input_st);
    
% concatenate ar_int across all runs vertically 
ar_int_all = [];
for run_ind = 1:numel(fd_run)
    ar_int_curr = input_st.(fd_run{run_ind}).(ar_int_str);
    ar_int_all = [ar_int_all; ar_int_curr];
end

% in the array obtained above, select rows with the same ied type and
% channel name, concatenate their entries for inst. ied rate and feature
% (reho/alff)
ar_ds_ied_feat_comb = {};   % initialize array
row_num = 1;   % initialize row number for storing info

% each pair of ied type and channel name can only be used once,
% initialize array, ied_type_ch_used, to record pairs used in the
% selection process
ied_type_ch_used = {};
     
% scan through each row of ar_int_all array, obtain ied type
% from col. 2 and channle name from col. 3
for i = 1:size(ar_int_all, 1)
    ied_type = ar_int_all{i, 2};   % ied type from col. 2
    ch_name = ar_int_all{i, 3};   % channel name from col. 3

    % initialize variable to indicate if the pair of ied type and
    % channel name has been used 
    ied_type_ch_found = 0;

    % search if current pair of ied type and channel name exists in
    % ied_type_ch_used array
    for j = 1:size(ied_type_ch_used, 1)   
        if isequal(ied_type_ch_used{j, 1}, ied_type) && strcmp(ied_type_ch_used{j, 2}, ch_name)
            ied_type_ch_found = 1;   % if so, change value to 1
        end
    end

    % if the pair of ied type and channel name has not been used,
    % proceed to combine inst. ied rate and reho together, otherwise,
    % proceed to next row in ar_int_all (next pair of ied type and channel
    % name)
    if isequal(ied_type_ch_found, 0)

        % initialize arrays
        ds_ied_comb = [];   % ied_during scan in each segment across sessions
        feat_comb = [];   % feature (reho/alff) in each segment across sessions

        % find row indices with matching channel name in ar_int_all
        ind_req_ch = find(strcmp(ar_int_all(:, 3), ch_name));

        for item = 1:numel(ind_req_ch)   % for each row index identified
            % check if ied type of current row (col. 2) matches with curr.
            % ied type interested, if so, add info to array
            if isequal(ar_int_all{ind_req_ch(item), 2}, ied_type)
                ds_ied_comb = [ds_ied_comb ar_int_all{ind_req_ch(item), end-1}];
                feat_comb = [feat_comb ar_int_all{ind_req_ch(item), end}];
            end
        end

        % store info obtained across sessions to array initialized with
        % format given below, 
        % 1st col. = subject number
        % 2nd col. = ied type
        % 3rd col. = channel name
        % 4th col. = inst. ied type in every segment combined from all
        % sessions (A)
        % 5th col. = feature (reho/alff) in every segment combined from all sessions (B)
        % 6th col. = spearman's corr. coeff. btw. (A) and (B), called it (C)
        % 7th col. = coeff. of variation of (A)
        % 8th col. = fisher z-score of (C)
        ar_ds_ied_feat_comb(row_num, 1:3) = ar_int_all(i, 1:3);
        ar_ds_ied_feat_comb{row_num, 4} = ds_ied_comb;
        ar_ds_ied_feat_comb{row_num, 5} = cell2mat(feat_comb);
        ar_ds_ied_feat_comb{row_num, 6} = corr(ds_ied_comb', cell2mat(feat_comb)', 'Type', 'Spearman');
        ar_ds_ied_feat_comb{row_num, 7} = std(ds_ied_comb)/mean(ds_ied_comb);
        ar_ds_ied_feat_comb{row_num, 8} = atanh(ar_ds_ied_feat_comb{row_num, 6});

        % increment row number by 1 for next pair of ied type and channel name
        row_num = row_num + 1;  

        % update array with current pair of ied type and channel name so it
        % will not be used again 
        ied_type_ch_used = [ied_type_ch_used; {ied_type ch_name}];
    end   % if isequal(ied_type_ch_found, 0)

end   % end for i = 1:size(ar_int_all, 1)

% update output array (concatenation of ied rate and feature in every
% segment across sessions for each pair of ied type and channel
output_ar = ar_ds_ied_feat_comb;

end   % end function output_ar = combine_ar(input_st, ar_int_str)
