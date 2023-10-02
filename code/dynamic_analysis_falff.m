clear all 
close all

% this script does the following dynamic analysis:
% (i) get fractional ALFF in every segment of time series of each channel in every session
% (ii) get IED rate in every segment of time series measured during scan

% Note: not all channels have (ii) registerd

% quantities calculated are saved in a struct with format given below, 
% grand --- subject number --- session number --- quantities
% e.g. to access falff calculated at session 7 for subject 14, 
% grand.sub14.run7.falff

% Oct. 12, 2022: updated list of "ied onset" and "ied propagating" channels
% is used and can be found here:
% /work/levan_lab/mtang/fmri_project/ied_onset_ied_propa_ch.txt

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

tic 

% enter path to directory where info. of all subjects are stored
% subnum_dir = 'C:\Users\siumichael.tang\Downloads\fmri_project'; 
subnum_dir = '/work/levan_lab/mtang/fmri_project';

% enter window size in seconds (must be divisible by tr, i.e. 1.5 s)
window_size = 120;

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% (A): identify list of subject numbers with path input

% use dir to get a struct for all folders/files found in subject direct
folder_list = dir(subnum_dir);

% access name of every folder/file identified, then search for pattern of
% 'sub\d*' in it, if matched, extract the subject number
sub_list = [];   % initialize array for storing subject numbers
for item = 1:numel(folder_list)   % every item in folders/files found
    
    % search filename with matching pattern 'sub\d*', then use strrep to
    % replace 'sub' by '' to obtain subject number
    matched_str = regexp(folder_list(item).name, '^sub\d*$', 'match');
    matched_sub_num = strrep(matched_str, 'sub', '');
    
    % if pattern is found, add subject number to sub_list
    if ~isempty(matched_str)
        sub_list = [sub_list matched_sub_num];
    end
end

% END (A): identify list of subject numbers with path input
%---------------------------------------------------------------------------


% (B): compute (i) ReHo and (ii) IED rate measured during long-term monitoring,
% and (iii) IED rate measured over entire scan, and get spearman's rho btw.
% (i) and (ii), and (i) and (iii), for every session of each subject

% initialize structure, grand, to store relevant info. of every subject
grand = struct;

% loop through each subject number identified and conduct analysis
for item = 1:numel(sub_list)

% get current subject number
subnum = sub_list{item};

% sprintf(['begin working on sub', subnum])

%---------------------------------------------------------------------------
% (B1): format paths of required files for current subject

% enter path to directory where all input files are located
directname = [subnum_dir, filesep, 'sub', subnum];
% directname = ['C:\Users\siumichael.tang\Downloads\fmri_project\', 'sub', subnum];

% format filenames of processed fmri images, explicit mask, electrode,
% motion parameters
filename_swraimg = ['swra*.nii'];   % processed func. images
filename_expmask = 'wEPI_bet_mask.nii';   % explicit mask
filename_elect =  [subnum, '_*Koordinaten*.xlsx'];   % file containing mni coord. of all electrode pairs
filename_event_onset = ['find_onset_times_sub', subnum, '.mat'];   % file of event onsets   
filename_motion = ['rp_*.txt'];   % motion parameters generated by SPM in the realignment step

% enter path where ouput struct. is stored at
fname_op = [subnum_dir, filesep, 'matrices' filesep 'dynamic_analysis_falff'];   % direct. of output matrix
filename_op = 'dynamic_analysis_falff.mat';   % filename of output file

% END (B1): format paths of required files for current subject
%---------------------------------------------------------------------------
% (B2): obtain full paths of required files for current subject

% use get_path function to obtain full path of input files
[~, swra_file_path, ~] = get_path(directname, filename_swraimg);   % swra func. images
[~, expmask_file_path, ~] = get_path(directname, filename_expmask);   % normalized mask images
[~, elect_file_path, ~] = get_path(directname, filename_elect);   % file containing info. about electrodes coordinates
[~, event_onsets_file_path, ~] = get_path([directname, filesep, 'matrices'], filename_event_onset);   % file of event onsets
[~, motion_file_path, ~] = get_path(directname, filename_motion);   % motion parameters

% END (B2): obtain full paths of required files for current subject
%---------------------------------------------------------------------------

% (B3): conduct analysis for current subject

% createterms structure to store info. of every session of each subject
terms = [];

% loop through every session of current subject and calculate relevant
% quantities

for run_ind = 1:numel(swra_file_path)   % for every session

    %---------------------------------------------------------------------------
    
    % PART (0): read input files for current session of current subject 
   
    swra_img = niftiread(swra_file_path{run_ind});   % swra func. images
    swra_info = niftiinfo(swra_file_path{run_ind});   % info about this particular swra func. images
    expmask = single(niftiread(expmask_file_path{:}));   % convert image data type to single
    electar = readcell(elect_file_path{:});   % cell array containing info. about electrodes
    event_onsets_st = load(event_onsets_file_path{:});   % struct. containing event onsets

    % access info for current session
    field_runs = fieldnames(event_onsets_st.terms);   % fieldname under event onsets structure
    event_onsets_curr_sess = event_onsets_st.terms.(field_runs{run_ind});   % array of event onsets for current session
    
    motion_file_path_curr_sess = motion_file_path{run_ind};   % motion parameter path for current session
    
    % END PART (0): read input files for current session of current subject 
    %---------------------------------------------------------------------------
    % Part (I): remove motion artifacts from processed images

    % convert explicit mask and processed image to double, 
    % apply element-wise multiplciation
    exp_swra_img = double(expmask).*double(swra_img);
                
    % call filter_motion.m to remove motion artifacts from image
    motion_filtered_img = filter_motion(motion_file_path_curr_sess, exp_swra_img);

    % END Part (I): remove motion artifacts from processed images
    %---------------------------------------------------------------------------
    
    % Part (II): get fractional amplitude of low-frequency fluctuation (fALFF) map (every voxel)
    % in every segment 

    % repetition time, in seconds
    tr = swra_info.PixelDimensions(4); 

    % use get_falff to get map of normalized fractional alff in every segment 
    norm_falff = get_falff(motion_filtered_img, expmask, tr, window_size);

    % if empty normalized fractional alff map is obtained, time series of 
    % current session is not long enough for making psd estimates using 
    % pwelch method, so, pass control to next session, output empty array 
    % for current session
    if isempty(norm_falff)
        continue
    end

    % END Part (II): get fractional amplitude of low-frequency fluctuation (fALFF) in
    % every segment
    %---------------------------------------------------------------------------

    % Part (III): sort electrodes by their types

    % use sort_elect to sort electrode array by their types
    ele_sorted = sort_elect(electar);

    % END PART (III): sorte electrodes by their types
    %---------------------------------------------------------------------------

    % % Part (IV): find midpoint of each pair of electrode contacts in both mni
    % and image space

    % use get_midpt_elect to locate midpoint of all electrode channels
    % col. 1 = name of channel, 
    % col. 2 = midpt in mni space, 
    % col. 3 = midpt in image space
    midpt = get_midpt_elect(ele_sorted, swra_info);

    % cell array midpt groups electrodes of same type into one cell, 
    % convert midpt to single layer cell array, with same column spec. 
    midpt_all = {};   % initialize cell array
    rownum = 1;   % initiaize row number

    for i = 1:numel(midpt)   % for each cell (each electrode type)
        for j = 1:size(midpt{i}, 1)   % for each row in current cell
            midpt_all(rownum, :) = midpt{i}(j, :);   % assign all col. to array
            rownum = rownum + 1;   % increase row number by 1 for next row read
        end
    end

    % END Part (IV): find midpoint of each pair of electrode contacts in both mni
    % and image space
    %---------------------------------------------------------------------------

    % Part (V): get names and event types of ied onset and ied
    % propagation channels 

    % note: 
    % event type of 'NA' is assigned to no-ied channels as they have
    % no ieds registered during fmri scan

    % use get_ied_onset_ied_propa to get array containing names of 
    % ied onset and ied propagation channels, as well as their event types
    % ch_list_onset = ied onset channels
    % ch_list_propa = ied propagation channels
    % ch_list_ext = ch_list_onset + ch_list_propa, any channels with ieds
    % registered during fmri scan
    [ch_list_onset, ch_list_propa, ch_list_ext] = get_ied_onset_ied_propa_ch(subnum);

    % END Part (V): get names and coordinates of channels with ied rates
    % registered during fmri scan
    % %---------------------------------------------------------------------------

    % Part (VI): select channels in midpt_all with ieds registered during
    % scan

    % midpt_all lists all channels available, ch_list_ext lists channels 
    % that have ieds registered during scan, 
    % find intersection btw. the two arrays and select corresponding rows in 
    % midpt_all to obtain their midpoint in image space
    [~, ia, ~] = intersect(midpt_all(:, 1), ch_list_ext(:, end));

    % select rows in midpt_all for channels with ieds registered during
    % scan, store info in midpt_all_sel array
    midpt_all_sel = midpt_all(ia, :);

    % END Part (VI): select channels in midpt_all with ieds registered during scan
    % %---------------------------------------------------------------------------

    % % Part (VII): construct 3X3X3 boxes centered at midpoint of each electrode
    % % pair to store signals enclosed thereof
 
    % get basic info. about the processed image
    tr = swra_info.PixelDimensions(4);   % repetition time, in seconds
    fs = 1/tr; % sampling frequency, in Hz

    % % initialize array, sig_box, for storing the following for every type of elect.
    % % 1st column lists the names of channels,
    % % 2nd column lists mean normalized alff in each time segment of every channel
    sig_box_all = {};

    % use for loop to get signals within each box from motion-filtered image
    for i = 1:size(midpt_all_sel, 1)
        sig_box_all{i, 1} = midpt_all_sel{i, 1};   % get name tag of electrode pair
        mxyz = midpt_all_sel{i, end};   % get midpoint in image space

        % create regional mask from explicit mask using the coordinates of
        % the 3X3X3 box
        reg_mask = double(expmask(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
            mxyz(3)-1:mxyz(3)+1));

        % get number non-zero entries in explicit mask in current box
        num_nonzero = nnz(reg_mask);

        % get number of time segment at current midpoint, should be the
        % same for all channels
        num_seg = numel(norm_falff{mxyz(1), mxyz(2), mxyz(3)});

        % for every segment in time, scan across every voxel in current
        % 3X3X3 box, (i) sum all normalized frac. alff values thereof, 
        % (ii) then get mean normalized frac. alff in that box by dividing
        % sum by number of non-zero entries of that box in explicit mask
        for seg_ind = 1:num_seg   % for every segment in time
            seg_sum_norm_falff = 0;   % initialize sum of norm. frac. alff as 0 
            for x_vox = mxyz(1)-1:mxyz(1)+1   % scan across every voxel in curr. box
                for y_vox = mxyz(2)-1:mxyz(2)+1
                    for z_vox = mxyz(3)-1:mxyz(3)+1
                        % if current voxel is included in explicit mask
                        if expmask(x_vox, y_vox, z_vox) == 1                 
                            % compute sum of all norm. frac. alff values in curr. box
                            seg_sum_norm_falff = seg_sum_norm_falff + ...
                                norm_falff{x_vox, y_vox, z_vox}{seg_ind};
                        end
                    end
                end
            end

            % compute mean normalized frac. alff in curr. box by dividing
            % the sum by number of non-zero entries of curr. box. in explicit 
            % mask, assign calculated mean to curr. cell in 2nd col. of
            % sig_box_all
            sig_box_all{i, 2}{seg_ind} = seg_sum_norm_falff/num_nonzero;
        end   % end for each segment in time

    end   % for i = 1:size(midot_all_sel, 1), for each channel

    % END Part (V): construct 3X3X3 boxes centered at midpoint of each electrode
    % pair to store signals enclosed thereof
    %---------------------------------------------------------------------------

    % Part (VI): obtain instantaneous ied rate [/min] of every event type
    % in current session

    % obtain instantaneous ied rate [/min] of every event type in curr.
    % session
    inst_ied_st = get_instant_ied_rate(swra_img, tr, window_size, event_onsets_curr_sess);

    % obtain instantaneous ied rate [/min] from inst_ied_st
    inst_ied = inst_ied_st.inst_ied;

    % END Part (VI): obtain instantaneous ied rate [/min] of every event type
    % in current session
    %---------------------------------------------------------------------------

    % Part (VII): assign instantaneous ied rate (per minute) and mean norm.
    % fractional alff in every time segment to 
    % corresponding channel (ied onset, ied propagation, and others)

    % use match_ds_ied_falff, defined at the end of script, to assign during 
    % scan ied rate (per minute) and frac. alff measured over entire scan to every 
    % channel listed in each of the arrays obtained above (ied onset, ied propagation, 
    % and ext = onset + propa, as only channels with ieds registered during
    % scan is considered in this case)
    ch_onset_ds_ied_falff = match_ds_ied_falff(ch_list_onset, inst_ied, sig_box_all);
    ch_propa_ds_ied_falff = match_ds_ied_falff(ch_list_propa, inst_ied, sig_box_all);
    ch_ext_ds_ied_falff = match_ds_ied_falff(ch_list_ext, inst_ied, sig_box_all);

    % END Part (VII): assign instantaneous ied rate (per minute) and mean
    % norm. frac. alff in every time segment to 
    % corresponding channel (ied onset, ied propagation, and others)
    %---------------------------------------------------------------------------

    % Part (VIII): store quantities calcu. in output struct.

    opstruct = struct;
    opstruct.ele_sorted = ele_sorted;   % info. of electrode (sorted by name)
    opstruct.midpt_all = midpt_all;   % coordinates of midpoint between each pair of electrode contacts in single layer
    opstruct.sig_box_all = sig_box_all;   % signals within 3X3X3 box centered at each midpoint of all channels in single laye

    % inst. ied rate [/min] in each time segment for every event type
    opstruct.inst_ied = inst_ied;   

    % arrays listing channels in each category and their event type 
    opstruct.ch_list_onset = ch_list_onset;   % onset 
    opstruct.ch_list_propa = ch_list_propa;   % propa 
    opstruct.ch_list_ext = ch_list_ext;   % ext = onset + propa
    
    % ied rate and mean normalized frac. alff of every time segment of every channel
    % that has ied registed during scan, divided into categories
    opstruct.ch_onset_ds_ied_falff = ch_onset_ds_ied_falff;   % onset 
    opstruct.ch_propa_ds_ied_falff = ch_propa_ds_ied_falff;   % propa
    opstruct.ch_ext_ds_ied_falff = ch_ext_ds_ied_falff;   % ext = onset + propa

    % END Part (VIII): store quantities calcu. in output struct.
    %---------------------------------------------------------------------------
    %---------------------------------------------------------------------------

    % get indices of current run number embedded in processed image filename
    sind = regexp(lower(swra_file_path{run_ind}), 'run');   % start index
    eind = regexp(lower(swra_file_path{run_ind}), '.nii');   % end index
    sess_num = swra_file_path{run_ind}(sind+3:eind-1);   % get indices of run number

    % create fieldname representing current session number
    fieldname = sprintf('run%s', sess_num);

    % update terms struct with info. of current session
    terms.(fieldname) = opstruct;
   
    % print message to terminal
    sprintf(['done working on run', num2str(run_ind), ' of sub', subnum])

end   % end for run_ind = 1:numel(swra_file_path)   % for every session

%---------------------------------------------------------------------------

% format fieldname for current subject number
g_fdname = ['sub', subnum];   % format fieldname for current subject number
grand.(g_fdname) = terms;   % update grand, struct, with info of curr. sub.

% output structure created in current session
if op_results == 1
    save(fullfile(fname_op, filename_op), 'grand');
end

sprintf(['done working on sub', subnum])

end   % end for item = 1:numel(sub_list)

toc


% USER DEFINED FUNCTIONS FOR THIS SCRIPT:
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% FUNCTION OUTPUT_AR = MATCH_DS_IED_FALFF(INPUT_AR, IED_DS, SIG_BOX_ALL)
% this function requires the following inputs,
% INPUT_AR = array containing channels to which during scan ied and reho
% are assigned, with format given below
% 1st col. = subject number
% 2nd col. = event type
% 3rd col. = channel name
% IED_DS = array containing during scan ied rate per minute, with format
% given below, 
% 1st col. = event type
% 2nd col. = ied rate per minute
% SIG_BOX_ALL_AR = array containing reho of all channels available, with
% format given eblow, 
% 1st col. = channel id
% 2nd col. = channel name
% 3rd col. = feautre interested (e.g. ReHo, ALFF, etc.)
% this function generates the following output, 
% OUTPUT_AR = output array with during scan ied rate and frac. alff assigned to
% each channel listed in INPUT_AR, with format given below
% 1st col. = subject number
% 2nd col. = ied type 
% 3rd col. = channel name
% 4th col. = during scan ied rate per minute
% 5th col. = time series feature (e.g. ReHo, ALFF)

function output_ar = match_ds_ied_falff(input_ar, ied_ds, sig_box_all_ar)

% initialize output_ar
output_ar = {};

% if input_ar is empty, return control to line calling this function
if isempty(input_ar)
    return
end

% initialize arrary storing info required, format is given below
% 1st col. = subject number
% 2nd col. = event type
% 3rd col. = channel name
% 4th col. = during scan ied rate
ch_ds_ied = {};
row_num = 1;    % initialize row number 
ncol = size(input_ar, 2);   % get number of col. established in input_ar

for i = 1:size(input_ar, 1)   % for every channel in input_ar
    % obtain curr. event type, convert it from str to num
    cur_ev_type = str2double(input_ar{i, 2});
    for j = 1:size(ied_ds, 1)   % for every row in ied_during_scan array
        if ismember(cur_ev_type, ied_ds{j, 1})   % find row with matching event type
            % assign all info related to the channel to current row
            ch_ds_ied(row_num, 1:ncol) = input_ar(i, :);
            % assign during-scan ied rate (per minute) to rightmost
            % col. of current row
            ch_ds_ied{row_num, ncol+1} = ied_ds{j, end};
            row_num = row_num + 1;   % increment row number by 1 for next matched channel
        end
    end
end

% get number of col. established in ch_ds_ied array 
ncol = size(ch_ds_ied, 2);

% initialize ch_ds_ied_falff, with the same format as
% ch_ds_ied, but with an additional column storing frac. alff
% computed over the entire scan
ch_ds_ied_falff = {};
row_num = 1;   % initialize row number as 1

% assign frac. alff measured over entire scan to each channel that has
% matched ied type found, and matched channel name in sig_box_all_ar
for i = 1:size(ch_ds_ied, 1)   % for every row in ch_ds_ied
    ch_name_cur = ch_ds_ied{i, 3};   % get current channel name from col. 3
    for j = 1:size(sig_box_all_ar, 1)   % for every row in sig_box_all_ar
        if isequal(ch_name_cur, sig_box_all_ar{j, 1})   % find row with matching channel name
            % assign alff of matched channel to rightmost col. of current row
            ch_ds_ied_falff(row_num, 1:ncol) = ch_ds_ied(i, :);
            ch_ds_ied_falff{row_num, ncol+1} = sig_box_all_ar{j, end};
            row_num = row_num + 1;   % increment row number by 1 for next matched channel
        end
    end
end

% assign variable computed to output
output_ar = ch_ds_ied_falff;

end   % end function output_ar = match_ds_ied_falff(input_ar, ied_ds, sig_box_all_ar)

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
