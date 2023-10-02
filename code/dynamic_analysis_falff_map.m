clear all 
close all

% this script does the following to assist with making the plot of falff
% overlaid on anatomical image in dynamic analysis
% (i) get falff in every segment of time series of each channel in every session
% regardless of whether the channel has during-scan IED registered. This is
% done only for plotting purposes.

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
subnum_dir = '/work/levan_lab/mtang/fmri_project/';

% enter list of subjects interested (subjects and runs that had significant
% variations in IED rate during scan)
sub_list = {'18', '27', '33', '34'};

% enter window size in seconds (must be divisible by tr, i.e. 1.5 s)
window_size = 120;

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% % (A): identify list of subject numbers with path input
% 
% % use dir to get a struct for all folders/files found in subject direct
% folder_list = dir(subnum_dir);
% 
% % access name of every folder/file identified, then search for pattern of
% % 'sub\d*' in it, if matched, extract the subject number
% sub_list = [];   % initialize array for storing subject numbers
% for item = 1:numel(folder_list)   % every item in folders/files found
%     
%     % search filename with matching pattern 'sub\d*', then use strrep to
%     % replace 'sub' by '' to obtain subject number
%     matched_str = regexp(folder_list(item).name, '^sub\d*$', 'match');
%     matched_sub_num = strrep(matched_str, 'sub', '');
%     
%     % if pattern is found, add subject number to sub_list
%     if ~isempty(matched_str)
%         sub_list = [sub_list matched_sub_num];
%     end
% end

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
filename_op = 'dynamic_analysis_falff_map.mat';   % filename of output file

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
    
    % Part (II): get amplitude of low-frequency fluctuation (falff) map (every voxel)
    % in every segment 

    % repetition time, in seconds
    tr = swra_info.PixelDimensions(4); 

    % use get_falff to get map of normalized falff in every segment 
    norm_falff = get_falff(motion_filtered_img, expmask, tr, window_size);

    % if empty normalized falff map is obtained, time series of current session
    % is not long enough for making psd estimates using pwelch method, 
    % so, pass control to next session, output empty array for current
    % session
    if isempty(norm_falff)
        continue
    end

    % END Part (II): get amplitude of low-frequency fluctuation (falff) in
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

    % Part (V): construct ROI of 3X3X3 box centered at midpoint of each electrode
    % pair to compute falff enclosed for every channel in midpt_all array,
    % regardless if the channel has during-scan IED registered
 
    % get basic info. about the processed image
    tr = swra_info.PixelDimensions(4);   % repetition time, in seconds
    fs = 1/tr; % sampling frequency, in Hz

    % % initialize array, sig_box, for storing the following for every type of elect.
    % % 1st column lists the names of channels,
    % % 2nd column lists mean normalized falff in each time segment of every channel
    sig_box_all = {};

    % use for loop to get signals within each box from motion-filtered image
    for i = 1:size(midpt_all, 1)
        sig_box_all{i, 1} = midpt_all{i, 1};   % get name tag of electrode pair
        mxyz = midpt_all{i, end};   % get midpoint in image space

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
        % 3X3X3 box, (i) sum all normalized falff values thereof, 
        % (ii) then get mean normalized falff in that box by dividing
        % sum by number of non-zero entries of that box in explicit mask
        for seg_ind = 1:num_seg   % for every segment in time
            seg_sum_norm_falff = 0;   % initialize sum of norm. falff as 0 
            for x_vox = mxyz(1)-1:mxyz(1)+1   % scan across every voxel in curr. box
                for y_vox = mxyz(2)-1:mxyz(2)+1
                    for z_vox = mxyz(3)-1:mxyz(3)+1
                        % if current voxel is included in explicit mask
                        if expmask(x_vox, y_vox, z_vox) == 1                 
                            % compute sum of all norm. falff values in curr. box
                            seg_sum_norm_falff = seg_sum_norm_falff + ...
                                norm_falff{x_vox, y_vox, z_vox}{seg_ind};
                        end
                    end
                end
            end

            % compute mean normalized falff in curr. box by dividing the sum
            % by number of non-zero entries of curr. box. in explicit mask,
            % assign calculated mean to curr. cell in 3rd col. of ch_mnfalff
            sig_box_all{i, 2}{seg_ind} = seg_sum_norm_falff/num_nonzero;
        end   % end for each segment in time

    end   % for i = 1:size(midot_all_sel, 1), for each channel

    % END Part (V): construct ROI of 3X3X3 box centered at midpoint of each electrode
    % pair to compute falff enclosed for every channel in midpt_all array,
    % regardless if the channel has during-scan IED registered
    %---------------------------------------------------------------------------

    % Part (VI): store quantities calcu. in output struct.

    opstruct = struct;
    opstruct.ele_sorted = ele_sorted;   % info. of electrode (sorted by name)
    opstruct.midpt_all = midpt_all;   % coordinates of midpoint between each pair of electrode contacts in single layer
    opstruct.sig_box_all = sig_box_all;   % signals within 3X3X3 box centered at each midpoint of all channels in single laye

    % END Part (VI): store quantities calcu. in output struct.
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

end   % for run_ind = 1:numel(swra_file_path)   % for every session

%---------------------------------------------------------------------------

% format fieldname for current subject number
g_fdname = ['sub', subnum];   % format fieldname for current subject number
grand.(g_fdname) = terms;   % update grand, struct, with info of curr. sub.

% output structure created in current session
if op_results == 1
%     save(fullfile(fname_op, filename_op), 'grand');
    save(fullfile(fname_op, filename_op), 'grand', '-v7.3'); 
end

sprintf(['done working on sub', subnum])

end   % end for item = 1:numel(sub_list)

toc