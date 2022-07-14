clear all 
close all

% this script does the following static analysis:
% (i) get ALFF over entire scan for each channel in every session
% (ii) get IED rate measured during long-term monitoring 
% (iii) get IED rate measured during scan
% (iv) get Spearman's corr. coeff. btw. (i) and (ii) for every session of
% each subject
% (v) get Spearman's corr. coeff. btw. (i) and (iii) for every session of
% each subject
% Note: not all channels have (iii) registerd

% quantities calculated are saved in a struct with format given below, 
% grand --- subject number --- session number --- quantities
% e.g. to access alff calculated at session 7 for subject 14, 
% grand.sub14.run7.alff

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

tic 

% enter path to directory where info. of all subjects are stored
subnum_dir = '/work/levan_lab/mtang/fmri_project/' ;

% enter window size in seconds (must be divisible by tr, i.e. 1.5 s)
% this is required to obtain instantaneous spike rate, thoguh ALFF is
% calculated for the entire time series (not used)
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

% (B): compute (i) ALFF and (ii) IED rate measured during long-term monitoring,
% and (iii) IED rate measured over entire scan, and get spearman's rho btw.
% (i) and (ii), and (i) and (iii), for every session of each subject

% initialize structure, grand, to store relevant info. of every subject
grand = struct;

% loop through each subject number identified and conduct analysis
for item = 1:numel(sub_list)
    
% get current subject number
subnum = sub_list{item};

sprintf(['begin working on sub', subnum])

%---------------------------------------------------------------------------
% (B1): format paths of required files for current subject

% enter path to directory where all input files are located
directname = ['/work/levan_lab/mtang/fmri_project/', 'sub', subnum];
% directname = ['C:\Users\siumichael.tang\Downloads\fmri_project\', 'sub', subnum];

% format filenames of processed fmri images, explicit mask, electrode,
% motion parameters
filename_swraimg = ['swra*.nii'];   % processed func. images
filename_expmask = 'wEPI_bet_mask.nii';   % explicit mask
filename_elect =  [subnum, '_*Koordinaten*.xlsx'];   % file containing mni coord. of all electrode pairs
filename_event_onset = ['onset_times_events_sub', subnum, '.mat'];   % file of event onsets   
filename_motion = ['rp_*.txt'];   % motion parameters generated by SPM in the realignment step
filename_static_spikes = ['subject', subnum, '_rates.txt'];   % file of static spike rates (long-term monitoring)

% enter path where ouput struct. is stored at
fname_op = [subnum_dir, filesep, 'matrices' filesep 'static_analysis_alff'];   % direct. of output matrix
filename_op = 'static_analysis_alff.mat';   % filename of output file

% END (B1): format paths of required files for current subject
%---------------------------------------------------------------------------
% (B2): obtain full paths of required files for current subject

% use get_path function to obtain full path of input files
[~, swra_file_path, ~] = get_path(directname, filename_swraimg);   % swra func. images
[~, expmask_file_path, ~] = get_path(directname, filename_expmask);   % normalized mask images
[~, elect_file_path, ~] = get_path(directname, filename_elect);   % file containing info. about electrodes coordinates
[~, event_onsets_file_path, ~] = get_path([directname, filesep, 'matrices'], filename_event_onset);   % file of event onsets
[~, motion_file_path, ~] = get_path(directname, filename_motion);   % motion parameters
[~, static_spikes_file_path] = get_path(directname, filename_static_spikes);   % static spike rates (long-term monitoring)

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
    spikesar = readcell(static_spikes_file_path{:});   % cell array containing info. about static spike rates
    
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
    
    % Part (II): get amplitude of low-frequency fluctuation (ALFF) map (every voxel)

    % repetition time, in seconds
    tr = swra_info.PixelDimensions(4); 

    % use get_alff to get map of normalized alff
    norm_alff = get_alff(motion_filtered_img, expmask, tr, []);

    % if empty normalized alff map is obtained, time series of current session
    % is not long enough for making psd estimates using pwelch method, 
    % so, pass control to next session, output empty array for current
    % session
    if isempty(norm_alff)
        continue
    end

    % END Part (II): get amplitude of low-frequency fluctuation (ALFF) in each 
    % segment
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
    % % Part (V): construct 3X3X3 boxes centered at midpoint of each electrode
    % % pair to store signals enclosed thereof
 
    % get basic info. about the processed image
    tr = swra_info.PixelDimensions(4);   % repetition time, in seconds
    fs = 1/tr; % sampling frequency, in Hz

    % % initialize array, sig_box, for storing the following for every type of elect.
    % % 1st column lists the names of channels,
    % % 2nd column lists mean normalized alff of every channel
    sig_box = {};

    % use for loop to get signals within each box from motion-filtered image
    for i = 1:numel(midpt)   % for each electrode type
        for j = 1:size(midpt{i}, 1)   % for each electrode pair under current type
            
            sig_box{i}{j, 1} = midpt{i}{j, 1};   % get name tag of electrode pair
            mxyz = midpt{i}{j, 3};   % get midpoint in image space

            % create regional mask from explicit mask using the coordinates of
            % the 3X3X3 box
            reg_mask = double(expmask(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
                mxyz(3)-1:mxyz(3)+1));

            % get mean normalized alff in current 3X3X3 box, assign result to
            % 3rd column of current cell within sig_box
            if ~isempty(norm_alff)   % if norm_alff array is not empty
                
                % get number of non-zero entries in current regional mask
                num_nonzero = nnz(reg_mask);
            
                % get mean normalized alff in current 3X3X3 box
                sig_box{i}{j, 2} = sum(sum(sum(norm_alff(mxyz(1)-1:mxyz(1)+1, ...
                    mxyz(2)-1:mxyz(2)+1, mxyz(3)-1:mxyz(3)+1))))/num_nonzero;
            else
                sig_box{i}{j, 2} = [];   % assign empty value if norm_alff is empty
            end
        
        end
    end

    % END Part (V): construct 3X3X3 boxes centered at midpoint of each electrode
    % pair to store signals enclosed thereof
    %---------------------------------------------------------------------------
    % Part (VI): group signals of all channels in single layer

    % group all electrode types in one single layer of cell array
    % col. 1 = id assigned to electrode pair
    % col. 2 = name tag of electrode pair
    % col. 3 = mean normalized alff of every channel
    sig_box_all = {};   % initialize sig_box_all array
    rownum = 1;   % initialize row number as 1 

    for i = 1:numel(sig_box)   % for each electrode type in sig_box
        for j = 1:size(sig_box{i}, 1)   % for each pair of electrode contacts
            sig_box_all(rownum, 1) = {rownum};   % assign row number as id 
            sig_box_all(rownum, 2:2+size(sig_box{i}, 2)-1) = sig_box{i}(j, :);   % assign relevant info.
            rownum = rownum + 1;   % increment row number by 1
        end
    end

    % END Part (VI): group signals of all channels in single layer
    %---------------------------------------------------------------------------
    % Part (VII): match channels in spike rates array (IED rates measured 
    % during long-term monitoring) with that of midpt

    % use clean_match_spikes to get spikes struct.
    static_spikes = clean_match_spikes(spikesar, sig_box, sig_box_all);

    % get spike rates array in single layer with channels matched with sig_box_all
    static_spikes_all = static_spikes.spikes_all;
    
    % get ied rate measured during long-term monitoring for every channel
    % in single layer (3rd. col. from static_spikes_all array)
    ied_long_term = static_spikes_all(:, 1:3);

    % END Part (VII): match channels in spike rates array with that of signal box
    %---------------------------------------------------------------------------
    % Part (VIII): get IED rates (per minute) measured during the entirety of scan 

    % get length of time series and repetition time (tr) in seconds
    nt = size(swra_img, 4);   % length in time
    tr = swra_info.PixelDimensions(4);   % tr, in seconds
    
    % get number of event types recorded in current session
    num_events = numel(event_onsets_curr_sess.results.evti);
    
    % initialize array, ied_during_scan, with format given below
    % col. 1 = event type
    % col. 2 = number of ieds [per minute] over the entire scan 
    ied_during_scan = cell(num_events, 2);
    
    for ev_ind = 1:num_events   % for each event type
        ied_during_scan{ev_ind, 1} = event_onsets_curr_sess.results.event_def{ev_ind};   % get event type def.
        ied_during_scan{ev_ind, 2} = (numel(event_onsets_curr_sess.results.evti{ev_ind}))/(nt*tr)*60;   % get number of ieds per minute
    end
    
    % END Part (VIII): get IED rates (per minute) measured during the entirety of scan 
    %---------------------------------------------------------------------------
    % Part (IX): get names and coordinates of channels with ied rates
    % registered during fmri scan

    % get names of channels associated with event types assigned to current 
    % subject number, and match the names with those in midpt_all
    ch_name_matched = get_channel_assoc_event(subnum, midpt_all);

    % END Part (IX): get names and coordinates of channels with ied rates
    % registered during fmri scan
    % %---------------------------------------------------------------------------
    % Part (IX): assign ied rate (per minute) measured over entire scan 
    % and mean normalized alff to each channel that has ied rate registered during scan
    
    % initialize array, with same format as ch_name_matched, but extra col.
    % added to store info. about ied rate measured over entire scan and
    % reho
    % col. 1 = event type
    % col. 2 = channel name
    % col. 3 = midpoint coordinates in mni space
    % col. 4 = midpoint coordinates in image space
    % col. 5 = ied rate measured over entire scan
    % col. 6 = mean normalized alff
    ch_name_matched_ied_alff = ch_name_matched;
    
    % get number of col. established
    ncol = size(ch_name_matched_ied_alff, 2);
  
    % assign ied rate measured over entire scan to each channel available
    for row = 1:size(ch_name_matched_ied_alff, 1)   % for each row in ch_name_matched_ied_reho
        ev_type = ch_name_matched_ied_alff{row, 1};   % get event type of curr. channel
        for i = 1:size(ied_during_scan, 1)   % for every row in ied_during_scan array
            % get event type from col. 1, and check is there is match
            if ismember(ev_type, ied_during_scan{i, 1})
                % if so, assign ied rate (per min.) measured over entire
                % scan to curr. channel
                ch_name_matched_ied_alff{row, ncol+1} = ied_during_scan{i, end};
            end
        end
    end
    
    % assign mean normalized alff measured over entire scan to each channel available
    for row = 1:size(ch_name_matched_ied_alff, 1)   % for every row in ch_name_matched_ied_alff
        ch_name_cur = ch_name_matched_ied_alff{row, 2};   % get current channel name
        for i = 1:size(sig_box_all, 1)   % for every row in sig_box_all
            if isequal(ch_name_cur, sig_box_all{i, 2})   % find row with matching channel name
                % assign mean norm. alff of matched channel to rightmost col. of current row
                ch_name_matched_ied_alff{row, ncol+2} = sig_box_all{i, end};
            end
        end
    end
    
    % END Part (IX): assign ied rate (per minute) measured over entire scan 
    % and mean normalized alff to each channel that has ied rate registered during scan
    %---------------------------------------------------------------------------
    % Part (X): get spearman's corr. coeff. btw. (i) ALFF and (ii) IED
    % rate measured during long-term monitoring, and btw. (i) ALFF and
    % (iii) IED rate measured during scan (entire)
    
    % initialize arrays of interest
    rho_alff_ied_lt = [];   % rho btw. alff and ied measured during long-term monitoring
    pval_alff_ied_lt = [];   % assoc. p-value
    
    rho_alff_ied_ds = [];   % rho btw. alff and ied measured during entire scan
    pval_alff_ied_ds = [];   % assoc. p-value

    % get spearman's rho btw. alff and ied rate_{long term} for each session
    [rho_alff_ied_lt, pval_alff_ied_lt] = corr([sig_box_all{:, end}]', [ied_long_term{:, 3}]', ...
    'Type', 'Spearman');

    % get spearman's rho btw. alff and ied rate_{during scan} for each session
    [rho_alff_ied_ds, pval_alff_ied_ds] = corr([ch_name_matched_ied_alff{:, end-1}]', ...
        [ch_name_matched_ied_alff{:, end}]', 'Type', 'Spearman');

    % END Part (X): get spearman's corr. coeff. btw. (i) ALFF and (ii) IED
    % rate measured during long-term monitoring, and btw. (i) ALFF and
    % (iii) IED rate measured during scan (entire)
    %---------------------------------------------------------------------------
    
    % Part (XI): store quantities calcu. in output struct.

    opstruct = struct;
    opstruct.ele_sorted = ele_sorted;   % info. of electrode (sorted by name)
    opstruct.midpt = midpt;   % coordinates of midpoint between each pair of electrode contacts 
    opstruct.sig_box = sig_box;   % signals within 3X3X3 box centered at each midpoint
    opstruct.sig_box_all = sig_box_all;   % signals within 3X3X3 box centered at each midpoint of all channels in single layer
    opstruct.static_spikes_all = static_spikes_all;   % list of spike rates (long term) of all electrode pairs in single layer
    
    opstruct.ied_long_term = ied_long_term;   % ied rate measured during long-term monitoring for every channel
    opstruct.ied_during_scan = ied_during_scan;   % ied rate measured during entire scan for every channel
    
    % only portions of channels have ied rates registered during scan,
    % ch_name_matched_ied_alff holds info about ied rate measured during scan 
    % and alff of those channels
    opstruct.ch_name_matched_ied_alff = ch_name_matched_ied_alff;
    
    % so, we have for every channel, (i) mean norm. ALFF, (ii) IED rate measured
    % during long-term monitoring, and (iii) IED rate measured during
    % entire scan, then, we computed spearman's rho btw. (i) and (ii), and
    % (i) and (iii)
    opstruct.rho_alff_ied_lt = rho_alff_ied_lt;   % spearman's rho btw. (i) and (ii)
    opstruct.pval_alff_ied_lt = pval_alff_ied_lt;   % assoc. p-value
    
    opstruct.rho_alff_ied_ds = rho_alff_ied_ds;   % spearman's rho btw. (i) and (iii)
    opstruct.pval_alff_ied_ds = pval_alff_ied_ds;   % assoc. p-value

    % END Part (XI): store quantities calcu. in output struct.
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
    
end   % end for every session of current subject
%---------------------------------------------------------------------------

% format fieldname for current subject number
g_fdname = ['sub', subnum];   % format fieldname for current subject number
grand.(g_fdname) = terms;   % update grand, struct, with info of curr. sub.

% output table and structure created in current session
if op_results == 1
    save(fullfile(fname_op, filename_op), 'grand');
end

sprintf(['done working on sub', subnum])

end   % end for item = 1:numel(sub_list)
%---------------------------------------------------------------------------
toc
    