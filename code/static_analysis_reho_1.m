clear all 
close all

% this script does the following static analysis:
% (i) get ReHo over entire scan for each channel in every session
% (ii) get IED rate measured during long-term monitoring 
% (iii) get IED rate measured during scan
% (iv) get Spearman's corr. coeff. btw. (i) and (ii) for every session of
% each subject
% (v) get Spearman's corr. coeff. btw. (i) and (iii) for every session of
% each subject
% (vi) separate (iv) and (v) into categries defined by their channel type
% (e.g. promin
% Note: not all channels have (iii) registerd

% quantities calculated are saved in a struct with format given below, 
% grand --- subject number --- session number --- quantities
% e.g. to access reho calculated at session 7 for subject 14, 
% grand.sub14.run7.reho 

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
filename_static_spikes = ['subject', subnum, '_rates.txt'];   % file of static spike rates (long-term monitoring)

% enter path where ouput struct. is stored at
fname_op = [subnum_dir, filesep, 'matrices' filesep 'static_analysis_reho'];   % direct. of output matrix
filename_op = 'static_analysis_reho_1.mat';   % filename of output file

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
    % Part (II): sort electrodes by their types

    % use sort_elect to sort electrode array by their types
    ele_sorted = sort_elect(electar);

    % END PART (II): sorte electrodes by their types
    %---------------------------------------------------------------------------
    % % Part (III): find midpoint of each pair of electrode contacts in both mni
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

    % END Part (III): find midpoint of each pair of electrode contacts in both mni
    % and image space
    %---------------------------------------------------------------------------
    % % Part (IV): construct 3X3X3 boxes centered at midpoint of each electrode
    % % pair to store signals enclosed thereof
 
    % get basic info. about the processed image
    tr = swra_info.PixelDimensions(4);   % repetition time, in seconds
    fs = 1/tr; % sampling frequency, in Hz

    % % initialize array, sig_box, for storing the following for every type of elect.
    % % 1st column lists the names of channels,
    % % 2nd column lists reho (aka kendall_w coeff.) of every channel
    sig_box = {};

    % use for loop to get signals within each box from motion-filtered image
    for i = 1:numel(midpt)   % for each electrode type
        for j = 1:size(midpt{i}, 1)   % for each electrode pair under current type
            sig_box{i}{j, 1} = midpt{i}{j, 1};   % get name tag of electrode pair
            mxyz = midpt{i}{j, 3};   % get midpoint in image space
        
            % get signals from motion-filtered images
            signal =  double(motion_filtered_img(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
            mxyz(3)-1:mxyz(3)+1, :));  
        
            % create regional mask from explicit mask using the coordinates of
            % the 3X3X3 box
            reg_mask = double(expmask(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
                mxyz(3)-1:mxyz(3)+1));
        
            % use Chebyshev Type 1 bandpass filter to allow freq. compo. btw. 
            % [0.01 and 0.1] Hz to pass through, use cheby1 to get 
            % transfer function coefficients (normalized freqs. are used)
            [b,a] = cheby1(2,0.5,[0.01 0.1]/(fs/2));
         
            % apply bandpass filter (Chebyshev filter) to remove low-freq. drift 
            % and noises from signal
            signal_noise_filt = signal;   % initialize array
            for x_vox = 1:size(signal, 1)   % for each voxel in x-dir.
                for y_vox = 1:size(signal, 2)   % for each voxel in y-dir.
                    for z_vox = 1:size(signal, 3)   % for each voxel in z-dir.
                        % extract signal at voxel interested, squeeze the
                        % array to remove dim. of length 1, then use Chebyshev
                        % filter to remove noises
                        signal_noise_filt(x_vox, y_vox, z_vox, :) = ...
                            filtfilt(b, a, squeeze(signal(x_vox, y_vox, z_vox, :))); 
                    end
                end
            end
       
            % get Kendall's W coefficient of filtered signal, assign
            % quantities calcu. to 3rd column of current cell within sig_box
            sig_box{i}{j, 2} = kendall_w(signal_noise_filt, reg_mask);       
        end
    end

    % END Part (IV): construct 3X3X3 boxes centered at midpoint of each electrode
    % pair to store signals enclosed thereof
    %---------------------------------------------------------------------------
    % Part (V): group signals of all channels in single layer

    % group all electrode types in one single layer of cell array
    % col. 1 = id assigned to electrode pair
    % col. 2 = name tag of electrode pair
    % col. 3 = regional homogeneity (reho) of every channel
    sig_box_all = {};   % initialize sig_box_all array
    rownum = 1;   % initialize row number as 1 

    for i = 1:numel(sig_box)   % for each electrode type in sig_box
        for j = 1:size(sig_box{i}, 1)   % for each pair of electrode contacts
            sig_box_all(rownum, 1) = {rownum};   % assign row number as id 
            sig_box_all(rownum, 2:2+size(sig_box{i}, 2)-1) = sig_box{i}(j, :);   % assign relevant info.
            rownum = rownum + 1;   % increment row number by 1
        end
    end

    % END Part (V): group signals of all channels in single layer
    %---------------------------------------------------------------------------
    % Part (VI): match channels in spike rates array (IED rates measured 
    % during long-term monitoring) with that of sig_box_all

    % use clean_match_spikes to get spikes struct.
    static_spikes = clean_match_spikes(spikesar, sig_box, sig_box_all);

    % get spike rates array in single layer with channels matched with sig_box_all
    static_spikes_all = static_spikes.spikes_all;
    
    % get ied rate measured during long-term monitoring for every channel
    % in single layer (3rd. col. from static_spikes_all array),
    ied_long_term = static_spikes_all(:, 1:3);

    % END Part (VI): match channels in spike rates array with that of
    % signal box
    %---------------------------------------------------------------------------
    % Part (VII): get IED rates (per minute) measured during the entirety of scan 

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
    
    % END Part (VII): get IED rates (per minute) measured during the entirety of scan 
    %---------------------------------------------------------------------------
    % Part (VIII): get names and event types of ied onset and ied
    % propagation, and no-ied channels

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

    % get names of no ied channels (those with ieds registered during fmri
    % scan)
    [~, ia, ~] = intersect(sig_box_all(:, 2), ch_list_ext(:, end));   % indices required
    ch_noied_name = sig_box_all(setdiff(1:end, ia), 2);   % select names from sig_box_all

    % assemble ch_list_noied array in the same format as those listed above
    % using indices and names acquired, 'NA' is used for event type, as no
    % during scan ied is recorded
    ch_list_noied = [repmat({subnum}, [size(ch_noied_name, 1) 1]) ...
        repmat({'NA'}, [size(ch_noied_name, 1) 1]) ch_noied_name];

    % obtain array listing names of all channels available 
    % ch_list_all = ch_list_ext + ch_list_noied, any channels regardless if
    % ieds were registered during fmri scan, note there could be duplicate
    % entries in ch_list_ext, as propa channels in one ied type could be
    % ied onset channels in another ied type
    ch_all_name = [unique(ch_list_ext(:, end)); ch_noied_name];
    ch_list_all = [repmat({subnum}, [size(ch_all_name, 1) 1]) ...
        repmat({'NA'}, [size(ch_all_name, 1) 1]) ch_all_name]; 

    % note: ch_list_onset/propa lists channel names derived from eeg data
    % set, while sig_box_all lists channel names derived from fmri data
    % set,
    % ch_list_all may contain a different number of channels than that in
    % sig_box_all, as onset and propa channels are determined using eeg data 
    % set, which could contain more channels than the fmri data set. 
    % Nonetheless, ied rate and time series features are assigned to 
    % channel name in ch_list_onset/propa/noied/all only if there is a match btw. 
    % the array and sig_box_all, see below for details

    % END Part (VIII): get names and coordinates of channels with ied rates
    % registered during fmri scan
    % %---------------------------------------------------------------------------

    % Part (IX): assign ied rate (per minute) measured over entire scan 
    % and reho to every channel in each array (ied onset, propagation, no ied, all)

    % use match_lt_ied_reho, defined at the end of script, to assigan long term
    % monitoring ied rate (per minute) and reho measured over entire scan to 
    % every channel listed in each of the arrays obtained above (ied onset, 
    % ied propagation, no ied, all = onset + propa + no ied)
    ch_onset_lt_ied_reho = match_lt_ied_reho(ch_list_onset, ied_long_term, sig_box_all);
    ch_propa_lt_ied_reho = match_lt_ied_reho(ch_list_propa, ied_long_term, sig_box_all);
    ch_noied_lt_ied_reho = match_lt_ied_reho(ch_list_noied, ied_long_term, sig_box_all);
    ch_all_lt_ied_reho = match_lt_ied_reho(ch_list_all, ied_long_term, sig_box_all);

    % use match_ds_ied_reho, defined at the end of script, to assign during 
    % scan ied rate (per minute) and reho measured over entire scan to every 
    % channel listed in each of the arrays obtained above (ied onset, ied propagation, 
    % and ext = onset + propa, as only channels with ieds registered during
    % scan is considered in this case)
    ch_onset_ds_ied_reho = match_ds_ied_reho(ch_list_onset, ied_during_scan, sig_box_all);
    ch_propa_ds_ied_reho = match_ds_ied_reho(ch_list_propa, ied_during_scan, sig_box_all);
    ch_ext_ds_ied_reho = match_ds_ied_reho(ch_list_ext, ied_during_scan, sig_box_all);

    % END Part (IX): assign ied rate (per minute) measured over entire scan 
    % and reho to every channel in each array (ied onset, propagation, no ied, all)
    %---------------------------------------------------------------------------

    % Part (X): get spearman's corr. coeff. btw. (i) ReHo and (ii) IED
    % rate measured during long-term monitoring, and btw. (i) ReHo and
    % (iii) IED rate measured during scan (entire)
    
    % corr. btw. reho and ied measured during long-term monitoring and its 
    % p-value are classified into the following categories based on channel
    % type, onset, propa, no ied, and all,
    [rho_reho_ied_lt_onset, pval_reho_ied_lt_onset] = compute_corr(ch_onset_lt_ied_reho);
    [rho_reho_ied_lt_propa, pval_reho_ied_lt_propa] = compute_corr(ch_propa_lt_ied_reho);
    [rho_reho_ied_lt_noied, pval_reho_ied_lt_noied] = compute_corr(ch_noied_lt_ied_reho);
    [rho_reho_ied_lt_all, pval_reho_ied_lt_all] = compute_corr(ch_all_lt_ied_reho);

    % corr. btw. reho and ied measured during scan and its 
    % p-value are classified into the following categories based on channel
    % type, onset, propa, ext,
    % initialize arrays for each category
    [rho_reho_ied_ds_onset, pval_reho_ied_ds_onset] = compute_corr(ch_onset_ds_ied_reho);
    [rho_reho_ied_ds_propa, pval_reho_ied_ds_propa] = compute_corr(ch_propa_ds_ied_reho);
    [rho_reho_ied_ds_ext, pval_reho_ied_ds_ext] = compute_corr(ch_ext_ds_ied_reho);

    % END Part (X): get spearman's corr. coeff. btw. (i) ReHo and (ii) IED
    % rate measured during long-term monitoring, and btw. (i) ReHo and
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
    
    % list of channels derived from eeg data set based on their category
    opstruct.ch_list_onset = ch_list_onset;   % ied onset channels
    opstruct.ch_list_propa = ch_list_propa;   % ied propagation chanenls
    opstruct.ch_list_ext = ch_list_ext;   % ext = ied onset + propa
    opstruct.ch_list_noied = ch_list_noied;   % no ied channels (during fmri scan)
    opstruct.ch_list_all = ch_list_all;   % all channels (onset + propa + no ied)

    % arrays listing channels with reho and long term monitoring ied rate per min
    % assigned, divided into categories 
    opstruct.ch_onset_lt_ied_reho = ch_onset_lt_ied_reho;   % onset 
    opstruct.ch_propa_lt_ied_reho = ch_propa_lt_ied_reho;   % propa 
    opstruct.ch_noied_lt_ied_reho = ch_noied_lt_ied_reho;   % no ied 
    opstruct.ch_all_lt_ied_reho = ch_all_lt_ied_reho;   % all = onset+propa+no ied

    % arrays listing channels with reho and during scan ied rate per min
    % assigned, divided into categories 
    opstruct.ch_onset_ds_ied_reho = ch_onset_ds_ied_reho;   % onset 
    opstruct.ch_propa_ds_ied_reho = ch_propa_ds_ied_reho;   % propa 
    opstruct.ch_ext_ds_ied_reho = ch_ext_ds_ied_reho;   % no ied 
    
    % so, we have for every channel, (i) ReHo, (ii) IED rate measured
    % during long-term monitoring, and (iii) IED rate (per minute) measured 
    % during entire scan, then, we computed spearman's rho btw. (i) and (ii), 
    % and (i) and (iii), based on channel category

    % for spearman's rho btw. (i) and (ii) and its p-value, we have 
    opstruct.rho_reho_ied_lt_onset = rho_reho_ied_lt_onset;   % onset
    opstruct.pval_reho_ied_lt_onset = pval_reho_ied_lt_onset;

    opstruct.rho_reho_ied_lt_propa = rho_reho_ied_lt_propa;   % propa
    opstruct.pval_reho_ied_lt_propa = pval_reho_ied_lt_propa;

    opstruct.rho_reho_ied_lt_noied = rho_reho_ied_lt_noied;   % no ied
    opstruct.pval_reho_ied_lt_noied = pval_reho_ied_lt_noied;

    opstruct.rho_reho_ied_lt_all = rho_reho_ied_lt_all;   % all = onset+propa+no ied
    opstruct.pval_reho_ied_lt_all = pval_reho_ied_lt_all;
    
    % for spearman's rho btw. (i) and (iii) and its p-value, we have 
    opstruct.rho_reho_ied_ds_onset = rho_reho_ied_ds_onset;   % onset
    opstruct.pval_reho_ied_ds_onset = pval_reho_ied_ds_onset;

    opstruct.rho_reho_ied_ds_propa = rho_reho_ied_ds_propa;   % propa
    opstruct.pval_reho_ied_ds_propa = pval_reho_ied_ds_propa;

    opstruct.rho_reho_ied_ds_ext = rho_reho_ied_ds_ext;   % ext = onset+propa
    opstruct.pval_reho_ied_ds_ext = pval_reho_ied_ds_ext;


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
    
% USER DEFINED FUNCTIONS FOR THIS SCRIPT:
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% FUNCTION OUTPUT_AR = MATCH_LT_IED_REHO(INPUT_AR, IED_LT, SIG_BOX_ALL)
% this function requires the following inputs,
% INPUT_AR = array containing channels to which during scan ied and reho
% are assigned, with format given below
% 1st col. = subject number
% 2nd col. = event type
% 3rd col. = channel name
% IED_LT = array containing long term monitoring ied rate per minute, with format
% given below, 
% 1st col. = channel id
% 2nd col. = channel name
% 3rd col. = long term monitoring ied rate per minute
% SIG_BOX_ALL_AR = array containing reho of all channels available, with
% format given eblow, 
% 1st col. = channel id
% 2nd col. = channel name
% 3rd col. = feautre interested (e.g. ReHo, ALFF, etc.)
% this function generates the following output, 
% OUTPUT_AR = output array with during scan ied rate and reho assigned to
% each channel listed in INPUT_AR, with format given below
% 1st col. = subject number
% 2nd col. = ied type 
% 3rd col. = channel name
% 4th col. = long term monitoring ied rate per minute
% 5th col. = time series feature (e.g. ReHo, ALFF)

function output_ar = match_lt_ied_reho(input_ar, ied_lt, sig_box_all_ar)

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
% 4th col. = long term monitoring ied rate
ch_lt_ied = {};
row_num = 1;    % initialize row number 
ncol = size(input_ar, 2);   % get number of col. established in input_ar

% assign ied measured over long term monitoring to each channel
for i = 1:size(input_ar, 1)   % for every row in input array
    ch_name_cur = input_ar{i, 3};   % get current channel name from col. 3
    for j = 1:size(ied_lt, 1)   % for every row in ied_lt
        if isequal(ch_name_cur, ied_lt{j, 2})   % find row with matching channel name
            % assign all info related to the channel to current row
            ch_lt_ied(row_num, 1:ncol) = input_ar(i, :);
            % assign long-term monitoring ied rate (per minute) to rightmost
            % col. of current row
            ch_lt_ied{row_num, ncol+1} = ied_lt{j, end};
            row_num = row_num + 1;   % increment row number by 1 for next matched channel
        end
    end
end

% get number of col. established in ch_lt_ied array 
ncol = size(ch_lt_ied, 2);

% initialize ch_lt_ied_reho, with the same format as
% ch_lt_ied, but with an additional column storing reho
% computed over the entire scan
ch_lt_ied_reho = {};
row_num = 1;   % initialize row number as 1

% assign reho measured over entire scan to each channel that has
% matched ied type found, and matched channel name in sig_box_all_ar
for i = 1:size(ch_lt_ied, 1)   % for every row in ch_lt_ied
    ch_name_cur = ch_lt_ied{i, 3};   % get current channel name from col. 3
    for j = 1:size(sig_box_all_ar, 1)   % for every row in sig_box_all_ar
        if isequal(ch_name_cur, sig_box_all_ar{j, 2})   % find row with matching channel name
            % assign reho of matched channel to rightmost col. of current row
            ch_lt_ied_reho(row_num, 1:ncol) = ch_lt_ied(i, :);
            ch_lt_ied_reho{row_num, ncol+1} = sig_box_all_ar{j, end};
            row_num = row_num + 1;   % increment row number by 1 for next matched channel
        end
    end
end

% assign variable computed to output
output_ar = ch_lt_ied_reho;

end   % end function output_ar = match_lt_ied_reho(input_ar, ied_lt, sig_box_all_ar)

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% FUNCTION OUTPUT_AR = MATCH_DS_IED_REHO(INPUT_AR, IED_DS, SIG_BOX_ALL)
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
% OUTPUT_AR = output array with during scan ied rate and reho assigned to
% each channel listed in INPUT_AR, with format given below
% 1st col. = subject number
% 2nd col. = ied type 
% 3rd col. = channel name
% 4th col. = during scan ied rate per minute
% 5th col. = time series feature (e.g. ReHo, ALFF)

function output_ar = match_ds_ied_reho(input_ar, ied_ds, sig_box_all_ar)

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

% initialize ch_ds_ied_reho, with the same format as
% ch_ds_ied, but with an additional column storing reho
% computed over the entire scan
ch_ds_ied_reho = {};
row_num = 1;   % initialize row number as 1

% assign reho measured over entire scan to each channel that has
% matched ied type found, and matched channel name in sig_box_all_ar
for i = 1:size(ch_ds_ied, 1)   % for every row in ch_ds_ied
    ch_name_cur = ch_ds_ied{i, 3};   % get current channel name from col. 3
    for j = 1:size(sig_box_all_ar, 1)   % for every row in sig_box_all_ar
        if isequal(ch_name_cur, sig_box_all_ar{j, 2})   % find row with matching channel name
            % assign reho of matched channel to rightmost col. of current row
            ch_ds_ied_reho(row_num, 1:ncol) = ch_ds_ied(i, :);
            ch_ds_ied_reho{row_num, ncol+1} = sig_box_all_ar{j, end};
            row_num = row_num + 1;   % increment row number by 1 for next matched channel
        end
    end
end

% assign variable computed to output
output_ar = ch_ds_ied_reho;

end   % end function output_ar = match_ds_ied_reho(input_ar, ied_ds, sig_box_all_ar)

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% FUNCTION [OUTPUT_RHO, OUTPUT_PVAL] = COMPUTE_CORR(INPUT_AR)
% this function requires an input, 
% INPUT_AR = input array with format given below
% end-1 col. = ied rate (A)
% end col. = feature interested (ReHo, ALFF) (B)
% and generates the following outputs,
% output_rho = spearman's corr. coeff. btw. (A) and (B), and
% output_pval = p-value of spearman's corr. coeff.
function [output_rho, output_pval] = compute_corr(input_ar)

% initialize output variables
output_rho = [];   % spearman's rho
output_pval = [];   % p-value

% if input ar is empty, return control to line
% calling this function
if isempty(input_ar) 
    return
end

% if size of input array is less than 9 and rho is exactly -1, p-value 
% cannot be computed by matlab, a bug that needs to be fixed, manually set
% p-value as 0 and return both rho_val and p-value

% get rho, then get diff, and check if diff is <= tol 
rho_val = corr([input_ar{:, end-1}]', [input_ar{:, end}]', 'Type', 'Spearman');

tol = 1e-5;   % set tolerance
diff = rho_val - (-1);   % compute diff btw. rho_val and -1
if size(input_ar, 1) <= 9 && diff <= tol
    output_rho = rho_val;   % assign rho computed
    output_pval = 0;   % manually set p-value as 0
    return
end

% compute spearman's corr. coeff. btw. the rightmost two columns of input
% array
[output_rho, output_pval] = corr([input_ar{:, end-1}]', ...
        [input_ar{:, end}]', 'Type', 'Spearman');

end

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------