clear all 
close all

tic 
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '14';

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

% enter path where ouput struct. is stored at
fname_op = [directname, filesep, 'matrices' filesep 'alff_spikes'];   % direct. of output matrix
filename_op = 'alff_spikes.mat';   % filename of output file

% enter window size in seconds (must be divisible by tr, i.e. 1.5 s)
window_size_list = [60, 90, 120, 150, 180];

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% locate paths of func. images and their details 

% use get_path function to obtain full path of input files
[~, swra_file_path, ~] = get_path(directname, filename_swraimg);   % swra func. images
[~, expmask_file_path, ~] = get_path(directname, filename_expmask);   % normalized mask images
[~, elect_file_path, ~] = get_path(directname, filename_elect);   % file containing info. about electrodes coordinates
[~, event_onsets_file_path, ~] = get_path([directname, filesep, 'matrices'], filename_event_onset);   % file of event onsets
[~, motion_file_path, ~] = get_path(directname, filename_motion);   % motion parameters

% createterms structure to store output data calcu. by function defined
% below
terms = struct;

% read input files and use get_rho_spikes function to get required data
for run_ind = 1:numel(swra_file_path)   % for each row in path of processed func. images
   
for window_size_ind = 1:numel(window_size_list)   % for each window size input
% run_ind = 1;
% read input files
swra_img = niftiread(swra_file_path{run_ind});   % swra func. images
swra_info = niftiinfo(swra_file_path{run_ind});   % info about this particular swra func. images
expmask = single(niftiread(expmask_file_path{:}));   % convert image data type to single
electar = readcell(elect_file_path{:});   % cell array containing info. about electrodes
event_onsets = load(event_onsets_file_path{:});   % struct. containing event onsets

% assemble input array of function, get_rho_spikes 
input_array = struct;   % initialize struct.
input_array.swra_img = swra_img;   % processed func. images
input_array.swra_info = swra_info;   % headers of processed func. images
input_array.expmask = expmask;   % explicit mask
input_array.electar = electar;   % array of electrode coordinates

field_runs = fieldnames(event_onsets.terms);
input_array.event_onsets = event_onsets.terms.(field_runs{run_ind});   % array of event onsets
input_array.motion_file_path = motion_file_path{run_ind};   % path of motion parameter file

input_array.subnum = subnum;   % subject number of interest
input_array.window_size = window_size_list(window_size_ind);   % current window size 

% get indices of current run number embedded in processed image filename
sind = regexp(lower(swra_file_path{run_ind}), 'run');   % start index
eind = regexp(lower(swra_file_path{run_ind}), '.nii');   % end index
sess_num = swra_file_path{run_ind}(sind+3:eind-1);   % get indices of run number

% create fieldname representing current session number
fieldname = sprintf('run%s_winsize_%s', sess_num, num2str(window_size_list(window_size_ind)));

% store output calcu. under current fieldname of terms
terms.(fieldname) = get_mnalff_num_spikes(input_array);

% execute lines below if current field is not empty
if ~isempty(terms.(fieldname))
    
% access tables of current session created by function
% get_mnalff_num_spikes
table_ch_rho_ns_mnalff = terms.(fieldname).table_ch_rho_ns_mnalff;   

% format filenames and full paths of tables
filename_table = ['sub', subnum, '_table_ns_mnalff_run', sess_num, ...
    '_winsize_', num2str(window_size_list(window_size_ind)), '.csv'];
table_path = fullfile(fname_op, filename_table);   % path of table of clin. det. channels

%---------------------------------------------------------------------------
% output table and structure created in current session
if op_results == 1
    writetable(table_ch_rho_ns_mnalff, table_path);   % table of clin. det. ch.
%     save(fullfile(fname_op, filename_op), 'terms');   % struct.
    save(fullfile(fname_op, filename_op), 'terms', '-v7.3');
end
%---------------------------------------------------------------------------

end   % end if ~isempty(terms.(fieldname))

end   % for window_size_ind = 1:numel(window_size)   % for each window size input

end   % end for run_ind = 1:numel(swra_file_path)

toc
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% define function get_reho_alff_spikes requiring an input array, which
% has the fields defined below,
% SWRA_IMG is the data of processed image file
% SWRA_INFO is the info (headers) associated with the processed image file
% EXPMASK is the explicit mask image
% ELECTAR is the array containing info. about all electrodes used,
% EVENT_ONSETS is the onset times of each event type during fmri scan
% MOTION_FILE_PATH is the path of motion parameter file
% and outputs,
% OPSTRUCT, a structure that contains info. about coordinates of
% electrdoes, spearman corr. coeff. btw. number of spikes and mean
% normalized alff for each channel selected,
% see end of function for details of variables included in this struct.

function opstruct = get_mnalff_num_spikes(input_array)

% %---------------------------------------------------------------------------
% % Prelim: access variables defined in input array
swra_img = input_array.swra_img;   % processed func. images
swra_info = input_array.swra_info;   % headers of processed func. images
expmask = input_array.expmask;   % explicit mask
electar = input_array.electar;   % array of electrode coordinates
event_onsets = input_array.event_onsets;   % struct. of event onsets
motion_file_path = input_array.motion_file_path;   % path of motion parameter file
subnum = input_array.subnum;   % subject number of interest
window_size = input_array.window_size;   % current window size in (s)

% initialize output struct. as empty array
opstruct = [];

% check if length of time series is long enough, 256 is used as it's chosen
% as the width of window in welch's method of power spectral density (psd) 
% estimate. If so, change value of com_alff to 1 so that the script will
% calculate ALFF, the relevant corr. coeff., and p-values.
com_alff = 0;
if size(swra_img, 4) >= 256
    com_alff = 1;
end

% end Prelim: access variables defined in input array

%---------------------------------------------------------------------------

% execute the lines below only if com_alff == 1 (long enough time series)
if com_alff ~= 1
    disp('time series is not long enough for pwelch method')
else
      
% Part (I): remove motion artifacts from processed image (entire brain)

% convert explicit mask and processed image to double, 
% apply element-wise multiplciation
exp_swra_img = double(expmask).*double(swra_img);
                
% call filter_motion.m to remove motion artifacts from image
motion_filtered_img = filter_motion(motion_file_path, exp_swra_img);

% END Part (I): remove motion artifacts from processed images
%---------------------------------------------------------------------------

% Part (II): get spike rate in each segment 

% get basic info. about the processed image
tr = swra_info.PixelDimensions(4);   % repetition time, in seconds
fs = 1/tr; % sampling frequency, in Hz

% get dimensions of motion-filtered images
ni = size(motion_filtered_img, 1);   % dim. in x-dir.
nj = size(motion_filtered_img, 2);   % dim. in y-dir.
nk = size(motion_filtered_img, 3);   % dim. in z-dir.
nt = size(motion_filtered_img, 4);   % dim. in t-dir.

% create time axis (slice-time corrected) spanning the entire 
% time series of siganl, 
% 1st image is shifted to t=0 (s), while the last image is
% shifted to t_end-tr (s)
taxis = (0:nt - 1)*tr;

% % set window size, in seconds, must be disivible by TR
% window_size = 120;

% create arrays indicating time points of each segment, with
% 50% overlapping btw. neighboring segments
ti_list = 0;   % list of ti (initial time) 
tf_list = window_size;   % list of tf (final time)
while tf_list(end) + window_size/2 < taxis(end)   % while tf < taxis(end)
    ti_list = [ti_list ti_list(end) + window_size/2];   % update lists
    tf_list = [tf_list tf_list(end) + window_size/2];
end

% get number of segments expected to form from time series, 
num_seg = numel(tf_list);  
 
% get arrays of event def and event onsets (in seconds, relative to 
% the start of fmri scan) from structure
event_def = event_onsets.results.event_def;
evti = event_onsets.results.evti;

% initialize array for storing number of spikes in each segment, for each
% type of event. i.e. col. 1 = event type, col. 2 = num. of spikes assoc.
num_spikes = cell(length(event_def), 2);
% num_spikes = cell(1, numel(evti));

% for each type of event, loop through each segment in time, then get spike
% rate in that seg. as number of entries that falls into that segment 
for event_ind = 1:numel(evti)   % for each event type
evti_int = evti{event_ind};   % get list of onset times for current event type
    for seg_ind = 1:length(ti_list)   % for each segment in time
        seg_ti = ti_list(seg_ind);   % get ti (initial time) of curr. seg.
        seg_tf = tf_list(seg_ind);   % get tf (final time) of curr. seg.
        % get number of spikes in curr. seg. by counting number of entries
        % in the list of event onsets that falls into time interval of curr.
        % seg.
        num_spikes{event_ind, 1} = event_def{event_ind};
        num_spikes{event_ind, 2}{seg_ind} = numel(find(evti_int >= seg_ti & evti_int <= seg_tf));
    end
end

% END Part (II): get spike rate in each segment 
%---------------------------------------------------------------------------

% Part (III): get amplitude of low-frequency fluctuation (ALFF) in each segment 

% freq. range interested in Hz, witth freq. soln. within the order of
% magnitude of default (1/256*Ts), Ts = sampling period in seconds
f_range = 0.01:0.001:0.08;   % freq. range interested in Hz, 

% initialize arrays of interest
alff_3d = cell(ni, nj, nk);   % map of alff
% mean_alff = zeros(1, num_seg);   % global mean alff for every segment
% mean_alff = [];   % global mean alff
% norm_alff = [];   % normalized alff
       
% scan through each voxel and calculate alff, which is defined as mean
% of square root of the power spectrum density over the low freq. band (0.01 -
% 0.08) Hz
for i = 1:ni
    for j = 1:nj
        for k = 1:nk
            % if current voxel at explicit mask is not zero
            if expmask(i, j, k) ~= 0  
                
                for seg_ind = 1:length(ti_list)   % for every segment
                    % find indices of ti and tf of current segment in taxis
                    ti_ind = find(taxis == ti_list(seg_ind));   % index of ti of cur. seg.
                    tf_ind = find(taxis == tf_list(seg_ind));   % index of tf of cur. seg.
                    
                    % obtain segment of time series required,
                    % then use pwelch to get psd estimate
                    time_series = squeeze(motion_filtered_img(i, j, k, ti_ind:tf_ind));
                    pxx = pwelch(time_series, [], [], f_range, fs);
                    
                    % get alff, which is mean of sq. root of psd estimate
                    alff_3d{i,j,k}{seg_ind} = mean(sqrt(pxx));
                    
                end 
            else
                % assign empty matrix if vox. is not in explicit mask
                alff_3d{i, j, k}{seg_ind} = [];
            end
        end
    end
end

% get global average alff for every segment 
mean_alff = zeros(1, num_seg);   % initialize array for each seg.
for seg_ind = 1:length(ti_list)   % for every segment
    seg_sum = 0;   % initialize sum as 0
    for i = 1:ni   % for each voxel in x-dir.
        for j = 1:nj
            for k = 1:nk
                % check if curr. voxel is included in epxlicit mask, 
                % if so, add value of curr. seg. at voxel to sum
                if expmask(i, j, k) ~= 0  
                    seg_sum = seg_sum + alff_3d{i, j, k}{seg_ind};
                end
            end
        end
    end
    % after looping through all voxels, get mean of curr. seg. by dividing
    % sum by number of non-zero entries in explicit mask
    mean_alff(seg_ind) = seg_sum/(nnz(expmask));
end
            
% get normalized alff for every segment, which is defined as alff/its global mean 
norm_alff = cell(ni, nj, nk);   % initialize array 
for seg_ind = 1:length(ti_list)   % for each segment, then for each voxel
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                % check if curr. voxel is included in explicit mask, 
                % if so, calculate normalized alff for curr. seg. and voxel
                if expmask(i, j, k) ~= 0
                    norm_alff{i, j, k}{seg_ind} = alff_3d{i, j, k}{seg_ind}/mean_alff(seg_ind);
                else
                    % if not, assign empty value to curr. seg.
                    norm_alff{i, j, k}{seg_ind} = [];
                end
            end
        end
    end
end

% END Part (III): get amplitude of low-frequency fluctuation (ALFF) in each 
% segment
%---------------------------------------------------------------------------

% Part (IV): sort electrodes by their types

% use sort_elect to sort electrode array by their types
ele_sorted = sort_elect(electar);

% END PART (IV): sorte electrodes by their types
%---------------------------------------------------------------------------
% % Part (V): find midpoint of each pair of electrode contacts in both mni
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

% END Part (V): find midpoint of each pair of electrode contacts in both mni
% and image space
%---------------------------------------------------------------------------

% Part (VI): get names and coordinates of channels with high spike rates
% during fmri scan

% open input file, assign read permission, and obtain file identifier 
% each row of the input file contains subject number, event marking, the 
% corresponding event type, and its location
event_file_path = '/work/levan_lab/mtang/fmri_project/event_types_locations.txt';
file_id = fopen(event_file_path, 'r');

% % read header of input file by appyling the format spec. (%s) four times, 
% % then scan the remainder of file using format spec. (%s %s %s %s),
% which denotes the four columns of data of string type,
% with empty space separating each field in format spec., but delimiter of
% comma is used for extracting data in each line.
data_header = textscan(file_id, '%s', 4, 'Delimiter', ',');
data_cell = textscan(file_id, '%s %s %s %s', 'Delimiter', ',');   % comma-delimited file

% close file using the file identifier
fclose(file_id);

% use strcmp (string compare) to find rows that have subject number equal
% to the one that is being examined, then use find to get non-zero entries
% from output of strcmp, assign results to row_req
row_req = find(strcmp(data_cell{1}, subnum));

% use required row numbers found to obtain event types and channel names 
% associated with the current subject number interesed, 
% get event types and their corresponding locations (channels)

% initializr arrays, with each cell storing relevant info. about current
% event type. i.e., 1st cell in ev. location array stores locations
% corresponding to the 1st event type recorded in ev types array
ev_types = cell(length(row_req), 1);   % initialize event types array
ev_locations = cell(length(row_req), 1);   % intiialize ev. locations array

% scan across each row number found, select the respective entry of event
% type, and the associcated names of channels
% if more than one location is assigned for the event type, each location
% is separated from others by empty space char., use strsplit to obtain
% every location recorded.
for row_ind = 1:length(row_req)   % for each row number interested
    ev_types{row_ind} = str2num(data_cell{3}{row_req(row_ind)});   % get event type
    ev_str = strsplit(data_cell{4}{row_req(row_ind)}, ' ');   % get every location
    
    % sometimes, extra spaces are found in input file, clean out extra
    % empty spaces contained in location char.
    ev_str_cleaned = {};   % initialize array to store cleaned location names
    item_count = 1;   % initialize item count as 1
    for item = 1:numel(ev_str)   % for each location obtained
        if ~isempty(ev_str{item})   % if current cell is not empty
            % replace extra empty space ' '  by ''
            ev_str_cleaned{item_count} = strrep(ev_str{item}, ' ' , '');
            item_count = item_count + 1;   % increment item count by 1 for the next location found
        end
    end
    % assign cleaned event location cell array to current row number
    ev_locations{row_ind} = ev_str_cleaned;  
end

% initialize ch_name array, with the following format, 
% i^th cell in 1st layer denotes the i^th event type identified, 
% j^th cell in the 2nd layer denotes the j^th electrode name for the i^th
% event type
% k^th cell in the 3rd layer denotes the k^th channel for the j^th
% electrode name that corresponds to the i^th event type
% e.g. ch_name{1}{2}{3} denotes the 3rd channel name under the 
% 2nd electrode type that corresponds to the 1st event type
ch_name = cell(length(ev_locations), 1);

% scan across each event type interested, then scan through each location
% name recorded, use regexp to search for alphabetic and numeric patterns, 
% e.g. given LA1-5, further decompose the number sequence to intergal
% increments, so that we have LA1-2, LA2-3, LA3-4, LA4-5
for ev_ind = 1:length(ev_locations)
    for ch_ind = 1:length(ev_locations{ev_ind})
        ch_type = regexp(ev_locations{ev_ind}{ch_ind}, '[A-Z]*', 'match');   % alphabetic compo.
        ele_con_num = regexp(ev_locations{ev_ind}{ch_ind}, '[\d]*', 'match');   % numeric comp.
    
        % initialize count of number pairs for current electrode type
        pair_count = 1;  
        
        % sometimes, numbers are not recorded (i.e. LA instead LA1-6), if
        % check if numeric pattern is present, if so, form number sequence
        if ~isempty(ele_con_num)
            num_i = str2num(ele_con_num{1});   % intitial number of sequence
            num_f = str2num(ele_con_num{2});   % final number of sequence
            num_s = num_i;   % assign start number
            
            while num_s < num_f   % as long as start number is less than final number
                num_pair = num_s:num_s+1;   % form number pair
                ch_name{ev_ind}{ch_ind}{pair_count} = [ch_type{:}, num2str(num_pair(1)), '-', ...
                ch_type{:}, num2str(num_pair(2))];   % construct channel name
                num_s = num_s + 1;   % increment start number by 1 for next entry
                pair_count = pair_count + 1;   % increment pair count by 1 for next entry
            end
        else   % if number is not recorded, leave it as it is
            ch_name{ev_ind}{ch_ind}{pair_count} = ch_type{:};
            pair_count = pair_count + pair_count;
        end
    end
end

% convert ch_name to single layer, with event type stored under the 1st
% column, and the associated channel names stored under the 2nd column
ch_name_all = {};   % initialize array
ch_ind = 1;   % initiliaze channel index (row number)
for ev_ind = 1:length(ch_name)   % for every event type
    for ch_type = 1:length(ch_name{ev_ind})   % for every electrode under curr. event type
        for item = 1:length(ch_name{ev_ind}{ch_type})   % for every name under curr. electrode 
            ch_name_all{ch_ind, 1} = ev_types{ev_ind};   % get event type associated
            ch_name_all{ch_ind, 2} = ch_name{ev_ind}{ch_type}{item};   % get channel names associated
            ch_ind = ch_ind + 1;   % increment channel index by 1 
        end
    end
end

% match channel name obtained to those in midpt_all array, if matched,
% assign all relevant info., including midpt corrdinates in mni and image
% spaces.
% if contact number is not given, match all available contacts for the channel
% type. i.e. if given LMT, instead of LMT1-2, then match all available
% contacts in midpt_all, which gives LMT1-2, LMT2-3, and so on.

% initialize array with format given below,
% 1st col. denotes event type
% 2nd col. denotes the corresponding channel name 
% 3rd col. denotes midpt coordinates in mni space
% 4th col. denotes midpt coordinates in image space
ch_name_matched = {};   % initialize array
ch_matched_ind = 1;   % initialize matched channel index as 1

for i = 1:size(ch_name_all, 1)   % for every channel name found
    % if contact numbers are given, '-' is included, check if name string
    % contains this '-' pattern, if so, check if name strings in ch_name
    % and midpt_all agree, if so, add channel name to ch_name_matched array
    if contains(ch_name_all{i, 2}, '-')
        if ~isempty(find(strcmp(ch_name_all{i, 2}, midpt_all(:, 1))))
            ind_found = find(strcmp(ch_name_all{i, 2}, midpt_all(:, 1)));   % get index
            ch_name_matched{ch_matched_ind, 1} = ch_name_all{i, 1};   % event type
            ch_name_matched(ch_matched_ind, 2:4) = midpt_all(ind_found, :);   % get all info.
            ch_matched_ind = ch_matched_ind + 1;
        end
    % if contact number is not given, search if name exists in midpt_all, 
    % if so, assign all channels with name matched in midpt_all to array
    else 
        % format search pattern, use regexp to find rows with matched
        % pattern, then use find to get those rows with non-zero values
        search_pat = [ch_name_all{i, 2}, '[\d*]-', ch_name_all{i, 2}, '[\d*]']; 
        % if matched pattern is found, non-zero value is assigned, else,
        % empty value is assigned to the cell
        row_srch = regexp(midpt_all(:, 1), search_pat);  
        % assign value of 0 to all empty cells 
        for row_num = 1:numel(row_srch)
            if isempty(row_srch{row_num})
                row_srch{row_num} = 0;
            end
        end
        % use find to get indices of non-zero entries
        ind_req = find(cell2mat(row_srch));
        % for each rows with non-zero entries found, get corresponding
        % channel name from midpt_all array
        for item = 1:numel(ind_req)
            ch_name_matched{ch_matched_ind, 1} = ch_name_all{i, 1};   % event type
            ch_name_matched(ch_matched_ind, 2:4) = midpt_all(ind_req(item), :);   % chnannel name
            ch_matched_ind = ch_matched_ind + 1;
        end
    end
end

% END Part (VI): get names and coordinates of channels with high spike rates
% during fmri scan
%---------------------------------------------------------------------------

% Part(VII): compute mean normalized ALFF in 3X3X3 box of each channel
% that has high spike rates recorded during fmri scan

ch_mnalff = ch_name_matched(:, 1:2);

for i = 1:size(ch_name_matched, 1)   % for each channel
    mxyz = ch_name_matched{i, 4};
    
    % create regional mask from explicit mask using the coordinates of 
    % the 3X3X3 box
    reg_mask = double(expmask(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ... 
        mxyz(3)-1:mxyz(3)+1));
    
    % get number non-zero entries in explicit mask in current box
    num_nonzero = nnz(reg_mask);
    
    % for every segment in time, scan across every voxel in current
    % 3X3X3 box, (i) sum all normalized alff values thereof, 
    % (ii) then get mean normalized alff in that box by dividing
    % sum by number of non-zero entries of that box in explicit mask
    for seg_ind = 1:length(ti_list)   % for every segment in time
        seg_sum_norm_alff = 0;   % initialize sum of norm. alff as 0 
        for x_vox = mxyz(1)-1:mxyz(1)+1   % scan across every voxel in curr. box
            for y_vox = mxyz(2)-1:mxyz(2)+1
                for z_vox = mxyz(3)-1:mxyz(3)+1
                    % if current voxel is included in explicit mask
                    if expmask(x_vox, y_vox, z_vox) == 1                 
                        % compute sum of all norm. alff values in curr. box
                        seg_sum_norm_alff = seg_sum_norm_alff + ...
                            norm_alff{x_vox, y_vox, z_vox}{seg_ind};
                    end
                end
            end
        end
        % compute mean normalized alff in curr. box by dividing the sum
        % by number of non-zero entries of curr. box. in explicit mask
        % assign calculated mean to curr. cell in 3rd col. of ch_mnalff
        ch_mnalff{i, 3}{seg_ind} = seg_sum_norm_alff/num_nonzero;
    end   % end for each segment in time
end   % end for each channel


% Part(VII): compute mean normalized ALFF in 3X3X3 box of each channel
% that has high spike rates recorded during fmri scan
%---------------------------------------------------------------------------

% Part(VIII): get Spearman's corr. coeff. btw. spike rates and mean
% normalized alff of every channel selected

% initialize arrays for storing spearman's corr. coeff. and its p-value for
% each channel selected, format of array is given below, 
% 1st col. = event type
% 2nd col.= channel name
% 3rd col. = spearman corr. coeff. OR p-value, depending on array
ch_rho_sp = cell(size(ch_mnalff, 1), 3);   % spearman's corr. coeff. 
ch_pval_sp = cell(size(ch_mnalff, 1), 3);   % p-value
rnum = 1;   % initialize row number as 1

for row = 1:size(num_spikes, 1)   % for each event type 
    curr_ev = num_spikes{row, 1};   % get current event type
    for ch = 1:size(ch_mnalff, 1)   % for each channel 
        % check if current event type agress with that of current channel, 
        % if so, get spearman corr. coeff. and its p-value  
        if ismember(ch_mnalff{ch, 1}, curr_ev)  
            [rho_sp, pval_sp] = corr(cell2mat(num_spikes{row, 2})', ...
                cell2mat(ch_mnalff{ch, 3})', 'Type', 'Spearman');
            ch_rho_sp(rnum, 1:2) = ch_mnalff(ch, 1:2);   % assign event type and channel name
            ch_rho_sp{rnum, 3} = rho_sp;   % assign spearman corr. coeff.
            ch_pval_sp(rnum, 1:2) = ch_mnalff(ch, 1:2);   % assign event type and channel name
            ch_pval_sp{rnum, 3} = pval_sp;   % assign p-value
            rnum = rnum + 1;   % increment row number by 1 for next entry
        end
    end
end

% End Part(VIII): get Spearman's corr. coeff. btw. spike rates and mean
% normalized alff of every channel selected

%---------------------------------------------------------------------------
% Part (IX): convert cell array to table

% create cell array for table,
% col. 1 = event type
% col. 2 = channel name
% col. 3 = spearman rho btw. num. spikes and mean normalized alff of
% channel
% col. 4 = p-value associated,
% all numeric values are corrected to 3 decimal place
ch_rho_ns_mnalff = cell(size(ch_rho_sp, 1), 4);
for i = 1:size(ch_rho_sp, 1)
    ch_rho_ns_mnalff(i, 1:2) = ch_rho_sp(i, 1:2);   % event type and channel name
    ch_rho_ns_mnalff{i, 3} = round(ch_rho_sp{i, 3}, 3);   % spearman rho, corr. to 3 d.p.
    ch_rho_ns_mnalff{i, 4} = round(ch_pval_sp{i, 3}, 3);   % p-value, corr. to 3 d.p.
end

% convert cell array to table and add table headers
table_ch_rho_ns_mnalff = cell2table(ch_rho_ns_mnalff, ...
    'VariableNames',{'event type', 'channel name', ['corr. coeff._', num2str(window_size)], ...
    ['p-val._', num2str(window_size)]});

% End Part (IX): convert cell array to table
%---------------------------------------------------------------------------

% Part (X): store quantities calcu. in output struct.

opstruct = struct;
opstruct.ele_sorted = ele_sorted;   % info. of electrode (sorted by name)
opstruct.midpt = midpt;   % coordinates of midpoint between each pair of electrode contacts 
opstruct.midpt_all = midpt_all;   % midpt, but in single layer

opstruct.window_size_used = window_size;   % current window size (s) 
opstruct.ti_list = ti_list;   % list of initial times  for all segments constructed
opstruct.tf_list = tf_list;   % list of final times  for all segments constructed
opstruct.num_spikes = num_spikes;   % number of spikes in each seg., for each event type

opstruct.alff_3d = alff_3d;   % alff at each voxel
opstruct.mean_alff = mean_alff;   % global mean alff
opstruct.norm_alff = norm_alff;   % normalized alff at each voxel

opstruct.ch_mnalff = ch_mnalff;   % matched channel names with mean normalized alff in 3X3X3 box of each channel
opstruct.ch_rho_sp = ch_rho_sp;   % spearman rho btw. num. of spikes and mean normalized alff for each channel
opstruct.ch_pval_sp = ch_pval_sp;   % p-value of spearman rho

opstruct.ch_rho_ns_mnalff = ch_rho_ns_mnalff;   % array spearman rho btw. num. spikes and mnalff with p-values attached
opstruct.table_ch_rho_ns_mnalff = table_ch_rho_ns_mnalff;   % table listing calcu. quant.

% END Part (X): store quantities calcu. in output struct.
%---------------------------------------------------------------------------

end   % end if com_alff == 1

end   % end function opstruct = get_mnalff_num_spikes(input_array)