clear all 
close all

% This script computes Spearman's corr. coeff. btw. Reho and number of
% spikes recorded DURING fmri scan for every channel selected. Motion
% artifacts and noises are filtered before calculating Reho. Fixed window
% size is used for segmentation of time series. 

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
fname_op = [directname, filesep, 'matrices' filesep 'reho_spikes_fixed_winsize'];   % direct. of output matrix
filename_op = 'reho_spikes_fixed_windsize.mat';   % filename of output file

% enter window size in seconds (must be divisible by tr, i.e. 1.5 s)
window_size = 120;

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
input_array.window_size = window_size;   % window size entered

% get indices of current run number embedded in processed image filename
sind = regexp(lower(swra_file_path{run_ind}), 'run');   % start index
eind = regexp(lower(swra_file_path{run_ind}), '.nii');   % end index
sess_num = swra_file_path{run_ind}(sind+3:eind-1);   % get indices of run number

% create fieldname representing current session number
fieldname = sprintf('run%s', sess_num);

% store output calcu. under current fieldname of terms
terms.(fieldname) = get_reho_num_spikes_fixed_winsize(input_array);

% execute lines below if current field is not empty
if ~isempty(terms.(fieldname))
    
% access tables of current session created by function get_rho_spikes
table_ch_rho_sp_reho_spikes = terms.(fieldname).table_ch_rho_sp_reho_spikes;   % table of spearman's rho btw. ReHo and num. spikes

% format filenames and full paths of tables
filename_table = ['sub', subnum, '_table_rho_sp_reho_spikes_run_', sess_num, '.csv'];
table_path = fullfile(fname_op, filename_table);   % path of table of clin. det. channels

%---------------------------------------------------------------------------
% output table and structure created in current session
if op_results == 1
    writetable(table_ch_rho_sp_reho_spikes, table_path);   % table of clin. det. ch.
    save(fullfile(fname_op, filename_op), 'terms');   % struct.
end
%---------------------------------------------------------------------------

end   % end if ~isempty(terms.(fieldname))

end   % end for run_ind = 1:numel(swra_file_path)

toc
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
% (A): EXECUTE THE FOLLOWING AFTER CALLING FUNCTION ON EACH RUN

% access fieldnames of struct. terms
fds = fieldnames(terms);

% get non-empty fieldnames
fd_ind = [];   % initialize array for storing indices of non-empty fields
for sess_ind = 1:length(fds)   % for each field
    if ~isempty(terms.(fds{sess_ind}))   % check if field is empty
        fd_ind = [fd_ind sess_ind];   % if not, append index to array
    end
end
fdnames = fds(fd_ind);   % access all non-empty fields using indices acquired

% execute the lines below if there exists at least one non-empty field
if length(fdnames) >= 1

% transform array to single layer, concatenating num_spikes and
% sig_box arrays from diff. runs together
num_spikes_all = {};   % initialize array
sig_box_all = {};
for sess_ind = 1:length(fdnames)   % for each session
    % get num_spikes array for current session;   
    num_spikes_cur = terms.(fdnames{sess_ind}).num_spikes;
    sig_box_cur = terms.(fdnames{sess_ind}).sig_box;
    % concatenate arrays vertically
    num_spikes_all = [num_spikes_all; num_spikes_cur];
    sig_box_all = [sig_box_all; sig_box_cur];
end

%--------------------------
% (A1): get unique event types

% get unique event types for all sessions
unique_ev_type = {};   % initialize array

for item = 1:size(num_spikes_all, 1)   % for each item available
    count = 0;   % initialize count as 0 for current item
    for j = 1:size(unique_ev_type, 1)   % for each unique entry identified
        % if any entry in num_spikes_all exists in uni. array, increment 
        % count by 1
        if isequal(num_spikes_all(item), unique_ev_type(j))
            count = count + 1;
        end
    end
    % if count remains zero at the end of the scan, it means no entry in
    % uni. array agrees with any item in num_spikes_all, add current
    % item to uni. array
    if count == 0
        unique_ev_type = [unique_ev_type; num_spikes_all(item, 1)];
    end
end

% END (A1): get unique event types
%--------------------------

% (A2): find all sessions at which each event type appears

% some event types do not appear in every session, the following code
% identities sessions in which every event type is found

% initialize array, with intended format given below
% 1st col. = unique event type
% 2nd col. = session indices where unique ev. type appear
ev_sess_ind = unique_ev_type;  

for item = 1:numel(unique_ev_type)   % for each unique event type
    cur_uni_ev_type = unique_ev_type{item};   % get current unique event type
    ind_found = [];   % initialize array for storing sess. ind. found for this uni. ev. type
    count = 0;   % initialize count as zero
    for sess_ind = 1:length(fdnames)   % for each session
        cur_num_spikes = terms.(fdnames{sess_ind}).num_spikes;   % get current num. spikes array
%         count = 0;
        for i = 1:size(cur_num_spikes, 1)   % for every row in array
            % check if event type of current row matches with current
            % unique event type, if so, increase count by 1
            if isequal(cur_num_spikes{i, 1}, cur_uni_ev_type)
                count = count + 1;
            end
        end
        % if count is greater than 0, rows with matching event type are
        % found, add session index to array, ind_found
        if count > 0
            ind_found = [ind_found sess_ind];
        end
    end
    % store sess. ind. found to 2nd col. of array
    ev_sess_ind{item, 2} = ind_found;
end

% END (A2): find all sessions at which each event type appears
%--------------------------

%--------------------------
% (A3): combine number of spikes from different runs together
num_spikes_comb = unique_ev_type;
for item = 1:size(ev_sess_ind, 1)
    cur_ev_type = ev_sess_ind{item, 1};
    spikes_found = {};
    for sess_ind = 1:length(fdnames)   % for each session
        num_spikes_cur_sess = terms.(fdnames{sess_ind}).num_spikes;
        for i = 1:size(num_spikes_cur_sess, 1)
            if isequal(num_spikes_cur_sess{i, 1}, cur_ev_type)
                spikes_found = [spikes_found num_spikes_cur_sess{i, 2}];
            end
        end
    end
    num_spikes_comb{item, 2} = spikes_found;
end

% END (A3): combine number of spikes from different runs together
%--------------------------

% (A4): combine ReHo of channels selected from different runs together

% initialize array for storing ReHo from all sessions together, 
% 1st col. = unique event type
% 2nd col. = channel name
% 3rd col. = ReHo from all sessions associated
ch_reho_comb = {};
rnum = 1;   % initialize row number as 1

for item = 1:size(ev_sess_ind, 1)   % every unique event type
    cur_ev_type = ev_sess_ind{item, 1};   % get current unique event type
    sess_ind_req = ev_sess_ind{item, 2};   % get sess. ind. required
    
    % store channel names and their reho values for current event type
    ch_reho_req = {};   
    for sess_ind = sess_ind_req(1):sess_ind_req(end)   % every sess. req.
        % get sig_box array in current session
        sig_box_cur = terms.(fdnames{sess_ind}).sig_box;   
        
        % scan across every row in current sig_box array, select rows with
        % event type matched with current unique event type
        for row = 1:size(sig_box_cur, 1)   % every row in sig_box array
            if ismember(sig_box_cur{row, 1}, cur_ev_type)   % if ev. types agree
                ch_reho_req =  [ch_reho_req; sig_box_cur(row, :)];   % concatenate rows vertically
            end
        end
    end
    
    % from the list of channels (and their ReHo) obtained, get unique
    % channel names
    if ~isempty(ch_reho_req)   % if channel name is matched for cur. ev. type
    unique_ch = unique(ch_reho_req(:, 2));
   
    % for each unique channel name, combine their mnalff from all sessions
    % required
    for ch = 1:numel(unique_ch)
        % get row indices for current unique channel name
        ind_req = find(strcmp(ch_reho_req(:, 2), unique_ch{ch}));   
        cell_obt = ch_reho_req(ind_req, :);   % obtain cells 
        ch_reho_comb(rnum, 1:2) = cell_obt(1, 1:2);   % assign event type and ch. name
        ch_reho_comb{rnum, 3} = [cell_obt{:, 3}];   % combine reho from all sess. 
        rnum = rnum + 1;   % increment row number by 1 for next channel
    end
    
    end   % end if ~isempty(ch_mnalff_req)
    
end   % end for item = 1:size(ev_sess_ind, 1)

% END (A4): combine ReHo of channels selected from different runs together
%--------------------------

% (A5): get spearman's corr. coeff. btw. num. of spikes and ReHo from all runs combined,

% initialize array with format given below, 
% 1st col. = event type
% 2nd col.= channel name
% 3rd col. = spearman corr. coeff.
% 4th col. = p-value

ch_rho_sp_comb = ch_reho_comb(:, 1:2);   

% scan across every row of ch_reho_comb array
for i = 1:size(ch_reho_comb, 1)
    ch_ev = ch_reho_comb{i, 1};   % get cur. uni. ev. type
    ch_name = ch_reho_comb{i, 2};   % get cur. channel name
    
    % scan across every row in num_spikes_comb array
    for j = 1:size(num_spikes_comb, 1)
        % if event type in num_spikes_comb matches with cur. uni. ev. type
        % assign num. spikes to variable num_spikes_found, converted to
        % matrix
        if ismember(ch_ev, num_spikes_comb{j, 1})
            num_spikes_found = cell2mat(num_spikes_comb{j, 2})';
        end
    end
    
    % convert current cell of ReHo to matrix 
    ch_reho_found = cell2mat(ch_reho_comb{i, 3})';
    
    % calculate Spearman's rho and p-value btw. num. of spikes and ReHo
    % matrices
    [rho_sp_comb, pval_sp_comb] = corr(num_spikes_found, ch_reho_found, ...
        'Type', 'Spearman');
    
    % assign data calculated to corresponding columns, round them to 3 d.p.
    ch_rho_sp_comb{i, 3} = round(rho_sp_comb, 3);   % spearman's rho
    ch_rho_sp_comb{i, 4} = round(pval_sp_comb, 3);   % p-value
    
end

% END (A5): get spearman's corr. coeff. btw. num. of spikes and ReHo from all runs combined,
%--------------------------

% (A6): combine table arrays from all sessions together

% read table arrays (of all runs) from struct., combine them 
% format of table array is given below,
% 1st col. = event type
% 2nd col. = channel name
% 3rd col. = spearman rho
% 4th col. = p-value, all table arrays stored in struct have event types
% and channel names arranged in same order
rho_pval_comb_all = ch_rho_sp_comb(:, 1:2);   % initialize array 

for sess_ind = 1:length(fdnames)   % for each run
    ar_read = terms.(fdnames{sess_ind}).ch_rho_sp;   % read table array in current run
    for i = 1:size(rho_pval_comb_all, 1)    % scan every row in rho_pval_comb_all
        for j = 1:size(ar_read, 1)   % scan every row in table array read
            % compare event types and channel names btw. both arrays, if
            % they agree, add corresponding col. to rho_pval_comb_all
            if isequal(ar_read{j, 1}, rho_pval_comb_all{i, 1}) && strcmp(ar_read{j, 2}, rho_pval_comb_all{i, 2})
                rho_pval_comb_all{i, 3+(sess_ind-1)*2} = ar_read{j, 3};   % spearman rho
                rho_pval_comb_all{i, 4+(sess_ind-1)*2} = ar_read{j, 4};   % p-value
            end
        end
    end
end

% next, add spearman rho and p-value obtained from all num. spikes and ReHo 
% obtained from all runs to last 2 col. of array 
ncol = size(rho_pval_comb_all, 2);   % get num. of col. established
rho_pval_comb_all(:, ncol+1:ncol+2) = ch_rho_sp_comb(:, 3:4);

% initialize array to store headers required for table
headers_req = {};
runind = 1;   % initialize run index

% first two col. are event type and channel name
for j = 1:ncol-2   % for every remaining col. 
    if mod(j, 2) == 1   % if col. num. is odd
        headers_req{j} = ['corr. coeff._', num2str(window_size), ...
            '_run_ind_', num2str(runind)];   % construct details of corr. coeff.
    end
    if mod(j, 2) == 0   % if col. num is even
        headers_req{j} = ['pval_', num2str(window_size), ...
            '_run_ind_', num2str(runind)];   % construct details of p-value
        %   once both details about corr. coeff. and p-value are
        %   constructed, increment runind by 1 for next run
        runind = runind + 1;
    end
end
headers_req{ncol-2+1} = ['corr. coeff._', num2str(window_size), ...
            '_run_ind_', 'all'];   % construct details of corr. coeff.
headers_req{ncol-2+2} = ['pval_', num2str(window_size), ...
            '_run_ind_', 'all'];
        
% convert cell array to table and add table headers
table_rho_pval_comb = cell2table(rho_pval_comb_all, ...
    'VariableNames',{'event type', 'channel name', headers_req{:}});

% format filename for table
filename_table_comb = ['sub', subnum, '_table_ns_reho_run_all' ...
    '_winsize_', num2str(window_size), '.csv'];
table_comb_path = fullfile(fname_op, filename_table_comb);

% write resultant table to path 
if op_results == 1
    writetable(table_rho_pval_comb, table_comb_path);   % table of clin. det. ch.
end

% store rho_pval_comb array (spearman corr. coeff. btw. num. spikes and 
% reho for every channel selected) under every fieldname in struct.
for sess_ind = 1:length(fdnames)
    terms.(fdnames{sess_ind}).rho_pval_comb_all = rho_pval_comb_all;
end

% update terms struct. to path
if op_results == 1
    save(fullfile(fname_op, filename_op), 'terms', '-v7.3');
end

end   % end if length(fdnames) >= 1



%---------------------------------------------------------------------------
% FUNCTION FOR GETTING MEAN ALFF AND NUMBER OF SPIKES IN EACH SEGEMENT IN
% EVERY RUN
%---------------------------------------------------------------------------
% define function get_mnalff_num_spikes_fixed_winsize requiring an input array, which
% has the fields defined below,
% SWRA_IMG is the data of processed image file
% SWRA_INFO is the info (headers) associated with the processed image file
% EXPMASK is the explicit mask image
% ELECTAR is the array containing info. about all electrodes used,
% EVENT_ONSETS is the onset times of each event type during fmri scan
% MOTION_FILE_PATH is the path of motion parameter file
% and outputs,
% WINDOW_SIZE is the window size chosen in seconds
% OPSTRUCT, a structure that contains info. about coordinates of
% electrdoes, spearman corr. coeff. btw. number of spikes and ReHo for each 
% channel selected,
% see end of function for details of variables included in this struct.

function opstruct = get_reho_num_spikes_fixed_winsize(input_array)

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

% Part (V): get names and coordinates of channels with high spike rates
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

% match channel names obtained to those in midpt_all array, if matched,
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

% END Part (V): get names and coordinates of channels with high spike rates
% during fmri scan
%---------------------------------------------------------------------------

% Part (VI): get ReHo for every channel matched

% ReHo is calculated for all voxels enclosed in 3X3X3 box centered at
% midpoint of every channel matched

% get basic info. about the processed image
tr = swra_info.PixelDimensions(4);   % repetition time, in seconds
fs = 1/tr; % sampling frequency, in Hz
% 
% % initialize array, sig_box, for storing the following for every channel matched
% % 1st column lists event type assoc. with each channel
% % 2nd column lists the names of channels,
% % 3rd column lists reho (aka kendall_w coeff.) of every channel
sig_box = {};

for item = 1:size(ch_name_matched, 1)   % every channel matched
    
    % assign event type and assoc. channel name to first two columns
    sig_box(item, 1:2) = ch_name_matched(item, 1:2);
    
    % get coordinates of midpt in image space
    mxyz = ch_name_matched{item, 4};
    
    % using coordinates of mid. pt. in image space to obtain signal
    % enclosed in 3X3X3 box centered thereof
    signal = double(motion_filtered_img(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
            mxyz(3)-1:mxyz(3)+1, :));  
    
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
    
    % create regional mask from explicit mask using the coordinates of
    % the 3X3X3 box
    reg_mask = double(expmask(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
            mxyz(3)-1:mxyz(3)+1));
      
    % after motion artifacts and noises are removed, signal is then
    % decomposed into segments of window size specificed by user. Each
    % segment is then used to compute Kendall's W coeff., (ReHo)
    
    for seg_ind = 1:length(ti_list)   % for every segment
        % find indices of ti and tf of current segment in taxis
        ti_ind = find(taxis == ti_list(seg_ind));   % index of ti of cur. seg.
        tf_ind = find(taxis == tf_list(seg_ind));   % index of tf of cur. seg.
            
        % extract portion of time series using indices obtained above
        time_series_ext = signal_noise_filt(:, :, :, ti_ind:tf_ind);
            
        % for current segment, get Kendall's W coefficient of filtered signal, 
        % assign quantities calcu. to 2nd column of current cell within sig_box
        sig_box{item, 3}{seg_ind} = kendall_w(time_series_ext, reg_mask);
    
    end   % end for seg_ind = 1:length(ti_list)   % for every segment

    
end   % end for item = 1:size(ch_name_matched, 1)   % every channel matched

% END Part (VI): get ReHo for every channel matched
%---------------------------------------------------------------------------

% Part (VII): get Spearman's corr. coeff. btw. ReHo and number of spikes
% recorded during fmri scan

 % initialize array, with intended format given below,
 % 1st col. = event type
 % 2nd col. = channel name
 % 3rd col. = spearman's rho btw. ReHo and num. of spikes
 % 4th col. = p-value
ch_rho_sp = cell(size(sig_box, 1), 3);  
rnum = 1;

for row = 1:size(num_spikes, 1)   % each event type 
    cur_ev = num_spikes{row, 1};   % get current event type
    for ch = 1:size(sig_box, 1)   % every channel in sig_box
        % check if event type of channel agrees with current event type, 
        % if so, compute spearman's corr. coeff. and p-value
        if ismember(sig_box{ch, 1}, cur_ev)   
            [rho_sp, pval_sp] = corr(cell2mat(num_spikes{row, 2})', ...
                cell2mat(sig_box{ch, 3})', 'Type', 'Spearman');
            ch_rho_sp(rnum, 1:2) = sig_box(ch, 1:2);   % assign event type and channel name
            ch_rho_sp{rnum, 3} = round(rho_sp, 3);   % assign spearman's rho, corr. to 3 d.p.
            ch_rho_sp{rnum, 4} = round(pval_sp, 3);   % assign p-value, corr. to 3 d.p.
            rnum = rnum + 1;   % increment row number by 1 for next channel
        end
    end
end

% convert cell array to table and add table headers
table_ch_rho_sp_reho_spikes = cell2table(ch_rho_sp, ...
    'VariableNames',{'event type', 'channel name', ['corr. coeff._', num2str(window_size)], ...
    ['p-val._', num2str(window_size)]});

% END Part (VII): get Spearman's corr. coeff. btw. ReHo and number of spikes
% recorded during fmri scan
%---------------------------------------------------------------------------

% Part (VIII): store quantities calcu. in output struct.

opstruct = struct;
opstruct.ele_sorted = ele_sorted;   % info. of electrode (sorted by name)
opstruct.midpt = midpt;   % coordinates of midpoint between each pair of electrode contacts 
opstruct.midpt_all = midpt_all;   % midpt, but in single layer

opstruct.window_size_used = window_size;   % current window size (s) 
opstruct.ti_list = ti_list;   % list of initial times  for all segments constructed
opstruct.tf_list = tf_list;   % list of final times  for all segments constructed
opstruct.num_spikes = num_spikes;   % number of spikes in each seg., for each event type
opstruct.sig_box = sig_box;   % array storing ReHo computed for selected channels

opstruct.ch_rho_sp = ch_rho_sp;   % array listing Spearman's rho and p-values btw. ReHo and num. spikes
opstruct.table_ch_rho_sp_reho_spikes  = table_ch_rho_sp_reho_spikes;   % table created using the above array

% END Part (VIII): store quantities calcu. in output struct.
%---------------------------------------------------------------------------

end   % end function opstruct = get_reho_num_spikes_fixed_winsize(input_array)
