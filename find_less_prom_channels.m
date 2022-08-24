close all
clear all

% Given a list of prominent channels for every event type, this script 
% identifies 'less prominent channels' and 'other channels' by first (i) 
% computing the mean-square eeg time series of each type of prominent
% channels, (ii) then aligning actual IED peak at every event onset with
% the marked onset time to obtain the waveform of that IED peak, then (iii)
% averaging waveforms across all event onsets to get the average waveform
% representing each event type. 
% Afterwards, using avg. waveform of each event type as a reference to
% compute cross-covariance (xcov) btw. ref. and every eeg channel to identify the
% lag required such that max xcov is attained at every event onset. Then, 
% align IED peak with the marked event onset using the lag acquired, and
% average across all event onsets to obtain average waveform of every eeg 
% channel with respect to each event type.

% Plotting avg. waveforms of eeg channels together with the
% reference to identify IED-peak-bearing time series for determining if
% channels are classified as 'less prom. ch.' or 'other ch.'
% 
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '27';

% enter directory where all input files are located at
directname =  '/work/levan_lab/mtang/fmri_project/';
% directname =  'C:\Users\siumichael.tang\Downloads\fmri_project\';

% format path to eeg data 
% directname_eeg = [directname, 'sub', subnum, '/eeg_data'];
directname_eeg = [directname, 'sub', subnum, filesep, 'eeg_data'];
filename_eeg = 'filt_eeg_ts.mat';

% format path to event onsets data
% directname_event_onset = [directname, 'sub', subnum, '/matrices'];
directname_event_onset = [directname, 'sub', subnum, filesep, 'matrices'];
filename_event_onset = ['onset_times_events_sub', subnum, '.mat'];

% enter path where ouput struct. is stored at
fname_op = [directname, filesep, 'sub', subnum, filesep, 'eeg_data'];   % direct. of output matrix
filename_op = 'find_less_prom_channels.mat';   % filename of output file

% enter if output is saved to path indicated (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% get paths of input files
% [~, stat_reho_file_path, ~] = get_path(directname_stat_reho, filename_stat_reho);   % stat. analysis, reho
[~, eeg_file_path, ~] = get_path(directname_eeg, filename_eeg);   % filtered eeg data
[~, event_onset_file_path, ~] = get_path(directname_event_onset, filename_event_onset);   % event onset times

%----------------------------------------------------------------------------

% PRELIMINARY: 

% initialize output structure 
opst = struct;

% load structure
% st_reho = load(stat_reho_file_path{:});
st_eeg = load(eeg_file_path{:});
st_event_onsets = load(event_onset_file_path{:});

% obtain fieldnames under each data set, eeg and event onsets, respectively
% then access data for current session
fdnames_eeg = fieldnames(st_eeg.filt_eeg_st);   % fieldnames of st_eeg
fdnames_run_events = fieldnames(st_event_onsets.terms);   % fieldnames of event_onsets

% END PRELIMINARY: 
%----------------------------------------------------------------------------

for run_ind = 1:numel(fdnames_eeg)   % for each session
    
% access data for current session
st_eeg_cur = st_eeg.filt_eeg_st.(fdnames_eeg{run_ind});   % filt. eeg data for current session
event_def = st_event_onsets.terms.(fdnames_run_events{run_ind}).results.event_def;   % event def. of curr. session
evti = st_event_onsets.terms.(fdnames_run_events{run_ind}).results.evti;   % event onsets (in seconds) of curr. sess.

% %----------------------------------------------------------------------------
% % PART (I): ARRANGE EEG TIME SERIES OF EVERY CHANNEL PROVIDED BY INPUT 
% EEG DATASET IN SINGLE LAYER
% 
% initialize array, eeg_ch, with format given below,
% 1st col. = channel name
% 2nd col. = eeg time series
eeg_ch = {}; 
row_num = 1;   % initialize row number
for i = 1:size(st_eeg_cur.filt_eeg, 1)   % for every channel listed in eeg data file
    ch_name = st_eeg_cur.filt_eeg{i, 1};   % assmeble channel name from col. 1

    eeg_ch{row_num, 1} = ch_name;
    eeg_ch{row_num, 2} = st_eeg_cur.filt_eeg{i, end};   % assign full eeg time series
    row_num = row_num + 1;   % increment row numbr by 1 for next channel matched
end

% get fs in Hz of current run from eeg struct.
fs = st_eeg_cur.fs;

% END PART (I): ARRANGE EEG TIME SERIES OF EVERY CHANNEL PROVIDED BY INPUT 
% EEG DATASET IN SINGLE LAYER
% %----------------------------------------------------------------------------
% PART (II): get names, event types, and eeg time series of prominent channels, 
% defined as those with sig. ied rates registered during fmri scan 

% get names of channels associated with event types assigned to current 
% subject number
ch_list = get_channel_assoc_event(subnum, []);

% from list of names obtained, get the corresponding eeg time series from
% eeg_ch

% 1st col. denotes event type
% 2nd col. denotes the corresponding channel name 
% 3rd col. denotes eeg time series (filtered)
prom_ch = {};   % initialize array
ch_matched_ind = 1;   % initialize matched channel index as 1

for i = 1:size(ch_list, 1)   % for every channel name found
    % if contact numbers are given, '-' is included, check if name string
    % contains this '-' pattern, if so, check if name strings in ch_list
    % and eeg_ch agree, if so, add channel name to prom_ch array
    if contains(ch_list{i, 2}, '-')
        if ~isempty(find(strcmp(ch_list{i, 2}, eeg_ch(:, 1))))
            ind_found = find(strcmp(ch_list{i, 2}, eeg_ch(:, 1)));   % get index
            prom_ch{ch_matched_ind, 1} = ch_list{i, 1};   % event type
            prom_ch(ch_matched_ind, 2:3) = eeg_ch(ind_found, :);   % get all info.
            ch_matched_ind = ch_matched_ind + 1;
        end
    % if contact number is not given, search if name exists in eeg_ch, 
    % if so, assign all channels with name matched in midpt_all to array
    else 
        % format search pattern, use regexp to find rows with matched
        % pattern, then use find to get those rows with non-zero values
        search_pat = [ch_list{i, 2}, '[\d*]-', ch_list{i, 2}, '[\d*]']; 
        % if matched pattern is found, non-zero value is assigned, else,
        % empty value is assigned to the cell
        row_srch = regexp(eeg_ch(:, 1), search_pat);  
        % assign value of 0 to all empty cells 
        for row_num = 1:numel(row_srch)
            if isempty(row_srch{row_num})
                row_srch{row_num} = 0;
            end
        end
        % use find to get indices of non-zero entries
        ind_req = find(cell2mat(row_srch));
        % for each rows with non-zero entries found, get corresponding
        % channel name from eeg_ch array
        for item = 1:numel(ind_req)
            prom_ch{ch_matched_ind, 1} = ch_list{i, 1};   % event type
            prom_ch(ch_matched_ind, 2:3) = eeg_ch(ind_req(item), :);   % chnannel name
            ch_matched_ind = ch_matched_ind + 1;
        end
    end
end

% END PART (II): get names, event types, and eeg time series of prominent channels, 
% defined as those with sig. ied rates registered during fmri scan 
% %----------------------------------------------------------------------------
% PART (III): sort prominent channels by event types, average squared of all 
% time series of each type of prominent channels, and arrange event
% onsets (evti) by order of prom_ch

% get unique event types from col. 1 of prom_ch
uni_ev_type = unique(cell2mat(prom_ch(:, 1)));

% sort prom_ch by event type, format of prom_ch_type is given below, 
% in each cell (unique event type), 
% 1st col. = event type
% 2nd col. = channel name
% 3rd col. = eeg time series
prom_ch_type = {};   % initialize array
for item = 1:numel(uni_ev_type)   % for every unique event type
    prom_ch_type{item} = {};   % initialize cell for curr. event type
    for i = 1:size(prom_ch, 1)   % for each row in prom_ch
        if isequal(prom_ch{i, 1}, uni_ev_type(item))   % if event type agrees
            prom_ch_type{item} = [prom_ch_type{item}; prom_ch(i, :)];
        end
    end
end

% get average of sqaure (mean-square) of all time series of each type of prominent
% channels, initialize avg_sq_prom_ch_type to store the avg. of each type
% in each type
avg_sq_prom_ch_type = {};
for item = 1:numel(uni_ev_type)
    avg_sq_prom_ch_type{item} = mean([cell2mat(prom_ch_type{item}(:, end))].^2, 1);
end

% arrange evti by order of prom_ch_type, format of evti_prom_ch is given
% below, 
% in each cell (unique event type), 
% event onsets in seconds
evti_prom_ch = {};   % initialize array
for item = 1:numel(prom_ch_type)   % for each event type 
    cur_ev_type = prom_ch_type{item}{1, 1};   % get curr. event type
    for entry = 1:numel(event_def)   % for every cell in event_def
        if ismember(cur_ev_type, event_def{entry})   % if event type agrees
            evti_prom_ch{item} = evti{entry};   % add event onsets (evti) to current cell
        end
    end
end

% END PART (III): sort prominent channels by event types
% %----------------------------------------------------------------------------

% PART (IV): adjust part of eeg time series at each event onset for every
% prominent channel

% note: as the order at which event type represented by each cell 
% in evti_prom_ch and prom_ch_type is the same, proceed with the following
% to identify the amount of shifting required at each event onset of every
% channel 

% note: as marked event onset may differ from the actual IED peak
% registered during scan, shift time series at every marked event onset to
% match with the actual IED peak for every prominent channel

% initialize eeg_ts_shifted_ev_type, with format given below
% 1st layer, each cell denotes an unique event type
% 2nd layer, shifted portion of eeg time series of every prominent channel
eeg_ts_shifted_ev_type = {};

% as channels in prom_ch_type are grouped by unique event type, 
for ch_type_ind = 1:numel(evti_prom_ch)   % for every unique event type

eeg_ts_shifted_ev_type{ch_type_ind} = {};   % initialize cell of curr. event type

% for each event onset (each occurrence) of curr. even type
for ev_ind = 1:numel(evti_prom_ch{ch_type_ind})   
    
    % create time series (in seconds) using length of eeg time series of 
    % current event type
    t_axis = [1:1:length(prom_ch_type{ch_type_ind}{1, end})]/fs;

    % obtain array of event onsets (in seconds) of curr. event type
    evti_req = evti_prom_ch{ch_type_ind};
    spac = 1000;   % set spacing allot in data points

    % get index in t_axis representing the curr. event onset
    diff = abs(t_axis - evti_req(ev_ind));   % get diff. btw. t_axis and curr. event onset
    ind_req = find(diff == min(diff));   % use find to get index that as the min. diff. 
    
    % check if ind_req +/- spac is within limits of time series, if so, 
    % execute the following code, if not, continue to
    % next iteration (ev_ind)
    if ~(ind_req - spac >= 1 && ind_req + spac <= length(t_axis))
        sprintf('window out of range at %d', ev_ind)
        continue;   % continue to next iteration (event onset)
    end
        
    % get mean-square eeg time series of each type of prominent channels, 
    % then extract portion of time series required (only section covering 
    % curr. event onset is needed) using index found and
    % spacing specified above
    eeg_ts = avg_sq_prom_ch_type{ch_type_ind};   % get avg. of sq. of all eeg time series of curr. type
    eeg_ts_req = eeg_ts(ind_req-spac:ind_req+spac);   % extract portion of eeg time series 

    % eeg_ts denotes mean-square time series of each type of prominent channels,
    % get indices of local max., '1' is returned for entries that are local
    % max, while '0' is returned for other entries, for every prom. ch.
    TF = islocalmax(eeg_ts);   % indices of local max. for every prom. ch.

    % find closest non-zero index to ind_req in each direction (left and right),
    
    % initialize eeg_ts_shifted_req with each row storing shifted time series 
    % at curr. event onset for every prominent channel. 
    eeg_ts_shifted_req = {};
    
    % determine number of entries to the right of curr. marked event
    % onset required to get local max. 
    right_ind_req = 0;   % initialize var. 
    found = 0;
    while isequal(found, 0)   % as long as found is 0 
        % check if TF at curr. index is 1, if not, increment 
        % right_ind_req by 1
        if ~isequal(TF(ind_req + right_ind_req), 1)
            right_ind_req = right_ind_req + 1;
        else    % if found, termnate while loop
            found = 1;   
        end
    end
    
    % determine number of entries to the left of curr. marked event
    % onset to get local max, similar to above
    left_ind_req = 0;
    found = 0; 
    while isequal(found, 0)
        if ~isequal(TF(ind_req - left_ind_req), 1)
            left_ind_req = left_ind_req + 1;
        else
            found = 1;
        end
    end

    % select the min. btw. number of entries to the left and right of 
    % curr. marked event onset, set index of min. found as shift_mode
    [~, shift_mode] = min([left_ind_req, right_ind_req]);
    
    % using shift_mode from above, if local max is to the left of curr.
    % marked event onset (shift_mode = 1), shift curr. portion of time series to the
    % right so that actual IED peak align with marked onset.
    % Similarly, if local max is to the right of curr. marked
    % event onset (shift_mode = 2), shift curr. portion of time series to the left.
    if isequal(shift_mode, 1)   % local max is on the left, shift to right
        shift_req = left_ind_req;   % append zeros to the beginning of time series
        eeg_ts_shifted_req = [zeros(1, shift_req) eeg_ts_req(1:end - shift_req)];
    else
        shift_req = right_ind_req;   % local max is on the right, shift to left, append zeros to the end of time series
        eeg_ts_shifted_req = [eeg_ts_req(1 + shift_req:end) zeros(1, shift_req)];
    end
 
    % assign shifted (aligned) portion of time series of curr. event onset
    % to curr. cell (event type)
    eeg_ts_shifted_ev_type{ch_type_ind}{ev_ind} = eeg_ts_shifted_req;

end   % end for ev_ind = 1:numel(evti_prom_ch{ch_type_ind}) 

end   % end for ch_type_ind = 1:numel(evti_prom_ch)

% after shifting time series at every onset to align actual IED peak 
% with the marked event onset, arrange all shifted time series at 
% all onsets of each event type in singe layer 

% initialize array to store mean waveform (mean eeg time series) for each
% event type
mean_eeg_ts = {};

for ch_type_ind = 1:numel(eeg_ts_shifted_ev_type)
% mean_eeg_ts = 0;
eeg_ts_single_layer = [];
for i = 1:numel(eeg_ts_shifted_ev_type{ch_type_ind})
    curr_ts = eeg_ts_shifted_ev_type{ch_type_ind}{i};
    eeg_ts_single_layer = [eeg_ts_single_layer; curr_ts];
end

% then, get mean along dim. 1, aross all channels at every time point
mean_eeg_ts{ch_type_ind} = mean(eeg_ts_single_layer, 1);

end   % end for ch_type_ind = 1:numel(evti_prom_ch)

% END PART (IV): adjust part of eeg time series at each event onset for every
% prominent channel
% %----------------------------------------------------------------------------

% PART (V): get avg. waveform of all channels with respect to each event
% type

% note: after obtaining mean time series (waveform) (mean ts) of every type of
% prominent channels, compute cross-covariance (xcov) btw. ref (mean ts) and every
% channel, and locate the amount of shifting required to attain max.
% xcov

% initialize avg_waveform_ch with format given below, 
% 1st layer: each cell denotes info obtained for each event type
% 2nd layer: col. 1 = channel name
% 2nd layer: col. 2 = mean of eeg time series that max xcov is attained at
% every event onset
avg_waveform_ch = {};

for ch_type_ind = 1:numel(evti_prom_ch)   % for every event type (prom channel type)
    
    % construct time axis (t_axis) in unit of seconds for current event
    % type
    t_axis = [1:1:length(avg_sq_prom_ch_type{ch_type_ind})]/fs;
    
    % initialize cell and row number for curr. event type
    avg_waveform_ch{ch_type_ind} = {};
    row_num = 1;
    for ch_ind = 1:size(eeg_ch, 1)   % for every eeg channel provided

        % check name of curr. eeg channel from col. 1 to see if it is a
        % member of current type of prominent channels (col. 2), if not,
        % proceed 
        if ~ismember(eeg_ch{ch_ind, 1}, prom_ch_type{ch_type_ind}(:, 2))
            
            % initialize array to store shifted time series at every event
            % onset (occurrence)
            eeg_ts_targ_shifted = {};  
            
            % get list of event onsets in seconds from evti_prom_ch
            evti_req = evti_prom_ch{ch_type_ind};
            
            for ev_ind = 1:numel(evti_req)   % for each marked event onset
                
                % assign spacing in data points on either side of marked
                % event onset for searching 
                spac = 1000;   
                
                % get index in t_axis corresponding to the time of marked
                % event onset
                diff = abs(t_axis - evti_req(ev_ind));
                ind_req = find(diff == min(diff));
                
                % check if ind_req +/- spac is within limits of time series, if so, 
                % execute the following code, if not, continue to
                % next iteration (ev_ind)
                if ~(ind_req - spac >= 1 && ind_req + spac <= length(t_axis))
                    continue;   % continue to next iteration (event onset)
                end
    
                % use indices obtained to get portion of time series
                % surrounding the currnet marked event onset, 
                % *** sqaure each time point of the series, then store info to
                % eeg_ts_targ, 
                % then get cross-covariance (xcov) btw. curr. mean eeg ts
                % and target
                eeg_ts_targ = [eeg_ch{ch_ind, end}(ind_req-spac:ind_req+spac)].^2;
                [c_val, lags] = xcov(mean_eeg_ts{ch_type_ind}, eeg_ts_targ);
                
                % get max c_val and use its index to determine lag required
                c_val_ind_req = find(c_val == max(c_val));   % find index of max c_val
                lag_req = lags(c_val_ind_req);   % use index to get lag required
                
                % if lag_req < 0, it means target eeg time series needs to
                % be shifted to left to attain max xcov, and vice versa
                % so, if lag_req < 0, remove beginning of targ. time series
                % and append zeros to right of targ. time series. 
                % if lag_req > 0, append zeros at the beginning of targ.
                % time series and trim end of it
                if lag_req <= 0
                    eeg_ts_targ_shifted{ev_ind, 1} = [eeg_ts_targ(abs(lag_req) + 1:end) zeros(1, abs(lag_req))];
                else
                    eeg_ts_targ_shifted{ev_ind, 1} = [zeros(1, lag_req) eeg_ts_targ(1:end - lag_req)];
                end
                
                
            end   % end for ev_ind = 1:numel(evti_req)
            
            % after computing shifted time series at every event onset,
            % compute the mean across every shifted time series at every
            % time point (along dim. 1), assign info to col. 2 of
            % avg_wave_form
            avg_waveform_ch{ch_type_ind}{row_num, 1} = eeg_ch{ch_ind, 1};  % eeg channel name
            avg_waveform_ch{ch_type_ind}{row_num, 2} = mean(cell2mat(eeg_ts_targ_shifted), 1);   % mean time series 
            row_num = row_num + 1;
            
        end   % end if ~ismember(eeg_ch{ch_ind, 1}, prom_ch_type{ch_type_ind}(:, 2))
    end   % end for ch_ind = 1:size(eeg_ch, 1)   % for every eeg channel provided
end   % end for ch_type_ind = 1:numel(evti_prom_ch)   % for every event type (prom channel type)

% END PART (V): get avg. waveform of all channels with respect to each event
% type
% %----------------------------------------------------------------------------

% PART (VI): assign required var. to structure

opst.(fdnames_eeg{run_ind}).prom_ch = prom_ch;   % prominent channels
opst.(fdnames_eeg{run_ind}).prom_ch_type = prom_ch_type;   % prominent channels by event type
opst.(fdnames_eeg{run_ind}).mean_eeg_ts = mean_eeg_ts;   % mean time series of prominent channels by event type
opst.(fdnames_eeg{run_ind}).avg_waveform_ch = avg_waveform_ch;   % avg. waveform of target channels

% print message to terminal
sprintf('done working on run%d of sub%d', run_ind, str2num(subnum))

% END PART (XI): assign required var. to structure
% %----------------------------------------------------------------------------

% PART (XII): save output struct to path

if op_results == 1
    save(fullfile(fname_op, filename_op), 'opst');
end

% END PART (XII): save output struct to path
% %----------------------------------------------------------------------------


end   % end for run_ind = 1:numel(fdnames_eeg)   % for each session
