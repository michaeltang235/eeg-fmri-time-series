% REMARKS:
% (i) 'prominent' channels = those that have IEDs registered during fmri scan,
% their names are provided in .txt (e.g. event_type_locations.txt) file 

% (ii) seed time series = mean time series of all prominent channels of each
% event type

% (iii) 'less prominent' channels = channels with their names not included in .txt
% file, but have corr. coeff. >= thresh with seed time series
% 
% (iv) 'other' channels = those that have corr. coeff < thresh

% note: only channels not in the list of prominent channels can be 
% 'less-prominent channels' and 'other channels', regardless of event type.


% This file contains functions of GET_LESS_PROM_CH and GET_MEAN_SUBSEQ_TIME_SERIES

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% FUNCTION OP_ST = GET_LESS_PROM_CH(PROM_CH_ARRAY, SIG_BOX_ALL,
% EEG_FILE_PATH, WINDOW_SIZE, THRESH)
% PROM_CH_ARRAY is an array of prominent channels (event type, names, ...)
% SIG_BOX_ALL is an array with fmri features calculated for every channel available
% EEG_FILE_PATH is path to eeg .mat file
% WINDOW_SIZE is length of window size in seconds require when computing
% mean subsequence time series 
% THRESH is threshold required to be classified as 'less-prominent' channels
% with respect to each type of prominent chanels, after computing corr.
% coeff. btw. time series of seed and every channel

function op_st = get_less_prom_ch(prom_ch_array, sig_box_all, eeg_file_path, window_size, thresh)

% initialize output struct, which contains the following
% op_st.prom_ch = prominent channels grouped by their event type
% op_st.less_prom_ch = less prominent channels relative to prominent channels
% op_st.other_ch = other channels (remaining channels from the input
% channel list)
op_st = struct;

% load eeg structure from path given
st_eeg = load(eeg_file_path);

% get fs in Hz from eeg struct.
fs = st_eeg.fs;

%----------------------------------------------------------------------------
% PART (I): GET EEG TIME SERIES FOR EACH CHANNEL IN SIG_BOX_ALL (THOSE WITH
% FMRI FEATURES CALCULATED)

% sometimes, different channels are found in eeg data and sig_box_all, only
% select channels that are found in sig_box_all (channels that have 
% fmri features available) for computation of mean time series, 
% and the following analysis. i.e. channels in eeg_ch is a subset of
% those in sig_box_all.
% Channels only available in sig_box_all, but
% not in eeg data means that they do not have IEDs recorded during scan, 
% hence, they are classified as 'other channels' with zero IED rate'.

% build eeg time series array, with format given below,
% col. 1 = channel name 
% col. 2 = eeg time series (full)
% col. 3 = mean of subsequence time series based on window size input
eeg_ch = {};   % initialize array
row_num = 1;   % initialize row number
for i = 1:size(st_eeg.montage_names, 1)   % for every channel listed in eeg data file
    ch_name = [st_eeg.montage_names{i, 1}, '-', st_eeg.montage_names{i, 2}];   % assmeble channel name

    % check if channel name is also in sig_box_all, which contains fmri
    % features for all channels available, if so, assign name, time series
    % and mean time series to channel
    if ismember(ch_name, sig_box_all(:, 2))
        eeg_ch{row_num, 1} = ch_name;
        eeg_ch{row_num, 2} = st_eeg.data_gradient(:, i);   % assign full eeg time series
        % call get_mean_subseq_time_series fuction, defined below, to get
        % mean of subsequence time series
        mean_subseq_ts = get_mean_subseq_time_series(st_eeg.data_gradient(:, i), window_size, fs);
        eeg_ch{row_num, 3} = mean_subseq_ts';
        row_num = row_num + 1;   % increment row numbr by 1 for next channel matched
    end
end

% END PART (I): GET EEG TIME SERIES FOR EACH CHANNEL IN SIG_BOX_ALL (THOSE WITH
% FMRI FEATURES CALCULATED)
%----------------------------------------------------------------------------
% PART (II): GET EEG TIME SERIES FOR EACH PROMINENT CHANNEL (THOSE HAVING
% IEDs REGISTERED DURING SCAN)

% get eeg time series for every channel in ch_name_matched_array

% initialize prom. channel array, with format given below
% col. 1 = event type
% col. 2 = prominent channel name
% col. 3 = mean subsequence eeg time series 
prom_ch = {};   % intiailize array
for i = 1:size(prom_ch_array, 1)   % for every row in ch_name_matched_array
    cur_ch_name = prom_ch_array{i, 2};   % get curr. ch. name
    % find row index in col. 1 of eeg_ch with matching name
    [~, ind_found] = ismember(cur_ch_name, eeg_ch(:, 1));   
    prom_ch{i, 1} = prom_ch_array{i, 1};   % assign event type
    prom_ch{i, 2} = prom_ch_array{i, 2};   % assign channel name
    prom_ch{i, 3} = eeg_ch{ind_found, end};   % assign eeg time series
end

% END PART (II): GET EEG TIME SERIES FOR EACH PROMINENT CHANNEL (THOSE HAVING
% IEDs REGISTERED DURING SCAN)
%----------------------------------------------------------------------------
% PART (III): GET MEAN TIME SERIES FOR EVERY SEED, DEFINED AS GROUP OF
% CHANNELS HAVING THE SAME EVENT TYPE

% get unique event types listed in col. 1 of prom. ch. array
uni_ev_type = unique(cell2mat(prom_ch(:, 1)));

% define seed as the average time series of each set of prominent channels 
% initialize array, with format given below
% 1st col. = event type
% 2nd col. = average of (mean subsequence time series) of all prominent
% channels of current event type

seed_ts = cell(numel(uni_ev_type), 2);   % initialize array
for ev_type = 1:numel(uni_ev_type)   % for every unique event type
    % find row indices in col. 1 of prom_ch with matching event types
    ind_found = find(cell2mat(prom_ch(:, 1)) == uni_ev_type(ev_type));
    seed_ts{ev_type, 1} = uni_ev_type(ev_type);   % assign event type
    seed_ts{ev_type, 2} = mean(cell2mat(prom_ch(ind_found, end)), 1);   % get mean eeg time series
end

% END PART (III): GET MEAN TIME SERIES FOR EVERY SEED, DEFINED AS GROUP OF
% CHANNELS HAVING THE SAME EVENT TYPE

%----------------------------------------------------------------------------
% PART (IV): GET CORR. COEFF. OF EEG TIME SERIES BTW. EACH SEED AND EVERY CHANNEL

% once time series of every seed (event type) is defined, compute
% correlation coefficient btw. every channel in eeg_ch array and each seed,
% except channels that are defined as prominent channels

% initialize array, with format given below
% 1st col. = event type
% 2nd col. = cell having the following format, 
% ---------- 1st col. = channel name
% ---------- 2nd col. = corr. of time series btw. channel and seed

rho_ts = cell(numel(uni_ev_type), 1);   
for ev_type = 1:numel(uni_ev_type)  % for each event type
    
    % get curr. event type from uni_ev_type array
    rho_ts{ev_type, 1} = uni_ev_type(ev_type);  

    ar_req = {};   % initialize array 
    row_num = 1;   % initialize row number
    for i = 1:size(eeg_ch, 1)   % for every channel in eeg_ch array
        
        cur_ch = eeg_ch{i, 1};   % get curr. ch. name in eeg_ch
        
        % check is curr. channel is a member of prominent channels, if so,
        % proceed with calculating corr. coeff. btw. channel and seed
        if ~ismember(cur_ch, prom_ch(:, 2))
            ar_req{row_num, 1} = eeg_ch{i, 1};
            
            % compute corr. coef. btw. curr. eeg ch. and seed
            ar_req{row_num, 2} = corr([seed_ts{ev_type, end}]', ...
                [eeg_ch{i, end}]', 'Type', 'Pearson');
            row_num = row_num + 1;   % increment row number by 1 for next eeg ch.
        end
    end
    % after looping through all channels in eeg_ch, assign ar_req to
    % 2nd col. of curr. row in rho_ts, listing corr. of curr. seed with
    % other channels
    rho_ts{ev_type, 2} = ar_req;
end

% END PART (IV): GET CORR. COEFF. OF EEG TIME SERIES BTW. EACH SEED AND EVERY CHANNEL
%----------------------------------------------------------------------------
% PART (V): GET LIST OF LESS-PROMINENT CHANNELS FOR EACH TYPE OF PROMINENT
% CHANNELS

% define less prominent channels as those having corr. coef. greater than
% or equal to threshold, relative to each type of prominent channels
% initialize less_prom_ch array with format below,
% 1st col. = event type
% 2nd col. = cell containing list of less prominent channels relative to
% curr. event type (prominent channels)

less_prom_ch = {};   % initialize array
ch_num = 1;   % initialize ch number (row number in array)
for ev_type = 1:size(rho_ts, 1)   % for every event type 
    
    ar_filt = {};   % initialize ar_filt array for curr. ev. type
    row_num = 1;   % initialize row number
    
    % search for channels with rho >= thresh from cell within 2nd col. of
    % rho_ts
    for i = 1:size(rho_ts{ev_type, 2}, 1)   % for every channel in cell (2nd. col. of rho_ts)
        % if rho (2nd col. of cell) of curr. ch. is >= thresh
        if rho_ts{ev_type, 2}{i, 2} >= thresh   
            ar_filt{row_num} = rho_ts{ev_type, 2}{i, 1};   % assign ch. name to ar_filt
            row_num = row_num + 1;   % increment row number by 1 for next qualifying ch.
        end
    end

    % assign names of less-prominent channels and their event type to
    % less_prom_ch array
    for item = 1:numel(ar_filt)   % for each channel in ar_filt
        less_prom_ch{ch_num, 1} = rho_ts{ev_type, 1};   % assign assoc. event type
        less_prom_ch{ch_num, 2} = ar_filt{item};   % assign channel name
        ch_num = ch_num + 1;   % increment ch. num. by 1 for next channel
    end
end

% END PART (V): GET LIST OF LESS-PROMINENT CHANNELS FOR EACH TYPE OF PROMINENT
% CHANNELS
%----------------------------------------------------------------------------
% PART (VI): GET LIST OF OTHER CHANNELS, DEFINED AS THOSE THAT HAVE CORR. COEFF.
% LESS THAN THRESHOLD 

% REMARKS:
% (i) PROM_CH and LESS_CH are subsets of SIG_BOX_ALL
% (ii) No duplicate channels are found in SIG_BOX_ALL 
% (iii) duplicate channels are found in PROM_CH, each channel can have more
% than one event type assigned

% OTHER CHANNELS are those with corr. coeff. less than threshdol, so they
% can be obtained by SIG_BOX_ALL - PROM_CH - LESS_PROM_CH

% get other channels relative to each type of prominent channels
deduct_prom_ch = setdiff(sig_box_all(:, 2), prom_ch(:, 2));
deduct_less_prom_ch = setdiff(deduct_prom_ch, less_prom_ch(:, 2));

% assign deduct_less_prom_ch to other_ch array
other_ch = deduct_less_prom_ch;

% END PART (VI): GET LIST OF OTHER CHANNELS, DEFINED AS THOSE THAT HAVE CORR. COEFF.
% LESS THAN THRESHOLD 
%----------------------------------------------------------------------------
% PART (VII): ASSIGN VARIABLES TO OUTPUT ARRAYS 

% update output structure, op_st
op_st.prom_ch = prom_ch;   % prominent channels
op_st.less_prom_ch = less_prom_ch;   % less-prominent channels relative to each type of prominent channels
op_st.other_ch = other_ch;   % other channels

end   % end function less_prom_ch = get_less_prom_ch(prom_ch_array, eeg_file_path, window_size)

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% This function takes input time series and output mean subsequence time
% series based on window size and fs
function mean_time_series = get_mean_subseq_time_series(input_time_series, window_size, fs)

% initialize subsequence based on input length (window_size * fs)
subseq = zeros(window_size*fs, 1);

% get number of time points input time series contains
nt = length(input_time_series);

% loop through each window (the i^th subsequence), append time series in
% window to subseq, then get element-wise mean afterwards
for i = 1:nt/(window_size*fs)
    
    % for current window, 
    ti_ind = 1 + (i-1)*window_size*fs;   % get index of initial time 
    tf_ind = window_size*fs + (i-1)*window_size*fs;   % get index of final time
    
    % extract segment of time series (subsequence) using indices found
    % above
    cur_time_seg = input_time_series(ti_ind:tf_ind);
    
    % concatenate (element-wise) curr. segment to subseq array
    subseq = subseq + cur_time_seg;
end

% get element-wise mean of subsequence
mean_time_series = subseq/(nt/(window_size*fs));

end