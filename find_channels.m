close all
clear all

subnum = '27';

% directname =  '/work/levan_lab/mtang/fmri_project/matrices';
directname =  'C:\Users\siumichael.tang\Downloads\fmri_project\';

% enter path where input struct. is stored at
% directname = [directname, filesep, 'matrices', filesep 'static_analysis_reho'];   % direct. of output matrix
directname_stat_reho = [directname, 'matrices\static_analysis_reho'];
filename_stat_reho = 'static_analysis_reho.mat';   % filename of output file

% directname_eeg = [directname, 'sub', subnum, '/eeg_data'];
directname_eeg = [directname, 'sub', subnum, '\eeg_data'];
filename_eeg = 'run*_gradient.mat';

% directname_event_onset = [directname, 'sub', subnum, '/matrices'];
% directname_event_onset = [directname, 'sub', subnum, '\matrices'];
% filename_event_onset = ['onset_times_events_sub', subnum, '.mat']; 

%----------------------------------------------------------------------------

% get paths of input files
[~, stat_reho_file_path, ~] = get_path(directname_stat_reho, filename_stat_reho);   % stat. analysis, reho
[~, eeg_file_path, ~] = get_path(directname_eeg, filename_eeg);   % eeg data
% [~, event_onset_file_path, ~] = get_path(directname_event_onset, filename_event_onset);   % event onset times

eeg_file_path = eeg_file_path{:};

% % load structure
st_reho = load(stat_reho_file_path{:});
% % st_events = load(event_onset_file_path{:});
% 
% % for run_ind = 1:numel(eeg_file_path)
% for run_ind = 1:1
% st_eeg = load(eeg_file_path{1});
% end
% 
fdnames_sub_list_reho = fieldnames(st_reho.grand);
sub_list_reho_ind = find(contains(fdnames_sub_list_reho, subnum));
fdnames_run_reho = fieldnames(st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}));
ch_name_matched_ied_reho = st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}).(fdnames_run_reho{1}).ch_name_matched_ied_reho;
% 
window_size = 2;   % in seconds
% fs = st_eeg.fs;

%----------------------
% ch_name_matched_array is a subset of sig_box_all, representing channels
% with non-zero IEDs registered during scan.

%----------------------
ch_name_matched_array = ch_name_matched_ied_reho;

sig_box_all = st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}).(fdnames_run_reho{1}).sig_box_all;

% less_prom_ch_test = get_less_prom_ch(ch_name_matched_ied_reho, eeg_file_path{:}, window_size, threshold);

% prominent channels = those that have IEDs registered during fmri scan,
% their names are provided in .txt (e.g. event_type_locations.txt) file 
% seed time series = mean time series of all prominent channels of each
% event type
% less prominent channels = channels with their names not included in .txt
% file, but have corr. coeff. >= thresh with seed time series
% other channels = those that have corr. coeff < thresh

% note: only channels not in the list of prominent channels can be 
% 'less-prominent channels' and 'other channels', regardless of event type.
% 
function less_prom_ch = get_less_prom_ch(ch_name_matched_array, eeg_file_path, window_size, threshold)

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
for i = 1:size(ch_name_matched_array, 1)   % for every row in ch_name_matched_array
    cur_ch_name = ch_name_matched_array{i, 2};   % get curr. ch. name
    % find row index in col. 1 of eeg_ch with matching name
    [~, ind_found] = ismember(cur_ch_name, eeg_ch(:, 1));   
    prom_ch{i, 1} = ch_name_matched_array{i, 1};   % assign event type
    prom_ch{i, 2} = ch_name_matched_array{i, 2};   % assign channel name
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

% input threshold for corr. coef. 
thresh = 0.4;

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

% end   % end function less_prom_ch = get_less_prom_ch(input_channel_array, eeg_file_path, window_size)

%-------------------
% TO BE ADDED TO STATIC_ANALYSIS_REHO.m BEOFRE CALLING THIS SCRIPT

% ind_swra_path = regexp(swra_file_path{run_ind}, 'run\d*', 'ignorecase');
% run_num = swra_file_path{run_ind}(ind_swra_path+3);
% 
% ind_eeg_path_req = 0;
% for item = 1:size(eeg_file_path, 1)
%     ind_eeg_path = regexp(eeg_file_path{item}, 'run\d*', 'ignorecase');
%     run_num_eeg_cur = eeg_file_path{item}(ind_eeg_path+3);
%     if run_num_eeg_cur == run_num
%         ind_eeg_path_req = item;
%     end
% end
%-------------------

% CHECK RUN NUMBER FROM .NII FILES BEFORE ASSIGNING THE PROPER EEG FILE AND
% EVENT ONSET MATRIX



% setdiff command

% rho_ts = cell(numel(prom_ch), 1);
% for ev_type = 1:numel(prom_ch)
%     rho_ts{ev_type, 1} = prom_ch{ev_type}{1, 1};
%     ar_req = {};
%     row_num = 1;
%     for i = 1:size(eeg_ch, 1)
%         if ~ismember(eeg_ch{i, 1}, prom_ch{ev_type}(:, 2))
%             rho_ts{ev_type}{row_num, 1} = prom_ch{ev_type}{1, 1};
%             rho_ts{ev_type}{row_num, 2} = eeg_ch{i, 1};
%             rho_ts{ev_type}{row_num, 3} = corr([seed_ts{ev_type, end}]', [eeg_ch{i, end}], 'Type', 'Pearson');
%             row_num = row_num + 1;
%         end
%     end
% end

% rho_mean_eeg_ts = mean_eeg_ts;
% for i = 1:size(data_gradient, 2)
%     rho_mean_eeg_ts(i, 1:2) = mean_eeg_ts(i, :);
%     [rho, pval] = corr(seed_eeg_ts', mean_eeg_ts{i, end}', 'Type', 'Pearson');
%     rho_mean_eeg_ts{i, 3} = corr(seed_eeg_ts', mean_eeg_ts{i, end}', 'Type', 'Pearson');
% end

% sig_box_all = st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}).(fdnames_run_reho{1}).sig_box_all;
% ch_list = sig_box_all(:, 2);


% data_gradient = st_eeg.data_gradient;   % get eeg time series of every channel
% fs = st_eeg.fs;   % sampling freq. in Hz
% montage_names = st_eeg.montage_names;   % get channel names
% 
% window_size = 2;   % in seconds
% 
% mean_eeg_ts = cell(size(data_gradient, 2), 2);
% for i = 1:size(data_gradient, 2)
%     mean_eeg_ts{i, 1} = strcat(montage_names{i, 1}, '-', montage_names{i, 2});
%     mean_eeg_ts{i, 2} = [get_mean_subseq_time_series(data_gradient(:, i), window_size, fs)]';
% end
% 
% seed_eeg_ts = mean(cell2mat(mean_eeg_ts(1:7, end)), 1);
% 
% rho_mean_eeg_ts = mean_eeg_ts;
% for i = 1:size(data_gradient, 2)
%     rho_mean_eeg_ts(i, 1:2) = mean_eeg_ts(i, :);
%     [rho, pval] = corr(seed_eeg_ts', mean_eeg_ts{i, end}', 'Type', 'Pearson');
%     rho_mean_eeg_ts{i, 3} = corr(seed_eeg_ts', mean_eeg_ts{i, end}', 'Type', 'Pearson');
% end
% 
% % % get list of
% ch_name_matched_ied_reho = st_reho.grand.sub27.run1_10min.ch_name_matched_ied_reho;
% 
% thresh  = 0.25;
% 
% ind_req = find(cell2mat(rho_mean_eeg_ts(:, end)) >= thresh);
% 
% less_prom_ch = rho_mean_eeg_ts(ind_req, :);
% 
% % t_axis = 1:1:length(seed_eeg_ts);
% % 
% % figure
% % plot(t_axis, seed_eeg_ts); hold on
% % plot(t_axis, rho_mean_eeg_ts{1, 2}); hold on
% % plot(t_axis, rho_mean_eeg_ts{4, 2}); hold on
% 
% % DONE METHOD 1: DECOMPOSE TIME SERIES INTO SEGMENT OF [WINDOW_SIZE], GET
% % AVG. FOR EACH CHANNEL, THEN GET MEAN OF ALL PROMINENT CHANNELS TO OBTAIN
% % TIME SERIES OF SEED, THEN COMPUTE CORRELATION OF OTHER CHANNELS WITH
% % SEED, SET THRESHOLD, AND IDENTIFY CHANNELS WITH CORR. HIGHER THAN OR
% % EQUAL TO THE THRESHOLD.
% 
% % METHOD 2: GET AVG. BY ONLY SELECTING SEGEMENTS CORRESPONDING TO
% % OCCURENCES OF SPIKES.
% 
% ev_onset_times = st_events.terms.run1.results.onsettimes{1}(:, 1);   % 1st col.
% 
% t_axis = [1:1:size(data_gradient, 1)]/fs;
% 
% spacing = (500/1000)*fs;
% 
% ti_ind = [];
% tf_ind = [];
% 
% for i = 1:numel(ev_onset_times)
%     diff = abs(t_axis - ev_onset_times(i));
%     ind_req = find(diff == min(diff), 1);
%     
%     ti_ind = [ti_ind, ind_req - spacing];
%     tf_ind = [tf_ind, ind_req + spacing];
% end
% 
% mean_eeg_ts2 = cell(size(data_gradient, 2), 2);
% for i = 1:size(data_gradient, 2)
%     mean_eeg_ts2{i, 1} = strcat(montage_names{i, 1}, '-', montage_names{i, 2});
%     mean_eeg_ts2{i, 2} = [get_mean2(data_gradient(:, i), ti_ind, tf_ind, spacing)]';
% end
% 
% seed_eeg_ts2 = mean(cell2mat(mean_eeg_ts2(1:7, end)), 1);
% 
% rho_mean_eeg_ts2 = mean_eeg_ts2;
% for i = 1:size(data_gradient, 2)
%     rho_mean_eeg_ts2(i, 1:2) = mean_eeg_ts2(i, :);
%     [rho, pval] = corr(seed_eeg_ts2', mean_eeg_ts2{i, end}', 'Type', 'Pearson');
%     rho_mean_eeg_ts2{i, 3} = corr(seed_eeg_ts2', mean_eeg_ts2{i, end}', 'Type', 'Pearson');
% end
% 
% ind_req2 = find(cell2mat(rho_mean_eeg_ts2(:, end)) >= thresh);
% 
% less_prom_ch2 = rho_mean_eeg_ts2(ind_req2, :);
% 
% % FINE TUNE BOTH METHODS, CHECK WHICH ONE IS BETTER
% 
% 
% 
% function mean_time_series2 = get_mean2(input_time_series, ti_ind, tf_ind, spacing)
% 
% sum_time_series2 = zeros(spacing*2+1, 1);
% 
% for i = 1:numel(ti_ind)-1
%     sum_time_series2 = [sum_time_series2 + input_time_series(ti_ind(i):tf_ind(i))];
% end
% 
% mean_time_series2 = sum_time_series2./(numel(ti_ind)-1);
% 
% end
% 
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
% 
% % plot(t_axis, ch_time_series);
