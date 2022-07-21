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

% load structure
st_reho = load(stat_reho_file_path{:});
% st_events = load(event_onset_file_path{:});

% for run_ind = 1:numel(eeg_file_path)
for run_ind = 1:1
st_eeg = load(eeg_file_path{run_ind});
end

fdnames_sub_list_reho = fieldnames(st_reho.grand);
sub_list_reho_ind = find(contains(fdnames_sub_list_reho, subnum));
fdnames_run_reho = fieldnames(st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}));
ch_name_matched_ied_reho = st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}).(fdnames_run_reho{1}).ch_name_matched_ied_reho;


uni_ev_type = unique(cell2mat(ch_name_matched_ied_reho(:, 1)));
prom_ch_list = {};
for item = 1:numel(uni_ev_type)
    ind_req = find(cell2mat(ch_name_matched_ied_reho(:, 1)) == uni_ev_type(item));
    prom_ch_list{item} = ch_name_matched_ied_reho(ind_req, :);
end

window_size = 2;   % in seconds
fs = st_eeg.fs;

eeg_ch = cell(size(st_eeg.montage_names, 1), 3);
for i = 1:size(st_eeg.montage_names, 1)
    eeg_ch{i, 1} = [st_eeg.montage_names{i, 1}, '-', st_eeg.montage_names{i, 2}];
    eeg_ch{i, 2} = st_eeg.data_gradient(:, i);
    eeg_ch{i, 3} = get_mean_subseq_time_series(st_eeg.data_gradient(:, i), window_size, fs);
end

prom_ch = cell(numel(prom_ch_list), 1);
for ev_type = 1:numel(prom_ch_list)
    prom_ch{ev_type} = prom_ch_list{ev_type}(:, 1:2);
    for row = 1:size(prom_ch_list{ev_type}, 1)
        prom_ch_cur = prom_ch_list{ev_type}{row, 2};
        ch_found = 0;
        for i = 1:size(eeg_ch, 1)
            if isequal(prom_ch_cur, eeg_ch{i, 1})
                prom_ch{ev_type}{row, 3} = eeg_ch{i, end}';
                ch_found = 1;
            end
            if ch_found == 0
                prom_ch{ev_type}{row, 3} = [];
            end
        end
    end
end

seed_ts = cell(numel(prom_ch), 1);
for ev_type = 1:numel(prom_ch)
    seed_ts{ev_type, 1} = prom_ch{ev_type}{1, 1};
    seed_ts{ev_type, 2} = mean(cell2mat(prom_ch{ev_type}(:, end)), 1);
end

rho_ts = cell(numel(prom_ch), 1);
for ev_type = 1:numel(prom_ch)
    rho_ts{ev_type, 1} = prom_ch{ev_type}{1, 1};
    ar_req = {};
    row_num = 1;
    for i = 1:size(eeg_ch, 1)
        if ~ismember(eeg_ch{i, 1}, prom_ch{ev_type}(:, 2)) && ...
                ~ismember(eeg_ch{i, 1}, ch_name_matched_ied_reho(:, 2))
            ar_req{row_num, 1} = eeg_ch{i, 1};
            ar_req{row_num, 2} = corr([seed_ts{ev_type, end}]', [eeg_ch{i, end}], 'Type', 'Pearson');
            row_num = row_num + 1;
        end
    end
    rho_ts{ev_type, 2} = ar_req;
end

thresh = 0.4;

less_prom_ch = cell(size(rho_ts, 1), 2);
for ev_type = 1:size(rho_ts, 1)
    less_prom_ch{ev_type, 1} = rho_ts{ev_type, 1};
    ar_filt = {};
    row_num = 1;
    for i = 1:size(rho_ts{ev_type, 2}, 1)
        if rho_ts{ev_type, 2}{i, 2} >= thresh
            ar_filt{row_num} = rho_ts{ev_type, 2}{i, 1};
            row_num = row_num + 1;
        end
    end
    less_prom_ch{ev_type, 2} = ar_filt';
end

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
