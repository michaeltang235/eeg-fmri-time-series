close all
clear all

subnum = '27';

% directname =  '/work/levan_lab/mtang/fmri_project/matrices';
directname =  'C:\Users\siumichael.tang\Downloads\fmri_project\';

% enter path where input struct. is stored at
% directname = [directname, filesep, 'matrices', filesep 'static_analysis_reho'];   % direct. of output matrix
directname_stat_reho = [directname, 'matrices\static_analysis_reho'];
filename_stat_reho = 'static_analysis_reho.mat';   % filename of output file

% directname_eeg = [directname, '/eeg_data'];
% directname_eeg = [directname, 'sub', subnum, '/eeg_data'];
directname_eeg = [directname, 'sub', subnum, '\eeg_data'];
filename_eeg = 'filt_eeg_ts.mat';

filename_eeg_unfilt = 'run*_gradient.mat';


% directname_event_onset = [directname, 'sub', subnum, '/matrices'];
directname_event_onset = [directname, 'sub', subnum, '\matrices'];
filename_event_onset = ['onset_times_events_sub', subnum, '.mat']; 

%----------------------------------------------------------------------------

% get paths of input files
[~, stat_reho_file_path, ~] = get_path(directname_stat_reho, filename_stat_reho);   % stat. analysis, reho
[~, eeg_file_path, ~] = get_path(directname_eeg, filename_eeg);   % eeg data
[~, event_onset_file_path, ~] = get_path(directname_event_onset, filename_event_onset);   % event onset times

[~, eeg_file_unfilt_path, ~] = get_path(directname_eeg, filename_eeg_unfilt);   % unfilt. eeg data

% set run index
run_ind = 1;

% % load structure
st_reho = load(stat_reho_file_path{:});
st_eeg = load(eeg_file_path{:});
st_event_onsets = load(event_onset_file_path{:});

st_eeg_unfilt = load(eeg_file_unfilt_path{run_ind});


fdnames_sub_list_reho = fieldnames(st_reho.grand);
sub_list_reho_ind = find(contains(fdnames_sub_list_reho, subnum));
fdnames_run_reho = fieldnames(st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}));
ch_name_matched_ied_reho = st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}).(fdnames_run_reho{run_ind}).ch_name_matched_ied_reho;

ch_name_matched_array = ch_name_matched_ied_reho;

sig_box_all = st_reho.grand.(fdnames_sub_list_reho{sub_list_reho_ind}).(fdnames_run_reho{run_ind}).sig_box_all;

fdnames_eeg = fieldnames(st_eeg.filt_eeg_st);
st_eeg_cur = st_eeg.filt_eeg_st.(fdnames_eeg{run_ind});

fdnames_run_events = fieldnames(st_event_onsets.terms);
event_def = st_event_onsets.terms.(fdnames_run_events{run_ind}).results.event_def;
evti = st_event_onsets.terms.(fdnames_run_events{run_ind}).results.evti;

window_size = 1;   % in seconds

% input threshold for corr. coef. 
thresh = 0.4;


% get fs in Hz of current run from eeg struct.
fs = st_eeg_cur.fs;

% %----------------------------------------------------------------------------
% % PART (I): SELECT CHANNELS IN EEG DATA SET THAT ARE ALSO FOUND IN
% % SIG_BOX_ALL (CHANNELS WITH FMRI FEATURES COMPUTED)
% 
% initialize array, eeg_ch, with format given below,
% 1st col. = channel name
% 2nd col. = eeg time series
eeg_ch = {}; 
row_num = 1;   % initialize row number
for i = 1:size(st_eeg_cur.filt_eeg, 1)   % for every channel listed in eeg data file
    ch_name = st_eeg_cur.filt_eeg{i, 1};   % assmeble channel name from col. 1

    % check if channel name is also in sig_box_all, which contains fmri
    % features for all channels available, if so, assign name and time series
    if ismember(ch_name, sig_box_all(:, 2))
        eeg_ch{row_num, 1} = ch_name;
        eeg_ch{row_num, 2} = st_eeg_cur.filt_eeg{i, end};   % assign full eeg time series
        row_num = row_num + 1;   % increment row numbr by 1 for next channel matched
    end
end

% END PART (I): SELECT CHANNELS IN EEG DATA SET THAT ARE ALSO FOUND IN
% SIG_BOX_ALL (CHANNELS WITH FMRI FEATURES COMPUTED)
%----------------------------------------------------------------------------
% PART (II): ASSIGN EEG TIME SERIES FOR EACH PROMINENT CHANNEL (THOSE WITH
% SIGNIFICANT IED RATES REGISTERED DURING FMRI SCAN)
 
% initialize array, prom_ch, with format given below
% 1st col. = event type
% 2nd col. = channel name
% 3rd col. = mean subsequence time series 
prom_ch = {};
row_num = 1;   % initialize row number
for i = 1:size(ch_name_matched_array, 1)   % for each channel in ch_name_matched_array
    
    % check if curr. channel name is included in eeg_ch array, 
    cur_ch_name = ch_name_matched_array{i, 2};   % get curr. channel name from col. 2
    [~, ind_found] = ismember(cur_ch_name, eeg_ch(:, 1));   % get index of matching row
    
    % if channel name is included in eeg_ch, add required info. to prom_ch array
    if ~isequal(ind_found, 0)
        prom_ch{row_num, 1} = ch_name_matched_array{i, 1};   % event type
        prom_ch{row_num, 2} = cur_ch_name;   % channel name
        prom_ch{row_num, 3} = eeg_ch{ind_found, end};   % eeg time series
        row_num = row_num + 1;   % incrememnt row number by 1 for next channel
    end
    
end

% END PART (II): ASSIGN EEG TIME SERIES FOR EACH PROMINENT CHANNEL (THOSE WITH
% SIGNIFICANT IED RATES REGISTERED DURING FMRI SCAN)
%----------------------------------------------------------------------------

% SORT PROM. CHANNELS BY EVENT TYPE

uni_ev_type = unique(cell2mat(prom_ch(:, 1)));

prom_ch_type = {};
for item = 1:numel(uni_ev_type)
    prom_ch_type{item} = {};
    for i = 1:size(prom_ch, 1)
        if isequal(prom_ch{i, 1}, uni_ev_type(item))
            prom_ch_type{item} = [prom_ch_type{item}; prom_ch(i, :)];
        end
    end
end

% %----------------------------------------------------------------------------

% ARRANGE EVTI BY ORDER OF PROM_CH
evti_prom_ch = {};
for item = 1:numel(prom_ch_type)
    cur_ev_type = prom_ch_type{item}{1, 1};
    for entry = 1:numel(event_def)
        if ismember(cur_ev_type, event_def{entry})
            evti_prom_ch{item} = evti{entry};
        end
    end
end

% %----------------------------------------------------------------------------

% UNFILTERED EEG TIME SERIES
eeg_ch_unfilt = {}; 
row_num = 1;   % initialize row number
for i = 1:size(st_eeg_unfilt.montage_names, 1)   % for every channel listed in eeg data file
    ch_name = [st_eeg_unfilt.montage_names{i, 1}, '-', st_eeg_unfilt.montage_names{i, 2}];   % assmeble channel name

    % check if channel name is also in sig_box_all, which contains fmri
    % features for all channels available, if so, assign name and time series
    if ismember(ch_name, sig_box_all(:, 2))
        eeg_ch_unfilt{row_num, 1} = ch_name;
        eeg_ch_unfilt{row_num, 2} = st_eeg_unfilt.data_gradient(:, i)';   % assign full eeg time series
        row_num = row_num + 1;   % increment row numbr by 1 for next channel matched
    end
end

prom_ch_unfilt = {};
row_num = 1;   % initialize row number
for i = 1:size(ch_name_matched_array, 1)   % for each channel in ch_name_matched_array
    
    % check if curr. channel name is included in eeg_ch array, 
    cur_ch_name = ch_name_matched_array{i, 2};   % get curr. channel name from col. 2
    [~, ind_found] = ismember(cur_ch_name, eeg_ch_unfilt(:, 1));   % get index of matching row
    
    % if channel name is included in eeg_ch, add required info. to prom_ch array
    if ~isequal(ind_found, 0)
        prom_ch_unfilt{row_num, 1} = ch_name_matched_array{i, 1};   % event type
        prom_ch_unfilt{row_num, 2} = cur_ch_name;   % channel name
        prom_ch_unfilt{row_num, 3} = eeg_ch_unfilt{ind_found, end};   % eeg time series
        row_num = row_num + 1;   % incrememnt row number by 1 for next channel
    end
    
end

prom_ch_type_unfilt = {};
for item = 1:numel(uni_ev_type)
    prom_ch_type_unfilt{item} = {};
    for i = 1:size(prom_ch_unfilt, 1)
        if isequal(prom_ch_unfilt{i, 1}, uni_ev_type(item))
            prom_ch_type_unfilt{item} = [prom_ch_type_unfilt{item}; prom_ch_unfilt(i, :)];
        end
    end
end

% %----------------------------------------------------------------------------

show_figure = 0;

eeg_ts_shifted_ev_type = {};

for ch_type_ind = 1:numel(evti_prom_ch)

eeg_ts_shifted_ev_type{ch_type_ind} = {};
    
for ev_ind = 1:numel(evti_prom_ch{ch_type_ind})
% for ev_ind = 4:4
% ev_ind = 5;

% if ev_ind == 5
%     show_figure = 1;
% else 
%     show_figure = 0;
% end

t_axis = [1:1:length(prom_ch_type{ch_type_ind}{1, end})]/fs;

evti_req = evti_prom_ch{ch_type_ind};
spac = 1000;

diff = abs(t_axis - evti_req(ev_ind));
ind_req = find(diff == min(diff));

eeg_ts = cell2mat(prom_ch_type{ch_type_ind}(:, end));
eeg_ts_req = eeg_ts(:, ind_req-spac:ind_req+spac);

TF = islocalmax(eeg_ts, 2);
TF_req = TF(:, ind_req-spac:ind_req+spac);

% find closest non-zero index to ind_req to determine amount of shifting
eeg_ts_shifted_req = {};
for item = 1:size(eeg_ts, 1)
    
right_ind_req = 0;
found = 0;
while isequal(found, 0)
    if ~isequal(TF(item, ind_req + right_ind_req), 1)
        right_ind_req = right_ind_req + 1;
    else 
        found = 1;
    end
end

left_ind_req = 0;
found = 0;
while isequal(found, 0)
    if ~isequal(TF(item, ind_req - left_ind_req), 1)
        left_ind_req = left_ind_req + 1;
    else 
        found = 1;
    end
end

[~, shift_mode] = min([left_ind_req, right_ind_req]);

shift_req = 0;
% eeg_ts_shifted = 0;
if isequal(shift_mode, 1)
    shift_req = left_ind_req;
    eeg_ts_shifted_req{item, 1} = [zeros(1, shift_req) eeg_ts_req(item, 1:end - shift_req)];
else
    shift_req = right_ind_req;
    eeg_ts_shifted_req{item, 1} = [eeg_ts_req(item, 1 + shift_req:end) zeros(1, shift_req)];
end

end

eeg_ts_shifted_ev_type{ch_type_ind}{ev_ind} = eeg_ts_shifted_req;

%-------------------------------------------------------------
% FOR PLOTTING:
eeg_ts_unfilt = cell2mat(prom_ch_type_unfilt{ch_type_ind}(:, end));
eeg_ts_req_unfilt = eeg_ts_unfilt(:, ind_req-spac:ind_req+spac);

t_axis_req = t_axis(ind_req-spac:ind_req+spac);

%-------------------------------------------------
% % VISUAL INSPECTION FOR EACH EVENT OCCURENCE:

% if show_figure == 1
% for item = 1:size(eeg_ts_req, 1)
% figure
% plot(t_axis_req, eeg_ts_req(item, :), 'b'); hold on
% plot(t_axis_req, eeg_ts_shifted_req{item}, 'r'); hold on
% xline(evti_req(ev_ind), 'LineWidth', 1.2); hold off
% legend('filtered', 'shifted');
% 
% end
% end
%----------------------


end   % end for ev_ind = 3:8

end   % end for ch_type_ind = 1:numel(evti_prom_ch)

mean_eeg_ts = {};

for ch_type_ind = 1:numel(eeg_ts_shifted_ev_type)
% mean_eeg_ts = 0;
eeg_ts_single_layer = [];
for i = 1:numel(eeg_ts_shifted_ev_type{ch_type_ind})
    curr_ts = cell2mat(eeg_ts_shifted_ev_type{ch_type_ind}{i});
    eeg_ts_single_layer = [eeg_ts_single_layer; curr_ts];
end

mean_eeg_ts{ch_type_ind} = mean(eeg_ts_single_layer, 1);

if show_figure == 1
figure
plot(1:1:length(mean_eeg_ts{ch_type_ind}), mean_eeg_ts{ch_type_ind}); hold on
xline(1001, 'LineWidth', 1, 'Color', 'k'); hold off
end

end   % end for ch_type_ind = 1:numel(evti_prom_ch)


%--------------------------------------------------------------------------
% AFTER OBTAINING AVERAGE WAVEFORM OF IED FOR EACH EVENT TYPE

avg_waveform_ch = {};

for ch_type_ind = 1:numel(evti_prom_ch)   
% for ch_type_ind = 1:1
    t_axis = [1:1:length(prom_ch_type{ch_type_ind}{1, end})]/fs;
    
    avg_waveform_ch{ch_type_ind} = {};
    row_num = 1;
    for ch_ind = 1:size(eeg_ch, 1)
%     for ch_ind = 8:8
        if ~ismember(eeg_ch{ch_ind, 1}, prom_ch_type{ch_type_ind}(:, 2))
            
            eeg_ts_targ_shifted = {};
            
            evti_req = evti_prom_ch{ch_type_ind};
            for ev_ind = 1:numel(evti_req)
                
                spac = 1000;

                diff = abs(t_axis - evti_req(ev_ind));
                ind_req = find(diff == min(diff));
                
                eeg_ts_targ = eeg_ch{ch_ind, end}(ind_req-spac:ind_req+spac);
                [c_val, lags] = xcov(mean_eeg_ts{ch_type_ind}, eeg_ts_targ);
                
                c_val_ind_req = find(c_val == max(c_val));
                lag_req = lags(c_val_ind_req);
                
                if lag_req <= 0
                    eeg_ts_targ_shifted{ev_ind, 1} = [eeg_ts_targ(abs(lag_req) + 1:end) zeros(1, abs(lag_req))];
                else
                    eeg_ts_targ_shifted{ev_ind, 1} = [zeros(1, lag_req) eeg_ts_targ(1:end - lag_req)];
                end
                
                
            end   % end for ev_ind = 1:numel(evti_req)
            
            avg_waveform_ch{ch_type_ind}{row_num, 1} = eeg_ch{ch_ind, 1};
            avg_waveform_ch{ch_type_ind}{row_num, 2} = mean(cell2mat(eeg_ts_targ_shifted), 1);
            row_num = row_num + 1;
            
        end
    end
end

corr_ch = avg_waveform_ch;
for ch_type_ind = 1:numel(evti_prom_ch)
    
    ref_ts = mean_eeg_ts{ch_type_ind};
    
    for i = 1:size(avg_waveform_ch{ch_type_ind}, 1)
        targ_ts = avg_waveform_ch{ch_type_ind}{i, end};
        [rho_ref_targ, pval_ref_targ] = corr(ref_ts', targ_ts', 'Type', 'Pearson');
        corr_ch{ch_type_ind}{i, 3} = rho_ref_targ;
        corr_ch{ch_type_ind}{i, 4} = pval_ref_targ;
    end
end


figure
plot(t_axis_req, mean_eeg_ts{2}, 'k'); hold on
plot(t_axis_req, avg_waveform_ch{2}{end, end}, 'b'); hold off
% % xline(evti_req(ev_ind), 'LineWidth', 1.2); hold off
legend('avg ied waveform', 'targ');

% figure
% plot(t_axis_req, mean_eeg_ts{3}, 'k'); hold on
% plot(t_axis_req, corr_ch{3}{11, 2}, 'b'); hold off
% % % xline(evti_req(ev_ind), 'LineWidth', 1.2); hold off
% legend('avg ied waveform', 'targ');

% lag_req = 25;
% eeg_ts_targ_shifted = [eeg_ts_targ(lag_req+1:end) zeros(1, lag_req)];

% figure
% plot(t_axis_req, mean_eeg_ts{1}, 'k'); hold on
% plot(t_axis_req, eeg_ts_targ, 'b'); hold on
% plot(t_axis_req, eeg_ts_targ_shifted, 'r'); hold off
% % xline(evti_req(ev_ind), 'LineWidth', 1.2); hold off
% legend('avg ied waveform', 'targ', 'targ shifted');


% for i = 1:size(eeg_ts_req, 1)
% figure
% plot(t_axis_req, eeg_ts_req(i, :), 'b'); hold on
% plot(t_axis_req, eeg_ts_req_unfilt(i, :), 'r'); hold on
% xline(evti_req(ev_ind), 'LineWidth', 1.2); hold off
% legend('filtered', 'unfiltered');
% 
% end

% figure 
% plot(t_axis_req, mean(eeg_ts_req, 1)); hold on
% xline(evti_req(ev_ind), 'LineWidth', 1.2); hold off
