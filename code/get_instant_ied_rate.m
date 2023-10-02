% This function gets instantaneous spike rate in every time segment in the
% unit of minutes for every event type

% FUNCTION INST_SPIKES = GET_INSTANT_IED_RATE(SWRA_IMG, SWRA_INFO,
% WINDOW_SIZE, EVENT_ONSETS)
% SWRA_IMG is processed fmri images
% TR is repetition time (seconds) of processed fmri images
% WINDOW_SIZE is window size in unit of seconds
% EVENT_ONSETS is event onsets structure
% INST_IED_ST is an output structure containing 
% - instantaneous spike rate [/min] in every time segment, 
% - time axis [/s] of input fmri images
% - list of initial times (ti) in seconds for every time segment created
% - list of final times (tf) in seconds for every time segment created

function inst_ied_st = get_instant_ied_rate(swra_img, tr, window_size, event_onsets)

% get dimensions of motion-filtered images
% ni = size(swra_img, 1);   % dim. in x-dir.
% nj = size(swra_img, 2);   % dim. in y-dir.
% nk = size(swra_img, 3);   % dim. in z-dir.
nt = size(swra_img, 4);   % dim. in t-dir.

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

% convert number of spikes in every segment to unit of minutes
inst_spikes = num_spikes(:, 1);   % initialize inst_spikes array
for item = 1:size(num_spikes, 1)   % for every event 
    % convert num. of spikes in every seg. to unit of minutes
    inst_spikes{item, 2} = cell2mat(num_spikes{item, 2})./window_size*60;
end

% store varaibles of interest in output structure, inst_ied_st
inst_ied_st = struct;
inst_ied_st.inst_ied = inst_spikes;   % instantaneous ied rate [/min]
inst_ied_st.taxis = taxis;   % time series of input fmri image, unit in seconds
inst_ied_st.ti_list = ti_list;   % initial time (s) for each segment constructed
inst_ied_st.tf_list = tf_list;   % final time (s) for each segment constructed


end   % end function num_spikes = get_instant_spike_rate(swra_img, swra_info, window_size, event_onsets)