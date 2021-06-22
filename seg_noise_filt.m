% This function takes an array of input signal and returns two outputs,
% (i) signal_seg, signal decomposed into segments, 
% (ii) signal_noise_filt, segments of signal with noise filtered

% FUNCTION [SIGNAL_SEG, SIGNAL_NOISE_FILT] = SEG_NOISE_FILT(SIGNAL_INPUT,
% TR, WINDOW_SIZE)
% SIGNAL_INPUT is the input signal that needs to be segmented and
% noise-filtered
% TR is repetition time (in seconds) of input signal
% WINDOW_SIZE is the length of each segment in seconds
% SIGNAL_SEG is an output array consisting of 2 layers, with 1st layer
% displaying time series in the region of interest, as defined in input
% siganl, and 2nd layer listing the segments of signal formed at each voxel
% within the region of interest
% SIGNAL_NOISE_FILT is similar to SIGNAL_SEG, but time series in 2nd layer
% are noise-removed using Chebyshev type 1 filter

function [signal_seg, signal_noise_filt] = seg_noise_filt(signal_input, tr, window_size)

% window_size = 60;   % in seconds
% tr = 1.5;   % repetition time, in seconds
fs = 1/tr; % sampling frequency, in Hz

% use Chebyshev Type 1 bandpass filter to allow freq. compo. btw. 
% [1/window size and 0.1] Hz to pass through, use cheby1 to get 
% transfer function coefficients (normalized freqs. are used)
[b,a] = cheby1(2,0.5,[1/window_size 0.1]/(fs/2));

% get dimensions of input signal arrray
ni = size(signal_input, 1);   % dim. in x-dir.
nj = size(signal_input, 2);   % dim. in y-dir.
nk = size(signal_input, 3);   % dim. in z-dir.
nt = size(signal_input, 4);   % dim. in t-dir.

% initialize output arrays 
signal_seg = cell(ni, nj, nk);   % segmented signal
signal_noise_filt = cell(ni, nj, nk);   % noise-filtered segments of signal

% loop through each time series in input signal array
for i = 1:ni   % for each voxel in x-dir.
    for j = 1:nj   % for each voxel in y-dir.
        for k = 1:nk   % for each voxel in z-dir.
            
            % select the time series at current voxel, then use 
            % squeeze to remove dim. of length 1
            signal_req = squeeze(signal_input(i, j, k, :));   
            
            % create time axis (slice-time corrected) spanning the entire 
            % time series of siganl interested
            % 1st image is shifted to t=0 (s), while the last image is
            % shifted to t_end-1 (s)
            taxis = (0:length(signal_req) - 1)*tr;
            
            % check if time axis is long enough to ouput meaningful
            % statistics, if not, return function
            % for window size of 60 s, time series has to be at least 10
            % mins long
            if taxis(end) < 10*60
                return
            end
            
            % create arrays indicating time points of each segment
            ti_list = 0;   % list of ti (initial time) 
            tf_list = window_size;   % list of tf (final time)
            while tf_list(end) + window_size/2 < taxis(end)   % while tf < taxis(end)
                ti_list = [ti_list ti_list(end) + window_size/2];   % update lists
                tf_list = [tf_list tf_list(end) + window_size/2];
            end
            
            % get number of segments expected to form from time series, 
            num_seg = numel(tf_list);  

            % initialize 2nd layers of arrays
            signal_seg{i, j, k} = cell(1, num_seg);   % segmented signal
            signal_noise_filt{i, j, k} = cell(1, num_seg);   % noise-filted seg. of signal
            
            % loop through each segment, find the indices of ti (intiial time)
            % tf (final time) in taxis, then extract points within 
            % the segment using the indices
            for seg_ind = 1:num_seg
                
                % find indices of ti and tf of current segment in taxis
                ti_ind = find(taxis == ti_list(seg_ind));   % index of ti of cur. seg.
                tf_ind = find(taxis == tf_list(seg_ind));   % index of tf of cur. seg.
                
                % extract points within the segment (ti:tf), assign data to
                % current segment of of current voxel in signal_seg array
                signal_seg{i, j, k}{seg_ind} = signal_req(ti_ind:tf_ind);
                
                % apply bandpass filter (filtfilt) to remove noise from
                % current segment of signal
                signal_noise_filt{i, j, k}{seg_ind} = filtfilt(b, a, signal_req(ti_ind:tf_ind));
                
            end   % end for seg = 1:num_seg
        
        end   % end for k = 1:nk
    end   % end for j = 1:nj
end   % end for i = 1:ni
            


end   % end function signal_filt = seg_nois_filt