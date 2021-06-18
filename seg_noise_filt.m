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

function [signal_seg, signal_noise_filt ] = seg_noise_filt(signal_input, tr, window_size)

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
            
            % create time axis spanning the entire time series of siganl interested
            taxis = (1:length(signal_req))*tr;
            
            % initialize indices for segmentation of time series at current
            % voxel
            ti_ind = 1;   % index of initial time
            seg_ind = 1;   % index of segment
            
            % get number of segments expected to form from time series
            num_seg = taxis(end)/window_size;  

            % initialize 2nd layers of arrays
            signal_seg{i, j, k} = cell(1, num_seg);   % segmented signal
            signal_noise_filt{i, j, k} = cell(1, num_seg);   % noise-filted seg. of signal
            
            % loop through each segment, find the index of final time (tf) in
            % taxis, then extract points within the segment using the indices
            for seg = 1:num_seg
                
                % find index of tf in taxis
                tf_ind = find(taxis == seg_ind*window_size); 
                
                % extract points within the segment (ti:tf), assign data to
                % current segment of signal_seg array
                signal_seg{i, j, k}{seg_ind} = signal_req(ti_ind:tf_ind);
                
                % apply bandpass filter (filtfilt) to remove noise from
                % current segment of signal
                signal_noise_filt{i, j, k}{seg_ind} = filtfilt(b, a, signal_req(ti_ind:tf_ind));
            
                % increment segment index by 1 for the next one
                seg_ind = seg_ind + 1;
            
                % assign index of tf to index of ti for the next segment
                ti_ind = tf_ind;
                
            end   % end for seg = 1:num_seg
        
        end   % end for k = 1:nk
    end   % end for j = 1:nj
end   % end for i = 1:ni
            


end   % end function signal_filt = seg_nois_filt