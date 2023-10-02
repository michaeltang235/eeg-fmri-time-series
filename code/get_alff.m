% FUNCTION NORM_ALFF = GET_ALFF(FMRI_IMG, EXP_MASK, TR, WINDOW_SIZE)
% FMRI_IMG is processed fmri images (4d array)
% EXP_MASK is explicit mask 
% TR is repetition time in seconds 
% WINDOW_SIZE is window size in seconds
% NORM_ALFF is normalized ALFF map (voxels NOT included in explicit mask is
% assigned with empty value)

% REMARKS: if window_size is empty, this function calculates normalized
% alff for the entire time series, otherwise, time series is decomposed
% into segments of length specified by window_size, and normalized alff is
% calculated for every segment.

function norm_alff = get_alff(fmri_img, exp_mask, tr, window_size)

% initialize output variable, normalized alff
norm_alff = [];

% get sampling frequency, in Hz
fs = 1/tr;

% freq. range interested in Hz, witth freq. soln. within the order of
% magnitude of default (1/256*Ts), Ts = sampling period in seconds
f_range = 0.01:0.001:0.08;   % freq. range interested in Hz, 

% get dimensions of fmri image
ni = size(fmri_img, 1);   % dim. in x
nj = size(fmri_img, 2);   % dim. in y
nk = size(fmri_img, 3);   % dim. in z
nt = size(fmri_img, 4);   % dim. in time

% check if length of time series is long enough, 256 is used as it's chosen
% as the width of window in welch's method of power spectral density (psd) 
% estimate. If not, return control to line calling this function
if nt < 256
    disp('time series is not long enough for pwelch method')
    return
end

%------------------------------------------------------------------
% (I): get normalized alff for entire time series (compute_mode == 1)

% only execute lines in this section if window size is empty 
if isempty(window_size)
    
    % *** only execute the following if nt >= 256, time series is long
    % enough for welch's method ***

    % set window size (num. of data points) for welch's method
    % fixed for all subjects regardless of length of time series, so same
    % number of sample are used for constructing peridograms
    win_len = 256;   % default value,

    % initialize arrays of interest
    alff_3d = zeros(ni, nj, nk);   % map of alff
    mean_alff = [];   % global mean alff
    norm_alff = zeros(ni, nj, nk);   % normalized alff

    % scan through each voxel and calculate alff, which is defined as mean
    % of square root of the power spectrum density over the low freq. band (0.01 -
    % 0.08) Hz
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                % if current voxel at explicit mask is not zero
                if exp_mask(i, j, k) ~= 0               
                    % extract time series at current voxel, call pwelch to get
                    % psd estimate, assign result to alff_3d array
                    time_series = squeeze(fmri_img(i, j, k, :));
                    pxx = pwelch(time_series, win_len, [], f_range, fs);
                    alff_3d(i, j, k) = mean(sqrt(pxx));     
                end
            end
        end
    end

    % get global average alff
    mean_alff = sum(sum(sum(alff_3d)))/(nnz(exp_mask));

    % get normalized alff, which is defined as alff/its global mean 
    norm_alff = alff_3d;   % initialize normalized alff as the original 3d map
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                norm_alff(i, j, k) = alff_3d(i, j, k)/mean_alff;
            end
        end
    end
    
end   % end if isempty(windows_size)

% END (I): get normalized alff for entire time series (compute_mode == 1)
%------------------------------------------------------------------
% (II): get normalized alff for every time segment in time series (compute_mode == 2)

% only execute lines in this section if window_size is not empty
if ~isempty(window_size)
    
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

    % initialize map of alff
    alff_3d = cell(ni, nj, nk);

    % scan through each voxel and calculate alff, which is defined as mean
    % of square root of the power spectrum density over the low freq. band (0.01 -
    % 0.08) Hz
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                for seg_ind = 1:length(ti_list)   % for every segment
                    % if current voxel at explicit mask is not zero
                    if exp_mask(i, j, k) ~= 0  
                        % find indices of ti and tf of current segment in taxis
                        ti_ind = find(taxis == ti_list(seg_ind));   % index of ti of cur. seg.
                        tf_ind = find(taxis == tf_list(seg_ind));   % index of tf of cur. seg.
                    
                        % obtain segment of time series required,
                        % then use pwelch to get psd estimate
                        time_series = squeeze(fmri_img(i, j, k, ti_ind:tf_ind));
                        pxx = pwelch(time_series, [], [], f_range, fs);
                    
                        % get alff, which is mean of sq. root of psd estimate
                        alff_3d{i,j,k}{seg_ind} = mean(sqrt(pxx));
                    else
                        % assign empty matrix if vox. is not in explicit mask
                        alff_3d{i, j, k}{seg_ind} = [];
                    end
                    
                end
            end
        end
    end   % end for i = 1:ni

    % get global average alff for every segment 
    mean_alff = zeros(1, num_seg);   % initialize array for each seg.
    for seg_ind = 1:length(ti_list)   % for every segment
        seg_sum = 0;   % initialize sum as 0
        for i = 1:ni   % for each voxel in x-dir.
            for j = 1:nj
                for k = 1:nk
                    % check if curr. voxel is included in epxlicit mask, 
                    % if so, add value of curr. seg. at voxel to sum
                    if exp_mask(i, j, k) ~= 0  
                        seg_sum = seg_sum + alff_3d{i, j, k}{seg_ind};
                    end
                end
            end
        end
        % after looping through all voxels, get mean of curr. seg. by dividing
        % sum by number of non-zero entries in explicit mask
        mean_alff(seg_ind) = seg_sum/(nnz(exp_mask));
    end
            
    % get normalized alff for every segment, which is defined as alff/its global mean 
    norm_alff = cell(ni, nj, nk);   % initialize array 
    for seg_ind = 1:length(ti_list)   % for each segment, then for each voxel
        for i = 1:ni
            for j = 1:nj
                for k = 1:nk
                    % check if curr. voxel is included in explicit mask, 
                    % if so, calculate normalized alff for curr. seg. and voxel
                    if exp_mask(i, j, k) ~= 0
                        norm_alff{i, j, k}{seg_ind} = alff_3d{i, j, k}{seg_ind}/mean_alff(seg_ind);
                    else
                        % if not, assign empty value to curr. seg.
                        norm_alff{i, j, k}{seg_ind} = [];
                    end
                end
            end
        end
    end
    
end   % end if ~isempty(windows_size)

% END (II): get normalized alff for every time segment in time series (compute_mode == 2)
%------------------------------------------------------------------

% display message on terminal
sprintf('normalized alff map calculated')

end   % end alff = get_alff(fmri_img, tr, window_size)
