% This function intakes (i) input signal and (ii) regional mask and outputs
% average signal across time in each segment.

% As regional mask gives info. about which voxels are included in our
% calculations, the avg. is computed by first summing all voxels at an instant 
% in time, then dividing the sum by the number of nonzero entries in the mask.
% In other words, this function creates a time series of avg. signal in
% every segment established.

% FUNCTION SIG_SEGMENT_AVG = SEGMENT_AVG(INPUT_SIGNAL, REG_MARK,
% NUM_NONZERO)
% INPUT_SIGANL is an array of input siganl of interest
% REG_MASK is an image of regional mask
% NUM_NONZERO is number of nonzero entries in the regional mark
% SIG_SEGMENT_AVG is an output array consisting of 2 layers, with signal 
% segments in the 1st layer, and the average signal of each segment in
% the 2nd layer.

function sig_segment_avg = segment_avg(input_signal, reg_mask)

% get number of nonzero entries in regional mask
num_nonzero = nnz(reg_mask);

% % get dimensions of input signal array
ni = size(input_signal, 1);   % dim. in x-dir.
nj = size(input_signal, 2);   % dim. in y-dir.
nk = size(input_signal, 3);   % dim. in z-dir.

% assume number of segments formed is the same in all voxels in input array
% get number of segments constructed at each voxel
num_seg = numel(input_signal{1});

% initialize array for storing avg. signal across time 
avg_signal = cell(1, num_seg);

% loop through each segment
for seg_ind = 1:num_seg   % for each segment
    nt = length(input_signal{1}{seg_ind});   % get length in time of current segment
    avg_signal{seg_ind} = zeros(1, nt);   % initialize array for storing avg. signal
    for t = 1:nt   % for each time point
        % initialize array for storing signals at all voxels at each time
        % point
        sig_cube = zeros(ni, nj, nk);   
        for i = 1:ni   % for each voxel in x-dir.
            for j = 1:nj   % for each voxel in y-dir.
                for k = 1:nk   % for each voxel in z-dir.
                    % assign the signal requested to sig_cube array
                    sig_cube(i, j, k) = input_signal{i, j, k}{seg_ind}(t);
                end
            end
        end
        % within current segment, at curren time point, get average signal 
        % by (i) summing the element-wise product btw. mask and signal cube 
        % (to remove signal outside region of interest), 
        % then (ii) dividing the sum by number of nonzero entries in the mask
        avg_signal{seg_ind}(t) = sum(sum(sum(reg_mask.*sig_cube)))/num_nonzero;
    end
end

% assign avg signal array to output variable after the loop
sig_segment_avg = avg_signal;

end   % end function sig_segment_avg = segment_avg(...)