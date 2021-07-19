% This function intakes input signal array and regional mask, and outputs
% Kendall's W coefficient of concordance (ReHo). Tied ranks within a time
% series are considered in the calculations.

% some information about Kendall's W is given on the links listed below,
% https://en.wikipedia.org/wiki/Kendall%27s_W
% https://www.real-statistics.com/reliability/interrater-reliability/kendalls-w/
% http://adn.biol.umontreal.ca/~numericalecology/Reprints/Legendre_Coefficient_of_concordance_2010.pdf

% FUNCTION KENDALL_W = GET_KENDALL_W(INPUT_SIGNAL, REG_MASK)
% INPUT SIGNAL is a 4-d signal array
% REG_MASK is regional mask related to the input signal array, only voxels
% included in the mask have value of 1. Both INPUT SIGNAL and REG_MASK must
% have same lengths in first 3 dimesions, i.e., they are representing the
% same cluster of voxels.
% KENDALL_W is the output computed, Kendall's W coefficient

function kendall_w = kendall_w(input_signal, reg_mask)

% initialize output variable, kendall_w, as empty array
kendall_w = [];

% get dimensions of input signal array
ni = size(input_signal, 1);   % dim. in x-dir.
nj = size(input_signal, 2);   % dim. in y-dir.
nk = size(input_signal, 3);   % dim. in z-dir.
nt = size(input_signal, 4);   % dim. in t-dir.

% reshape input 3d array into 2d (single layer), only including
% voxels within the mask, making it easier for further processing
        
% initialize 2d array with first dim. equal to number of non-zero
% entries in regional mask, and second dim. equal to number of time
% points in every siganl
sig_2d = zeros(nnz(reg_mask), nt);   % initialize 2d array
row_num = 1;   % initialize row number 

% loop through each voxel in input array, check if current voxel is
% included in regional mask, if so, store signal of that voxel in 2d array
for i = 1:ni
    for j = 1:nj
        for k = 1:nk
            % if current voxel is included in mask, that voxel should have
            % a value of 1 in the mask
            if reg_mask(i, j, k) > 0   
                sig_2d(row_num, :) = input_signal(i, j, k, :);
                row_num = row_num + 1;   % increment row number by 1
            end
        end
    end
end

% get ranks of all time points in every array in sig_2d, using
% tiedrank. If any ranks in the signal are tied, their average rank
% is computed.
rank_sig = zeros(size(sig_2d));   % initialize array
for judge = 1:size(sig_2d, 1)   % for every time series 
    rank_sig(judge, :) = tiedrank(sig_2d(judge, :));   % get ranks, including tied ranks
end

% Kendall's W coefficient = 12S/(m^2(n^3-n)-mT), where S = sum of (R_i
% - mean R)^2, m is the number of judges, n is the number of
% time points, and T is the correction factor for tied ranks
    
% get R_i, mean R_i, and the sum of squared difference of them, 
% according to the formula of Kendall's W coefficient
    
% for m judges and nt time points, we have, 
R_i = sum(rank_sig, 1);   % sum of all ratings at each time point, R_i
mean_R_i = mean(R_i);   % mean of all R_i
ssd = 0;   % initialize sum of sqaured differences
for item = 1:length(R_i)   % get sum of squared diff. at every time point
    ssd = ssd + (R_i(item) - mean_R_i)^2;
end

% sort ranks in every time series (every row) in ascending order
rank_sorted = sort(rank_sig, 2);
   
% initialize correction factor T for all judges, if no tied ranks are
% found, T is 0
T = 0;   

% for every time series, find tied ranks (if any) and 
% compute their correction factor (T_j)
for judge = 1:size(rank_sorted, 1)   % for every judge (every time series)
    
    T_j = 0;   % initialize T factor for current judge
    rank_sc = rank_sorted(judge, :);   % get ranked signal, sorted and current
    search_ind = 1;   % initialize search index as 1
        
    % scan through entries of current time series
    while search_ind <=  length(rank_sc)
        
        % find all entries in time series with value equal to the current
        % time point (group), if number of entries identified is greater than
        % 1, duplicate ranks (tied ranks) are found
        % get number of tied ranks for current group, compute t_g, then
        % add it to T_j (T factor for current judge), 
        % adjust search index for next entry in the sorted rank array
        if numel(find(rank_sc == rank_sc(search_ind))) > 1
            num_tied_rank = numel(find(rank_sc == rank_sc(search_ind)));
            t_g = num_tied_rank^3 - num_tied_rank;   % for current group
            T_j = T_j + t_g;
            search_ind = search_ind + num_tied_rank;   % adjust search index 
        else   % if no duplicate rank is found, adjust search index 
               % as the next entry in array
               search_ind = search_ind + 1;
        end
    end
    T = T + T_j;   % add T_j of current judge to T (all judges)
    
end   % end for judge = 1:size(rank_sorted, 1) 

num_judges = size(rank_sig, 1);   % get num. of judges
num_time_pts = size(rank_sig, 2);   % get num. of time points

% get Kendall's W coefficient using the formula
kendall_w = 12*ssd/((num_judges^2)*(num_time_pts^3 - num_time_pts) - num_judges*T);


end   % end function kendall_w = get_kendall_w(input_signal, reg_mask)