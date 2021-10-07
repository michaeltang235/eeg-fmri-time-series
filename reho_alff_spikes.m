clear all 
close all

tic 
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '14';

% enter path to directory where all input files are located
directname = ['/work/levan_lab/mtang/fmri_project/', 'sub', subnum];
% directname = ['C:\Users\siumichael.tang\Downloads\fmri_project\', 'sub', subnum];

% format filenames of processed fmri images, explicit mask, electrode,
% motion parameters
filename_swraimg = ['swra*.nii'];   % processed func. images
filename_expmask = 'wEPI_bet_mask.nii';   % explicit mask
filename_elect =  [subnum, '_*Koordinaten*.xlsx'];   % file containing mni coord. of all electrode pairs
filename_spikes = ['subject', subnum, '_rates.txt'];   % file of spike rates
filename_motion = ['rp_*.txt'];   % motion parameters generated by SPM in the realignment step

% enter path where ouput struct. is stored at
fname_op = [directname, filesep, 'matrices' filesep 'reho_alff_spikes'];   % direct. of output matrix
filename_op = 'reho_alff_spikes.mat';   % filename of output file

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% locate paths of func. images and their details 

% use get_path function to obtain full path of input files
[~, swra_file_path, ~] = get_path(directname, filename_swraimg);   % swra func. images
[~, expmask_file_path, ~] = get_path(directname, filename_expmask);   % normalized mask images
[~, elect_file_path, ~] = get_path(directname, filename_elect);   % file containing info. about electrodes coordinates
[~, spikes_file_path] = get_path(directname, filename_spikes);   % spike rates
[~, motion_file_path, ~] = get_path(directname, filename_motion);   % motion parameters

% createterms structure to store output data calcu. by function defined
% below
terms = struct;

% read input files and use get_rho_spikes function to get required data
for run_ind = 1:numel(swra_file_path)   % for each row in path of processed func. images
    
% run_ind = 1;
% read input files
swra_img = niftiread(swra_file_path{run_ind});   % swra func. images
swra_info = niftiinfo(swra_file_path{run_ind});   % info about this particular swra func. images
expmask = single(niftiread(expmask_file_path{:}));   % convert image data type to single
electar = readcell(elect_file_path{:});   % cell array containing info. about electrodes
spikesar = readcell(spikes_file_path{:});   % cell array containing info. about spike rates

% assemble input array of function, get_rho_spikes 
input_array = struct;   % initialize struct.
input_array.swra_img = swra_img;   % processed func. images
input_array.swra_info = swra_info;   % headers of processed func. images
input_array.expmask = expmask;   % explicit mask
input_array.electar = electar;   % array of electrode coordinates
input_array.spikesar = spikesar;   % array of spike rates
input_array.motion_file_path = motion_file_path{run_ind};   % path of motion parameter file

% get indices of current run number embedded in processed image filename
sind = regexp(lower(swra_file_path{run_ind}), 'run');   % start index
eind = regexp(lower(swra_file_path{run_ind}), '.nii');   % end index
sess_num = swra_file_path{run_ind}(sind+3:eind-1);   % get indices of run number

% create fieldname representing current session number
fieldname = sprintf('run%s', sess_num);

% store output calcu. under current fieldname of terms
terms.(fieldname) = get_reho_alff_spikes(input_array);

% execute lines below if current field is not empty
if ~isempty(terms.(fieldname))
    
% access tables of current session created by function get_rho_spikes
table_reho_alff = terms.(fieldname).table_reho_alff;   % table of reho and alff for each channel

% format filenames and full paths of tables
filename_table = ['sub', subnum, '_table_reho_alff_run_', sess_num, '.csv'];
table_path = fullfile(fname_op, filename_table);   % path of table of clin. det. channels

%---------------------------------------------------------------------------
% output table and structure created in current session
if op_results == 1
    writetable(table_reho_alff, table_path);   % table of clin. det. ch.
    save(fullfile(fname_op, filename_op), 'terms');   % struct.
end
%---------------------------------------------------------------------------

end   % end if ~isempty(terms.(fieldname))

end   % end for run_ind = 1:numel(swra_file_path)

toc
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% define function get_reho_alff_spikes requiring an input array, which
% has the fields defined below,
% SWRA_IMG is the data of processed image file
% SWRA_INFO is the info (headers) associated with the processed image file
% EXPMASK is the explicit mask image
% ELECTAR is the array containing info. about all electrodes used,
% CLIN_REF_CH_FILE_PATH is the path of the file of clinically determined
% channels,
% SUBNUM is the subject number in integer form
% MOTION_FILE_PATH is the path of motion parameter file
% and outputs,
% OPSTRUCT, a structure that contains info. about coordinates of
% electrdoes, signals within 3X3X3 boxes centered at midpoint between
% electrode pairs, corr. coeff. btw. ref. and target channels, etc 
% see end of function for details of variables included in this struct.

function opstruct = get_reho_alff_spikes(input_array)

% %---------------------------------------------------------------------------
% % Prelim: access variables defined in input array
swra_img = input_array.swra_img;   % processed func. images
swra_info = input_array.swra_info;   % headers of processed func. images
expmask = input_array.expmask;   % explicit mask
electar = input_array.electar;   % array of electrode coordinates
spikesar = input_array.spikesar;   % array of spike rates
motion_file_path = input_array.motion_file_path;   % path of motion parameter file

% initialize output struct. as empty array
opstruct = [];

% check if length of time series is long enough, 256 is used as it's chosen
% as the width of window in welch's method of power spectral density (psd) 
% estimate. If so, change value of com_alff to 1 so that the script will
% calculate ALFF, the relevant corr. coeff., and p-values.
com_alff = 0;
if size(swra_img, 4) >= 256
    com_alff = 1;
end

% % end Prelim: access variables defined in input array

%---------------------------------------------------------------------------

% Part (I): remove motion artifacts from processed image (entire brain)

% convert explicit mask and processed image to double, 
% apply element-wise multiplciation
exp_swra_img = double(expmask).*double(swra_img);
                
% call filter_motion.m to remove motion artifacts from image
motion_filtered_img = filter_motion(motion_file_path, exp_swra_img);

% END Part (I): remove motion artifacts from processed images
%---------------------------------------------------------------------------
% Part (II): calculate amplitude of low-frequency fluctuations (alff)

% get basic info. about the processed image
tr = swra_info.PixelDimensions(4);   % repetition time, in seconds
fs = 1/tr; % sampling frequency, in Hz

% set window size (num. of data points) for welch's method
% fixed for all subjects regardless of length of time series, so same
% number of sample are used for constructing peridograms
win_len = 256;   % default value,

% freq. range interested in Hz, witth freq. soln. within the order of
% magnitude of default (1/256*Ts), Ts = sampling period in seconds
f_range = 0.01:0.001:0.08;   % freq. range interested in Hz, 

% get dimensions of motion-filtered image
ni = size(motion_filtered_img, 1);   % dim. in x
nj = size(motion_filtered_img, 2);   % dim. in y
nk = size(motion_filtered_img, 3);   % dim. in z

% initialize arrays of interest
alff_3d = [];   % map of alff
mean_alff = [];   % global mean alff
norm_alff = [];   % normalized alff

if com_alff == 1
    
% initialize array for storing alff calculated at each voxel
alff_3d = zeros(ni, nj, nk);

% scan through each voxel and calculate alff, which is defined as mean
% of square root of the power spectrum density over the low freq. band (0.01 -
% 0.08) Hz
for i = 1:ni
    for j = 1:nj
        for k = 1:nk
            % if current voxel at explicit mask is not zero
            if expmask(i, j, k) ~= 0               
                % extract time series at current voxel, call pwelch to get
                % psd estimate, assign result to alff_3d array
                time_series = squeeze(motion_filtered_img(i, j, k, :));
                pxx = pwelch(time_series, win_len, [], f_range, fs);
                alff_3d(i, j, k) = mean(sqrt(pxx));     
            end
        end
    end
end

% get global average alff
mean_alff = sum(sum(sum(alff_3d)))/(nnz(expmask));

% get normalized alff, which is defined as alff/its global mean 
norm_alff = alff_3d;   % initialize normalized alff as the original 3d map
for i = 1:ni
    for j = 1:nj
        for k = 1:nk
            norm_alff(i, j, k) = alff_3d(i, j, k)/mean_alff;
        end
    end
end

% display message on terminal
sprintf('normalized alff map calculated')

end   % end if com_alff == 1

% END Part (II): calculate amplitude of low-frequency fluctuations (alff)
%---------------------------------------------------------------------------

% Part (III): sort electrodes by their types

% use sort_elect to sort electrode array by their types
ele_sorted = sort_elect(electar);

% END PART (III): sorte electrodes by their types
%---------------------------------------------------------------------------
% % Part (IV): find midpoint of each pair of electrode contacts in both mni
% and image space

% use get_midpt_elect to locate midpoint of all electrode channels
% col. 1 = name of channel, 
% col. 2 = midpt in mni space, 
% col. 3 = midpt in image space
midpt = get_midpt_elect(ele_sorted, swra_info);

% END Part (IV): find midpoint of each pair of electrode contacts in both mni
% and image space
%---------------------------------------------------------------------------
% % Part (V): construct 3X3X3 boxes centered at midpoint of each electrode
% % pair to store signals enclosed thereof
% 
% % initialize array, sig_box, for storing the following for every type of elect.
% % 1st column lists the names of channels,
% % 2nd column lists reho (aka kendall_w coeff.) of every channel
% % 3rd column lists mean alff of every channel
sig_box = {};

% use for loop to get signals within each box from motion-filtered image
for i = 1:numel(midpt)   % for each electrode type
    for j = 1:size(midpt{i}, 1)   % for each electrode pair under current type
        sig_box{i}{j, 1} = midpt{i}{j, 1};   % get name tag of electrode pair
        mxyz = midpt{i}{j, 3};   % get midpoint in image space
        
        % get signals from motion-filtered images
        signal =  double(motion_filtered_img(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
            mxyz(3)-1:mxyz(3)+1, :));  
        
        % create regional mask from explicit mask using the coordinates of
        % the 3X3X3 box
        reg_mask = double(expmask(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
            mxyz(3)-1:mxyz(3)+1));
        
        % use Chebyshev Type 1 bandpass filter to allow freq. compo. btw. 
        % [0.01 and 0.1] Hz to pass through, use cheby1 to get 
        % transfer function coefficients (normalized freqs. are used)
        [b,a] = cheby1(2,0.5,[0.01 0.1]/(fs/2));
         
        % apply bandpass filter (Chebyshev filter) to remove low-freq. drift 
        % and noises from signal
        signal_noise_filt = signal;   % initialize array
        for x_vox = 1:size(signal, 1)   % for each voxel in x-dir.
            for y_vox = 1:size(signal, 2)   % for each voxel in y-dir.
                for z_vox = 1:size(signal, 3)   % for each voxel in z-dir.
                    % extract signal at voxel interested, squeeze the
                    % array to remove dim. of length 1, then use Chebyshev
                    % filter to remove noises
                    signal_noise_filt(x_vox, y_vox, z_vox, :) = ...
                        filtfilt(b, a, squeeze(signal(x_vox, y_vox, z_vox, :))); 
                end
            end
        end
       
        % get Kendall's W coefficient of filtered signal, assign
        % quantities calcu. to 3rd column of current cell within sig_box
        sig_box{i}{j, 2} = kendall_w(signal_noise_filt, reg_mask);
        
        % get mean normalized alff in current 3X3X3 box, assign result to
        % 3rd column of current cell within sig_box
        if ~isempty(norm_alff)   % if norm_alff array is not empty
            sig_box{i}{j, 3} = mean(mean(mean(norm_alff(mxyz(1)-1:mxyz(1)+1, ...
             mxyz(2)-1:mxyz(2)+1, mxyz(3)-1:mxyz(3)+1))));
        else
            sig_box{i}{j, 3} = [];   % assign empty value if norm_alff is empty
        end
    end
end

% END Part (V): construct 3X3X3 boxes centered at midpoint of each electrode
% pair to store signals enclosed thereof
%---------------------------------------------------------------------------

% Part (VI): group signals of all channels in single layer

% group all electrode types in one single layer of cell array
% col. 1 = id assigned to electrode pair
% col. 2 = name tag of electrode pair
% col. 3 = regional homogeneity (reho) of every channel
% col. 4 = mean alff of every channel
sig_box_all = {};   % initialize sig_box_all array
rownum = 1;   % initialize row number as 1 

for i = 1:numel(sig_box)   % for each electrode type in sig_box
    for j = 1:size(sig_box{i}, 1)   % for each pair of electrode contacts
        sig_box_all(rownum, 1) = {rownum};   % assign row number as id 
        sig_box_all(rownum, 2:2+size(sig_box{i}, 2)-1) = sig_box{i}(j, :);   % assign relevant info.
        rownum = rownum + 1;   % increment row number by 1
    end
end

% END Part (VI): group signals of all channels in single layer
%---------------------------------------------------------------------------
% Part (VII): match channels in spike rates array with that of signal box

% use clean_match_spikes to get spikes struct.
spikes = clean_match_spikes(spikesar, sig_box, sig_box_all);

% get spike rates array in single layer with channels matched with sig_box_all
spikes_all = spikes.spikes_all;

% END Part (VII): match channels in spike rates array with that of signal box
%---------------------------------------------------------------------------

% Part (VIII): get spearman's corr. coeff. btw. (i) ReHo and spike rates,
% and (ii) ALFF and spike rates

% initialize arrays of interest
rho_reho_spikes = [];
pval_reho_sp = [];
rho_malff_spikes = [];
pval_malff_sp = [];

% get spearman's rho btw. reho and spike rates
[rho_reho_spikes, pval_reho_sp] = corr([sig_box_all{:, 3}]', [spikes_all{:, 3}]', ...
    'Type', 'Spearman');

% get spearman's rho btw. alff and spike rates
if ~isempty(norm_alff)
[rho_malff_spikes, pval_malff_sp] = corr([sig_box_all{:, 4}]', [spikes_all{:, 3}]', ...
    'Type', 'Spearman');
end

% END Part (VIII): get spearman's corr. coeff. btw. (i) ReHo and spike rates,
% and (ii) ALFF and spike rates
%---------------------------------------------------------------------------

% Part (IX): assemble array for storing all variables interested
reho_alff_array = {};

% round all entries to 3 d. p.
for row = 1:size(sig_box_all, 1)
    reho_alff_array{row, 1} = sig_box_all{row, 1};   % id of channel
    reho_alff_array{row, 2} = sig_box_all{row, 2};   % name of channel
    reho_alff_array{row, 3} = round(spikes_all{row, 3}, 3);   % spike rate [/min]
    reho_alff_array{row, 4} = round(sig_box_all{row, 3}, 3);   % reho
    reho_alff_array{row, 5} = round(sig_box_all{row, 4}, 3);   % mean alff
end

reho_alff_array{1, 6} = round(rho_reho_spikes, 3);   % corr. coeff. btw. reho and spike rates
reho_alff_array{1, 7} = round(pval_reho_sp, 3);   % p-val. of corr. coeff. btw. reho and spike rates
reho_alff_array{1, 8} = round(rho_malff_spikes, 3);    % corr. coeff. btw. malff and spike rates
reho_alff_array{1, 9} = round(pval_malff_sp, 3);   % p-val. of corr. coeff. btw. malff and spike rates

% END Part (IX): assemble array for storing all variables interested
%---------------------------------------------------------------------------

% Part (X): generate tables listing all variables calculated

% convert cell array to table and add table headers
table_reho_alff = cell2table(reho_alff_array, ...
    'VariableNames',{'id', 'channel name', 'spike rate [/min]', 'reho', 'malff', ...
    'corr reho spikes', 'pval_corr_rs', 'corr malff spikes', 'pval_corr_ap'});

% Part (X): generate tables listing all variables calculated
%---------------------------------------------------------------------------

% Part (XI): store quantities calcu. in output struct.

opstruct = struct;
opstruct.ele_sorted = ele_sorted;   % info. of electrode (sorted by name)
opstruct.midpt = midpt;   % coordinates of midpoint between each pair of electrode contacts 
opstruct.sig_box = sig_box;   % signals within 3X3X3 box centered at each midpoint
opstruct.sig_box_all = sig_box_all;   % signals within 3X3X3 box centered at each midpoint of all channels in single layer
opstruct.spikes_all = spikes_all;   % list of spike rates of all electrode pairs in single layer

% reho and alff in every 3X3X3 box are already stored within sig_box
% and sig_box_all arrays, var. below are for reference only
opstruct.alff_3d = alff_3d;   % alff at each voxel
opstruct.mean_alff = mean_alff;   % global mean alff
opstruct.norm_alff = norm_alff;   % normalized alff at each voxel

% store table created above in struct.
opstruct.table_reho_alff = table_reho_alff;   % table listing all var. int.

% END Part (XIII): store quantities calcu. in output struct.
%---------------------------------------------------------------------------

end   % end function opstruct = get_reho_alff_spikes(input_array)