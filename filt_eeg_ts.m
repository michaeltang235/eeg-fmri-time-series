close all
clear all

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '27';

% enter path to directory where all input files are located
directname_eeg = ['/work/levan_lab/mtang/fmri_project/', 'sub', subnum, '/eeg_data'];

% enter format of input eeg files
filename_eeg = 'run*_gradient.mat';

% enter frequency threshold (in Hz) for filtering
freq_thresh = 30;   % in Hz

% enter path where ouput struct. is stored at
fname_op = [directname_eeg];   % direct. of output matrix
filename_op = 'filt_eeg_ts.mat';   % filename of output file

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

%----------------------------------------------------------------------------

% get paths of input files
[~, eeg_file_path, ~] = get_path(directname_eeg, filename_eeg);   % eeg data

% initialize structure, filt_eeg_st, to store filtered eeg time series of
% every run
filt_eeg_st = struct;

% load eeg data from path for every run, then apply filtering, store output
% in filt_eeg_st

for run_ind = 1:numel(eeg_file_path)
    
% load structure of current run
st_eeg_cur = load(eeg_file_path{run_ind});

% call filt_eeg function defined below to filter unwanted noise
[filt_eeg_ch, fs_cur] = filt_eeg(st_eeg_cur, freq_thresh);

% get current session number as listed in current eeg file path
sind = regexp(lower(eeg_file_path{run_ind}), 'run');   % start index
eind = regexp(lower(eeg_file_path{run_ind}), '_gradient');   % end index
sess_num = eeg_file_path{run_ind}(sind+3:eind-1);   % get run number

% create fieldname representing current session number
fieldname = sprintf('run%s', sess_num);

% assign filtered eeg time series and sampling freq. to currnent session of
% filt_eeg_st structure
filt_eeg_st.(fieldname).filt_eeg = filt_eeg_ch;   % filtered eeg time series
filt_eeg_st.(fieldname).fs = fs_cur;   % sampling freq., in Hz

%---------------------------------------------------------------------------
% output structure created in current session
if op_results == 1
    save(fullfile(fname_op, filename_op), 'filt_eeg_st', '-v7.3');   % struct.
end
%---------------------------------------------------------------------------

end   % end for run_ind = 1:numel(eeg_file_path)

function [eeg_ch, fs] = filt_eeg(st_eeg, freq_thresh)

%----------------------------------------------------------------------------
% PART (I): ASSIGN EEG TIME SERIES TO ITS CHANNEL IN SINGLE-LAYER ARRAY 

% get sampling freq in Hz from input file
fs = st_eeg.fs;

% use Chebyshev Type 1 bandpass filter to allow high freq. compo.
% to pass through, use cheby1 to get 
% transfer function coefficients (normalized freqs. are used)
% [b,a] = cheby1(2, 0.5, freq_thresh/(fs/2), 'high');
% [b,a] = cheby1(2, 0.5, freq_thresh/(fs/2), 'low');
[b,a] = cheby1(2, 0.5, [0.01 freq_thresh]/(fs/2));
% [b,a] = cheby1(2, 0.5, [8 freq_thresh]/(fs/2));

% initialize array, eeg_ch, with format given below,
% 1st col. = channel name
% 2nd col. = filtered eeg time series
eeg_ch = {}; 
for i = 1:size(st_eeg.montage_names, 1)   % for every channel listed in eeg data file
    
    ch_name = [st_eeg.montage_names{i, 1}, '-', st_eeg.montage_names{i, 2}];   % assmeble channel name
    eeg_ts_int = st_eeg.data_gradient(:, i)';
    
    eeg_ch{i, 1} = ch_name;
    eeg_ch{i, 2} = filtfilt(b, a, eeg_ts_int);   % filtered eeg time series
%     eeg_ch{i, 2} = [filtfilt(b, a, eeg_ts_int)].^2;   % filtered eeg time series
end

% END PART (I): ASSIGN EEG TIME SERIES TO ITS CHANNEL IN SINGLE-LAYER ARRAY 
%----------------------------------------------------------------------------

end



