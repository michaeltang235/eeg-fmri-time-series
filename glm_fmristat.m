tic

close all
clear all

% This script uses fmristat toolbox to statistically analyze fmri data,
% including creating design matrix, fitting linear model, and calculating
% t-statistics.

% The following .m files in the toolbox are utilized
% (i): fmridesign.m for constructing design matrix of processed fmri data
% (ii): fmrilm.m for estimating beta parameters in design matrix
% (iii) multistat.m for calculating t-stats from beta images created in (ii)
%
% Remarks:
% (i): linear model created in each run is stored under a struct. named lm
% (ii): calling fmrilm.m will write images of beta and beta standard
% devations to path
% (iii): user can select if step of model estimation is run or not

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '27';

% enter path to directory where all input files are located
directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum, '_imthres0_exmask'];

% enter number of runs (num. of processed .nii file = % number of session(s))
numsess = 3;

% enter format of processed fmri image file(s)
proc_filename_format = 'swraRun%d_10min.nii';   % %d is added for string formatting later

% format directory and filename of onset times array and prelim. results
directnamemat = [directname, '/matrices'];
onset_times_filename = ['onset_times_events_', 'sub', subnum, '.mat'];   % onset times array for all runs
% prelim_filename = ['prelim_preproc_', 'sub', subnum, '.mat'];   % prelim. results

% enter if step of model estimation should be executed (takes a long time)
modest = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
% PART (I): PRELIMINARY

% creat an array of processed image filenames
proc_filename = cell(numsess, 1);
for i = 1:numsess
    proc_filename{i} = sprintf(proc_filename_format, i);
end

% read nifti image info. from file(s) 
imginfo = cell(numsess, 1);
for i = 1:numsess
    imginfo{i} = niftiinfo(fullfile(directname, proc_filename{i}));
end

% get number of processed images in different runs
numimage = cell(numsess, 1);
for i = 1:numsess
    numimage{i} = imginfo{i}.raw.dim(5);   
end

% while numslices and TR should be the same for all runs, numbers of images
% could differ among diff. runs, 
% check if numbers of images are the same in all runs
for i = 1:numsess-1 
    tmp = numimage{i};   % temporary var.
    tmp1 = numimage{i+1};   % temporary var1.
    if tmp ~= tmp1   % print warning message if num. of images  different
        sprintf('warning, numbers of images are different among runs')
    end
end

% get number of slices (length of z-direct.) and repetition time (s)
numslices = imginfo{1}.raw.dim(4);  
TR = imginfo{1}.raw.pixdim(5); 

% print processed image info. to console
for i = 1:numsess
sprintf('run%d, numimage = %d, numslices = %d, TR = %f', i, numimage{i}, numslices, TR)
end

% load onset times array
terms = load(fullfile(directnamemat, onset_times_filename));

% if no onset times are recorded for an event, an empty matrix is assigned
% select non-empty entries from term struct. array, get following variables
% (i): def. each event (str) identified in each run
% (ii) def. of each event (i.e. in matrix) identified in each run
% (iii): onset times of each evenet identified in each run
event_def_str = cell(1, numsess);   % initialize event_def_str cell array
event_def = cell(1, numsess);   % initialize event_def cell array
evti = cell(1, numsess);   % intialize evti cell array
for i = 1:numsess   
    % for specific run
    fieldname_run = sprintf('run%d', i);
    
    % get non-empty event def. and onset times (evti)
    event_def_str{i} = terms.terms.(fieldname_run).results.event_def_str;
    event_def{i} = terms.terms.(fieldname_run).results.event_def;
    evti{i} = terms.terms.(fieldname_run).results.evti;
end

% END PART (I): PRELIMINARY
%---------------------------------------------------------------------------
% PART (II): Model Specification

% create struct. lm for storing linear model created in each run
lm = struct;

% perform model specification and model estimation for each run
for runind = 1:numsess

% get frame times (image acqusition times in seconds)
frametimes = TR*(0:numimage{runind}-1);

% set relslicetimes (relative slice acquisition times). As all slices in 
% each processed image are assumed to be acquired at the same time, 
% the relative slice time is 0
relslicetimes = zeros(1, numslices);

% assemble event matrix with the following spec.
% column 1 = index of event type (1,2,3,...)
% col. 2 = event onset times 
% col. 3 = event durations
event_matrix = [];  % initialize empty event_matrix

% loop through each event type in evti array in current run, append 
% relevant info. to matrix vertically
for i = 1:numel(evti{runind})    % for each event type in current run
    event_matrix = [event_matrix; [i*ones(numel(evti{runind}{i}), 1) ...
        evti{runind}{i} zeros(numel(evti{runind}{i}), 1)]];
end

% get fmri design matrix (GLM) using fmrilm
X_lm = fmridesign(frametimes, relslicetimes, event_matrix);

% create fieldname for current run, which will be used in struct. lm
fieldname = sprintf('run%d', runind);

% store linear model created in current run in struct. lm
lm.(fieldname).X_lm = X_lm;

% END PART (II): Model Specification
%---------------------------------------------------------------------------
% PART (III): Model Estimation

% get full path to processed image file in current run
input_file = [directname, filesep, proc_filename{runind}];

% format full paths of output files in current run
output_file_base = cell(numel(event_def_str{runind}, 1));
for i = 1:numel(event_def_str{runind})
    output_file_base{i, 1} = [directname, filesep, event_def_str{runind}{i}];
end

% create contrast matrix for each condition in current run
cont_mat = eye(numel(evti{runind}));   % assign 1s along main diagonal 

% specify stats (str) interested from fmrilm 
which_stats = '_mag_t _mag_ef _mag_sd';

% estimate GLM defined earlier using fmrilm
if modest == 1
    fmrilm(input_file, char(output_file_base), X_lm, cont_mat, [], which_stats);
end

% print message to terminal
sprintf('linear model in run%d is estimated', runind)

end   % end for runind = 1:numsess

% [df, spatial_av] = fmrilm(input_file, char(output_file_base), X_lm, cont_mat, [], which_stats);
% fmrilm(input_file,char(output_file_base),X,contrast,[],'_mag_t _mag_ef _mag_sd');

% END PART (III): Model Estimation
%---------------------------------------------------------------------------
% PART (IV): Model Statistics

%------------------------------------
% PART (IV_A): create paths to required image files as inputs to multistat

% Generate combined activation maps from all the runs
% The *_ef and *_sd files are generated by the model estimation step
% (*_ef images are the beta images, and *_sd are the beta standard deviation images)

% get number of conditions established in all runs from term struct.
num_cond_tot = 0;
for i = 1:numsess
    fieldname = sprintf('run%d', i);
    num_cond_tot = num_cond_tot + terms.terms.(fieldname).results.num_cond;
end

% assemble arrays of input fmri effect files (*_ef and *_sd.nii)
% arrays are col. vector with each row describing path to each image file
beta_images_files = cell(num_cond_tot, 1);   % array of mag. of beta images
beta_sd_images_files = cell(num_cond_tot, 1);   % array of beta standard deviation images
rownum = 1;   % initialize row number as 1
for runind = 1:numsess   % for each run
    for j = 1:length(event_def_str{runind})   % for each event type defined in current run
        %   append strings together to form full paths of required images
        beta_images_files{rownum, 1} = [directname, filesep, event_def_str{runind}{j}, '_mag_ef', '.nii'];
        beta_sd_images_files{rownum, 1} = [directname, filesep, event_def_str{runind}{j}, '_mag_sd', '.nii'];
        rownum = rownum + 1;   % increment row number by 1
    end
end

% convert arrays from cell to character arrays for input of the function multistat
input_files_ef = char(beta_images_files);
input_files_sd = char(beta_sd_images_files);

% END PART (IV_A): create paths to required image files as inputs to multistat
%------------------------------------
% PART (IV_B): identify unique events in all runs and create associated
% contrast matrix for multistat

% create contrast matrix for multistat with each column representing a
% combination of event types identified in all runs

% identify unique events in all runs

% initialize array with the first entry of event_def in the first run
uni_event_def{1} = event_def{1}{1};

% scan through each entry of event_def (A) and compare with each existing entry
% in uni_event_def (B), if A is not in B, add A to B
for runind = 1:numel(event_def)   % for each run in event_def
    for i = 1:numel(event_def{runind})   % for each entry in that run in event_def
        count = 0;   % initialize count as 0
        for item = 1:numel(uni_event_def)   % for each item in uni_event_def
            % check if current entry in event_def equals to any item in
            % uni_event_def,
            if isequal(event_def{runind}{i}, uni_event_def{item})   
                count = count + 1;   % if so, increment count by 1
            end
        end
        % if count remains 0 after going through each item in
        % uni_event_def (i.e. not in it), add current event_def to array
        if count == 0   
            uni_event_def = [uni_event_def event_def{runind}{i}];
        end
    end
end

% create event_def_all storing all event types identified in all runs
event_def_all = {};
for i = 1:numel(event_def)
    event_def_all = [event_def_all event_def{i}];
end

% create contrast matrix: each row representing a combinatio of events
cont_mat_all_runs = [];

% use the unique event_def{i} for comparing with other cells in
% event_def_all, get indices of duplicate entries
for i = 1:numel(uni_event_def)   % for each item in uni_event_def
    indreq = [];   % indices required
    for j = 1:numel(event_def_all)   % for each entry in event_def_all
        if isequal(uni_event_def{i}, event_def_all{j})   % compare entries
            indreq = [indreq j];   % if yes, add index to indreq
        end
    end
    % after looping through all entries in event_def_all to search for
    % matching entries with values equal to the i^th element in uni_event_def
    arreq = zeros(1, num_cond_tot);   % create array of 0s, array required
    arreq(indreq) = 1;   % assign value of 1 at indices found earlier
    
    % concatenate cont_mat and arreq vertically 
    cont_mat_all_runs = [cont_mat_all_runs; arreq];
end

% END PART (IV_B): identify unique events in all runs and create associated
% contrast matrix for multistat
%------------------------------------
% PART (IV_C): assemble remaining input arg. and execute multistat 

% take tranpose of cont_mat_all_runs and assign value to X_multistat,
% a design matrix, which is one of the input arg. of multistat
% rows denote image files and columns denote combo. of events in all runs
% e.g.: 1st column combines files type(s)_1_11_run_1, type(s)_1_11_run_2 and type(s)_1_11_run_3
X_multistat = cont_mat_all_runs';

% create contrast_multistat matrix, also an input of multistat, rows are 
% contrasts for stat. images. e.g. for [1 0 ... 0], it picks out the first 
% column of X_multistat.
% num. of col. in X_multistat equals to num. of rows in contrast_multistat
contrast_multistat = eye(numsess);   % identity matrix 

% create contrast name strings for all runs
cont_name_str_all_runs = cell(length(uni_event_def), 1);   % initialize array and preallocating memory
for i = 1:numel(uni_event_def)   % for each unique event def. identified
    currlen = length(uni_event_def{i});   % get length of current item in uni_event_def
    placestr = repmat('%d_', 1, currlen);   % prepare place holder for formatting event_type  
    % assemble contrast name string for current type
    cont_name_str_all_runs{i, 1} = ['type_', sprintf(placestr, uni_event_def{i}), 'run_all'];
end

% create output filename base: filename base for output stat., one base for
% each contrast, will be padded with the type of stat. specified
contrast_name_base_multi = cell(length(uni_event_def), 1);   % initialize array
for i = 1:length(uni_event_def)   % for each unique event 
    contrast_name_base_multi{i,1} = [directname, filesep, cont_name_str_all_runs{i}];
end

% convert contrast_name_base_multi from cell to char array and use it as 
% input of multistat
output_file_base_multistat = char(contrast_name_base_multi);

% set fwhm_varatio to inf for fixed effects analysis
fwhm_varatio = Inf;

% call multistat function, specify to output t stat. images ('_t')
multistat(input_files_ef, input_files_sd, [], [], X_multistat, ...
    contrast_multistat, output_file_base_multistat, '_t', fwhm_varatio);

% END PART (IV_C): assemble remaining input arg. and execute multistat 
%------------------------------------

toc


