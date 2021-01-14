tic 

clear all 
close all
%---------------------------------------------------------------------------
% Based on template generated using "Batch" window in SPM, this script
% performs pre-processing steps on functional image dataset

% spm stores lines of command in a cell array named matlabbatch and execute the
% commands using a function called spm_jobman

% This script conducts pre-processing steps for functional images obtained
% and provides the option of running the pre-processing steps until
% Model Validation, the remaining steps can be done using another script
% "preproc_func_contrast.m"

% several options are provided in this script, user can either
% (i) set (run_spm = 0), check contrast matrix constructued and matlabbatch 
% array assembled in this analysis. This option does not call spm pre-proc. function!
% (ii) set (run_spm = 1), run spm analysis 
% (iii) set (run_spm = 1, run_contrast = 0), only run pre-proc. to the 
% step of model evalution, without running contrast manager and results
% report
% (iv) set (run_spm = 1, run_contrast = 1), run all pre-proc. steps, 
% including contrast maanger and results report.


% enter information about the image dataset and design matrix before
% running the script
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% enter subject number (str)
subnum = '41';

% enter path to directory where all input files are located
directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum];

% enter number of runs (num. of .nii file = % number of session(s))
numsess = 3;

% enter format of input file(s)
filename_run_format = 'Run%d.nii';   % %d is added for string formatting later

% enter directory and filename of onset times array and prelim. results
directnamemat = [directname, '/matrices'];
filenamemat = ['onset_times_events_', 'sub', subnum, '.mat'];   % onset times array for all runs
filename_prelim = ['prelim_preproc_', 'sub', subnum, '.mat'];   % prelim. results

% enter slice acqusiition order (odd = 1, even = 2)
sordernum = 1;

% enter final (normalized) voxel size in all 3 dim. (mm)
voxelsize = [3.75 3.75 3.75]; 

% enter info. about design matrix of 1st level analysis
unitdesign = 'secs';   % unit of design matrix

% sometimes, we are interested in checking the contrast matrix but not 
% running any pre-processing steps on image dataset, 
% this will only assesmble the matlabbatch arrays, without calling the spm
% function to run them
% enter if user wants to call spm function (1=yes, 0=no)
run_spm = 0;   

% enter if preliminary results (e.g. contrast matrix) 
% should be written to path (1 = yes, 0 = no)
op_prelim = 1;

% enter if the step of contrast manager should be run (1=yes, 0=no)
run_contrast = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
% PART (I): PRELIMINARY

%----------------------------------
% (IA and IB): input image filenames and details of images

% (IA): create array of input image filename(s)
% creat an array of image filename
filename_runs = cell(numsess, 1);
for i = 1:numsess
    filename_runs{i} = sprintf(filename_run_format, i);
end

% (IB): details of images

% read nifti image info. from file(s) 
imginfo = cell(numsess, 1);
for i = 1:numsess
    imginfo{i} = niftiinfo(fullfile(directname, filename_runs{i}));
end

% get number of images in different runs
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

% print image info. to console
for i = 1:numsess
sprintf('run%d, numimage = %d, numslices = %d, TR = %f', i, numimage{i}, numslices, TR)
end

% set slice acquisition order and reference slice
sliceorder = [];
refslice = 0;
if sordernum == 1   % odd slices first
    sliceorder = [1:2:numslices-1 2:2:numslices];   % odd then even
    refslice = 1;   % ref. slice 
end
if sordernum == 2   % even slices first
    sliceorder = [2:2:numslices 1:2:numslices-1];   % even then odd
    reslice = 2;   % ref. slice
end

% END (IA and IB): input image filenames and details of images
%----------------------------------
% (IC): definition of events and associated onset times

% load onset times array
terms = load(fullfile(directnamemat, filenamemat)); 

% if no onset times are recorded for an event, an empty matrix is assigned
% select non-empty entries from term struct. array, get following variables
% (i) number of conditions in each run and the total
% (ii) def. each event (str) identified in each run
% (iii) def. of each event (i.e. in matrix) identified in each run
% (iv) onset times of each evenet identified in each run
numcond = cell(1, numsess);   % initialize number of conditions cell array
numcond_est = 0;   % number of conditions established
event_def_str = cell(1, numsess);   % initialize event_def_str cell array
event_def = cell(1, numsess);   % initialize event_def cell array
evti = cell(1, numsess);   % intialize evti cell array
for i = 1:numsess   
    % for specific run
    fieldname_run = sprintf('run%d', i);
    
    % get non-empty event def. and onset times (evti)
    numcond{i} = terms.terms.(fieldname_run).results.num_cond;
    numcond_est = numcond_est + terms.terms.(fieldname_run).results.num_cond;
    event_def_str{i} = terms.terms.(fieldname_run).results.event_def_str;
    event_def{i} = terms.terms.(fieldname_run).results.event_def;
    evti{i} = terms.terms.(fieldname_run).results.evti;
end

% END (IC): definition of events and associated onset times
%----------------------------------
% (ID): create contrast matrix and name of each contrast

% create contrast matrix:

% create cond_mat for storing contrasts for each condition in each run
cont_mat = diag(ones(1, numcond_est));   % assign 1s along main diagonal 

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

% create condition matrix: 
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
    arreq = zeros(1, numcond_est);   % create array of 0s, array required
    arreq(indreq) = 1;   % assign value of 1 at indices found earlier
    
    % concatenate cont_mat and arreq vertically 
    cont_mat = [cont_mat; arreq];
end

%-----------------------
% create name of each contrast:

% create cont_name_str for storing string of each condition established,
% by going through each cell of event_def_str, take transpose, and add it to array
cont_name_str = cell(size(cont_mat, 1), 1);   % preallocating memory
rownum = 1;
for i = 1:numsess   % for each run
    rowspan = length(event_def_str{i}');   % get row span
    % fitting the i^th entry of event_def_str to cont_name_str
    cont_name_str(rownum:rownum + rowspan - 1) = event_def_str{i}';  
    % increment rownum by sum of current rownum and rowspan
    % for storing the next string in the row below
    rownum = rownum + rowspan;   
end

% then appened condition names to cont_name_str for each unique condition 
% applied to all runs
% e.g. type(s)_1_21_run123, denotes event types 1 and 21 in run 1, 2, and 3

for i = 1:numel(uni_event_def)
    rownum = numcond_est + i;   % for remaining rows, numcond_est is total num. of conditions
    currlen = length(uni_event_def{i});   % get length of current item in uni_event_def
    placestr = repmat('%d_', 1, currlen);   % prepare place holder for formatting event_type
    runstr = 'all';   % prepare place holder for run numbers
%     runstr = repmat('%d', 1, numsess);   % prepare place holder for run numbers
    runnum = 1:1:numsess;   % get total number of runs
    cont_name_str{rownum} = append('type(s)_', sprintf(placestr, uni_event_def{i}), ...
        'run_', sprintf(runstr, runnum));   % assemble name of contrasts
end

% END (ID): create contrast matrix and name of each contrast
%----------------------------------
% (IE): configure filenames each step in pre-processing of spm algorithm

% create arrays of filenames for steps of slice time correction and
% realignment
filename_path_slicetime = cell(numsess,1);   % slice time corr. requires info. about session
filename_path_realign = cell(numsess,1);  % realignment requires info. about session
for sess = 1:numsess   % for each session
    for i = 1:numimage{sess}   % for each image
    filename_path_slicetime{sess}{i,1} = append(fullfile(directname, ...
        filename_runs{sess}), sprintf(',%s', num2str(i)));
    filename_path_realign{sess}{i,1} = append(directname, '/a', ...
        filename_runs{sess}, sprintf(',%s', num2str(i)));   % assing prefix 'a'
    end
end

% create arrays of filenames for steps of normalization and smoothing

% only first image is selected for other images to realign to
filename_imagetoallignto = append(directname, '/', 'ra', filename_runs{1}, ',1');

% get path for pre-processed files (i.e. append 'ra' to name for normalization 
% and 'wra' to name for smoothing) 
filename_path_normalize = cell(sum([numimage{:}]),1);   % normalization, initialize arrays
filename_path_smooth = cell(sum([numimage{:}]), 1);   % smoothing
rowind= 1;   % initialize row index as 1

for sess = 1:numsess   % for each session
    for i = 1:numimage{sess}   % for each image in each session
        % format path for normalization
        filename_path_normalize{rowind} = append(directname, '/', 'ra', ...
            filename_runs{sess}, sprintf(',%s', num2str(i)));   
        % format path for smoothing
        filename_path_smooth{rowind} = append(directname, '/', 'wra', ...
            filename_runs{sess}, sprintf(',%s', num2str(i)));  
        rowind = rowind + 1;   % increment row index by 1
    end
end

% get path for pre-processed files (i.e. append 'swra' to name) for the
% step of model specification
filename_path_model = cell(numsess,1);
for sess = 1:numsess
    for i = 1:numimage{sess}
    filename_path_model{sess}{i,1} = append(directname, '/swra', ...
        filename_runs{sess}, sprintf(',%s', num2str(i)));
    end
end

% set up filename of spm.mat file for steps of model evaluation, 
% contrast manager, and results reporting
filename_spm = 'SPM.mat';   % name of spm.mat file (the GLM)

% END (IE): configure filenames each step in pre-processing of spm algorithm
%----------------------------------

% (IF): write quantities calculated in preliminary to file

% create prelim struct
prelim = struct;

% assign related quantities to struct.
prelim.imginfo = imginfo;   % image info from nifti headers
prelim.numimage = numimage;   % number of images in every run
prelim.numslices = numslices;   % number of slices (z-dim.)
prelim.TR = TR;   % repetition time (s)
prelim.numcond = numcond;   % number of conditions in every run 
prelim.event_def_str = event_def_str;   % event def. str in every run
prelim.event_def = event_def;   % event definition found in each run
prelim.evti = evti;   % onset time (ti) of each event identified in each run
prelim.uni_event_def = uni_event_def;   % unique events identified in all runs
prelim.cont_mat = cont_mat;   % contrast matrix (for contrast manager)
prelim.cont_name_str = cont_name_str;   % string of each contrast defined

% write prelim struct. to path
if op_prelim == 1
save(fullfile(directnamemat, filename_prelim),'prelim');
end

% END PART (I): PRELIMINARY
%----------------------------------------------------------------
%----------------------------------------------------------------

% PART (II): SPM script

% initialize spm module without graphical interface
spm('defaults','fmri');
spm_jobman('initcfg');

% %-------------------------------------------------
% % Slice time correction:
matlabbatch{1}.spm.temporal.st.scans = filename_path_slicetime';
matlabbatch{1}.spm.temporal.st.nslices = numslices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR - TR/numslices;
matlabbatch{1}.spm.temporal.st.so = sliceorder;
matlabbatch{1}.spm.temporal.st.refslice = refslice;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

% end Slice time correction
% %-------------------------------------------------
% Realign: Estimate and Reslice:

matlabbatch{2}.spm.spatial.realign.estwrite.data = filename_path_realign;

matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = voxelsize(1);   % separation (mm) between points in sampled image;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

% end Realign: Estimate and Reslice
% %-------------------------------------------------
% % Normalization: Estimate and write:

% % only first image is selected for other images to realign to
matlabbatch{3}.spm.spatial.normalise.estwrite.subj.vol = {filename_imagetoallignto};

matlabbatch{3}.spm.spatial.normalise.estwrite.subj.resample = filename_path_normalize;

matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/Users/michaeltang/Downloads/spm12/tpm/TPM.nii'};
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.vox = voxelsize;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

% % end Normalization: Estimate and write:
% %-------------------------------------------------
% % Smoothing:

matlabbatch{4}.spm.spatial.smooth.data = filename_path_smooth;
matlabbatch{4}.spm.spatial.smooth.fwhm = 2*voxelsize;   % 2 x size of voxel
matlabbatch{4}.spm.spatial.smooth.dtype = 0;
matlabbatch{4}.spm.spatial.smooth.im = 0;
matlabbatch{4}.spm.spatial.smooth.prefix = 's';

% % end Smoothing
% % %-------------------------------------------------
% % FMRI model specification:

% % assign variables to matlabbatch array
matlabbatch{5}.spm.stats.fmri_spec.dir = {directname};
matlabbatch{5}.spm.stats.fmri_spec.timing.units = unitdesign;
matlabbatch{5}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{5}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{5}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

% %----------------------------
% % loop through each session and each event
for sessind = 1: numsess
    for eventind = 1:numcond{sessind}
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).scans = filename_path_model{sessind};
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).cond(eventind).name = event_def_str{sessind}{eventind};
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).cond(eventind).onset = evti{sessind}{eventind};
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).cond(eventind).duration = 0;
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).cond(eventind).tmod = 0;
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).cond(eventind).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).cond(eventind).orth = 1;
    end
    
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).multi = {''};
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).regress = struct('name', {}, 'val', {});
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).multi_reg = {''};
    matlabbatch{5}.spm.stats.fmri_spec.sess(sessind).hpf = 128;
end
% %----------------------------

matlabbatch{5}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{5}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{5}.spm.stats.fmri_spec.volt = 1;
matlabbatch{5}.spm.stats.fmri_spec.global = 'None';
matlabbatch{5}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{5}.spm.stats.fmri_spec.mask = {''};
matlabbatch{5}.spm.stats.fmri_spec.cvi = 'AR(1)';
% 
% % end FMRI model specification:
% %-------------------------------------------------
% % Model estimation:

matlabbatch{6}.spm.stats.fmri_est.spmmat = {fullfile(directname, filename_spm)};
matlabbatch{6}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{6}.spm.stats.fmri_est.method.Classical = 1;

% % end Model estimation
% % -------------------------------------------------
% % Contrast manager: (only run this step if runcontrast == 1)

if run_contrast == 1

filename_spm = 'SPM.mat';   % name of spm.mat file (the GLM)
matlabbatch{7}.spm.stats.con.spmmat = {fullfile(directname, filename_spm)};

% use loop to assign required variables (i.e. name, contrast) to matlabbatch array
for connum = 1:size(cont_mat, 1)   % through each conditrast
    matlabbatch{7}.spm.stats.con.consess{connum}.tcon.name = cont_name_str{connum};
    matlabbatch{7}.spm.stats.con.consess{connum}.tcon.weights = cont_mat(connum,:);
    matlabbatch{7}.spm.stats.con.consess{connum}.tcon.sessrep = 'none';
end

matlabbatch{7}.spm.stats.con.delete = 1;

% % end Contrast manager
% %-------------------------------------------------
% % Results report

matlabbatch{8}.spm.stats.results.spmmat = {fullfile(directname, filename_spm)};

% use loop to assign required variables to matlabbatch array, applying
% contrast to each run
for i = 1:size(cont_mat, 1)
matlabbatch{8}.spm.stats.results.conspec(i).titlestr = '';
matlabbatch{8}.spm.stats.results.conspec(i).contrasts = i;
matlabbatch{8}.spm.stats.results.conspec(i).threshdesc = 'FWE';
matlabbatch{8}.spm.stats.results.conspec(i).thresh = 0.05;
matlabbatch{8}.spm.stats.results.conspec(i).extent = 0;
matlabbatch{8}.spm.stats.results.conspec(i).conjunction = 1;
matlabbatch{8}.spm.stats.results.conspec(i).mask.none = 1;
end

matlabbatch{8}.spm.stats.results.units = 1;
matlabbatch{8}.spm.stats.results.export{1}.ps = true;

end   % end if runcontrast == 1

% % end Result report
% % -------------------------------------------------

% % END PART (II): SPM script
% %-----------------------------------------------------------------------------

% calling spm function with run_spm is set to 1
if run_spm == 1
% % call spm_jobman to run the script above
spm_jobman('run', matlabbatch);
end


toc