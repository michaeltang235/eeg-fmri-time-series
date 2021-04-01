tic 

clear all 
close all
%---------------------------------------------------------------------------
% Based on template generated using "Batch" window in SPM, this script
% performs pre-processing steps on functional image dataset

% spm stores lines of command in a cell array named matlabbatch and execute the
% commands using a function called spm_jobman

% This script conducts pre-processing steps listed below for functional images 
% obtained, 
% (i) slice time correction,
% (ii) realignment
% (iii) normalization
% (iv) smoothing
% leaving the remaining steps startng from model specification to another
% script which uses fmristat

% several options are provided in this script, user can either
% (i) set (run_spm = 0), check contrast matrix constructued and matlabbatch 
% array assembled in this analysis. This option does not call spm pre-proc. function!
% (ii) set (run_spm = 1), run spm analysis 

% Remarks:
% (i) no excess functional images (.nii) should be stored in working
% directory. i.e. only include images that are used in this script 

% Updates:
% Feb. 4, 2021, to see the effect of explicit mask on the pre-processing 
% algorithm, run spm with (i) explicit mask imposed and (ii) no ex. mask 
% imposed, with implicit mask set to 0 and see what happens

% Mar. 30, 2021, (i) this script has been modified to only include portions
% that run steps from slice-timing correction to smoothing, leaving the
% steps starting from model specification to another script which uses fmristat

% (ii) modified the way paths of filenames are formatted, that is, file 
% separator is determined automatically by the platform in which this script 
% is run on

% (iii) requires users to specify path to spm package, so that this script can be
% executed in different platforms

% enter information about the image dataset and design matrix before
% running the script
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% enter subject number (str)
subnum = '14';

% enter path to directory where all input files are located
directname = ['C:\Users\siumichael.tang\Downloads\fmri_project\', 'sub', subnum];
% directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum, '_imthres0_exmask'];

% enter number of runs (num. of .nii file = % number of session(s))
numsess = 2;

% enter format of filename of input file(s)
filename_format = 'Run*.nii';   % wildcard char is added for searching

% enter slice acqusiition order (odd = 1, even = 2)
sordernum = 1;

% enter final (normalized) voxel size in all 3 dim. (mm)
voxelsize = [3.75 3.75 3.75]; 

% enter info. about design matrix of 1st level analysis
unitdesign = 'secs';   % unit of design matrix

% enter full path to spm package
spm_path = 'C:\Users\siumichael.tang\Downloads\spm12';

% sometimes, we are interested in checking the contrast matrix but not 
% running any pre-processing steps on image dataset, 
% this will only assesmble the matlabbatch arrays, without calling the spm
% function to run them
% enter if user wants to call spm function (1=yes, 0=no)
run_spm = 1;   

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
% PART (I): PRELIMINARY

%----------------------------------
% (IA): locate paths of func. images and their details 

% get number and full path of func. images using get_path function
[numfile, file_path, file_name_runs] = get_path(directname, filename_format);

% check if number of file(s) found equals to number of sessions entered, if
% not, print warning message
if numfile ~= numsess
    sprintf(['number of file(s) found is different than number of sessions', ...
        ' entered, please check'])
end

% read nifti image info. from file(s) 
imginfo = cell(numfile, 1);
for i = 1:numfile
    imginfo{i} = niftiinfo(file_path{i});
end

% get number of images in different runs
numimage = cell(numfile, 1);
for i = 1:numfile
    numimage{i} = imginfo{i}.raw.dim(5);   
end

% while numslices and TR should be the same for all runs, numbers of images
% could differ among diff. runs, 
% check if numbers of images are the same in all runs
for i = 1:numfile-1 
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

% END (IA): locate paths of func. images and their details 
%----------------------------------
% (IB): configure filenames each step in pre-processing of spm algorithm

% create arrays of filenames for steps of slice time correction and
% realignment
filename_path_slicetime = cell(numfile,1);   % slice time corr. requires info. about session
filename_path_realign = cell(numfile,1);  % realignment requires info. about session
for sess = 1:numfile   % for each session
    for i = 1:numimage{sess}   % for each image
    filename_path_slicetime{sess}{i,1} = [directname, ...
        filesep, file_name_runs{sess}, sprintf(',%s', num2str(i))];
    filename_path_realign{sess}{i,1} = [directname, ...
        filesep, 'a', file_name_runs{sess}, sprintf(',%s', num2str(i))];   % assing prefix 'a'
    end
end

% create arrays of filenames for steps of normalization and smoothing
% only first image from first session is selected for other images to realign to
filename_imagetoalignto = [directname, filesep, 'ra', file_name_runs{1}, ',1'];

% get path for pre-processed files (i.e. append 'ra' to name for normalization 
% and 'wra' to name for smoothing) 
filename_path_normalize = cell(sum([numimage{:}]),1);   % normalization, initialize arrays
filename_path_smooth = cell(sum([numimage{:}]), 1);   % smoothing
rowind= 1;   % initialize row index as 1

for sess = 1:numsess   % for each session
    for i = 1:numimage{sess}   % for each image in each session
        % format path for normalization
        filename_path_normalize{rowind} = [directname, filesep 'ra', ...
            file_name_runs{sess}, sprintf(',%s', num2str(i))];   
        % format path for smoothing
        filename_path_smooth{rowind} = [directname, filesep, 'wra', ...
            file_name_runs{sess}, sprintf(',%s', num2str(i))];  
        rowind = rowind + 1;   % increment row index by 1
    end
end

% format path to tissue probability map (tpm image) stored in SPM package
tpm_path = [spm_path, filesep, 'tpm', filesep, 'TPM.nii'];

% END (IB): configure filenames each step in pre-processing of spm algorithm
%----------------------------------

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
matlabbatch{3}.spm.spatial.normalise.estwrite.subj.vol = {filename_imagetoalignto};

matlabbatch{3}.spm.spatial.normalise.estwrite.subj.resample = filename_path_normalize;

matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.tpm = tpm_path;
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
% % % FMRI model specification:
% % % assign variables to matlabbatch array
% matlabbatch{1}.spm.stats.fmri_spec.dir = {directory_path_spmmat};
% matlabbatch{1}.spm.stats.fmri_spec.timing.units = unitdesign;
% matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
% matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = numslices;
% matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = refslice;
% 
% % %----------------------------
% % % loop through each session and each event
% for sessind = 1: numsess
%     for eventind = 1:numcond{sessind}
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).scans = filename_path_model{sessind};
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).cond(eventind).name = event_def_str{sessind}{eventind};
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).cond(eventind).onset = evti{sessind}{eventind};
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).cond(eventind).duration = 0;
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).cond(eventind).tmod = 0;
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).cond(eventind).pmod = struct('name', {}, 'param', {}, 'poly', {});
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).cond(eventind).orth = 1;
%     end
%     
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).multi = {''};
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).regress = struct('name', {}, 'val', {});
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).multi_reg = {''};
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sessind).hpf = 128;
% end
% % %----------------------------
% 
% matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
% matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
% matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
% matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
% matlabbatch{1}.spm.stats.fmri_spec.mthresh = imthresh;
% matlabbatch{1}.spm.stats.fmri_spec.mask = {filename_path_expmask};
% matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
% % 
% % % end FMRI model specification:
% %-------------------------------------------------
% % % Model estimation:
% % matlabbatch{2}.spm.stats.fmri_est.spmmat = {filename_path_modelest};
% % matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
% % matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
% % % end Model estimation
% % % -------------------------------------------------
% % % Contrast manager: (only run this step if runcontrast == 1)
% if run_contrast == 1
%     
% matlabbatch{3}.spm.stats.con.spmmat = {filename_path_contman};
% 
% % use loop to assign required variables (i.e. name, contrast) to matlabbatch array
% for connum = 1:size(cont_mat, 1)   % through each conditrast
%     matlabbatch{3}.spm.stats.con.consess{connum}.tcon.name = cont_name_str{connum};
%     matlabbatch{3}.spm.stats.con.consess{connum}.tcon.weights = cont_mat(connum,:);
%     matlabbatch{3}.spm.stats.con.consess{connum}.tcon.sessrep = 'none';
% end
% 
% matlabbatch{3}.spm.stats.con.delete = 1;
% 
% % % end Contrast manager
% % %-------------------------------------------------
% % % Results report
% 
% matlabbatch{4}.spm.stats.results.spmmat = {filename_path_resultsrep};
% 
% % use loop to assign required variables to matlabbatch array, applying
% % contrast to each run
% for i = 1:size(cont_mat, 1)
% matlabbatch{4}.spm.stats.results.conspec(i).titlestr = '';
% matlabbatch{4}.spm.stats.results.conspec(i).contrasts = i;
% matlabbatch{4}.spm.stats.results.conspec(i).threshdesc = 'FWE';
% matlabbatch{4}.spm.stats.results.conspec(i).thresh = 0.05;
% matlabbatch{4}.spm.stats.results.conspec(i).extent = 0;
% matlabbatch{4}.spm.stats.results.conspec(i).conjunction = 1;
% matlabbatch{4}.spm.stats.results.conspec(i).mask.none = 1;
% end
% 
% matlabbatch{4}.spm.stats.results.units = 1;
% matlabbatch{4}.spm.stats.results.export{1}.ps = true;
% 
% end   % end if runcontrast == 1
% 
% % % end Result report
% % % -------------------------------------------------

% % END PART (II): SPM script
% %-----------------------------------------------------------------------------

% calling spm function with run_spm is set to 1
if run_spm == 1
% % call spm_jobman to run the script above
spm_jobman('run', matlabbatch);
end


toc