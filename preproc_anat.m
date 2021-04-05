clear all 
close all
%---------------------------------------------------------------------------
% This script normalizes 3d anatomical image into standardized (MNI) space,
% using the preprocessing step of normalization provided by SPM

% enter info. about the image file before running the script
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str, e.g. sub31, sub27, ...)
subnum = '14';

% enter path to directory where all input files are located
directname = ['C:\Users\siumichael.tang\Downloads\fmri_project\', 'sub', subnum];
% directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum, '_imthres0_exmask'];
filename = '3danat.nii';   % name of anatomical 3d image

% enter final voxel size in (mm)
voxelsize = [1 1 1];

% enter full path to spm package
spm_path = 'C:\Users\siumichael.tang\Downloads\spm12';

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Preliminary:

% get full path of image file
filename_path = fullfile(directname, append(filename, ',1'));

% format path to tissue probability map (tpm image) stored in SPM package
tpm_path = [spm_path, filesep, 'tpm', filesep, 'TPM.nii'];

% % initialize spm module without graphical interface
spm('defaults','fmri');
spm_jobman('initcfg');

% END Preliminary
%---------------------------------------------------------------------------

% FROM SPM:
% Normalization: Estimate and write: only step required for anat. image
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {filename_path};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {filename_path};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {tpm_path};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = voxelsize;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

% end Normalization: Estimate and write:
%---------------------------------------------------------------------------

% call spm_jobman to run the script above
spm_jobman('run', matlabbatch);