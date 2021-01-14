clear all 
close all
%---------------------------------------------------------------------------
% This script normalizes 3d anatomical image into standardized (MNI) space,
% using the preprocessing step of normalization provided by SPM

% enter info. about the image file before running the script
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter directory and filename of input file(s)
directname = '/Users/michaeltang/Downloads/fmri_project/sub41';
filename = '3Danat.nii';   % name of anatomical 3d image

% enter final voxel size in (mm)
voxelsize = [1 1 1];

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Preliminary:

% get full path of image file
filename_path = fullfile(directname, append(filename, ',1'));

% initialize spm module without graphical interface
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
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/Users/michaeltang/Downloads/spm12/tpm/TPM.nii'};
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