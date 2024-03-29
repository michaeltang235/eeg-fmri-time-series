clear all 
close all
%---------------------------------------------------------------------------
% This script normalizes 3d mask image into standardized (MNI) space,
% using the preprocessing step of normalization provided by SPM

% enter info. about the image file before running the script
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str, e.g. sub31, sub27, ...)
subnum = '14';

% enter path to directory where all input files are located
directname = ['/work/levan_lab/mtang/fmri_project/', 'sub', subnum];
% directname = ['C:\Users\siumichael.tang\Downloads\fmri_project\', 'sub', subnum];
% directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum, '_imthres0_exmask'];

% format filename of input mask file
filename_mask = 'EPI_bet_mask.nii';   % name of mask image

% enter final voxel size in (mm)
voxelsize = [3.75 3.75 3.75];

% enter full path to spm package
spm_path = '/home/siumichael.tang/spm12';
% spm_path = 'C:\Users\siumichael.tang\Downloads\spm12';

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Preliminary:

% get path of image that input mask is aligned to 
% (the file required to generate the deformation field),
% use EPI.nii stored under 'spm12/toolbox/OldNorm' as 'image to align to'
file_path_imgtoalignto = [spm_path, filesep, 'toolbox', filesep, ...
    'OldNorm', filesep, 'EPI.nii'];

% get full path of mask image
[~, file_path_mask, ~] = get_path(directname, filename_mask);

% format path to tissue probability map (tpm image) stored in SPM package
tpm_path = [spm_path, filesep, 'tpm', filesep, 'TPM.nii'];

% add path of spm to current working directory
addpath(spm_path);

% % initialize spm module without graphical interface
spm('defaults','fmri');
spm_jobman('initcfg');

% END Preliminary
%---------------------------------------------------------------------------

% FROM SPM:
% Normalization: Estimate and write: only step required for mask image
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {file_path_imgtoalignto};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {file_path_mask{:}};
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
