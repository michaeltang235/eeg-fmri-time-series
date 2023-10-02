function ijk = mni2ijk(mnicoord, spmTmat)

% method 3: get (i,j,k) from (x,y,z) by using a matrix that is constructed 
% when general affine transformation took place from (i,j,k) to (x,y,z)

% get the required parameters for the transformation matrix from spmTmat
srow_x = spmTmat.raw.srow_x;
srow_y = spmTmat.raw.srow_y;
srow_z = spmTmat.raw.srow_z;

% assemble the trans. matrix with quantities obtained
srowmat = [srow_x(1) srow_x(2) srow_x(3);...
            srow_y(1) srow_y(2) srow_y(3);...
            srow_z(1) srow_z(2) srow_z(3)];

% from nifti documentation, (x,y,z) is related to (i,j,k) by a matrix
% equation in the form of Ax=B.

% by comparison, matrix A is the trans. matrix, while matrix B is the
% remaining terms in the equation given.
matb3 = transpose(mnicoord) - [srow_x(4); srow_y(4); srow_z(4)];

% solve the linear system using linsolve and store solution to ijk3
ijk3 = linsolve(srowmat, matb3);   % method 3 is used

% assign ijk3 to output variable ijk
ijk = ijk3;

end

%--------------------------------------------------------------------------------
% References:

% nifti documentation on coordinate transformation:
% https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html

% nifti documentation on rotation matrix:
% https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/quatern.html/view?searchterm=nifti1_io.c

% info on affine transformation (12 degrees of freedom):
%https://afni.nimh.nih.gov/sscc/staff/rwcox/ISMRM_2006/Syllabus%202006%20-%203340/files/B_05.pdf
