% FUNCTION RESIDUALS = FILTER_MOTION(MOTION_FILE_PATH, SIGNAL)
% MOTION_FILE_PATH is path of motion parameters file
% SIGNAL is array of signal required to be motion-filtered, dim. of signal
% array is ni X nj X nk X nt, where ni, nj, nk denote dim. in x, y, and z
% directions, and nt denotes dim. in time domain.
% and outputs
% RESIDUALS, the time series of signal with motion artifacts removed, with
% dim. the same as SIGNAL

function residuals = filter_motion(motion_file_path, signal)

%---------------------------------------------------------------------------
% Part (I): read motion parameters from input file

% open motion file with read permission assigned only, 
file_id = fopen(motion_file_path, 'r');

% read data from the file opened, then convert cells to mat
data_cell = textscan(file_id, '%f %f %f %f %f %f');
mo_para = cell2mat(data_cell);

% close motion file opened
fclose(file_id);

% END Part (I): read motion parameters from input file
%---------------------------------------------------------------------------
% Part (II): generate 24 regressors + constant term from those 6 motion time series

% The 24 regressors are the 6 motion time series, 
% those same 6 time series shifted by one scan, 
% and the squares of the previous 12 regressors

% create 6 time series shifted by one scan, horiz. conca. with
% motion para. read to form reg_shifted 
mo_shifted = [repmat([0], [1, 6]); mo_para(1:end -1, :)];   % motion para. shifted by one scan
reg_shifted = [mo_para mo_shifted];   % combine motion para. read with shifted para.

% square the 12 regressors and add them all toegther to form reg_sq
mo_sq = reg_shifted.^2;   % square of the previous 12 motion regressors
reg_sq = [reg_shifted mo_sq];   % add mo_sq to list of regressors

% create constant term to form all regressors required
mo_constant = ones(size(reg_sq, 1), 1);   % constant term
reg_mo = [reg_sq mo_constant];   % assemble all motion regressors

% get dimensions of input signals
ni = size(signal, 1);   % dim. in x
nj = size(signal, 2);   % dim. in y
nk = size(signal, 3);   % dim. in z
nt = size(signal, 4);   % dim. in time 

% initialize residuals array for storing motion-filtered signals
residuals = zeros(ni, ni, nk, nt);

% loop through input signal array, perform motion filtering in each voxel
for i = 1:ni   % for each voxel in x-dim.
    for j = 1:nj   % for each voxel in y-dim.
        for k = 1:nk   % for each voxel in z-dim.
            % extract time series required, squeeze array to remove
            % dim. of length 1
            timeseries = squeeze(signal(i, j, k, :));
            
            % solve system of eq. (Ax = b) for beta, where A is the matrix
            % of motion regressors and b is the time series
            beta = reg_mo\timeseries;
            
            % comput the residuals by subtracting fitted regressors
            % i.e. time series - motion regressor * beta (matrix multip.)
            residuals(i,j,k, :) = timeseries - reg_mo*beta;
        end
    end
end

% Part (II): generate 24 regressors + constant term from those 6 motion time series
%---------------------------------------------------------------------------

end   % end function residuals = filter_motion(...)









% end
