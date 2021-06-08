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

% format filenames of processed fmri images, explicit mask, and electrode
filename_swraimg = ['swra*.nii'];   % processed func. images
filename_expmask = 'wEPI_bet_mask.nii';   % explicit mask
filename_elect =  [subnum, '_*Koordinaten*.xlsx'];   % file containing mni coord. of all electrode pairs
filename_spikes = ['subject', subnum, '_rates.txt'];   % file of spike rates

% format path at which clinically determined reference channels located
clin_ref_ch_directname = '/work/levan_lab/mtang/fmri_project';   % directory where the file is
clin_ref_ch_filename = 'ref_channel_clin_det.txt';   % filename 

% enter path where ouput struct. is stored at
fname_op = [directname, filesep, 'matrices' filesep 'rho_spikes'];   % direct. of output matrix
filename_op = 'rho_spikes.mat';   % filename of output file

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%----------------------------------
% locate paths of func. images and their details 

% use get_path function to obtain full path of input files
[~, swra_file_path, ~] = get_path(directname, filename_swraimg);   % swra func. images
[~, expmask_file_path, ~] = get_path(directname, filename_expmask);   % normalized mask images
[~, elect_file_path, ~] = get_path(directname, filename_elect);   % file containing info. about electrodes coordinates
[~, spikes_file_path] = get_path(directname, filename_spikes);   % spike rates

 % get full path of file of clinically determined reference channels
clin_ref_ch_file_path = fullfile(clin_ref_ch_directname, clin_ref_ch_filename);  

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
input_array.clin_ref_ch_file_path = clin_ref_ch_file_path;   % path to file of clin. det. channels
input_array.subnum = str2num(subnum);   % subject number in format of integer

% get indices of current run number embedded in processed image filename
sind = regexp(lower(swra_file_path{run_ind}), 'run');   % start index
eind = regexp(lower(swra_file_path{run_ind}), '.nii');   % end index
sess_num = swra_file_path{run_ind}(sind+3:eind-1);   % get indices of run number

% create fieldname representing current session number
fieldname = sprintf('run%s', sess_num);

% store output calcu. under current fieldname of terms
terms.(fieldname) = get_rho_spikes(input_array);

% access tables of current session created by function get_rho_spikes
table_rho_clin_ch = terms.(fieldname).table_rho_clin_ch;   % table of clinically det. channels
table_rho_emp_ch = terms.(fieldname).table_rho_emp_ch;   % table of empirically det. channels

% format filenames and full paths of tables
filename_clin_table = sprintf('table_rho_clin_ch_run_%s.csv', sess_num);
filename_emp_table = sprintf('table_rho_emp_ch_run_%s.csv', sess_num);
clin_table_path = fullfile(fname_op, filename_clin_table);   % path of table of clin. det. channels
emp_table_path = fullfile(fname_op, filename_emp_table);   % path of table of emp. det. channels

% output tables created in current session
if op_results == 1
    writetable(table_rho_clin_ch, clin_table_path);   % table of clin. det. ch.
    writetable(table_rho_emp_ch, emp_table_path);   % table of emp. det. ch.
end

end   % end for run_ind = 1:numel(swra_file_path)

%---------------------------------------------------------------------------
% output struct. calcu. to path after looping through all sessions 
if op_results == 1
    save(fullfile(fname_op, filename_op), 'terms');
end

%---------------------------------------------------------------------------

toc

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% define function get_rho_spikes requiring an input array, which
% has the fields defined below,
% SWRA_IMG is the data of processed image file
% SWRA_INFO is the info (headers) associated with the processed image file
% EXPMASK is the explicit mask image
% ELECTAR is the array containing info. about all electrodes used,
% CLIN_REF_CH_FILE_PATH is the path of the file of clinically determined
% channels,
% SUBNUM is the subject number in integer form
% and outputs,
% OPSTRUCT, a structure that contains info. about coordinates of
% electrdoes, signals within 3X3X3 boxes centered at midpoint between
% electrode pairs, corr. coeff. btw. ref. and target channels, etc 
% see end of function for details of variables included in this struct.

function opstruct = get_rho_spikes(input_array)

%---------------------------------------------------------------------------
% Prelim: access variables defined in input array
swra_img = input_array.swra_img;   % processed func. images
swra_info = input_array.swra_info;   % headers of processed func. images
expmask = input_array.expmask;   % explicit mask
electar = input_array.electar;   % array of electrode coordinates
spikesar = input_array.spikesar;   % array of spike rates
clin_ref_ch_file_path = input_array.clin_ref_ch_file_path;   % path of clin. det. channels
subnum = input_array.subnum;   % subject number in format of integer

% end Prelim: access variables defined in input array
%---------------------------------------------------------------------------
% Part (I): sort electrodes by their types

% use sort_elect to sort electrode array by their types
ele_sorted = sort_elect(electar);

% % get number of electrode types identified
% ele_count = numel(ele_sorted);

% END PART (I): sorte electrodes by their types
%---------------------------------------------------------------------------
% % Part (II): find midpoint of each pair of electrode contacts in both mni
% and image space

% use get_midpt_elect to locate midpoint of all electrode channels
% col. 1 = name of channel, 
% col. 2 = midpt in mni space, 
% col. 3 = midpt in image space
midpt = get_midpt_elect(ele_sorted, swra_info);

% END Part (II): find midpoint of each pair of electrode contacts in both mni
% and image space
%---------------------------------------------------------------------------

% Part (III): construct 3X3X3 boxes centered at midpoint of each electrode
% pair to store signals enclosed thereof

% initialize cell array named sig_box to store signals of voxels within
% 3X3X3 boxes, 
% col. 1: name tags of electrode pairs
% col. 2: signal within the 3X3X3 box of each electrode pair
% col. 3: mean signal within the 3X3X3 box of each electrode pair
sig_box = cell(1, numel(midpt));

% use for loop to get signals within each box from processed (swra) func. images
for i = 1:numel(midpt)   % for each electrode type
    for j = 1:size(midpt{i}, 1)   % for each electrode pair under current type
        sig_box{i}{j, 1} = midpt{i}{j, 1};   % get name tag of electrode pair
        mxyz = midpt{i}{j, 3};   % get midpoint in image space
        
        % get signals from swra images, convert from int to double
        signal =  double(swra_img(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
            mxyz(3)-1:mxyz(3)+1, :));   
        sig_box{i}{j, 2} = signal;  % assign signal located to col. 2 of cell array
        
        % create regional mask from explicit mask using the coordinates of
        % the 3X3X3 box
        reg_mask = double(expmask(mxyz(1)-1:mxyz(1)+1, mxyz(2)-1:mxyz(2)+1, ...
            mxyz(3)-1:mxyz(3)+1));
        num_nonzero = nnz(reg_mask);   % get number of non-zero entries in regional mask

        % get avg. signal within the current box
        avg_signal = zeros(1, size(signal, 4));   % initialize avg. signal array
        for t = 1:size(signal, 4)   % for each time step
            % get sum of singals in box at current time, then get average
            % signal by dividing the sum by number of non-zero entries in
            % in the regional mask
            % could have used mean(nonzeros), but this method is adopted just 
            % in case if the signal is zero but the mask isn't at a voxel
            avg_signal(t) = sum(sum(sum(reg_mask.*(signal(:, :, :, t)))))/num_nonzero;
        end     
        
        % assign avg. signal calculated to col. 3 of cell array
        sig_box{i}{j, 3} = avg_signal;
    end
end

% END Part (III): construct 3X3X3 boxes centered at midpoint of each electrode
% pair to store signal enclosed thereof
%---------------------------------------------------------------------------
% Part (IV): group signals of all channels in single layer

% group all electrode types in one single layer of cell array
% col. 1 = id assigned to electrode pair
% col. 2 = name tag of electrode pair
% col. 3 = signals within 3X3X3 box constructued
% col. 4 = mean signal within the 3X3X3 box
sig_box_all = {};   % initialize sig_box_all array
rownum = 1;   % initialize row number as 1 

for i = 1:numel(sig_box)   % for each electrode type in sig_box
    for j = 1:size(sig_box{i}, 1)   % for each pair of electrode contacts
        sig_box_all(rownum, 1) = {rownum};   % assign row number as id 
        sig_box_all(rownum, 2:2+size(sig_box{i}, 2)-1) = sig_box{i}(j, :);   % assign relevant info.
        rownum = rownum + 1;   % increment row number by 1
    end
end

% END Part (IV): group signals of all channels in single layer
%---------------------------------------------------------------------------
% Part (V): match channels in spike rates array with that of signal box

% use clean_match_spikes to get spikes struct.
spikes = clean_match_spikes(spikesar, sig_box, sig_box_all);

% get spike rates array in single layer with channels matched with sig_box_all
spikes_all = spikes.spikes_all;

% END Part (V): match channels in spike rates array with that of signal box
%---------------------------------------------------------------------------
% Part (VI): get correlation coefficient between electrode channels

% with ids assigned to pairs of electrode contacts, get correlation
% coefficients between each pair, store data in corr_array, with
% col. 1 = id of reference electrode pair, 
% col. 2 = id of target electrode pair, 
% col. 3 = name of reference electrode pair
% col. 4 = name of tagret electrode pair
% col. 5 = pearson's correlation coefficient between ref. and target pairs

% get total number of pairs of ele. contacts
num_contact_pair = size(sig_box_all, 1);

% initialize corr_array 
corr_array = cell(num_contact_pair^2, 5);
row = 1;   % initialize row number as 1

for i = 1:size(sig_box_all, 1)   % for each row in sig_box_all
    for j = 1:size(sig_box_all, 1)   % again, for each row in sig_box_all
        corr_array{row, 1} = sig_box_all{i, 1};   % id of ref. ele. pair
        corr_array{row, 2} = sig_box_all{j, 1};   % id of target ele. pair
        corr_array{row, 3} = sig_box_all{i, 2};   % name of ref. ele. pair
        corr_array{row, 4} = sig_box_all{j, 2};   % name of target ele. pair
        ref_sig = [sig_box_all{i, 4}]';   % ref. signal
        tar_sig = [sig_box_all{j, 4}]';   % target signal
        % corr. coef. btw. ref. and targ. sig.
        corr_array{row, 5} = corr(ref_sig, tar_sig, 'Type', 'Pearson');  
        row = row + 1;   % increment row number by 1
    end
end

% END Part (VI): get correlation coefficient between electrode channels
%---------------------------------------------------------------------------
% Part (VII): build correlation matrix (with channels on each axis)

% get number of channels 
num_ch = size(sig_box_all, 1);

% initialize corr_mat (correlation matrix)
corr_mat = zeros(num_ch);

% with ids assigned, scan across rows in corr_array, extract relevant info.
for row = 1:size(corr_array, 1)
    corr_mat(corr_array{row, 1}, corr_array{row, 2}) = corr_array{row, 5};
end

% set up list of thresholds for filtering entries 
corr_mat_thresh = [0.6:0.1:0.9];

% filter entries in corr. matrix by each threshold defined above, 
% store filtered matrix in struct. named cor_mat_filt
corr_mat_filt = struct;

% loop through each entry in the list of thresh. defined for corr. matrix, 
% replace entries below thresh. with value of 0, then store filtered corr.
% matrix in cor_mat_filt struct. under the fieldname thresh_X, where X
% stands for current threshold.
for item = 1:length(corr_mat_thresh)   % for each item in the list of thresholds
    threshold = corr_mat_thresh(item);   % get current threshold
    corr_mat_current = corr_mat;   % set corr_mat_current as corr_mat initially
    for i = 1:size(corr_mat_current, 1)   % for each row in corr. mat.
        for j = 1:size(corr_mat_current, 2)   % for each col. in corr. mat.
            if abs(corr_mat_current(i, j)) < threshold   % if abs. value of entry is less than thresh.
                corr_mat_current(i, j) = 0;   % set it as zero
            end
        end
    end
    fldname = sprintf('thresh_%.1f', threshold);   % build fieldname using current thresh.
    fldname_str = strrep(fldname, '.', '_');   % replace '.' by '' to avoid confusion 
    corr_mat_filt.(fldname_str) = corr_mat_current;   % store current filtered corr. mat. 
end

% END Part (VII): build correlation matrix (with channels on each axis)
%---------------------------------------------------------------------------
% Part (VIII): grouped corr. coeff. and spike rates by ref. elect. channel

% exclude entries of corr. coeff. btw. 
% the same ref. and targ. electrode (e.g. {1, 1}, {2, 2}, ...)

% get list of ids and channel names
ids_all = cell2mat(sig_box_all(:, 1));

% sort corr. coeff. by electrode channels
cor_ref = {};   % initialize array
for ref_id = 1:length(ids_all)   % for each channel
    cor_ref{ref_id} = {};   % initialize cell for current channel
    for i = 1:size(corr_array, 1)   % for each row in corr_array
        % if id of ref. ch. (col. 1) agrees with current id, 
        % and if id of target ch. (col. 2) isn't the same as ref. ch.
        % assign relevant info. to cell of current ref. channel
        if (corr_array{i, 1} == ids_all(ref_id)) && (corr_array{i, 1} ~= corr_array{i, 2})
            cor_ref{ref_id} = [cor_ref{ref_id}; corr_array(i, :)];   % concatenate current cell
        end
    end
end

% get spike rates for all target elect. under each ref. elect. channel 
% in cor_ref,
% names of channels in both cor_ref and spikes_ref are in the same order
spikes_ref = {};    % initialize array
for ref_id = 1:length(ids_all)   % for each ref. elect. channel
    spikes_ref{ref_id} = {};   % initialize cell for current channel
    for i = 1:size(cor_ref{ref_id}, 1)   % for each row in current channel of cor_ref
        for j = 1:size(spikes_all, 1)   % for each row in spikes_all array
            % compare channel names in the two arrays, if name in current
            % row in cor_ref agrees with that in spikes_all, assign
            % relevant info. to current cell of spikes_ref by vertical
            % concatenation
            if strcmp(cor_ref{ref_id}{i, 4}, spikes_all{j, 2})
                spikes_ref{ref_id} = [spikes_ref{ref_id}; spikes_all(j, :)];  
                break;   % break loop if spike rates are found for current targ. ch.
            end
        end
    end
end

% Part (VIII): grouped corr. coeff. and spike rates by ref. elect. channel
%---------------------------------------------------------------------------
% Part (IX): get correlations between corr. coeff. and spikes rates, then
% filter entries 

% for each ref. channel, we have obtained in the following
% (i) corr. coeff. btw. ref. and target channels 
% (ii) spike rates of all target channels, 
% then get corr. coeff. between (i) and (ii)

% get corr. coeff. (pearson and spearman) between corr coeff. btw. channels and spikes rates
rho_cs_pe = {};   % pearson corr.
pval_cs_pe = {};   % p-value of pearson corr.
rho_cs_sp = {};   % spearman corr.
pval_cs_sp = {};   % p-value of spearman corr.

for ref_id = 1:numel(cor_ref)   % for each ref. channel
    % get pearson and spearman corr. coeff. and their p-values btw. 
    % (i) and (ii)
    [rho_cs_pe{ref_id}, pval_cs_pe{ref_id}] = corr([cor_ref{ref_id}{:, 5}]', ...
        [spikes_ref{ref_id}{:, 3}]', 'Type', 'Pearson');
    [rho_cs_sp{ref_id}, pval_cs_sp{ref_id}] = corr([cor_ref{ref_id}{:, 5}]', ...
        [spikes_ref{ref_id}{:, 3}]', 'Type', 'Spearman');
end

% get entries in rho_cs_sp (Spearman) with values greater than or equal 
% to the thresh. 
rho_cs_sp_filt = {};   % initialize filtered spearman rho array
pval_cs_sp_filt = {};   % initialize filtered p-values array
thresh_rho_cs_sp = 0.3;   % set thresh. for filtering entries 
row_ind = 1;   % intialize row index as 1
for ref_id = 1:numel(rho_cs_sp)   % for each entry in rho_cs_sp
    if abs(rho_cs_sp{ref_id}) >= thresh_rho_cs_sp   % if current entry exceeds thresh.
        rho_cs_sp_filt{row_ind, 1} = ids_all(ref_id);   % get id of ref. ch.
        rho_cs_sp_filt{row_ind, 2} = rho_cs_sp{ref_id};   % get entry in rho_cs_sp
        
        pval_cs_sp_filt{row_ind, 1} = ids_all(ref_id);   % get id of ref. ch.
        pval_cs_sp_filt{row_ind, 2} = pval_cs_sp{ref_id};   % get entry in pval_cs_sp
        row_ind = row_ind + 1;   % increment row index by 1 
    end
end

% END Part (IX): get correlations between corr. coeff. and spikes rates, then
% filter entries
%---------------------------------------------------------------------------
% Part (X): get spearman's corr. coeff. (both same and other types) btw.
% filtered ref. ch. and their corresponding target ch., and spike rates

% get ids for each electrode type
ids_type = cell(1, numel(sig_box));   % initialize cell array
row_i = 1;   % initialize row_id (initial row number)
for i = 1:numel(sig_box)   % for each type in sig_box
    rowspan = size(sig_box{i}, 1);   % get rowspan of current electrode type
    % get ids of current electrode type
    ids_type{i} = [sig_box_all{row_i:(row_i + rowspan - 1), 1}];   
    row_i = row_i + rowspan;   % adjust row_i for next electrode type
end

% initialize structure for storing data of each ref. channel with spearman
% corr. coeff. exceeding the thresh.
ch_filt = struct;

% for each ref. channel with spearman corr. coeff. exceeding the thresh., 
% get target spike rates and corr. coeff. with all target channels in 
% current electrode type, 
% do the same for electrodes of other types
% store output in struct.

for item = 1:size(rho_cs_sp_filt, 1)   % for each ref. channel in filtered array
    
    ref_id = rho_cs_sp_filt{item, 1};   % get current ch. id.
    
    % get current electrode type of ref_id
    for i = 1:numel(ids_type)   % for each id type
        if ismember(ref_id, ids_type{i})   % if ref_id belongs to the type
            type = i;   % get type of ref_id
        end
    end
       
    spikes_same_type = {};   % initialize arrays
    spikes_other_type = {};
    cor_same_type = {};
    cor_other_type = {};
    row_ind_same_type = 1;   % initialize row number for arrays of same and other types
    row_ind_other_type = 1;
    
    % for each row in spikes_ref (spike rates sorted by ref. ch.)
    for i = 1:size(spikes_ref{ref_id}, 1)  
        % if id of targ. ch. is a memeber of current elect. type 
        if ismember(spikes_ref{ref_id}{i, 1}, ids_type{type})
            
            % store ids and spike rates of same type to array
            spikes_same_type{row_ind_same_type, 1} = spikes_ref{ref_id}{i, 1};   % id of targ. ch.
            spikes_same_type{row_ind_same_type, 2} = spikes_ref{ref_id}{i, 3};   % spike rate of targ. ch.
            
            % store corr. coeff. of same type to array
            cor_same_type{row_ind_same_type, 1} = cor_ref{ref_id}{i, 2};   % id of targ. ch.
            cor_same_type{row_ind_same_type, 2} = cor_ref{ref_id}{i, 5};   % corr. coeff. of targ. ch.
            row_ind_same_type = row_ind_same_type + 1;   % increment row index by 1
        else
            % else (ids do not belong the type of ref. ch.), 
            % store equi. info. of other types to array
            spikes_other_type{row_ind_other_type, 1} = spikes_ref{ref_id}{i, 1};   % id of targ. ch.
            spikes_other_type{row_ind_other_type, 2} = spikes_ref{ref_id}{i, 3};   % spike rate of targ. ch.
            
            cor_other_type{row_ind_other_type, 1} = cor_ref{ref_id}{i, 2};   % id of targ. ch.
            cor_other_type{row_ind_other_type, 2} = cor_ref{ref_id}{i, 5};   % corr. coeff. of targ. ch.
            row_ind_other_type = row_ind_other_type + 1;   % increment row index by 1
        end
    end
    
    % get spearman corr. coeff. between cor. and spikes that are of same type,
    % and of different types
    [cor_sp_same_type, pval_sp_same_type] = corr(cell2mat(cor_same_type(:, 2)), ...
        cell2mat(spikes_same_type(:, 2)), 'Type', 'Spearman');
    [cor_sp_other_type, pval_sp_other_type] = corr(cell2mat(cor_other_type(:, 2)), ...
        cell2mat(spikes_other_type(:, 2)), 'Type', 'Spearman');
    
    if isnan(cor_sp_same_type)
        cor_sp_same_type = 0;
        pval_sp_same_type = 0;
    end
    if isnan(cor_sp_other_type)
        cor_sp_other_type = 0;
        pval_sp_other_type = 0;
    end
    
    % sort cor_sp_other_type and spikes_sp_other_type by types
    cor_ot_sorted = {};   % initialize arrays
    spikes_ot_sorted = {};
    entry_id = 1;   % intialize cell index
    row_req = 1;   % initialize row number
    for i = 1:numel(ids_type)   % for each type of electrodes
        for j = 1:size(cor_other_type, 1)   % for each row in cor_other_type
            % if id of current row belongs to curren type, extra into
            if ismember(cor_other_type{j, 1}, ids_type{i})
                cor_ot_sorted{entry_id}{row_req, 1} = cor_other_type{j, 1};   % id of targ. ch.
                cor_ot_sorted{entry_id}{row_req, 2} = cor_other_type{j, 2};   % corr. coeff. btw. ref. and targ. ch.
                spikes_ot_sorted{entry_id}{row_req, 1} = spikes_other_type{j, 1};   % id of targ. ch.
                spikes_ot_sorted{entry_id}{row_req, 2} = spikes_other_type{j, 2};   % spike rates of targ. ch.
                row_req = row_req + 1;   % increment row number by 1
            end
        end
        if row_req > 1   % if current type is found
        entry_id = entry_id + 1;
        end
        row_req = 1;   % reset row_req to 1 for next type
    end
    
    % get name of ref. and target channel types
    % get end index in ch. name before any numeric and hyphen
    alph_eind_ref_ch = regexp(sig_box_all{ref_id, 2}, '\d*-');   
    ch_same_type_name = sig_box_all{ref_id, 2}(1:alph_eind_ref_ch - 1);   % select required places in str
    
    ch_other_type_name = cell(numel(spikes_ot_sorted), 1);   % initialize array for storing targ. ch. name 
    for i = 1:numel(spikes_ot_sorted)   % for each type of targ. channel
        id_targ = spikes_ot_sorted{i}{1, 1};   % get id of first entry in current type
        str_ext = sig_box_all{id_targ, 2};   % get name tag of targ. ch. using id
        alph_eind = regexp(str_ext, '\d*-');   % get end index in name tag before any numeric and hyphen
        ch_other_type_name{i} = str_ext(1:alph_eind - 1);   % select required places in str
    end
        
    % store quantities in ch_filt struct.
    fldname = sprintf('ref_id_%d', ref_id);   % get current fieldname 
    ch_filt.(fldname).ch_same_type_name = ch_same_type_name;   % name of ref. ch. type
    ch_filt.(fldname).ch_other_type_name = ch_other_type_name;   % names of targ. ch. types
    ch_filt.(fldname).spikes_same_type = spikes_same_type;   % spike rates of same type (A)
    ch_filt.(fldname).cor_same_type = cor_same_type;   % corr. coff. btw. ref. and targ. ch. of same type (B)
    ch_filt.(fldname).spikes_other_type = spikes_other_type;   % spike rates of other types (C)
    ch_filt.(fldname).cor_other_type = cor_other_type;   % corr. coff. btw. ref. and targ. ch. of other types (D)
    ch_filt.(fldname).cor_sp_same_type = cor_sp_same_type;   % spear corr. btw. (A) and (B)
    ch_filt.(fldname).pval_sp_same_type = pval_sp_same_type;   % p-value of spear. corr. btw. (A) and (B)
    ch_filt.(fldname).cor_sp_other_type = cor_sp_other_type;   % spear corr. btw. (C) and (D)
    ch_filt.(fldname).pval_sp_other_type = pval_sp_other_type;   % p-value of spear. corr. btw. (C) and (D)
    ch_filt.(fldname).cor_ot_sorted = cor_ot_sorted;   % spear corr. btw. (C) and (D), sorted by types
    ch_filt.(fldname).spikes_ot_sorted = spikes_ot_sorted;   % (C), sorted by types
    
end   % end for item = 1:size(rho_cs_sp_filt, 1)

% END Part (X): get spearman's corr. coeff. (both same and other types) btw.
% filtered ref. ch. and their corresponding target ch., and spike rates
%---------------------------------------------------------------------------

% Part (XI): get ranks of spearman's corr. coeff. and spike rates of
% filtered ref. channel

for item  = 1:size(rho_cs_sp_filt, 1)   % for each ref. ch. in filtered rho_cs_sp
    
    % get id of current ref. channel
    ref_id = rho_cs_sp_filt{item, 1};   % get id of current ref. ch.
    
    % get current electrode type of ref_id
    for i = 1:numel(ids_type)   % for each id type
        if ismember(ref_id, ids_type{i})   % if ref_id belongs to the type
            type = i;   % get type of ref_id
        end
    end
    
    % select spike rates and corr. based on ref_id
    spikes_col1 = spikes_ref{ref_id}(:,1);   % get targ. ch. ids
    spikes_col2 = spikes_ref{ref_id}(:,3);   % get targ. ch. spike rates
    spikes_selected = horzcat(spikes_col1, spikes_col2);   % concatenate them horizontally
    
    cor_col1 = cor_ref{ref_id}(:,2);   % get targ. ch. ids
    cor_col2 = cor_ref{ref_id}(:,5);   % get corr. coeff. btw. ref. and targ. ch.
    cor_selected = horzcat(cor_col1, cor_col2);   % concatenate them horizontally
    
    % sort both cell arrays (spikes/cor_selected) in ascending order, then get
    % corresponding indices
    [~, ind_spikes] = sort([spikes_selected{:, 2}]);   % sort spike rates
    [~, ind_cor] = sort([cor_selected{:, 2}]);   % sort corr. coeff. btw. ref. and targ. ch.
     
    % assign ranks to spikes and cor using indices obtained,
    % e.g. ind_spikes = [3, 1, 2], denotes that the first entry in sorted
    % array (ascending order) corresponds to the 3rd entry in the original
    % array
    rank_spikes = zeros(1, length(ind_spikes));   % initialize arrays
    rank_cor = zeros(1, length(ind_cor));
    for i = 1:length(ind_spikes)   % for each entry in index arrays
        rank_spikes(ind_spikes(i)) = i;   % assign ranks of spikes to entries in original array
        rank_cor(ind_cor(i)) = i;   % assign ranks of corr. coeff. to entries in original array
    end

    % get ranks of same type and other types for both arrays
    rank_spikes_same_type = {};   % initialize arrays
    rank_spikes_other_type = {};
    rank_cor_same_type = {};
    rank_cor_other_type = {};
    entry_same_type = 1;
    entry_other_type = 1;

    for i = 1:length(rank_spikes)   % for each element in rank arrays
        % if id of current spike rates belong to current electrode type
        if ismember(spikes_selected{i, 1}, ids_type{type})   
            rank_spikes_same_type{entry_same_type, 1} = spikes_selected{i, 1};   % assign id of spike rates
            rank_spikes_same_type{entry_same_type, 2} = rank_spikes(i);   % assign rank of spike rates
            rank_cor_same_type{entry_same_type, 1} = cor_selected{i, 1};   % assign id of corr. coeff.
            rank_cor_same_type{entry_same_type, 2} = rank_cor(i);   % assign rank of corr. coeff.
            entry_same_type = entry_same_type + 1;   % increment row number by 1
        else   % if id of current spike rates does not belong to current electrode type
            rank_spikes_other_type{entry_other_type, 1} = spikes_selected{i, 1};
            rank_spikes_other_type{entry_other_type, 2} = rank_spikes(i);
            rank_cor_other_type{entry_other_type, 1} = cor_selected{i, 1};
            rank_cor_other_type{entry_other_type, 2} = rank_cor(i);
            entry_other_type = entry_other_type + 1;
        end
    end
    
    % sort ranks of spike rates of other types by types
    rank_spikes_ot_sorted = {};   % initialize arrays
    rank_cor_ot_sorted = {};
    entry_id = 1;   % initialize cell index
    row_req = 1;   % intialize row number
    for i = 1:numel(ids_type)   % for each type of electrodes
        for j = 1:size(rank_cor_other_type, 1)   % for each row in rank_cor_other_type
            % if id of current row belongs to current ele. type
            if ismember(rank_cor_other_type{j, 1}, ids_type{i})   
                rank_cor_ot_sorted{entry_id}{row_req, 1} = rank_cor_other_type{j, 1};   % assign id
                rank_cor_ot_sorted{entry_id}{row_req, 2} = rank_cor_other_type{j, 2};   % assign rank of corr. coeff.
                rank_spikes_ot_sorted{entry_id}{row_req, 1} = rank_spikes_other_type{j, 1};   % assign id
                rank_spikes_ot_sorted{entry_id}{row_req, 2} = rank_spikes_other_type{j, 2};   % assign rank of spike rates
                row_req = row_req + 1;   % increment row number by 1
            end
        end
        if row_req > 1   % if current type is found
            entry_id = entry_id + 1;   % increment cell index by 1 for next type
        end
        row_req = 1;   % reset row_req to 1 for next type
    end
   
    % store quantities in ch_filt struct.
    fldname = sprintf('ref_id_%d', ref_id);   % get current fieldname 
    ch_filt.(fldname).rank_spikes_same_type = rank_spikes_same_type;   % ranks of spike rates of same type
    ch_filt.(fldname).rank_cor_same_type = rank_cor_same_type;   % ranks of corr. coeff. of same type
    ch_filt.(fldname).rank_spikes_other_type = rank_spikes_other_type;   % ranks of spike rates of other types
    ch_filt.(fldname).rank_cor_other_type = rank_cor_other_type;   % ranks of corr. coeff. of other types
    ch_filt.(fldname).rank_cor_ot_sorted = rank_cor_ot_sorted;   % ranks of corr. of other types, sorted by types
    ch_filt.(fldname).rank_spikes_ot_sorted = rank_spikes_ot_sorted;   % ranks of corr. of other types, sorted by types
    
end   % end for item = 1:size(rho_cs_sp_filt, 1)

% END Part (XI): get ranks of spearman's corr. coeff. and spike rates of
% filtered ref. channel
%---------------------------------------------------------------------------

% Part (XII): create tables listing corr. coeff. of ref. channels

% use get_clin_ref_ch to get clinically determined ref. channels and their ids 
ch_name_clinical = get_clin_ref_ch(subnum, clin_ref_ch_file_path, sig_box_all);

% initialize tables and assign values to them depending on ch_name_clinical
table_rho_clin_ch = {};   % table of clinically determined channels
table_rho_emp_ch = {};   % table of empirically determined channels

% create variable to indicate if table of corr. coeff. of clinically
% determined reference channels is to be updated (0 = no, 1 = yes)
make_clin_table = 0;   % initially assumed to be 0 

% get if ch_name_clinical array is empty, if so, disp message, otherwise,
% execute the lines below to create table of corr. coeff. of clinically
% determined reference channels
if isempty(ch_name_clinical)
    disp('warning, no channel names recorded in clinical file, no table is created')
else
    make_clin_table = 1;   % if not empty, change value of make_clin_table 
end

% get corr. coeff., p-value, and spike rates of each ref. channel listed 
% in ch_name_clinical, if array is not empty
if make_clin_table == 1

% initialize array, with format as follows, 
% 1st layer denotes type of electrode, 2nd layer denotes details of
% channels of current electrodet type, with 
% col. 1 and 2 listing id and name of channel, 
% col. 3 and 4 listing spearman corr. coeff. and p-value of channel, 
% col. 5 listing spike rates of channel
rho_clin_ch = ch_name_clinical;   
for i = 1:numel(ch_name_clinical)   % for each type of electrode
    for j = 1:size(ch_name_clinical{i}, 1)   % for each ch. in current type
        clin_ref_ch_id = ch_name_clinical{i}{j, 1};   % get channel id
        % use channel id to get corr. coeff. and p-val. in rho_cs_sp array
        rho_clin_ch{i}{j, 3} = rho_cs_sp{clin_ref_ch_id};   % spearman's rho
        rho_clin_ch{i}{j, 4} = pval_cs_sp{clin_ref_ch_id};   % p-value
        rho_clin_ch{i}{j, 5} = spikes_all{clin_ref_ch_id, 3};   % spike rate of ch.
    end
end

% convert rho_clin_ch into single layer array, round numeric entries to 3 d.p.
rho_clin_ch_all = {};   % initialize array
row_num_table_clin = 1;   % initialize row number 
for i = 1:numel(rho_clin_ch)   % for each type of electrode det. clinically
    for j = 1:size(rho_clin_ch{i}, 1)   % for each channel under cur. type
        rho_clin_ch_all{row_num_table_clin, 1} = rho_clin_ch{i}{j, 1};   % id of channel
        rho_clin_ch_all{row_num_table_clin, 2} = rho_clin_ch{i}{j, 2};   % name of channel
        rho_clin_ch_all{row_num_table_clin, 3} = round(rho_clin_ch{i}{j, 3}, 3);   % corr. coeff. 
        rho_clin_ch_all{row_num_table_clin, 4} = round(rho_clin_ch{i}{j, 4}, 3);   % p-value
        rho_clin_ch_all{row_num_table_clin, 5} = round(rho_clin_ch{i}{j, 5}, 3);   % spike rate
        row_num_table_clin = row_num_table_clin + 1;   % increment row number by 1
    end
end

% convert cell array to table and add table headers
table_rho_clin_ch = cell2table(rho_clin_ch_all, ...
    'VariableNames',{'id' 'channel name' 'corr. coeff.' 'p-value' 'spike rate'});

end   % end if make_clin_table == 1

% similarly, get corr. coeff. of ref. channel determined empirically 
% (those exceeding threshold)
rho_emp_ch_all = cell(size(rho_cs_sp_filt, 1), 5);   % initialize array
for i = 1:size(rho_emp_ch_all, 1)   % for each row in array
    emp_ch_id = rho_cs_sp_filt{i, 1};   % get id of channel from rho_cs_sp_filt  
    rho_emp_ch_all{i, 1} = rho_cs_sp_filt{i, 1};   % assign id to array
    % use id to get channel name from sig_box_all
    rho_emp_ch_all{i, 2} = sig_box_all{emp_ch_id, 2};
    % assign corr. coeff., p-val, spike rate (corr. to 3 d. p.)
    rho_emp_ch_all{i, 3} = round(rho_cs_sp_filt{i, 2}, 3);   % corr. coeff.
    rho_emp_ch_all{i, 4} = round(pval_cs_sp_filt{i, 2}, 3);   % p-val.
    rho_emp_ch_all{i, 5} = round(spikes_all{clin_ref_ch_id, 3}, 3);   % spike rate
end

% convert cell array to table and add table headers
table_rho_emp_ch = cell2table(rho_emp_ch_all, ...
    'VariableNames',{'id' 'channel name' 'corr. coeff.' 'p-value' 'spike rate'});

% END Part (XII): create tables listing corr. coeff. of ref. channels
%---------------------------------------------------------------------------

% Part (XIII): store quantities calcu. in output struct.

opstruct = struct;
opstruct.ele_sorted = ele_sorted;   % info. of electrode (sorted by name)
opstruct.midpt = midpt;   % coordinates of midpoint between each pair of electrode contacts 
opstruct.sig_box = sig_box;   % signals within 3X3X3 box centered at each midpoint
opstruct.sig_box_all = sig_box_all;   % signals within 3X3X3 box centered at each midpoint of all channels in single layer
opstruct.spikes_all = spikes_all;   % list of spike rates of all electrode pairs in single layer

opstruct.corr_array = corr_array;   % corr. coeff. (pearson's) btw. channels
opstruct.corr_mat = corr_mat;   % corr. matrix btw. channels
opstruct.corr_mat_filt = corr_mat_filt;   % (struct.) filtered corr. matrix with entries exceeding thresh.

opstruct.ids_all = ids_all;   % ids of all channels
opstruct.ids_type = ids_type;   % ids of all channels, sorted by types

opstruct.cor_ref = cor_ref;   % corr. coeff. (pearson's) btw. ref. (1st layer) and targ. (2nd layer) channels
opstruct.spikes_ref = spikes_ref;   % spike rates btw. ref. (1st layer) and targ. (2nd layer) channels

% with corr. coeff. (Y) and spike rates (X) btw. ref. and targ. ch.
% determined,
opstruct.rho_cs_pe = rho_cs_pe;   % pearson's corr. coeff. btw. X and Y
opstruct.pval_cs_pe = pval_cs_pe;   % p-value of pearson's corr. coeff.
opstruct.rho_cs_sp = rho_cs_sp;   % spearman's corr. coeff. btw. X and Y
opstruct.pval_cs_sp = pval_cs_sp;   % p-value of spearman's corr. coeff.
opstruct.rho_cs_sp_filt = rho_cs_sp_filt;   % for ref. ch. exceeding thresh., spearman's corr. coeff. btw. X and Y
opstruct.pval_cs_sp_filt = pval_cs_sp_filt;   % for ref. ch. exceeding thresh., p-values of spear. corr. coeff.

% for each ref. channel with spearman's corr. coeff. btw. X and Y exceeding
% the theshold, another struct. is created to store corr. coeff. and spike
% rates btw. ref. and targ. channels, with data sorted by electrode types, 
% ranks of corr. coeff. and spike rates are also given (see above for
% details)

% struct. for each ref. ch. with spearman's corr. coeff. btw. X and Y exceeding threshold
opstruct.ch_filt = ch_filt;   

% store tables created above in struct.
opstruct.table_rho_clin_ch = table_rho_clin_ch;   % table of clinically det. channels
opstruct.table_rho_emp_ch = table_rho_emp_ch;   % table of empirically det. channels

% END Part (XIII): store quantities calcu. in output struct.
%---------------------------------------------------------------------------



end   % end function
