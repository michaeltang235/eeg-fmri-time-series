clear all 
close all

% This script requires filenames of t-stats generated by fmristat, 
% coordinates of electrode contacts, their corresponding spike rates, 
% and explicit mask image as inputs, then return a struct. named terms, 
% which contains avg. non-zero t-values at volume enclosed by
% each pair of electrode contacts, and the associated spike rates.

% Remarks: file containing info. about electrode contacts must follow the
% format that electrode name and its x, y, and z coordinates are stored 
% in columns 2 to 5.

% Updates: Apr. 6, 2021, updated code for cleaning cell array read from 
% input file of electrode contacts, in that only include cells with both
% name tags and their coordinates present in electrode_cleaned array

tic 
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '14';

% enter path to directory where all input files are located
directname = ['C:\Users\siumichael.tang\Downloads\fmri_project\', 'sub', subnum];
% directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum, '_imthres0_exmask'];

% format filename to all t-statistics images generated by fmristat
% and that of explicit mask image
filename_tstat = ['*run_all_t.nii'];
filename_expmask = 'wEPI_bet_mask.nii';

% format path to directory and name(s) of input file(s)
% fname = [directname, filesep, 'matrices'];  % directory of file(s) below
filename_elect =  [subnum, '_*_Koordinaten.xlsx'];   % file containing mni coord. of all electrode pairs
filename_spikes = ['subject', subnum, '_rates.txt'];   % file of spike rates

% enter path where ouput file is stored at
fname_op = [directname, filesep, 'matrices'];   % direct. of output matrix
filename_op = 'tval_spikes_fmristat.mat';   % filename of output file

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% get number of contrasts interested and path of file containing info. about
% t-stats images, explicit mask, electrode contacts, and spike rates, repsectively
[numcont, tstat_file_path] = get_path(directname, filename_tstat);   % t-stat. images
[~, expmask_file_path] = get_path(directname, filename_expmask);   % explicit mask
[~, elect_file_path] = get_path(directname, filename_elect);   % electrode contacts
[~, spikes_file_path] = get_path(directname, filename_spikes);   % spike rates

% get avg. non-zero t-values and spike rates using function defined below,
% store output of function to cont cell array. 
% also, get names of contrasts
cont = cell(numcont, 1);   % initialize cont cell array
cont_name_str = cell(numcont, 1);   % initialize contrast seq. cell array
for i = 1:numcont
    [~, tstat_filename, ~] = fileparts(tstat_file_path{i});
    if ~isempty(regexp(tstat_filename, '_mag_t'))   % if t stat. images are from only one run
        eind = regexp(tstat_filename, '_mag_t');
    else
        eind = regexp(tstat_filename, '_t');   % if t stat. images are from multiple runs
    end        
    cont_name_str{i} = tstat_filename(1:eind - 1);   % get def. of contrast (str) from filename
    cont{i} = tvalues_spikes(elect_file_path, spikes_file_path, tstat_file_path{i}, expmask_file_path);
end

% store data of each contrast to terms struct.
for i = 1:numcont
    terms.contrast_seq = cont_name_str;   % describes the sequence of contrast fieldnames are referring to
    fieldname = sprintf('cont%d', i); 
    terms.(fieldname) = cont{i};   % quantities cal. by func. tvalues_spikes
end

%---------------------------------------------------------------------------
% output avg. t-values and spike rates to path specified
if op_results == 1
    save(fullfile(fname_op, filename_op), 'terms');
end

toc 

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% define function tvalues_spikes that requires the following inoputs
% (i) ELECT_FILE_PATH = path of file containing info. about electrodes
% (ii) SPIKESFILE_PATH = path of file containing info. about spike rates
% (iii) TSTAT_PART_FILE_PATH = path of file containing info. about
% t-statistics generated by fmristat
% (iv) EXPMASK_FILE_PATH = path of explicit mask 
% and output
% OPSTRUCT, a struct. that contains info. about coordinates of each
% type of electrode contacts, the associated t-values, avg. non-zero
% t-values, related spike rates, etc.

function opstruct = tvalues_spikes(elect_file_path, spikes_file_path, tstat_part_file_path, expmask_file_path)

% read input files
electar = readcell(elect_file_path{:});   % cell array containing info. about electrodes
spikesar = readcell(spikes_file_path{:});   % cell array containing info. about spike rates
tstat_info = niftiinfo(tstat_part_file_path);   % info about this particular t-stat. image
tstat_img = niftiread(tstat_part_file_path);   % particular t-stat. image of single data type
expmask = single(niftiread(expmask_file_path{:}));   % convert image data type to single

% apply explicit mask to t-stat. image 
tstat_img_mask = expmask.*tstat_img;

%---------------------------------------------------------------------------
% Part (I): SORT electrodes by their names

% file containing info. about electrodes may contains rows with value of 0
% assigned to missing name tags and coordinates,
% remove such rows and return resultant array to electar_cleaned

% initialize empty array named electar_cleaned
electar_cleaned = {};

% use loop to go through each row of electar array , if name tag (col. 2) 
% is NOT missing, assign info. of that row to electar_cleaned array
rownum = 1;   % initialize row number as 1
for i = 1:size(electar, 1)   % for each row in electar
    % if name tag is found or if coordinates are not missing
    % ismissing finds missing entries in array, have to separtely evaluate
    % col. 2 and col. 3:5, as they are of different data types, 
    % any means if any of the entries is non-zero or logical 1
    if any(~ismissing(electar{i, 2})) && any(~ismissing([electar{i, 3:5}]))
        electar_cleaned(rownum, :) = electar(i,:);   % assign relevant info. to array
        rownum = rownum + 1;   % increment row number by 1
    end
end

% initialize cell array named ele_sorted with first entry having info. about
% the first electrode in the cleaned electrode array (electar_cleaned)
ele_sorted{1} = electar_cleaned(2, 2:5);

% initialize electrode count
ele_count = 1;

% each electrode name contains a compo. of string and a compo. of number
% e.g. 'RPF1' has string compo. of 'RPF' and numeric compo. of '1'
% use for loop to go through each entry in the name column, 
% find max. index of string in name tag and number associated with name,
% compare name in current name tag with the next one, 
% if names match and number of the next one is larger than the current
% one, store the info. within current electrode, 
% else create another cell with new name tag

for i = 2:size(electar_cleaned, 1)-1
    ele_name = electar_cleaned{i, 2};   % name of current electrode 
    ele_name1 = electar_cleaned{i+1, 2};   % name of next electrode 
    
    if length(ele_name) <= length(ele_name1)
        maxind = max(regexp(ele_name, '[A-Z]'));   % find max. index of string in name
        if strcmp(ele_name(1:maxind), ele_name1(1:maxind)) && ...
                str2double(ele_name(maxind+1:end)) < str2double(ele_name1(maxind+1:end))   % if names match and num. agree
            ele_sorted{ele_count} = [ele_sorted{ele_count}; electar_cleaned(i+1,2:5)];   % add info. of next one to current type
        else
            ele_count = ele_count + 1;   % increment electrode count by 1
            ele_sorted{ele_count} = electar_cleaned(i+1,2:5);   % add info. to new cell
        end
    else
        ele_count = ele_count + 1;   % increment electrode count by 1
        ele_sorted{ele_count} = electar_cleaned(i+1,2:5);   % add info. to new cell
    end
end

% print msg. after sorting
sprintf('electrodes sorted, %d of them found', ele_count)

% END Part (I): Sort electrodes
%--------------------------------------------------------------------------- 
    
% Part (II): find midpoint of each pair of electrode contacts

% initialize midpt cell array
midpt = {};

% use loop to go through each entry in ele_sorted
% get midpoint of each pair of ele. contacts and store value to array
for i = 1:ele_count
    for j = 1:size(ele_sorted{i}, 1) - 1
        midpt{i}{j, 1} = [ele_sorted{i}{j, 1} '-' ele_sorted{i}{j+1, 1}];   % string showing the ele. pair
        mx = 0.5*(ele_sorted{i}{j, 2} + ele_sorted{i}{j+1, 2});   % midpt. in x
        my = 0.5*(ele_sorted{i}{j, 3} + ele_sorted{i}{j+1, 3});   % midpt. in y
        mz = 0.5*(ele_sorted{i}{j, 4} + ele_sorted{i}{j+1, 4});   % midpt. in z
        midpt{i}{j, 2} = [mx, my, mz];   % assign value to array
    end
end
    
% END Part (II): find midpoint
%--------------------------------------------------------------------------- 
% Part (III): convert coordinates from mni space to voxel space

for i = 1:ele_count   % each electrode
    for j = 1:size(midpt{i}, 1)   % each pair of electrode contacts
        midpt{i}{j, 3} = round(mni2ijk(midpt{i}{j, 2}, tstat_info));   % convert from mni to voxel space
    end
end

% END Part (III): convert coordinates from mni space to voxel space
%--------------------------------------------------------------------------- 

% Part (IV): get average T-value of 3X3X3 box centered at midpoint of each
% pair of electrode contacts

% initialize tbox array for storing avg. t-value of each pair of electrode contacts
tbox = {};

% get name of each pair of electrode contacts and get coordinates of
% midpoint, get t-values of each entry in the resultant 3X3X3 box, 
% then, get average of non-zero t-values in the box
for i = 1:ele_count
    for j = 1:size(midpt{i}, 1)
        tbox{i}{j, 1} = midpt{i}{j,1};   % get name of each pair of electrode contacts
        mxyz = midpt{i}{j, 3};   % get coordinates of midpoint 
        tbox{i}{j, 2} = tstat_img_mask(mxyz(1)-1:mxyz+1, mxyz(2)-1:mxyz(2)+1, mxyz(3)-1:mxyz(3)+1);   % get t-values
        tbox{i}{j, 3} = mean(nonzeros(tbox{i}{j, 2}));   % get mean of non-zero t-values in each box
    end
end 

% END Part (IV): get average T-values of each 3X3X3 box
%--------------------------------------------------------------------------- 


% Part (V): sort spike rates by channel (pairs of electrode contacts)

% Part (VA): clean spike rates (.txt) file

% sometimes the file contains extra rows that aren't showing any electrode
% pairs, i.e. rows showing the name tags (e.g. LO1, L02, ...) instead of
% (L01 - L02), 
% select only rows that contain info. about electrode pairs (e.g. L01-L02)

% get number of rows in spikerates array
srni = size(spikesar, 1);   % length in i dimension, spike rates

% create empty array named spikerates_cleaned
spikesar_cleaned = {};

% assign the first row (header info., i.e. Channel, spike/min, ...)
% of input file to spikerates_cleaned
spikesar_cleaned(1,:) = spikesar(1,:);

% set row index as 1 for iteration in loop
rowind = 1;

% search str '-' in each row of input spike rate file, if found, assign all
% columns of that row to spikesar_cleaned array, then increment rownum by
% 1 and continue the loop
for i = 1:srni   % for each row
    if ~isempty(regexp(spikesar{i,1}, '-'))   % if str '-' is found in certain name tag (first col.)
        spikesar_cleaned(rowind,:) = spikesar(i, :);   % assign all info. to array
        rowind = rowind + 1;   % increment rownum by 1
    end
end

% PART (VB): match spike rates file with electrode pair names in tbox
% define search str as each name tag in tbox, find corresponding row number
% in spikesar_cleaned array, 
% if found, assign related info. to array,
% if not, assign value of 0 to that specific name tag
% note: assume tbox has NO MISSING electrode contacts

% initialize spike_sorted cell array
spike_matched = {};

for i = 1:ele_count   % for each type of electrode
    for j = 1:length(tbox{i}(:,1))   % for each pair of electrode contacts
        search_str = tbox{i}{j,1};   % search str accor. to name tag in tbox
        % get index (row number) of name tag in spikesar_cleaned
        if ~isempty(find(strcmp(search_str, spikesar_cleaned(:,1))))
            indreq = find(strcmp(search_str, spikesar_cleaned(:,1)));   % row index required
            % assign values to spike_matched
            spike_matched{i}(j,:) = spikesar_cleaned(indreq,:); 
        else
            % if specific electrode contacts (name tag) isn't found in
            % spikesar_cleaned array, assign value of 0 to that name
            spike_matched{i}(j,:) = {tbox{i}{j,1} [0] [0]};   
            sprintf('i=%d, j=%d, not found, but added 0 to name tag', i, j)            
        end
    end
end

% END Part (V): sort spike rates by channel (pairs of electrode contacts)
%--------------------------------------------------------------------------- 

% Part (VI): group data (t-values and spike rates) from all types of electrodes

% initialize arrays required
tbox_all = {};   % array listing all t-values (3X3X3 and avg.) found
spike_all = {};   % array listing all spike rates

% loop through each cell in arrays
for i = 1:ele_count
    tbox_all = [tbox_all; tbox{i}];
    spike_all = [spike_all; spike_matched{i}];
end

% check if order of name tags in tbox_all matches with spike_all 

% initialize error count as 0
errcount = 0;

% go through each entry of the two arrays, check if name tags (electorde
% pairs) match. For correlation analysis, the order at which the variables
% present in the arrays matters
for i = 1:size(tbox_all, 1)
    if strcmp(tbox_all{i,1}, spike_all{i,1}) ~= 1
        errcount = errcount + 1;
    end
end

% if error is found, print warning message
if errcount ~= 0
    sprintf(['Warning! Order at which eletrode pairs assigned in spike_matched', ...
        ' does not match with that in tbox'])
end

% convert avg. t-values and spike rates to numeric arrays
avgtvalues = cell2mat(tbox_all(:,3));
spikes = cell2mat(spike_all(:,2:3)); 

% get propotions of spikes with HFOs, store values to spikes matrix
% i.e. spikes+hfo/spikes
for i = 1:size(spikes,1)
    spikes(i, 3) = spikes(i, 2)/spikes(i, 1);
end

% if spike/min = 0, then spike+hfo/spike = inf,
% convert all entries of Inf to NaN, store values to spikes_re (spikes
% refined), keep original dimensions 
spikes_re = spikes;

for j = 1:size(spikes, 2)
    if find(isinf(spikes(:,j))) > 0        
        spikes_re(isinf(spikes(:,j)),j) = NaN;
    end
end

% END Part (VI): group data (t-values and spike rates) from all types of electrodes
%--------------------------------------------------------------------------- 

% Part (VII): construct structure to output data for plotting
opstruct.ele_sorted = ele_sorted;   % info. of electrode (sorted by name)
opstruct.tbox = tbox;   % cell array containing t-values (all and non-zero mean) of each electrode pair (with name tag)
opstruct.tbox_all = tbox_all;   % list of avg. t-values of all electode pairs
opstruct.avgtvalues = avgtvalues;   % avg. t-value of each electrode pair
opstruct.spike_matched = spike_matched;   % matched spike rates with name tag of tbox
opstruct.spike_all = spike_all;   % list of spike rates of all electrode pairs
opstruct.spikes_re = spikes_re;   % refined matrix of spike rates, inf replaced by NaN

opstruct.midpt = midpt;   % coordinates of midpoint between each pair of electrode contacts 

% END Part (VII): construct structre to ouput data
%--------------------------------------------------------------------------- 

end   % end function
%---------------------------------------------------------------------------





